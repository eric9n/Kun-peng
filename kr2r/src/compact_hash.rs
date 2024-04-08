use byteorder::{ByteOrder, LittleEndian};
use memmap2::{Mmap, MmapMut, MmapOptions};
use std::cmp::Ordering as CmpOrdering;
use std::fs::OpenOptions;
use std::io::{Error, ErrorKind, Result};
use std::marker::PhantomData;
use std::path::Path;

/// 1101010101 => left: 11010, right: 10101;
pub trait Compact: Default + PartialEq + Clone + Copy + Eq + Sized + Send + Sync + Debug {
    fn compacted(hash_key: u64, value_bits: usize) -> Self;
    fn hash_value(hash_key: u64, value_bits: usize, value: Self) -> Self;

    fn left(&self, value_bits: usize) -> Self;
    fn right(&self, value_mask: usize) -> Self;
    fn combined(left: Self, right: Self, value_bits: usize) -> Self;
    fn to_u32(&self) -> u32;
    fn from_u32(value: u32) -> Self;
}

impl Compact for u32 {
    fn hash_value(hash_key: u64, value_bits: usize, value: u32) -> u32 {
        Self::compacted(hash_key, value_bits) << value_bits | value
    }
    fn compacted(value: u64, value_bits: usize) -> u32 {
        (value >> (32 + value_bits)) as u32
    }

    fn left(&self, value_bits: usize) -> u32 {
        *self >> value_bits
    }

    fn right(&self, value_mask: usize) -> u32 {
        *self & value_mask as u32
    }
    fn combined(left: Self, right: Self, value_bits: usize) -> Self {
        left << value_bits | right
    }

    fn to_u32(&self) -> u32 {
        *self
    }
    fn from_u32(value: u32) -> Self {
        value
    }
}

impl Compact for u64 {
    fn hash_value(hash_key: u64, value_bits: usize, value: u64) -> u64 {
        Self::compacted(hash_key, value_bits) << (32 + value_bits) | value
    }
    fn compacted(value: u64, value_bits: usize) -> u64 {
        (value >> (32 + value_bits)) as u64
    }

    fn left(&self, value_bits: usize) -> u64 {
        *self >> (32 + value_bits)
    }

    fn right(&self, value_mask: usize) -> u64 {
        let mask: u64 = ((value_mask as u64) << 32) | 0xFFFFFFFF;
        mask & *self
    }

    fn combined(left: Self, right: Self, value_bits: usize) -> Self {
        left << (32 + value_bits) | right
    }

    fn to_u32(&self) -> u32 {
        *self as u32
    }
    fn from_u32(value: u32) -> Self {
        value as u64
    }
}

#[repr(C)]
#[derive(PartialEq, Clone, Copy, Eq, Debug)]
pub struct Row {
    pub value: u32,
    pub seq_id: u32,
    pub kmer_id: u32,
}

impl Row {
    pub fn new(value: u32, seq_id: u32, kmer_id: u32) -> Self {
        Self {
            value,
            seq_id,
            kmer_id,
        }
    }
    #[inline]
    pub fn as_slice(&self, row_size: usize) -> &[u8] {
        let slot_ptr = self as *const Self as *const u8;
        unsafe { std::slice::from_raw_parts(slot_ptr, row_size) }
    }
}

// 实现 PartialOrd，只比较 index 字段
impl PartialOrd for Row {
    fn partial_cmp(&self, other: &Self) -> Option<CmpOrdering> {
        self.kmer_id.partial_cmp(&other.kmer_id)
    }
}

// 实现 Ord，只比较 index 字段
impl Ord for Row {
    fn cmp(&self, other: &Self) -> CmpOrdering {
        self.kmer_id.cmp(&other.kmer_id)
    }
}

#[repr(C)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Slot<B>
where
    B: Compact,
{
    pub idx: usize,
    pub value: B,
}

impl<B> Slot<B>
where
    B: Compact,
{
    pub fn new(idx: usize, value: B) -> Self {
        Self { idx, value }
    }

    pub fn update_right(&mut self, right: B, value_bits: usize) {
        let left = self.value.left(value_bits);
        self.value = B::combined(left, right, value_bits);
    }

    #[inline]
    pub fn as_slice(&self, slot_size: usize) -> &[u8] {
        let slot_ptr = self as *const Self as *const u8;
        unsafe { std::slice::from_raw_parts(slot_ptr, slot_size) }
    }

    pub fn get_seq_id(&self) -> B {
        self.value.right(0xFFFFFFFF)
    }

    pub fn to_b(&self, left: B) -> B {
        B::combined(left, self.value.right(0), 0)
    }
}

// 实现 PartialOrd，只比较 index 字段
impl<B> PartialOrd for Slot<B>
where
    B: Compact,
{
    fn partial_cmp(&self, other: &Self) -> Option<CmpOrdering> {
        self.idx.partial_cmp(&other.idx)
    }
}

// 实现 Ord，只比较 index 字段
impl<B> Ord for Slot<B>
where
    B: Compact,
{
    fn cmp(&self, other: &Self) -> CmpOrdering {
        self.idx.cmp(&other.idx)
    }
}

/// 与Slot的区别,只是idx的类型
#[repr(C)]
pub struct Cell<B>
where
    B: Compact,
{
    pub idx: u32,
    value: B,
}

impl<B> Cell<B>
where
    B: Compact,
{
    pub fn new(idx: u32, value: B) -> Self {
        Self { idx, value }
    }

    pub fn as_slice(&self, cell_size: usize) -> &[u8] {
        let cell_ptr = self as *const Self as *const u8;
        unsafe { std::slice::from_raw_parts(cell_ptr, cell_size) }
    }

    pub fn as_slot(&self) -> Slot<B> {
        Slot::new(self.idx as usize, self.value)
    }
}

pub struct Page<B>
where
    B: Compact,
{
    pub index: usize,
    pub size: usize,
    pub data: Vec<B>,
}

impl<B> Page<B>
where
    B: Compact,
{
    pub fn new(index: usize, size: usize, data: Vec<B>) -> Self {
        Self { index, size, data }
    }

    pub fn start(&self) -> usize {
        self.index * self.size
    }

    pub fn end(&self, capacity: usize) -> usize {
        std::cmp::min((self.index + 1) * self.size, capacity)
    }
}

pub struct PagePtr<'a, B>
where
    B: Compact + 'a,
{
    #[allow(dead_code)]
    mmap: Mmap,
    pub index: usize,
    pub size: usize,
    pub data: &'a [B],
}

impl<'a, B> PagePtr<'a, B>
where
    B: Compact + 'a,
{
    pub fn new(mmap: Mmap, index: usize, size: usize, data: &'a [B]) -> Self {
        Self {
            mmap,
            index,
            size,
            data,
        }
    }
}

use std::fmt::{self, Debug};

#[derive(Clone, Copy)]
pub struct HashConfig<B>
where
    B: Compact,
{
    // value_mask = ((1 << value_bits) - 1);
    pub value_mask: usize,
    // 值的位数
    pub value_bits: usize,
    // 哈希表的容量
    pub capacity: usize,
    // 哈希表中当前存储的元素数量。
    pub size: usize,
    // 分区数
    pub partition: usize,
    // 分块大小
    pub hash_size: usize,
    _phantom: PhantomData<B>,
}

// 为HashConfig手动实现Debug trait
impl<B> fmt::Debug for HashConfig<B>
where
    B: Compact,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CompactHashTableConfig")
            .field("value_mask", &self.value_mask)
            .field("value_bits", &self.value_bits)
            .field("capacity", &self.capacity)
            .field("size", &self.size)
            // 注意，我们没有包括_phantom字段
            .finish()
    }
}

impl<B> HashConfig<B>
where
    B: Compact,
{
    // 使用常量替代硬编码的数字，增加代码可读性
    const PARTITION_OFFSET: usize = 0;
    const HASH_SIZE_OFFSET: usize = 8;
    const CAPACITY_OFFSET: usize = 0;
    const SIZE_OFFSET: usize = 8;
    const VALUE_BITS_OFFSET: usize = 24;
    const U64_SIZE: usize = 8;

    pub fn new(
        capacity: usize,
        value_bits: usize,
        size: usize,
        partition: usize,
        hash_size: usize,
    ) -> Self {
        let value_mask = (1 << value_bits) - 1;
        Self {
            capacity,
            value_bits,
            value_mask,
            size,
            partition,
            hash_size,
            _phantom: PhantomData,
        }
    }

    pub fn from_hash_header<P: AsRef<Path>>(filename: P) -> Result<Self> {
        let file = OpenOptions::new().read(true).open(&filename)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };

        let offset = 2 * Self::U64_SIZE;

        let partition = LittleEndian::read_u64(
            &mmap[Self::PARTITION_OFFSET..Self::PARTITION_OFFSET + Self::U64_SIZE],
        ) as usize;
        let hash_size = LittleEndian::read_u64(
            &mmap[Self::HASH_SIZE_OFFSET..Self::HASH_SIZE_OFFSET + Self::U64_SIZE],
        ) as usize;
        let capacity = LittleEndian::read_u64(
            &mmap[Self::CAPACITY_OFFSET + offset..Self::CAPACITY_OFFSET + offset + Self::U64_SIZE],
        ) as usize;
        let size = LittleEndian::read_u64(
            &mmap[Self::SIZE_OFFSET + offset..Self::SIZE_OFFSET + offset + Self::U64_SIZE],
        ) as usize;
        let value_bits = LittleEndian::read_u64(
            &mmap[Self::VALUE_BITS_OFFSET + offset
                ..Self::VALUE_BITS_OFFSET + offset + Self::U64_SIZE],
        ) as usize;

        Ok(Self::new(capacity, value_bits, size, partition, hash_size))
    }

    pub fn from<P: AsRef<Path>>(filename: P) -> Result<Self> {
        let file = OpenOptions::new().read(true).open(&filename)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        Ok(Self::from_mmap(&mmap))
    }

    pub fn from_mmap(mmap: &Mmap) -> Self {
        let capacity = LittleEndian::read_u64(
            &mmap[Self::CAPACITY_OFFSET..Self::CAPACITY_OFFSET + Self::U64_SIZE],
        ) as usize;
        let size =
            LittleEndian::read_u64(&mmap[Self::SIZE_OFFSET..Self::SIZE_OFFSET + Self::U64_SIZE])
                as usize;
        let value_bits = LittleEndian::read_u64(
            &mmap[Self::VALUE_BITS_OFFSET..Self::VALUE_BITS_OFFSET + Self::U64_SIZE],
        ) as usize;

        Self::new(capacity, value_bits, size, 0, 0)
    }

    pub fn index(&self, hash_key: u64) -> usize {
        hash_key as usize % self.capacity
    }

    pub fn slot(&self, hash_key: u64, taxid: B) -> Slot<B> {
        let idx = self.index(hash_key);
        Slot::<B>::new(idx, B::hash_value(hash_key, self.value_bits, taxid))
    }

    pub fn slot_u64(&self, hash_key: u64, seq_id: u64) -> Slot<u64> {
        let idx = self.index(hash_key);
        Slot::<u64>::new(idx, u64::hash_value(hash_key, self.value_bits, seq_id))
    }
}

pub trait K2Compact<B>: std::marker::Sync + Send
where
    B: Compact,
{
    fn get_idx_mask(&self) -> usize;
    fn get_idx_bits(&self) -> usize;
    fn get_value_mask(&self) -> usize;
    fn get_value_bits(&self) -> usize;
    fn get_from_page(&self, idx: usize, value: u64) -> B;
}

#[allow(unused)]
pub struct CHPage<'a, B>
where
    B: Compact,
{
    // 哈希表的容量
    pub config: HashConfig<B>,
    pub page: Page<B>,
    pub next_page: PagePtr<'a, B>,
}

fn read_page_from_file<P: AsRef<Path>, B: Compact>(filename: P) -> Result<Page<B>> {
    let file = OpenOptions::new().read(true).open(&filename)?;
    let mmap = unsafe { MmapOptions::new().populate().map(&file)? };
    let index = LittleEndian::read_u64(&mmap[0..8]) as usize;
    let capacity = LittleEndian::read_u64(&mmap[8..16]) as usize;

    let data = unsafe { std::slice::from_raw_parts(mmap.as_ptr().add(16) as *const B, capacity) };

    // 初始化Vec<B>，预分配足够容量
    let mut page_data: Vec<B> = Vec::with_capacity(capacity);
    // 为Vec<B>安全地设置长度
    unsafe {
        page_data.set_len(capacity);
    }
    // 使用copy_from_slice进行数据复制
    unsafe {
        let page_data_slice = std::slice::from_raw_parts_mut(page_data.as_mut_ptr(), capacity);
        page_data_slice.copy_from_slice(data);
    }

    // let page_data = data.to_vec();
    Ok(Page::new(index, capacity, page_data))
}

fn read_pageptr_from_file<'a, P: AsRef<Path>, B: Compact>(filename: P) -> Result<PagePtr<'a, B>> {
    let file = OpenOptions::new().read(true).open(&filename)?;
    let mmap = unsafe { MmapOptions::new().populate().map(&file)? };
    let index = LittleEndian::read_u64(&mmap[0..8]) as usize;
    let capacity = LittleEndian::read_u64(&mmap[8..16]) as usize;

    let page_data =
        unsafe { std::slice::from_raw_parts(mmap.as_ptr().add(16) as *const B, capacity) };

    Ok(PagePtr::new(mmap, index, capacity, page_data))
}

impl<'a, B> CHPage<'a, B>
where
    B: Compact + 'a,
{
    pub fn from<P: AsRef<Path> + Debug>(
        config: HashConfig<B>,
        chunk_file1: P,
        chunk_file2: P,
    ) -> Result<CHPage<'a, B>> {
        // let file2 = OpenOptions::new().read(true).open(&chunk_file2)?;
        // let mmap = unsafe { MmapOptions::new().map(&file2)? };
        // let index2 = LittleEndian::read_u64(&mmap[0..8]) as usize;
        // let capacity = LittleEndian::read_u64(&mmap[8..16]) as usize;

        // let next_page =
        //     unsafe { std::slice::from_raw_parts(mmap.as_ptr().add(16) as *const B, capacity) };

        let page = read_page_from_file(chunk_file1)?;
        let next_page = read_pageptr_from_file(chunk_file2)?;

        // let file1 = OpenOptions::new().read(true).open(&chunk_file1)?;
        // let mmap1 = unsafe { MmapOptions::new().map(&file1)? };

        // let index1 = LittleEndian::read_u64(&mmap1[0..8]) as usize;
        // let capacity = LittleEndian::read_u64(&mmap1[8..16]) as usize;

        // let table =
        //     unsafe { std::slice::from_raw_parts(mmap.as_ptr().add(16) as *const B, capacity) };
        // let page_data = table.to_vec();
        // let page = Page::<B>::new(index1, capacity, page_data);

        let chtm = CHPage {
            config,
            next_page,
            page,
        };
        Ok(chtm)
    }

    pub fn get_from_next_page(&self, index: usize, compacted_key: B) -> B {
        let value_mask = self.config.value_mask;
        let mut idx = index;

        loop {
            if let Some(cell) = self.next_page.data.get(idx) {
                if cell.right(value_mask) == B::default()
                    || cell.left(self.config.value_bits) == compacted_key
                {
                    return cell.right(value_mask);
                }

                idx = idx + 1;
                if idx >= self.next_page.size {
                    break;
                }
            } else {
                // 如果get(idx)失败，返回默认值
                return B::default();
            }
        }
        B::default()
    }
}

impl<'a, B> K2Compact<B> for CHPage<'a, B>
where
    B: Compact,
{
    fn get_idx_mask(&self) -> usize {
        let idx_bits = ((self.config.hash_size as f64).log2().ceil() as usize).max(1);
        (1 << idx_bits) - 1
    }

    fn get_idx_bits(&self) -> usize {
        ((self.config.hash_size as f64).log2().ceil() as usize).max(1)
    }

    fn get_value_mask(&self) -> usize {
        self.config.value_mask
    }

    fn get_value_bits(&self) -> usize {
        self.config.value_bits
    }

    fn get_from_page(&self, indx: usize, value: u64) -> B {
        let compacted_key = B::from_u32(value.left(self.config.value_bits) as u32);
        let value_mask = self.config.value_mask;
        let mut idx = indx;
        let first_idx = idx;

        loop {
            if let Some(cell) = self.page.data.get(idx) {
                if cell.right(value_mask) == B::default()
                    || cell.left(self.config.value_bits) == compacted_key
                {
                    return cell.right(value_mask);
                }

                idx = idx + 1;
                if idx >= self.page.size {
                    // 需要确定在table中的位置, page index 从0开始
                    let index = idx % self.page.size;
                    return self.get_from_next_page(index, compacted_key);
                }
                if idx == first_idx {
                    break;
                }
            } else {
                // 如果get(idx)失败，返回默认值
                return B::default();
            }
        }
        B::default()
    }
}

#[allow(unused)]
pub struct CHTable<'a, B>
where
    B: Compact + 'a,
{
    // memmap
    mmap: Mmap,
    // 哈希表的容量
    pub config: HashConfig<B>,
    pub table: &'a [B],
    pub page: Page<B>,
}

impl<'a, B> CHTable<'a, B>
where
    B: Compact + 'a,
{
    pub fn from<P: AsRef<Path>>(
        filename: P,
        page_index: usize,
        page_size: usize,
    ) -> Result<CHTable<'a, B>> {
        let file = OpenOptions::new().read(true).open(&filename)?;

        let mmap = unsafe { MmapOptions::new().populate().map(&file)? };
        let config = HashConfig::from_mmap(&mmap);
        let table = unsafe {
            std::slice::from_raw_parts(mmap.as_ptr().add(32) as *const B, config.capacity)
        };

        let start_index = page_index * page_size;
        let end_index = std::cmp::min((page_index + 1) * page_size, config.capacity);
        if start_index > config.capacity {
            return Err(Error::new(ErrorKind::Other, "out of capacity"));
        }

        let page_data = table[start_index..end_index].to_vec();
        let page = Page::<B>::new(page_index, page_size, page_data);

        let chtm = CHTable {
            config,
            table,
            mmap,
            page,
        };
        Ok(chtm)
    }

    pub fn get_none_counts(&self) -> usize {
        self.table
            .iter()
            .filter(|&&item| item == B::default())
            .count()
    }

    pub fn get_from_table(&self, index: usize, compacted_key: B) -> B {
        let value_mask = self.config.value_mask;
        let mut idx = index;
        let first_idx = idx;

        loop {
            if let Some(cell) = self.table.get(idx) {
                if cell.right(value_mask) == B::default()
                    || cell.left(self.config.value_bits) == compacted_key
                {
                    return cell.right(value_mask);
                }

                idx = (idx + 1) % self.config.capacity;
                if idx == first_idx {
                    break;
                }
            } else {
                // 如果get(idx)失败，返回默认值
                return B::default();
            }
        }
        B::default()
    }

    pub fn get(&self, hash_key: u64) -> B {
        let compacted_key = B::compacted(hash_key, self.config.value_bits);
        let value_mask = self.config.value_mask;
        let mut idx = self.config.index(hash_key);
        let first_idx = idx;

        loop {
            if let Some(cell) = self.table.get(idx) {
                if cell.right(value_mask) == B::default()
                    || cell.left(self.config.value_bits) == compacted_key
                {
                    return cell.right(value_mask);
                }
            } else {
                // 如果get(idx)失败，返回默认值
                return B::default();
            }

            idx = (idx + 1) % self.config.capacity;
            if idx == first_idx {
                break;
            }
        }
        B::default()
    }
}

impl<'a, B> K2Compact<B> for CHTable<'a, B>
where
    B: Compact + 'a,
{
    fn get_idx_mask(&self) -> usize {
        let idx_bits = ((self.config.hash_size as f64).log2().ceil() as usize).max(1);
        (1 << idx_bits) - 1
    }

    fn get_idx_bits(&self) -> usize {
        ((self.config.hash_size as f64).log2().ceil() as usize).max(1)
    }

    fn get_value_mask(&self) -> usize {
        self.config.value_mask
    }

    fn get_value_bits(&self) -> usize {
        self.config.value_bits
    }

    fn get_from_page(&self, indx: usize, value: u64) -> B {
        let compacted_key = B::from_u32(value.left(self.config.value_bits) as u32);
        let value_mask = self.config.value_mask;
        let mut idx = indx;
        let first_idx = idx;

        loop {
            if let Some(cell) = self.page.data.get(idx) {
                if cell.right(value_mask) == B::default()
                    || cell.left(self.config.value_bits) == compacted_key
                {
                    return cell.right(value_mask);
                }

                idx = idx + 1;
                if idx >= self.page.size {
                    // 需要确定在table中的位置, page index 从0开始
                    let index = self.page.size * self.page.index + idx;
                    return self.get_from_table(index, compacted_key);
                }
                if idx == first_idx {
                    break;
                }
            } else {
                // 如果get(idx)失败，返回默认值
                return B::default();
            }
        }
        B::default()
    }
}

#[allow(unused)]
pub struct CHTableMut<'a, B>
where
    B: Compact + 'a,
{
    // memmap
    mmap: MmapMut,
    // 哈希表的容量
    pub config: HashConfig<B>,
    pub table: &'a mut [B],
    pub page: Page<B>,
}

impl<'a, B> CHTableMut<'a, B>
where
    B: Compact + 'a,
{
    pub fn new<P: AsRef<Path>>(
        hash_file: P,
        config: HashConfig<B>,
        page_index: usize,
        page_size: usize,
    ) -> Result<CHTableMut<'a, B>> {
        let key_bits = 32 - config.value_bits;

        let file_len = 32 + std::mem::size_of::<B>() * config.capacity;
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .open(&hash_file)?;
        file.set_len(file_len as u64)?;

        let mut mut_mmap = unsafe { MmapOptions::new().len(file_len).map_mut(&file)? };

        mut_mmap[0..8].copy_from_slice(&config.capacity.to_le_bytes());
        // 不能直接写入size,有可能文件已经存在,原来的size被覆盖
        // mut_mmap[8..16].copy_from_slice(&config.size.to_le_bytes());
        mut_mmap[16..24].copy_from_slice(&key_bits.to_le_bytes());
        mut_mmap[24..32].copy_from_slice(&config.value_bits.to_le_bytes());

        let table = unsafe {
            std::slice::from_raw_parts_mut(mut_mmap.as_mut_ptr().add(32) as *mut B, config.capacity)
        };
        let start_index = page_index * page_size;
        let end_index = std::cmp::min((page_index + 1) * page_size, config.capacity);
        if start_index > config.capacity {
            return Err(Error::new(ErrorKind::Other, "out of capacity"));
        }

        let page_data: Vec<B> = table[start_index..end_index].to_vec();
        let page = Page::<B>::new(page_index, page_size, page_data);
        let chtm = Self {
            mmap: mut_mmap,
            config,
            page,
            table,
        };
        Ok(chtm)
    }
}

impl<'a, B> CHTableMut<'a, B>
where
    B: Compact + 'a,
{
    pub fn set_table_cell(&mut self, index: usize, value: B) -> Option<Slot<B>> {
        let mut idx = index;
        let first_idx = idx;
        let value_bits = self.config.value_bits; // 局部变量存储配置
        let value_mask = self.config.value_mask;
        let capacity = self.config.capacity;

        loop {
            if let Some(cell) = self.table.get_mut(idx) {
                if cell.right(value_mask) == B::default() {
                    *cell = value;
                    break;
                }

                if cell.left(value_bits) == value.left(value_bits) {
                    return Some(Slot::<B>::new(idx, cell.clone()));
                }

                idx = (idx + 1) % capacity; // 循环索引，避免超出范围
                if idx == first_idx {
                    break;
                }
            } else {
                return None;
            }
        }

        None
    }

    pub fn set_page_cell(&mut self, item: Slot<B>) -> Option<(usize, Slot<B>)> {
        let mut idx = item.idx;
        let first_idx = idx;
        let value_bits = self.config.value_bits; // 局部变量存储配置
        let value_mask = self.config.value_mask;

        loop {
            if let Some(cell) = self.page.data.get_mut(idx) {
                if cell.right(value_mask) == B::default() {
                    *cell = item.value;
                    break;
                }

                if cell.left(value_bits) == item.value.left(value_bits) {
                    return Some((0, Slot::<B>::new(idx, cell.clone())));
                }

                idx = idx + 1;
                if idx >= self.page.size {
                    // 需要确定在table中的位置
                    let index = self.page.size * self.page.index + idx;
                    match self.set_table_cell(index, item.value) {
                        None => return None,
                        Some(s) => return Some((1, s)),
                    }
                }
                if idx == first_idx {
                    break;
                }
            } else {
                return None;
            }
        }

        None
    }

    pub fn copy_from_page(&mut self) {
        self.table[self.page.start()..self.page.end(self.config.capacity)]
            .copy_from_slice(&self.page.data);
    }

    // 直接更新
    pub fn update_cell(&mut self, flag: &usize, item: Slot<B>) {
        // 需要确定在table中的位置
        match flag {
            0 => {
                self.page.data[item.idx] = item.value;
            }
            _ => {
                self.table[item.idx] = item.value;
            }
        }
    }

    pub fn update_size(&self, size: usize) {
        let size_bytes = size.to_le_bytes(); // 假设使用小端序
        unsafe {
            let size_ptr = self.mmap.as_ptr().add(8) as *mut u8;
            std::ptr::copy_nonoverlapping(size_bytes.as_ptr(), size_ptr, size_bytes.len());
        }
        self.mmap.flush().expect("Failed to flush mmap");
    }

    /// 在原始size的基础上再加上size
    pub fn add_size(&mut self, size: usize) {
        let old_size = LittleEndian::read_u64(&self.mmap[8..16]) as usize;
        let new_size = old_size + size;
        self.update_size(new_size);
    }
}
