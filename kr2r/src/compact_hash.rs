use byteorder::{ByteOrder, LittleEndian};
use memmap2::{Mmap, MmapOptions};
use std::cmp::Ordering as CmpOrdering;
use std::fs::OpenOptions;
use std::io::{Read, Result};
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
}

impl Slot<u64> {
    pub fn get_seq_id(&self) -> u64 {
        self.value.right(0xFFFFFFFF)
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

pub struct Page {
    pub index: usize,
    pub size: usize,
    pub data: Vec<u32>,
}

impl Page {
    pub fn new(index: usize, size: usize, data: Vec<u32>) -> Self {
        Self { index, size, data }
    }

    pub fn start(&self) -> usize {
        self.index * self.size
    }

    pub fn end(&self, capacity: usize) -> usize {
        std::cmp::min((self.index + 1) * self.size, capacity)
    }
}

pub struct PagePtr<'a> {
    #[allow(dead_code)]
    mmap: Option<Mmap>,
    pub index: usize,
    pub size: usize,
    pub data: &'a [u32],
}

impl<'a> PagePtr<'a> {
    pub fn new(mmap: Option<Mmap>, index: usize, size: usize, data: &'a [u32]) -> Self {
        Self {
            mmap,
            index,
            size,
            data,
        }
    }
}

impl<'a> Default for PagePtr<'a> {
    fn default() -> Self {
        // 创建一个包含0的静态数组
        static DEFAULT_DATA: [u32; 1] = [0];
        Self {
            mmap: None,
            index: 1,
            size: 1,
            data: &DEFAULT_DATA,
        }
    }
}

use std::fmt::{self, Debug};

#[derive(Clone, Copy)]
pub struct HashConfig {
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
    pub hash_capacity: usize,
}

// 为HashConfig手动实现Debug trait
impl fmt::Debug for HashConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CompactHashTableConfig")
            .field("value_mask", &self.value_mask)
            .field("value_bits", &self.value_bits)
            .field("capacity", &self.capacity)
            .field("size", &self.size)
            .field("hash_capacity", &self.hash_capacity)
            // 注意，我们没有包括_phantom字段
            .finish()
    }
}

impl HashConfig {
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
        hash_capacity: usize,
    ) -> Self {
        let value_mask = (1 << value_bits) - 1;
        Self {
            capacity,
            value_bits,
            value_mask,
            size,
            partition,
            hash_capacity,
        }
    }

    pub fn from_hash_header<P: AsRef<Path>>(filename: P) -> Result<Self> {
        let file = OpenOptions::new().read(true).open(&filename)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };

        let offset = 2 * Self::U64_SIZE;

        let partition = LittleEndian::read_u64(
            &mmap[Self::PARTITION_OFFSET..Self::PARTITION_OFFSET + Self::U64_SIZE],
        ) as usize;
        let hash_capacity = LittleEndian::read_u64(
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

        Ok(Self::new(
            capacity,
            value_bits,
            size,
            partition,
            hash_capacity,
        ))
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

    pub fn slot(&self, hash_key: u64, taxid: u32) -> Slot<u32> {
        let idx = self.index(hash_key);
        Slot::<u32>::new(idx, u32::hash_value(hash_key, self.value_bits, taxid))
    }

    pub fn slot_u64(&self, hash_key: u64, seq_id: u64) -> Slot<u64> {
        let idx = self.index(hash_key);
        Slot::<u64>::new(idx, u64::hash_value(hash_key, self.value_bits, seq_id))
    }
}

pub trait K2Compact: std::marker::Sync + Send {
    fn get_idx_mask(&self) -> usize;
    fn get_idx_bits(&self) -> usize;
    fn get_value_mask(&self) -> usize;
    fn get_value_bits(&self) -> usize;
    fn get_from_page(&self, idx: usize, value: u64, next: bool) -> u32;
}

#[allow(unused)]
pub struct CHPage<'a> {
    // 哈希表的容量
    pub config: HashConfig,
    pub page: Page,
    pub next_page: PagePtr<'a>,
}

fn _read_page_from_file<P: AsRef<Path>>(filename: P) -> Result<Page> {
    let file = OpenOptions::new().read(true).open(&filename)?;
    let mmap = unsafe { MmapOptions::new().populate().map(&file)? };
    let index = LittleEndian::read_u64(&mmap[0..8]) as usize;
    let capacity = LittleEndian::read_u64(&mmap[8..16]) as usize;

    let data = unsafe { std::slice::from_raw_parts(mmap.as_ptr().add(16) as *const u32, capacity) };

    // 初始化Vec<B>，预分配足够容量
    let mut page_data: Vec<u32> = Vec::with_capacity(capacity);
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

fn read_page_from_file<P: AsRef<Path>>(filename: P) -> Result<Page> {
    let mut file = std::fs::File::open(filename)?;

    // 读取索引和容量
    let mut buffer = [0u8; 16];
    file.read_exact(&mut buffer)?;

    let index = LittleEndian::read_u64(&buffer[0..8]) as usize;
    let capacity = LittleEndian::read_u64(&buffer[8..16]) as usize;

    // 读取数据部分
    let mut data = vec![0u32; capacity];
    let data_bytes = unsafe {
        std::slice::from_raw_parts_mut(
            data.as_mut_ptr() as *mut u8,
            capacity * std::mem::size_of::<u32>(),
        )
    };
    file.read_exact(data_bytes)?;

    Ok(Page::new(index, capacity, data))
}

fn _read_pageptr_from_file<'a, P: AsRef<Path>>(filename: P) -> Result<PagePtr<'a>> {
    let file = OpenOptions::new().read(true).open(&filename)?;
    let mmap = unsafe { MmapOptions::new().populate().map(&file)? };
    let index = LittleEndian::read_u64(&mmap[0..8]) as usize;
    let capacity = LittleEndian::read_u64(&mmap[8..16]) as usize;

    let raw_page_data =
        unsafe { std::slice::from_raw_parts(mmap.as_ptr().add(16) as *const u32, capacity) };

    let first_zero_end = raw_page_data
        .iter()
        .position(|&x| x == 0)
        .map(|pos| pos + 1)
        .unwrap_or(capacity);
    let page_data = &raw_page_data[..first_zero_end];

    // let page_data =
    //     unsafe { std::slice::from_raw_parts(mmap.as_ptr().add(16) as *const u32, capacity) };

    Ok(PagePtr::new(Some(mmap), index, first_zero_end, page_data))
}

fn read_pageptr_from_file_chunk<'a, P: AsRef<Path>>(filename: P) -> Result<PagePtr<'a>> {
    let file = OpenOptions::new().read(true).open(&filename)?;
    let mmap = unsafe { MmapOptions::new().populate().map(&file)? };
    let index = LittleEndian::read_u64(&mmap[0..8]) as usize;
    let capacity = LittleEndian::read_u64(&mmap[8..16]) as usize;

    let mut first_zero_end = capacity;
    let chunk_size = 1024; // 定义每次读取的块大小
    let mut found_zero = false;

    for i in (0..capacity).step_by(chunk_size) {
        let end = usize::min(i + chunk_size, capacity);
        let chunk = unsafe {
            std::slice::from_raw_parts(mmap.as_ptr().add(16 + i * 4) as *const u32, end - i)
        };
        if let Some(pos) = chunk.iter().position(|&x| x == 0) {
            first_zero_end = i + pos + 1;
            found_zero = true;
            break;
        }
    }

    if !found_zero {
        first_zero_end = capacity;
    }

    let page_data =
        unsafe { std::slice::from_raw_parts(mmap.as_ptr().add(16) as *const u32, first_zero_end) };

    Ok(PagePtr::new(Some(mmap), index, first_zero_end, page_data))
}

impl<'a> CHPage<'a> {
    pub fn from<P: AsRef<Path> + Debug>(
        config: HashConfig,
        chunk_file1: P,
        chunk_file2: P,
    ) -> Result<CHPage<'a>> {
        let page = read_page_from_file(chunk_file1)?;
        let next_page = if page.data.last().map_or(false, |&x| x == 0) {
            read_pageptr_from_file_chunk(chunk_file2)?
        } else {
            PagePtr::default()
        };

        let chtm = CHPage {
            config,
            next_page,
            page,
        };
        Ok(chtm)
    }

    pub fn get_from_next_page(&self, index: usize, compacted_key: u32) -> u32 {
        let value_mask = self.config.value_mask;
        let mut idx = index;

        loop {
            if let Some(cell) = self.next_page.data.get(idx) {
                if cell.right(value_mask) == u32::default()
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
                return u32::default();
            }
        }
        u32::default()
    }
}

impl<'a> K2Compact for CHPage<'a> {
    fn get_idx_mask(&self) -> usize {
        let idx_bits = ((self.config.hash_capacity as f64).log2().ceil() as usize).max(1);
        (1 << idx_bits) - 1
    }

    fn get_idx_bits(&self) -> usize {
        ((self.config.hash_capacity as f64).log2().ceil() as usize).max(1)
    }

    fn get_value_mask(&self) -> usize {
        self.config.value_mask
    }

    fn get_value_bits(&self) -> usize {
        self.config.value_bits
    }

    fn get_from_page(&self, indx: usize, value: u64, next: bool) -> u32 {
        let compacted_key = value.left(self.config.value_bits) as u32;
        let value_mask = self.config.value_mask;
        let mut idx = indx;
        let first_idx = idx;

        loop {
            if let Some(cell) = self.page.data.get(idx) {
                if cell.right(value_mask) == u32::default()
                    || cell.left(self.config.value_bits) == compacted_key
                {
                    return cell.right(value_mask);
                }

                if next {
                    idx = idx + 1;
                    if idx >= self.page.size {
                        // 需要确定在table中的位置, page index 从0开始
                        let index = idx % self.page.size;
                        return self.get_from_next_page(index, compacted_key);
                    }
                } else {
                    idx = (idx + 1) % self.page.size;
                }

                if idx == first_idx {
                    break;
                }
            } else {
                // 如果get(idx)失败，返回默认值
                return u32::default();
            }
        }
        u32::default()
    }
}
