use byteorder::{ByteOrder, LittleEndian};
use memmap2::{Mmap, MmapMut, MmapOptions};
use std::cmp::Ordering as CmpOrdering;
use std::fs::OpenOptions;
use std::io::{Error, ErrorKind, Result};
use std::marker::PhantomData;
use std::path::Path;

pub trait CompactValue: Default + PartialEq + Clone + Copy + PartialEq + Eq + Sized {
    fn compacted(hash_key: u64, value_bits: usize) -> Self;
    fn hash_value(hash_key: u64, value_bits: usize, value: Self) -> Self;
}

/// 1101010101 => left: 11010, right: 10101;
pub trait BitN: CompactValue {
    fn left(&self, value_bits: usize) -> Self;
    fn right(&self, value_mask: u32) -> Self;
    fn combined(left: Self, right: Self, value_bits: usize) -> Self;
}

impl CompactValue for u32 {
    fn hash_value(hash_key: u64, value_bits: usize, value: u32) -> u32 {
        Self::compacted(hash_key, value_bits) << value_bits | value
    }
    fn compacted(value: u64, value_bits: usize) -> u32 {
        (value >> (32 + value_bits)) as u32
    }
}

impl BitN for u32 {
    fn left(&self, value_bits: usize) -> u32 {
        *self >> value_bits
    }

    fn right(&self, value_mask: u32) -> u32 {
        // 为 u32 类型实现 right 方法
        *self & value_mask
    }
    fn combined(left: Self, right: Self, value_bits: usize) -> Self {
        left << value_bits | right
    }
}

#[allow(unused_variables)]
impl CompactValue for bool {
    fn hash_value(hash_key: u64, value_bits: usize, value: bool) -> bool {
        value
    }
    fn compacted(value: u64, value_bits: usize) -> bool {
        true
    }
}

#[allow(unused_variables)]
impl BitN for bool {
    fn left(&self, value_bits: usize) -> bool {
        true
    }

    fn right(&self, value_mask: u32) -> bool {
        // 为 u32 类型实现 right 方法
        *self
    }
    fn combined(left: Self, right: Self, value_bits: usize) -> Self {
        right
    }
}

#[repr(C)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Slot<B>
where
    B: BitN,
{
    pub idx: usize,
    pub value: B,
}

impl<B> Slot<B>
where
    B: BitN,
{
    pub fn new(idx: usize, value: B) -> Self {
        Self { idx, value }
    }

    pub fn update_right(&mut self, right: B, value_bits: usize) {
        let left = self.value.left(value_bits);
        self.value = B::combined(left, right, value_bits);
    }
}

// 实现 PartialOrd，只比较 index 字段
impl<B> PartialOrd for Slot<B>
where
    B: BitN,
{
    fn partial_cmp(&self, other: &Self) -> Option<CmpOrdering> {
        self.idx.partial_cmp(&other.idx)
    }
}

// 实现 Ord，只比较 index 字段
impl<B> Ord for Slot<B>
where
    B: BitN,
{
    fn cmp(&self, other: &Self) -> CmpOrdering {
        self.idx.cmp(&other.idx)
    }
}

/// 与Slot的区别,只是idx的类型
#[repr(C)]
pub struct Cell<B>
where
    B: BitN,
{
    pub idx: u32,
    value: B,
}

impl<B> Cell<B>
where
    B: BitN,
{
    pub fn new(idx: u32, value: B) -> Self {
        Self { idx, value }
    }

    pub fn as_slice(&self) -> &[u8] {
        let cell_ptr = self as *const Self as *const u8;
        let cell_size = std::mem::size_of::<Self>();
        let cell_bytes = unsafe {
            // 将Cell实例的内存表示转换为字节切片
            std::slice::from_raw_parts(cell_ptr, cell_size)
        };
        cell_bytes
    }

    pub fn as_slot(&self) -> Slot<B> {
        Slot::new(self.idx as usize, self.value)
    }
}

pub struct Page<B>
where
    B: BitN,
{
    pub index: usize,
    pub size: usize,
    pub data: Vec<B>,
}

impl<B> Page<B>
where
    B: BitN,
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

#[derive(Debug, Clone, Copy)]
pub struct HashConfig<B>
where
    B: BitN,
{
    // value_mask = ((1 << value_bits) - 1);
    pub value_mask: u32,
    // 值的位数
    pub value_bits: usize,
    // 哈希表的容量
    pub capacity: usize,
    // 哈希表中当前存储的元素数量。
    pub size: usize,
    _phantom: PhantomData<B>,
}

impl<B> HashConfig<B>
where
    B: BitN,
{
    // 使用常量替代硬编码的数字，增加代码可读性
    const CAPACITY_OFFSET: usize = 0;
    const SIZE_OFFSET: usize = 8;
    const VALUE_BITS_OFFSET: usize = 24;
    const U64_SIZE: usize = 8;

    pub fn new(capacity: usize, value_bits: usize, size: usize) -> Self {
        let value_mask = (1 << value_bits) - 1;
        Self {
            capacity,
            value_bits,
            value_mask,
            size,
            _phantom: PhantomData,
        }
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

        Self::new(capacity, value_bits, size)
    }

    pub fn index(&self, hash_key: u64) -> usize {
        hash_key as usize % self.capacity
    }

    pub fn slot(&self, hash_key: u64, taxid: B) -> Slot<B> {
        let idx = self.index(hash_key);
        Slot::<B>::new(idx, B::hash_value(hash_key, self.value_bits, taxid))
    }
}

// macro_rules! define_page {
//     ($struct_name:ident, $type:ty) => {
//         pub struct $struct_name<B>
//         where
//             B: BitN<$type>,
//         {
//             pub index: usize,
//             pub size: usize,
//             pub data: Vec<B>,
//         }

//         impl<B: BitN<$type>> $struct_name<B> {
//             pub fn new(index: usize, size: usize, data: Vec<B>) -> Self {
//                 Self { index, size, data }
//             }
//         }
//     };
// }

// define_page!(BoolPage, bool);
// define_page!(U32Page, u32);

#[allow(unused)]
pub struct CHTable<'a, B>
where
    B: BitN + 'a,
{
    // memmap
    mmap: Mmap,
    // 哈希表的容量
    pub config: HashConfig<B>,
    pub table: &'a [B],
}

impl<'a, B> CHTable<'a, B>
where
    B: BitN + 'a,
{
    pub fn from<P: AsRef<Path>>(filename: P) -> Result<CHTable<'a, B>> {
        let file = OpenOptions::new().read(true).open(&filename)?;

        let mmap = unsafe { MmapOptions::new().populate().map(&file)? };
        let config = HashConfig::from_mmap(&mmap);
        let table = unsafe {
            std::slice::from_raw_parts(mmap.as_ptr().add(32) as *const B, config.capacity)
        };

        let chtm = CHTable {
            config,
            table,
            mmap,
        };
        Ok(chtm)
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

#[allow(unused)]
pub struct CHTableMut<'a, B>
where
    B: BitN + 'a,
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
    B: BitN + 'a,
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
    B: BitN + 'a,
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

    pub fn set_page_cell(&mut self, item: Slot<B>) -> Option<Slot<B>> {
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
                    return Some(Slot::<B>::new(idx, cell.clone()));
                }

                idx = idx + 1;
                if idx >= self.page.size {
                    return self.set_table_cell(idx, item.value);
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
    pub fn update_cell(&mut self, item: Slot<B>) {
        if item.idx < self.page.size {
            self.page.data[item.idx] = item.value;
        } else {
            self.table[item.idx] = item.value;
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
