use byteorder::{ByteOrder, LittleEndian, ReadBytesExt};
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
        f.debug_struct("HashConfig")
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
        let mut file = OpenOptions::new().read(true).open(&filename)?;
        let partition = file.read_u64::<LittleEndian>()? as usize;
        let hash_capacity = file.read_u64::<LittleEndian>()? as usize;

        let capacity = file.read_u64::<LittleEndian>()? as usize;
        let size = file.read_u64::<LittleEndian>()? as usize;
        let _ = file.read_u64::<LittleEndian>()? as usize;
        let value_bits = file.read_u64::<LittleEndian>()? as usize;

        Ok(Self::new(
            capacity,
            value_bits,
            size,
            partition,
            hash_capacity,
        ))
    }

    pub fn get_idx_mask(&self) -> usize {
        let idx_bits = ((self.hash_capacity as f64).log2().ceil() as usize).max(1);
        (1 << idx_bits) - 1
    }

    pub fn get_idx_bits(&self) -> usize {
        ((self.hash_capacity as f64).log2().ceil() as usize).max(1)
    }

    pub fn get_value_mask(&self) -> usize {
        self.value_mask
    }

    pub fn get_value_bits(&self) -> usize {
        self.value_bits
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

fn read_first_block_from_file<P: AsRef<Path>>(filename: P) -> Result<Page> {
    let mut file = std::fs::File::open(filename)?;

    // Read the index and capacity
    let mut buffer = [0u8; 16];
    file.read_exact(&mut buffer)?;
    let index = LittleEndian::read_u64(&buffer[0..8]) as usize;
    let capacity = LittleEndian::read_u64(&buffer[8..16]) as usize;

    let mut first_zero_end = capacity;
    let chunk_size = 1024; // Define the chunk size for reading
    let mut found_zero = false;
    let mut data = vec![0u32; capacity];
    let mut read_pos = 0;

    while read_pos < capacity {
        let end = usize::min(read_pos + chunk_size, capacity);
        let bytes_to_read = (end - read_pos) * std::mem::size_of::<u32>();
        let mut chunk = vec![0u8; bytes_to_read];
        file.read_exact(&mut chunk)?;
        let chunk_u32 =
            unsafe { std::slice::from_raw_parts(chunk.as_ptr() as *const u32, end - read_pos) };
        data[read_pos..end].copy_from_slice(chunk_u32);

        if let Some(pos) = chunk_u32.iter().position(|&x| x == 0) {
            first_zero_end = read_pos + pos + 1;
            found_zero = true;
            break;
        }
        read_pos = end;
    }

    if !found_zero {
        first_zero_end = capacity;
        eprintln!("Warning: No zero value found in the data, using full capacity.");
    }

    data.truncate(first_zero_end);

    Ok(Page::new(index, first_zero_end, data))
}

#[derive(Clone)]
pub struct Page {
    pub index: usize,
    pub size: usize,
    pub data: Vec<u32>,
}

impl Default for Page {
    fn default() -> Self {
        Self {
            index: 1,
            size: 1,
            data: vec![0],
        }
    }
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

    pub fn merge(&mut self, other: Self) {
        self.size = self.size + other.size;
        self.data.extend_from_slice(&other.data);
    }

    pub fn find_index(
        &self,
        index: usize,
        value: u64,
        value_bits: usize,
        value_mask: usize,
    ) -> u32 {
        let compacted_key = value.left(value_bits) as u32;
        let mut idx = index;
        if idx > self.size {
            return u32::default();
        }

        loop {
            if let Some(cell) = self.data.get(idx) {
                if cell.right(value_mask) == u32::default()
                    || cell.left(value_bits) == compacted_key
                {
                    return cell.right(value_mask);
                }

                idx = idx + 1;
                if idx >= self.size {
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

#[allow(unused)]
pub struct CHTable {
    pub config: HashConfig,
    pub pages: Vec<Page>,
}

impl CHTable {
    pub fn from_pair<P: AsRef<Path> + Debug>(
        config: HashConfig,
        chunk_file1: P,
        chunk_file2: P,
    ) -> Result<CHTable> {
        let mut page = read_page_from_file(chunk_file1)?;
        let next_page = if page.data.last().map_or(false, |&x| x == 0) {
            read_first_block_from_file(chunk_file2)?
        } else {
            Page::default()
        };
        page.merge(next_page);
        let count = page.index;
        let mut pages = vec![Page::default(); count - 1];
        pages.push(page);
        let chtm: CHTable = CHTable { config, pages };
        Ok(chtm)
    }

    pub fn from_hash_files<P: AsRef<Path> + Debug>(
        config: HashConfig,
        hash_files: Vec<P>,
    ) -> Result<CHTable> {
        let mut pages = vec![Page::default(); hash_files.len() + 1];
        for hash_file in hash_files {
            let mut page = read_page_from_file(&hash_file)?;
            let next_page = if page.data.last().map_or(false, |&x| x == 0) {
                read_first_block_from_file(&hash_file)?
            } else {
                Page::default()
            };
            page.merge(next_page);
            if let Some(elem) = pages.get_mut(page.index) {
                *elem = page;
            }
        }

        let chtm = CHTable { config, pages };
        Ok(chtm)
    }

    pub fn from<P: AsRef<Path> + Debug>(config: HashConfig, chunk_file1: P) -> Result<CHTable> {
        let mut page = read_page_from_file(&chunk_file1)?;
        let next_page = if page.data.last().map_or(false, |&x| x == 0) {
            read_first_block_from_file(&chunk_file1)?
        } else {
            Page::default()
        };
        page.merge(next_page);
        let count = page.index;
        let mut pages = vec![Page::default(); count - 1];
        pages.push(page);
        let chtm = CHTable { config, pages };
        Ok(chtm)
    }

    pub fn get_from_page(&self, indx: usize, value: u64, page_index: usize) -> u32 {
        if let Some(page) = self.pages.get(page_index) {
            page.find_index(indx, value, self.config.value_bits, self.config.value_mask)
        } else {
            0
        }
    }
}
