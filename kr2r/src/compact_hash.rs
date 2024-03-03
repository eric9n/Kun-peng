use byteorder::{ByteOrder, LittleEndian};
use memmap2::{MmapMut, MmapOptions};
use rayon::prelude::*;
use std::cmp::Ordering as CmpOrdering;
use std::collections::HashMap;
use std::fmt;
use std::fs::OpenOptions;
use std::io::{Error, ErrorKind, Result};
use std::path::Path;
use std::sync::atomic::{AtomicU32, Ordering};

pub trait Compact {
    fn compacted(&self, value_bits: usize) -> u32;
    fn index(&self, capacity: usize) -> usize;
}

impl Compact for u64 {
    fn compacted(&self, value_bits: usize) -> u32 {
        (*self >> (32 + value_bits)) as u32
    }

    fn index(&self, capacity: usize) -> usize {
        *self as usize % capacity
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Cell {
    /// let hash_key: u64 = 123;
    /// let compacted_key: u32 = hash_key >> (32 + value_bits);
    pub compacted_key: u32,
    pub taxid: u32,
}

impl Cell {
    pub fn new(compacted_key: u32, taxid: u32) -> Self {
        Self {
            compacted_key,
            taxid,
        }
    }

    pub fn from_u32(cell_value: u32, value_bits: usize, value_mask: u32) -> Self {
        Self {
            compacted_key: cell_value.compacted_key(value_bits),
            taxid: cell_value.taxid(value_mask),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CellIndex {
    pub cell: Cell,
    pub index: usize,
}

impl CellIndex {
    pub fn new(index: usize, compacted_key: u32, taxid: u32) -> Self {
        Self {
            index,
            cell: Cell::new(compacted_key, taxid),
        }
    }
}

// 实现 PartialOrd，只比较 index 字段
impl PartialOrd for CellIndex {
    fn partial_cmp(&self, other: &Self) -> Option<CmpOrdering> {
        self.index.partial_cmp(&other.index)
    }
}

// 实现 Ord，只比较 index 字段
impl Ord for CellIndex {
    fn cmp(&self, other: &Self) -> CmpOrdering {
        self.index.cmp(&other.index)
    }
}

pub trait CompactHash {
    // 定义 trait 方法
    /// compacted_key
    fn compacted_key(&self, value_bits: usize) -> u32;
    fn taxid(&self, value_mask: u32) -> u32;
    fn populate(&self, cell: Cell, value_bits: usize);
}

// 实现 CompactHash trait 对于 u32 类型
impl CompactHash for u32 {
    #[inline]
    fn compacted_key(&self, value_bits: usize) -> u32 {
        // 为 u32 类型实现 hashed_key 方法
        *self >> value_bits
    }

    #[inline]
    fn taxid(&self, value_mask: u32) -> u32 {
        // 为 u32 类型实现 value 方法
        *self & value_mask
    }

    #[inline]
    #[allow(unused_variables)]
    fn populate(&self, cell: Cell, value_bits: usize) {
        // *self = (cell.compacted_key << value_bits) | cell.taxid;
    }
}

impl CompactHash for AtomicU32 {
    #[inline]
    fn compacted_key(&self, value_bits: usize) -> u32 {
        let value = self.load(Ordering::SeqCst);
        value >> value_bits
    }

    #[inline]
    fn taxid(&self, value_mask: u32) -> u32 {
        let value = self.load(Ordering::SeqCst);
        value & value_mask
    }

    #[inline]
    fn populate(&self, cell: Cell, value_bits: usize) {
        let new_value = (cell.compacted_key << value_bits) | cell.taxid;
        self.store(new_value, Ordering::SeqCst);
    }
}

// CompactHashTable 结构体
/// 实现一个紧凑的、固定大小的概率哈希表。
///
/// 支持的操作：
///
/// - 搜索（search）
/// - 比较并设置（compare-and-set）
///
/// 不支持的操作：
///
/// - 值为 0
///
/// 键（Keys）必须是 64 位无符号整数，小于 2 的 63 次方。
/// 值（Values）必须是 32 位无符号整数，大于 0。
///
/// 单元格（Cells）是 32 位无符号整数，截断的哈希键存储在最高有效位，
/// 截断的值存储在最低有效位。压缩级别在表创建时通过 HashTable 类的
/// WriteCompactedTable() 方法设置。
// #[derive(Debug)]
pub struct CompactHashTable<'a, T>
where
    T: CompactHash,
{
    // memmap
    #[allow(unused)]
    pub mmap: MmapMut,
    // 哈希表的容量。
    pub capacity: usize,
    // 哈希表中当前存储的元素数量。
    #[allow(unused)]
    size: usize,
    // 值的位数。
    pub value_bits: usize,
    // value_mask = ((1 << value_bits) - 1);
    pub value_mask: u32,
    // 存储哈希单元的向量。
    pub table: &'a [T],
}

impl<'a, T> fmt::Debug for CompactHashTable<'a, T>
where
    T: CompactHash,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CompactHashTable")
            .field("capacity", &self.capacity)
            .field("size", &self.size)
            .field("value_bits", &self.value_bits)
            .field("value_mask", &self.value_mask)
            .finish()
    }
}

impl<'a, T> CompactHashTable<'a, T>
where
    T: CompactHash,
{
    /// 找到hash_key所在的槽位
    #[inline]
    pub fn find_cell(&self, hash_key: u64) -> Option<CellIndex> {
        let compacted_key = hash_key.compacted(self.value_bits);
        let mut idx = hash_key.index(self.capacity);
        let first_idx = idx;
        let step = 1;
        while let Some(cell) = self.table.get(idx) {
            if cell.taxid(self.value_mask) == 0
                || cell.compacted_key(self.value_bits) == compacted_key
            {
                return Some(CellIndex::new(
                    idx,
                    compacted_key,
                    cell.taxid(self.value_mask),
                ));
            }
            idx = (idx + step) % self.capacity;
            if idx == first_idx {
                break;
            }
        }
        None
    }

    #[inline]
    pub fn get(&self, hash_key: u64) -> u32 {
        let compacted_key = hash_key.compacted(self.value_bits);
        let mut idx = hash_key.index(self.capacity);
        let first_idx = idx;
        let step = 1;
        while let Some(cell) = self.table.get(idx) {
            if cell.taxid(self.value_mask) == 0 {
                return 0;
            }
            if cell.compacted_key(self.value_bits) == compacted_key {
                return cell.taxid(self.value_mask);
            }
            idx = (idx + step) % self.capacity;
            if idx == first_idx {
                break;
            }
        }
        0
    }
}

impl<'a> CompactHashTable<'a, u32> {
    pub fn from<P: AsRef<Path>>(filename: P) -> Result<CompactHashTable<'a, u32>> {
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .open(&filename)?;

        let mmap = unsafe { MmapOptions::new().populate().map_mut(&file)? };

        let capacity = LittleEndian::read_u64(&mmap[0..8]) as usize;
        let size = LittleEndian::read_u64(&mmap[8..16]) as usize;
        // let _ = LittleEndian::read_u64(&mmap[16..24]) as usize;
        let value_bits = LittleEndian::read_u64(&mmap[24..32]) as usize;

        let table =
            unsafe { std::slice::from_raw_parts(mmap.as_ptr().add(32) as *const u32, capacity) };
        let chtm = CompactHashTable {
            value_mask: (1 << value_bits) - 1,
            capacity,
            size,
            value_bits,
            table,
            mmap,
        };
        Ok(chtm)
    }

    pub fn get_value_counts(&self) -> HashMap<u64, u64> {
        self.table
            .par_iter() // 使用并行迭代器
            .fold(HashMap::new, |mut acc, cell| {
                let val = cell.taxid(self.value_mask) as u64; // 使用正确的 value_mask
                if val != 0 {
                    *acc.entry(val).or_insert(0) += 1;
                }
                acc
            })
            .reduce(HashMap::new, |mut acc, m| {
                for (key, val) in m {
                    *acc.entry(key).or_insert(0) += val;
                }
                acc
            })
    }

    #[inline]
    pub fn find_index(&self, hash_key: u64) -> Option<usize> {
        if let Some(ci) = self.find_cell(hash_key) {
            if ci.cell.taxid != 0 {
                return Some(ci.index);
            }
        }
        None
    }
}

impl<'a> CompactHashTable<'a, AtomicU32> {
    pub fn new<P: AsRef<Path>>(
        hash_file: P,
        capacity: usize,
        value_bits: usize,
    ) -> Result<CompactHashTable<'a, AtomicU32>> {
        let key_bits = 32 - value_bits;
        if key_bits == 0 || value_bits == 0 {
            return Err(Error::new(
                ErrorKind::InvalidInput,
                "Key bits and value bits cannot be zero",
            ));
        }

        let file_len = 32 + 4 * capacity;
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .open(&hash_file)?;
        file.set_len(file_len as u64)?;

        let mut mut_mmap = unsafe { MmapOptions::new().len(file_len).populate().map_mut(&file)? };

        let size: usize = 0;

        mut_mmap[0..8].copy_from_slice(&capacity.to_le_bytes());
        mut_mmap[8..16].copy_from_slice(&size.to_le_bytes());
        mut_mmap[16..24].copy_from_slice(&key_bits.to_le_bytes());
        mut_mmap[24..32].copy_from_slice(&value_bits.to_le_bytes());

        let table = unsafe {
            std::slice::from_raw_parts_mut(
                mut_mmap.as_mut_ptr().add(32) as *mut AtomicU32,
                capacity,
            )
        };

        let chtm = Self {
            value_mask: (1 << value_bits) - 1,
            capacity,
            size,
            value_bits,
            table,
            mmap: mut_mmap,
        };
        Ok(chtm)
    }

    pub fn update_size(&self, size: usize) {
        let size_bytes = size.to_le_bytes(); // 假设使用小端序
        unsafe {
            let size_ptr = self.mmap.as_ptr().add(8) as *mut u8;
            std::ptr::copy_nonoverlapping(size_bytes.as_ptr(), size_ptr, size_bytes.len());
        }
        self.mmap.flush().expect("Failed to flush mmap");
    }

    // 直接更新
    pub fn update_cell(&self, ci: CellIndex) {
        let cell = &self.table[ci.index];
        cell.populate(ci.cell, self.value_bits);
    }

    /// 设置单元格内容
    /// 如果已经存在 compacted_key, 但是taxid不一致,就返回文件中的内容,重新lca之后再更新
    pub fn set_cell(&self, ci: CellIndex) -> Option<CellIndex> {
        let mut idx = ci.index;
        let first_idx = idx;
        let step = 1;
        while let Some(cell1) = self.table.get(idx) {
            if cell1.taxid(self.value_mask) == 0 {
                cell1.populate(ci.cell, self.value_bits);
                break;
            }
            if cell1.compacted_key(self.value_bits) == ci.cell.compacted_key {
                return Some(CellIndex::new(
                    idx,
                    ci.cell.compacted_key,
                    cell1.taxid(self.value_mask),
                ));
            }

            idx = (idx + step) % self.capacity;
            if idx == first_idx {
                break;
            }
        }
        None
    }
}
