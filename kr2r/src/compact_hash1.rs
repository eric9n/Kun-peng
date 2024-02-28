use byteorder::{ByteOrder, LittleEndian};
use memmap2::{Mmap, MmapMut, MmapOptions};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fmt;
use std::fs::{File, OpenOptions};
use std::io::{Error, ErrorKind, Result};
use std::path::Path;
use std::sync::atomic::{AtomicU32, Ordering};

#[allow(unused_variables)]
#[inline]
fn second_hash(first_hash: u64) -> usize {
    #[cfg(feature = "double_hashing")]
    {
        (first_hash >> 8) | 1 as usize // 使用双重哈希
    }

    #[cfg(not(feature = "double_hashing"))]
    {
        1 // 使用线性探测
    }
}

// value in the low bits
// hash of the key in the high bits
// CompactHashCell 结构体，用于存储哈希表的单个单元。
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct CompactHashCell(pub u32);

impl CompactHashCell {
    #[inline]
    pub fn hashed_key(&self, value_bits: usize) -> u64 {
        self.0 as u64 >> value_bits
    }

    /// value_mask = ((1 << value_bits) - 1);
    #[inline]
    pub fn value(&self, value_mask: u32) -> u32 {
        self.0 & value_mask
    }

    #[inline]
    pub fn populate(&mut self, compacted_key: u64, val: u32, value_bits: usize) {
        // 直接使用右移操作来判断
        if val >> value_bits != 0 {
            panic!(
                "Value length of {} is too small for value of {}",
                value_bits, val
            );
        }
        self.0 = (compacted_key << value_bits) as u32;
        self.0 |= val;
    }
}

impl Default for CompactHashCell {
    fn default() -> Self {
        CompactHashCell(0)
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
pub struct CompactHashTable<'a> {
    // meta_table: MetaTable<'a>,
    // memmap
    #[allow(unused)]
    mmap: Mmap,
    // 哈希表的容量。
    capacity: usize,
    // 哈希表中当前存储的元素数量。
    #[allow(unused)]
    size: usize,
    // 键的位数。
    // key_bits: usize,
    // 值的位数。
    value_bits: usize,
    // value_mask = ((1 << value_bits) - 1);
    value_mask: u32,
    // 存储哈希单元的向量。
    pub table: &'a [CompactHashCell],
}

impl<'a> fmt::Debug for CompactHashTable<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CompactHashTable")
            .field("capacity", &self.capacity)
            .field("size", &self.size)
            .field("value_bits", &self.value_bits)
            .field("value_mask", &self.value_mask)
            .finish()
    }
}

impl<'a> CompactHashTable<'a> {
    pub fn from<P: AsRef<Path>>(filename: P) -> Result<CompactHashTable<'a>> {
        let file = File::open(&filename)?;
        let mmap = unsafe { Mmap::map(&file)? };

        let capacity = LittleEndian::read_u64(&mmap[0..8]) as usize;
        let size = LittleEndian::read_u64(&mmap[8..16]) as usize;
        // let _ = LittleEndian::read_u64(&mmap[16..24]) as usize;
        let value_bits = LittleEndian::read_u64(&mmap[24..32]) as usize;

        let table = unsafe {
            std::slice::from_raw_parts(mmap.as_ptr().add(32) as *const CompactHashCell, capacity)
        };
        let chtm = CompactHashTable {
            value_mask: (1 << value_bits) - 1,
            capacity,
            size,
            // key_bits,
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
                let val = cell.value(self.value_mask) as u64; // 使用正确的 value_mask
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

    // 私有通用函数，用于查找元素
    fn find_cell(&self, hash_key: u64) -> Option<(usize, &CompactHashCell)> {
        let compacted_key = hash_key >> (32 + self.value_bits);
        let mut idx = hash_key as usize % self.capacity;
        let first_idx = idx;
        let mut step = 0;
        while let Some(cell) = self.table.get(idx as usize) {
            if cell.value(self.value_mask) == 0 {
                break;
            }
            if cell.hashed_key(self.value_bits) == compacted_key {
                return Some((idx as usize, cell));
            }
            if step == 0 {
                step = second_hash(hash_key);
            }
            idx = (idx + step) % self.capacity;
            if idx == first_idx {
                break;
            }
        }
        None
    }

    // 使用 find_cell 函数的公开方法
    pub fn get(&self, hash_key: u64) -> u32 {
        self.find_cell(hash_key)
            .map(|(_, cell)| cell.value(self.value_mask))
            .unwrap_or(0)
    }

    pub fn find_index(&self, hash_key: u64) -> Option<usize> {
        self.find_cell(hash_key).map(|(idx, _)| idx)
    }
}

#[derive(Debug)]
#[repr(C)]
pub struct CHCell(pub AtomicU32);

impl CHCell {
    pub fn from_key(hash_key: u64, value: u32, value_bits: usize) -> Self {
        let compacted_key = hash_key >> (64 - value_bits as u64);
        let new_value = ((compacted_key << value_bits) as u32) | value;
        CHCell(AtomicU32::new(new_value))
    }

    // 提供一个方法来安全地修改值
    #[inline]
    pub fn update(&self, value: u32) {
        self.0.store(value, Ordering::SeqCst);
    }

    // 如果需要，也可以提供获取当前值的方法
    #[inline]
    pub fn get(&self) -> u32 {
        self.0.load(Ordering::SeqCst)
    }

    #[inline]
    pub fn hashed_key(&self, value_bits: usize) -> u64 {
        let value = self.0.load(Ordering::SeqCst) as u64;
        value >> value_bits
    }

    /// value_mask = ((1 << value_bits) - 1);
    #[inline]
    pub fn value(&self, value_mask: u32) -> u32 {
        let value = self.0.load(Ordering::SeqCst);
        value & value_mask
    }

    #[inline]
    pub fn populate(&self, cell: Cell, value_bits: usize) {
        let new_value = ((cell.compacted_key << value_bits) as u32) | cell.taxid;
        self.0.store(new_value, Ordering::SeqCst);
    }
}

impl Default for CHCell {
    fn default() -> Self {
        CHCell(AtomicU32::new(0))
    }
}

use std::cmp::Ordering as CmpOrdering;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(C)]
pub struct Cell {
    pub index: usize,
    pub compacted_key: u64,
    pub taxid: u32,
}

impl Cell {
    pub fn new(index: usize, compacted_key: u64, taxid: u32) -> Self {
        Self {
            index,
            compacted_key,
            taxid,
        }
    }
}

// 实现 PartialOrd，只比较 index 字段
impl PartialOrd for Cell {
    fn partial_cmp(&self, other: &Self) -> Option<CmpOrdering> {
        self.index.partial_cmp(&other.index)
    }
}

// 实现 Ord，只比较 index 字段
impl Ord for Cell {
    fn cmp(&self, other: &Self) -> CmpOrdering {
        self.index.cmp(&other.index)
    }
}

#[derive(Debug)]
pub struct CompactHashTableMut<'a> {
    // meta_table: MetaTable<'a>,
    // memmap
    pub mmap: MmapMut,
    // 哈希表的容量。
    capacity: usize,
    // 哈希表中当前存储的元素数量。
    pub size: usize,
    // 键的位数。
    // key_bits: usize,
    // 值的位数。
    value_bits: usize,
    // value_mask = ((1 << value_bits) - 1);
    value_mask: u32,
    // 存储哈希单元的向量。
    pub table: &'a mut [CHCell],
}

impl<'a> CompactHashTableMut<'a> {
    pub fn new<P: AsRef<Path>>(
        hash_file: P,
        capacity: usize,
        value_bits: usize,
    ) -> Result<CompactHashTableMut<'a>> {
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
            std::slice::from_raw_parts_mut(mut_mmap.as_mut_ptr().add(32) as *mut CHCell, capacity)
        };

        // let size_ptr: *mut usize = unsafe { mut_mmap.as_mut_ptr().add(8) as *mut usize };
        // let size_ref: &mut usize = unsafe { &mut *size_ptr };

        let chtm = Self {
            value_mask: (1 << value_bits) - 1,
            capacity,
            size,
            // key_bits,
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

    /// 找到hash_key所在的槽位
    pub fn get_table_index(&self, hash_key: u64) -> Option<Cell> {
        let compacted_key = hash_key >> (32 + self.value_bits);
        let mut idx = (hash_key as usize) % self.capacity;
        let first_idx = idx;
        let step = 1;

        loop {
            let cell = &self.table[idx];
            if cell.value(self.value_mask) == 0 {
                return Some(Cell::new(idx, compacted_key, 0));
            }
            if cell.hashed_key(self.value_bits) == compacted_key {
                return Some(Cell::new(idx, compacted_key, cell.value(self.value_mask)));
            }

            idx = (idx + step) % self.capacity;
            if idx == first_idx {
                break; // 遍历完所有槽位，未找到匹配项
            }
        }

        None
    }

    /// 设置单元格内容
    pub fn set_cell(&self, index: usize, value: u32) {
        let cell = &self.table[index];
        cell.update(value);
    }

    pub fn compare_and_set(&self, cell: Cell, old_value: &mut u32) -> bool {
        if cell.taxid == 0 {
            return false;
        }

        let mut idx = cell.index;
        let first_idx = idx;
        let step = 1;
        loop {
            let cell1 = &self.table[idx as usize];
            if cell1.value(self.value_mask) == 0
                || cell1.hashed_key(self.value_bits) == cell.compacted_key
            {
                if *old_value == cell1.value(self.value_mask) {
                    cell1.populate(cell, self.value_bits);
                    if *old_value == 0 {
                        // self.size += 1;
                    }
                    return true;
                } else {
                    *old_value = cell1.value(self.value_mask);
                    return false;
                }
            }
            idx = (idx + step) % self.capacity;
            if idx == first_idx {
                break; // 退出循环
            }
        }

        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // use std::fs;

    // #[test]
    // fn test_new() {
    //     let cht = CHTMeta::new(1024, 16, 16);
    //     assert_eq!(cht.capacity, 1024);
    //     assert_eq!(cht.key_bits, 16);
    //     assert_eq!(cht.value_bits, 16);
    // }

    // #[test]
    // fn test_write_and_load_table() {
    //     // 创建临时文件
    //     let tmp_file = "temp_cht.bin";

    //     // 创建一个示例哈希表并写入文件
    //     let cht = CHTMeta::new(2, 16, 16);
    //     let cells = vec![CompactHashCell(12345), CompactHashCell(67890)];
    //     cht.write_table(tmp_file, &cells).unwrap();

    //     // 从文件中读取哈希表
    //     let (loaded_cht, loaded_cells) = CHTMeta::load_table(tmp_file).unwrap();

    //     assert_eq!(loaded_cht.capacity, 2);
    //     assert_eq!(loaded_cells.len(), 2);
    //     assert_eq!(loaded_cells[0].0, 12345);
    //     assert_eq!(loaded_cells[1].0, 67890);

    //     let file = File::open(tmp_file).unwrap();
    //     let mmap = unsafe { Mmap::map(&file).unwrap() };

    //     let (loaded_cht, loaded_cells) = CHTMeta::load_table_mmap(&mmap).unwrap();

    //     assert_eq!(loaded_cht.capacity, 2);
    //     assert_eq!(loaded_cells.len(), 2);
    //     assert_eq!(loaded_cells[0].0, 12345);
    //     assert_eq!(loaded_cells[1].0, 67890);

    //     // 清理临时文件
    //     fs::remove_file(tmp_file).unwrap();
    // }

    #[test]
    fn test_compact_hash_cell() {
        let mut cell = CompactHashCell(0u32);
        cell.populate(123u64, 456u32, 16);

        let value_mask = (1 << 16) - 1;
        assert_eq!(cell.hashed_key(16), 123);
        assert_eq!(cell.value(value_mask), 456);
    }
}
