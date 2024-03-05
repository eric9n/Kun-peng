use crate::compact_hash::{Cell, CellIndex, CompactHash};
use byteorder::{ByteOrder, LittleEndian};
use memmap2::{MmapMut, MmapOptions};
use std::fmt;
use std::fs::OpenOptions;
use std::io::Result;
use std::path::Path;

// pub const U32MAX: usize = u32::MAX as usize;

pub struct CompactHashTable<'a> {
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
    pub table: &'a mut [u32],

    pub page: Vec<u32>,
    pub page_index: usize,
    pub page_size: usize,
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
    pub fn new<P: AsRef<Path>>(
        hash_file: P,
        capacity: usize,
        value_bits: usize,
        page_index: usize,
        page_size: usize,
    ) -> Result<CompactHashTable<'a>> {
        let key_bits = 32 - value_bits;

        let file_len = 32 + 4 * capacity;
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .open(&hash_file)?;
        file.set_len(file_len as u64)?;

        let mut mut_mmap = unsafe { MmapOptions::new().len(file_len).map_mut(&file)? };

        let size: usize = 0;

        mut_mmap[0..8].copy_from_slice(&capacity.to_le_bytes());
        mut_mmap[8..16].copy_from_slice(&size.to_le_bytes());
        mut_mmap[16..24].copy_from_slice(&key_bits.to_le_bytes());
        mut_mmap[24..32].copy_from_slice(&value_bits.to_le_bytes());

        let table = unsafe {
            std::slice::from_raw_parts_mut(mut_mmap.as_mut_ptr().add(32) as *mut u32, capacity)
        };

        let start_index = page_index * page_size;
        let end_index = std::cmp::min((page_index + 1) * page_size, capacity);

        let page: Vec<u32> = if start_index < capacity {
            table[start_index..end_index].to_vec()
        } else {
            Vec::new() // start_index 超出范围，返回空 Vec
        };

        let chtm = Self {
            value_mask: (1 << value_bits) - 1,
            capacity,
            size,
            value_bits,
            table,
            mmap: mut_mmap,
            page,
            page_index,
            page_size,
        };
        Ok(chtm)
    }
}

impl<'a> CompactHashTable<'a> {
    pub fn set_page_cell(&mut self, ci: CellIndex) -> Option<CellIndex> {
        let mut idx = ci.index;
        let first_idx = idx;

        // 因为每页的大小是u32::MAX,如果初始idx大于这个数,就返回None,虽然实际运行中的数据以避免了这种情况.

        while let Some(cell) = self.page.get_mut(idx) {
            if cell.taxid(self.value_mask) == 0 {
                let val = (ci.cell.compacted_key << self.value_bits) | ci.cell.taxid;
                *cell = val;
                break;
            }
            if cell.compacted_key(self.value_bits) == ci.cell.compacted_key {
                return Some(CellIndex::new(
                    idx,
                    ci.cell.compacted_key,
                    cell.taxid(self.value_mask),
                ));
            }

            idx = idx + 1;

            if idx == self.page_size {
                return self.set_table_cell(idx, ci.cell);
            }
            if idx == first_idx {
                break;
            }
        }
        None
    }

    fn set_table_cell(&self, index: usize, value: Cell) -> Option<CellIndex> {
        let mut idx = index;
        let first_idx = idx;
        while let Some(cell1) = self.table.get(idx) {
            if cell1.taxid(self.value_mask) == 0 {
                cell1.populate(value, self.value_bits);
                break;
            }
            if cell1.compacted_key(self.value_bits) == value.compacted_key {
                return Some(CellIndex::new(
                    idx,
                    value.compacted_key,
                    cell1.taxid(self.value_mask),
                ));
            }

            idx = (idx + 1) % self.capacity;
            if idx == first_idx {
                break;
            }
        }
        None
    }

    pub fn copy_from_page(&mut self) {
        let start_index = self.page_index * self.page_size;
        let end_index = std::cmp::min((self.page_index + 1) * self.page_size, self.capacity);
        self.table[start_index..end_index].copy_from_slice(&self.page);
    }

    // 直接更新
    pub fn update_cell(&mut self, ci: CellIndex) {
        if ci.index < self.page_size {
            self.page[ci.index] = (ci.cell.compacted_key << self.value_bits) | ci.cell.taxid;
        } else {
            self.table[ci.index] = (ci.cell.compacted_key << self.value_bits) | ci.cell.taxid;
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
