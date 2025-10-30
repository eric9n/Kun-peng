use byteorder::{ByteOrder, LittleEndian, ReadBytesExt, WriteBytesExt};
#[cfg(target_endian = "little")]
use bytemuck::cast_slice_mut;
use std::cmp::Ordering as CmpOrdering;
use std::fmt::{self, Debug};
use std::fs::File;
use std::fs::OpenOptions;
use std::io::{BufWriter, Read, Result, Write};
use std::path::Path;

/// Trait for compact hash operations
pub trait Compact: Default + PartialEq + Clone + Copy + Eq + Sized + Send + Sync + Debug {
    /// Creates a compacted value from a hash key
    ///
    /// # Examples
    ///
    /// ```
    /// use kun_peng::compact_hash::Compact;
    ///
    /// let compacted = u32::compacted(0x1234567890ABCDEF, 16);
    /// assert_eq!(compacted, 0x1234);
    /// ```
    fn compacted(hash_key: u64, value_bits: usize) -> Self;

    /// Creates a hash value from a hash key and a value
    ///
    /// # Examples
    ///
    /// ```
    /// use kun_peng::compact_hash::Compact;
    ///
    /// let hash_value = u32::hash_value(0x1234567890ABCDEF, 16, 0xABCD);
    /// assert_eq!(hash_value, 0x1234ABCD);
    /// ```
    fn hash_value(hash_key: u64, value_bits: usize, value: Self) -> Self;

    /// Returns the left part of the value
    ///
    /// # Examples
    ///
    /// ```
    /// use kun_peng::compact_hash::Compact;
    ///
    /// let value: u32 = 0x1234ABCD;
    /// assert_eq!(value.left(16), 0x1234);
    /// ```
    fn left(&self, value_bits: usize) -> Self;

    /// Returns the right part of the value
    ///
    /// # Examples
    ///
    /// ```
    /// use kun_peng::compact_hash::Compact;
    ///
    /// let value: u32 = 0x1234ABCD;
    /// assert_eq!(value.right(0xFFFF), 0xABCD);
    /// ```
    fn right(&self, value_mask: usize) -> Self;

    /// Combines left and right parts into a single value
    ///
    /// # Examples
    ///
    /// ```
    /// use kun_peng::compact_hash::Compact;
    ///
    /// let combined = u32::combined(0x1234, 0xABCD, 16);
    /// assert_eq!(combined, 0x1234ABCD);
    /// ```
    fn combined(left: Self, right: Self, value_bits: usize) -> Self;

    /// Converts the value to u32
    ///
    /// # Examples
    ///
    /// ```
    /// use kun_peng::compact_hash::Compact;
    ///
    /// let value: u32 = 0x1234ABCD;
    /// assert_eq!(value.to_u32(), 0x1234ABCD);
    /// ```
    fn to_u32(&self) -> u32;

    /// Creates a value from u32
    ///
    /// # Examples
    ///
    /// ```
    /// use kun_peng::compact_hash::Compact;
    ///
    /// let value = u32::from_u32(0x1234ABCD);
    /// assert_eq!(value, 0x1234ABCD);
    /// ```
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
    /// Creates a new Row
    ///
    /// # Examples
    ///
    /// ```
    /// use kun_peng::compact_hash::Row;
    ///
    /// let row = Row::new(0x1234, 1, 2);
    /// assert_eq!(row.value, 0x1234);
    /// assert_eq!(row.seq_id, 1);
    /// assert_eq!(row.kmer_id, 2);
    /// ```
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

// Implement PartialOrd, comparing only the kmer_id field
impl PartialOrd for Row {
    fn partial_cmp(&self, other: &Self) -> Option<CmpOrdering> {
        self.kmer_id.partial_cmp(&other.kmer_id)
    }
}

// Implement Ord, comparing only the kmer_id field
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
    /// Creates a new Slot
    ///
    /// # Examples
    ///
    /// ```
    /// use kun_peng::compact_hash::{Slot, Compact};
    ///
    /// let slot = Slot::<u32>::new(1, 0x1234);
    /// assert_eq!(slot.idx, 1);
    /// assert_eq!(slot.value, 0x1234);
    /// ```
    pub fn new(idx: usize, value: B) -> Self {
        Self { idx, value }
    }

    #[inline]
    pub fn as_slice(&self, slot_size: usize) -> &[u8] {
        let slot_ptr = self as *const Self as *const u8;
        unsafe { std::slice::from_raw_parts(slot_ptr, slot_size) }
    }
}

impl Slot<u64> {
    /// Returns the sequence ID (lower 32 bits) for a u64 Slot
    ///
    /// # Examples
    ///
    /// ```
    /// use kun_peng::compact_hash::{Slot, Compact};
    ///
    /// let slot = Slot::<u64>::new(1, 0x1234567890ABCDEF);
    /// assert_eq!(slot.get_seq_id(), 0x90ABCDEF);
    /// ```
    pub fn get_seq_id(&self) -> u32 {
        self.value.right(0) as u32
    }
}

// Implement PartialOrd, comparing only the idx field
impl<B> PartialOrd for Slot<B>
where
    B: Compact,
{
    fn partial_cmp(&self, other: &Self) -> Option<CmpOrdering> {
        self.idx.partial_cmp(&other.idx)
    }
}

// Implement Ord, comparing only the idx field
impl<B> Ord for Slot<B>
where
    B: Compact,
{
    fn cmp(&self, other: &Self) -> CmpOrdering {
        self.idx.cmp(&other.idx)
    }
}

#[derive(Clone, Copy)]
pub struct HashConfig {
    // value_mask = ((1 << value_bits) - 1);
    pub value_mask: usize,
    // Number of bits for the value
    pub value_bits: usize,
    // Capacity of the hash table
    pub capacity: usize,
    // Current number of elements stored in the hash table
    pub size: usize,
    // Number of partitions
    pub partition: usize,
    // Size of each hash chunk
    pub hash_capacity: usize,
    // Database version (0 is converted from Kraken 2 database)
    pub version: usize,
}

// Manually implement Debug trait for HashConfig
impl fmt::Debug for HashConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("HashConfig")
            .field("version", &self.version)
            .field("partition", &self.partition)
            .field("hash_capacity", &self.hash_capacity)
            .field("capacity", &self.capacity)
            .field("size", &self.size)
            .field("value_bits", &self.value_bits)
            .field("value_mask", &self.value_mask)
            .finish()
    }
}

impl HashConfig {
    /// Creates a new HashConfig
    ///
    /// # Examples
    ///
    /// ```
    /// use kun_peng::compact_hash::HashConfig;
    ///
    /// let config = HashConfig::new(1, 1000, 16, 500, 10, 100);
    /// assert_eq!(config.capacity, 1000);
    /// assert_eq!(config.value_bits, 16);
    /// ```
    pub fn new(
        version: usize,
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
            version,
        }
    }

    pub fn write_to_file<P: AsRef<Path>>(&self, file_path: P) -> Result<()> {
        // Open the file for writing
        let file = File::create(file_path)?;
        let mut writer = BufWriter::new(file);
        writer.write_u64::<LittleEndian>(self.version as u64)?;
        writer.write_u64::<LittleEndian>(self.partition as u64)?;
        writer.write_u64::<LittleEndian>(self.hash_capacity as u64)?;
        writer.write_u64::<LittleEndian>(self.capacity as u64)?;
        writer.write_u64::<LittleEndian>(self.size as u64)?;
        writer.write_u64::<LittleEndian>(self.value_bits as u64)?;
        writer.flush()?;
        Ok(())
    }

    pub fn from_kraken2_header<P: AsRef<Path>>(filename: P) -> Result<Self> {
        let mut file = OpenOptions::new().read(true).open(&filename)?;
        let capacity = file.read_u64::<LittleEndian>()? as usize;
        let size = file.read_u64::<LittleEndian>()? as usize;
        let _ = file.read_u64::<LittleEndian>()? as usize;
        let value_bits = file.read_u64::<LittleEndian>()? as usize;
        Ok(Self::new(0, capacity, value_bits, size, 0, 0))
    }

    pub fn from_hash_header<P: AsRef<Path>>(filename: P) -> Result<Self> {
        let mut file = OpenOptions::new().read(true).open(&filename)?;
        let version = file.read_u64::<LittleEndian>()? as usize;
        let partition = file.read_u64::<LittleEndian>()? as usize;
        let hash_capacity = file.read_u64::<LittleEndian>()? as usize;
        let capacity = file.read_u64::<LittleEndian>()? as usize;
        let size = file.read_u64::<LittleEndian>()? as usize;
        let value_bits = file.read_u64::<LittleEndian>()? as usize;

        Ok(Self::new(
            version,
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

    pub fn compact(&self, hash_key: u64) -> (usize, u32) {
        (self.index(hash_key), hash_key.left(self.value_bits) as u32)
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

fn read_first_block_from_file<P: AsRef<Path>>(filename: P) -> Result<Page> {
    let mut file = std::fs::File::open(filename)?;

    let mut header = [0u8; 16];
    file.read_exact(&mut header)?;
    let index    = LittleEndian::read_u64(&header[0..8])  as usize;
    let capacity = LittleEndian::read_u64(&header[8..16]) as usize;

    let chunk_elems = 4096;
    let mut data = vec![0u32; capacity];

    let mut read_pos = 0usize;
    let mut found_zero = false;
    let mut first_zero_end = capacity;

    #[cfg(target_endian = "little")]
    {
        while read_pos < capacity {
            let to_read_elems = (capacity - read_pos).min(chunk_elems);
            let bytes = cast_slice_mut::<u32, u8>(&mut data[read_pos..read_pos+to_read_elems]);
            file.read_exact(bytes)?;

            if let Some(pos) = data[read_pos..read_pos+to_read_elems].iter().position(|&x| x == 0) {
                first_zero_end = read_pos + pos + 1;
                found_zero = true;
                break;
            }
            read_pos += to_read_elems;
        }
    }

    #[cfg(not(target_endian = "little"))]
    {
        // 回退到方案 A 的安全端序路径
        let mut tmp = vec![0u8; chunk_elems * 4];
        while read_pos < capacity {
            let to_read_elems = (capacity - read_pos).min(chunk_elems);
            let to_read_bytes = to_read_elems * 4;
            file.read_exact(&mut tmp[..to_read_bytes])?;
            LittleEndian::read_u32_into(&tmp[..to_read_bytes], &mut data[read_pos..read_pos+to_read_elems]);
            if let Some(pos) = data[read_pos..read_pos+to_read_elems].iter().position(|&x| x == 0) {
                first_zero_end = read_pos + pos + 1;
                found_zero = true;
                break;
            }
            read_pos += to_read_elems;
        }
    }

    if !found_zero {
        first_zero_end = capacity;
        eprintln!("Warning: No zero value found in the data, using full capacity.");
    }

    data.truncate(first_zero_end);

    Ok(Page::new(index, first_zero_end, data))
}

fn read_page_metadata(file: &mut File) -> Result<(usize, usize)> {
    let index = file.read_u64::<LittleEndian>()? as usize;
    let capacity = file.read_u64::<LittleEndian>()? as usize;
    Ok((index, capacity))
}

// fn read_page_data(file: &mut File, data: &mut [u32]) -> Result<()> {
//     let data_bytes = unsafe {
//         std::slice::from_raw_parts_mut(
//             data.as_mut_ptr() as *mut u8,
//             data.len() * std::mem::size_of::<u32>(),
//         )
//     };
//     file.read_exact(data_bytes)?;
//     Ok(())
// }

fn read_page_data(file: &mut File, data: &mut [u32]) -> Result<()> {
    #[cfg(target_endian = "little")]
    {
        // 零拷贝，安全，无端序转换（文件即小端）
        let bytes = cast_slice_mut::<u32, u8>(data);
        file.read_exact(bytes)?;
        return Ok(());
    }

    #[cfg(not(target_endian = "little"))]
    {
        // 非小端平台回退到端序转换路径（方案 A）
        const CHUNK: usize = 4096;
        let mut buf = vec![0u8; CHUNK * 4];
        let mut filled = 0usize;
        while filled < data.len() {
            let n = (data.len() - filled).min(CHUNK);
            let bytes = &mut buf[..n * 4];
            file.read_exact(bytes)?;
            LittleEndian::read_u32_into(bytes, &mut data[filled..filled + n]);
            filled += n;
        }
        return Ok(());
    }
}

fn read_page_from_file<P: AsRef<Path>>(filename: P) -> Result<Page> {
    let mut file = std::fs::File::open(filename)?;
    let (index, capacity) = read_page_metadata(&mut file)?;
    let mut data = vec![0u32; capacity];
    read_page_data(&mut file, &mut data)?;

    Ok(Page::new(index, capacity, data))
}

fn read_large_page_from_file<P: AsRef<Path>>(large_page: &mut Page, filename: P) -> Result<()> {
    let mut file = File::open(filename)?;

    let (index, capacity) = read_page_metadata(&mut file)?;

    let current_len = large_page.data.capacity();

    if capacity > current_len {
        // If the capacity in the file is greater than the current page's capacity, extend the memory
        large_page.data.resize(capacity, 0);
    } else if capacity < current_len {
        // If the capacity in the file is smaller than the current page's capacity, truncate and zero out the excess
        large_page.data.truncate(capacity);
        large_page.data.shrink_to_fit(); // Free up excess memory
    }

    read_page_data(&mut file, &mut large_page.data)?;

    large_page.index = index;
    large_page.size = capacity;

    Ok(())
}

pub fn read_next_page<P: AsRef<Path> + Debug>(
    large_page: &mut Page,
    hash_sorted_files: &Vec<P>,
    page_index: usize,
    config: HashConfig,
) -> Result<()> {
    let mut hash_file = &hash_sorted_files[page_index];
    let parition = config.partition;
    read_large_page_from_file(large_page, hash_file)?;

    let next_page = if large_page.data.last().map_or(false, |&x| x != 0) {
        if config.version < 1 {
            hash_file = &hash_sorted_files[(page_index + 1) % parition]
        }
        read_first_block_from_file(&hash_file)?
    } else {
        Page::default()
    };
    large_page.merge(next_page);

    Ok(())
}

#[derive(Clone)]
pub struct Page {
    pub index: usize,
    pub size: usize,
    pub data: Vec<u32>,
}

impl Default for Page {
    fn default() -> Self {
        Self::with_capacity(1, 1)
    }
}

impl Page {
    pub fn with_capacity(index: usize, capacity: usize) -> Self {
        Self::new(index, capacity, vec![0; capacity + 1024])
    }

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
        let new_size = self.size + other.size;
        if self.data.capacity() < new_size {
            self.data.reserve(new_size - self.data.len());
        }
        self.data.extend_from_slice(&other.data[..other.size]);
        self.size = new_size;
    }

    pub fn find_index(
        &self,
        index: usize,
        compacted_key: u32,
        value_bits: usize,
        value_mask: usize,
    ) -> u32 {
        let mut idx = index;
        if idx >= self.size {
            return 0;
        }

        loop {
            if let Some(cell) = self.data.get(idx) {
                if cell.right(value_mask) == 0 || cell.left(value_bits) == compacted_key {
                    return cell.right(value_mask);
                }

                idx = idx + 1;
                if idx >= self.size {
                    break;
                }
            } else {
                return 0;
            }
        }
        0
    }
}

#[allow(unused)]
pub struct CHTable {
    pub config: HashConfig,
    pub pages: Vec<Page>,
}

impl CHTable {
    pub fn from_hash_files<P: AsRef<Path> + Debug>(
        config: HashConfig,
        hash_sorted_files: &Vec<P>,
    ) -> Result<CHTable> {
        let end = hash_sorted_files.len();
        Self::from_range(config, hash_sorted_files, 0, end)
    }

    pub fn from_range<P: AsRef<Path> + Debug>(
        config: HashConfig,
        hash_sorted_files: &Vec<P>,
        start: usize,
        end: usize,
    ) -> Result<CHTable> {
        let mut pages = vec![Page::default(); start];
        let parition = hash_sorted_files.len();
        for i in start..end {
            let mut hash_file = &hash_sorted_files[i];
            let mut page = read_page_from_file(&hash_file)?;
            let next_page = if page.data.last().map_or(false, |&x| x != 0) {
                if config.version < 1 {
                    hash_file = &hash_sorted_files[(i + 1) % parition]
                }
                read_first_block_from_file(&hash_file)?
            } else {
                Page::default()
            };
            page.merge(next_page);
            pages.push(page);
        }

        let chtm = CHTable { config, pages };
        Ok(chtm)
    }

    pub fn get_from_page(&self, indx: usize, compacted: u32, page_index: usize) -> u32 {
        if let Some(page) = self.pages.get(page_index) {
            page.find_index(
                indx,
                compacted,
                self.config.value_bits,
                self.config.value_mask,
            )
        } else {
            0
        }
    }
}
