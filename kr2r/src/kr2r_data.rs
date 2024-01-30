use crate::{Meros, CURRENT_REVCOM_VERSION};
use std::fs::File;
use std::io::{Read, Result};
use std::path::Path;

/// 判断u64的值是否为0，并将其转换为Option<u64>类型
pub fn u64_to_option(value: u64) -> Option<u64> {
    Option::from(value).filter(|&x| x != 0)
}

#[repr(C)]
pub struct IndexOptions {
    pub k: usize,
    pub l: usize,
    pub spaced_seed_mask: u64,
    pub toggle_mask: u64,
    pub minimum_acceptable_hash_value: u64,
    pub dna_db: bool,
    pub revcom_version: i32, // 如果等于 0，就报错
    pub db_version: i32,     // 为未来的数据库结构变化预留
    pub db_type: i32,        // 为未来使用其他数据结构预留
}

impl IndexOptions {
    pub fn read_index_options<P: AsRef<Path>>(file_path: P) -> Result<Self> {
        let mut file = File::open(file_path)?;
        let mut buffer = vec![0; std::mem::size_of::<Self>()];
        file.read_exact(&mut buffer)?;

        let idx_opts = unsafe {
            // 确保这种转换是安全的，这依赖于数据的确切布局和来源
            std::ptr::read(buffer.as_ptr() as *const Self)
        };

        if idx_opts.revcom_version != CURRENT_REVCOM_VERSION as i32 {
            // 如果版本为0，直接触发 panic
            panic!("Unsupported version (revcom_version == 0)");
        }

        Ok(idx_opts)
    }

    pub fn as_meros(&self) -> Meros {
        Meros::new(
            self.k,
            self.l,
            u64_to_option(self.spaced_seed_mask),
            u64_to_option(self.toggle_mask),
            u64_to_option(self.minimum_acceptable_hash_value),
        )
    }
}
