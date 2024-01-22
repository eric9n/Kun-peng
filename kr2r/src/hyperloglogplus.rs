use crate::ertl;
use std::collections::HashSet;
use std::error::Error;

trait HyperLogLogPlusMinus {
    fn new(precision: u8, bit_mixer: fn(u64) -> u64) -> Self
    where
        Self: Sized;

    fn insert(&mut self, item: u64) -> Result<(), Box<dyn Error>>;

    fn insert_all(&mut self, items: &[u64]) {
        for &item in items {
            self.insert(item);
        }
    }

    fn cardinality(&mut self) -> u64;
}

struct HLLSparse {
    n_observed: u64,
    p: u8,
    sparse_list: HashSet<u32>,
    bit_mixer: fn(u64) -> u64,
    // number of registers
    m: usize,
    // sparse versions of p and m
    // precision when using a sparse representation
    // which encodes the rank + index in 32 bits:
    //  25 bits for index +
    //   6 bits for rank +
    //   1 flag bit indicating if bits p..pPrime are 0
    p_prime: u8,
    m_prime: u32,
}

impl HyperLogLogPlusMinus for HLLSparse {
    fn new(precision: u8, bit_mixer: fn(u64) -> u64) -> Self {
        let m = 1usize << precision;
        let mut set = HashSet::new();
        set.reserve(m / 4);

        Self {
            n_observed: 0,
            p: precision,
            sparse_list: set,
            bit_mixer,
            m,
            p_prime: 25,
            m_prime: 1 << 25,
        }
    }

    fn insert(&mut self, item: u64) -> Result<(), Box<dyn Error>> {
        // 自增观察到的项目数量
        self.n_observed += 1;

        // 计算哈希值
        let hash_value = (self.bit_mixer)(item);

        // 如果稀疏列表变得太大，切换到常规表示
        if self.sparse_list.len() + 1 > self.m / 4 {
            // 实现 switchToNormalRepresentation 方法
            return Err("Switching to non-sparse representation".into());
        }

        // 稀疏模式下添加哈希值到稀疏列表
        let encoded_hash_value = ertl::encode_hash_in_32bit(hash_value, self.p_prime, self.p);
        self.sparse_list.insert(encoded_hash_value);

        Ok(())
    }

    fn cardinality(&mut self) -> u64 {
        let q = 64 - self.p_prime;
        let m = self.m_prime;
        let c = ertl::sparse_register_histogram(&self.sparse_list, self.p_prime, self.p, q);
        0
    }
}

struct HLLDense {
    p: u8,
    m_registers: Vec<u8>,
    bit_mixer: fn(u64) -> u64,
    n_observed: u64,
}

impl HyperLogLogPlusMinus for HLLDense {
    fn new(precision: u8, bit_mixer: fn(u64) -> u64) -> Self {
        let m = 1usize << precision;
        Self {
            p: precision,
            m_registers: vec![0; m],
            bit_mixer,
            n_observed: 0,
        }
    }

    /// Inserts a new item into the HyperLogLog data structure.
    ///
    /// # Arguments
    /// * `item` - The item to be inserted.
    fn insert(&mut self, item: u64) -> Result<(), Box<dyn Error>> {
        self.n_observed += 1;

        let hash_value = (self.bit_mixer)(item);

        let idx = ertl::get_index_u64(hash_value, self.p) as usize;
        let rank = ertl::get_rank_u64(hash_value, self.p);

        if rank > self.m_registers[idx] {
            self.m_registers[idx] = rank;
        }

        Ok(())
    }

    fn cardinality(&mut self) -> u64 {
        0
    }
}
