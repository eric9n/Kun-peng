use seahash::SeaHasher;
use std::hash::{BuildHasher, Hasher};

// 定义 KeyValueStore trait
trait KeyValueStore {
    fn get(&self, key: u64) -> u32;
}

// 声明常量
const M1: u64 = 0xff51afd7ed558ccd;
const M2: u64 = 0xc4ceb9fe1a85ec53;

///
/// # Examples
///
/// ```
/// # use kr2r::fmix64;
/// let key: u64 = 123;
/// let hash = fmix64(key);
/// assert_eq!(hash, 9208534749291869864);
/// ```
#[inline]
pub fn fmix64(key: u64) -> u64 {
    let mut k = key;
    k ^= k >> 33;
    k = k.wrapping_mul(M1);
    k ^= k >> 33;
    k = k.wrapping_mul(M2);
    k ^= k >> 33;
    k
}

const C1: u64 = 0x87c37b91114253d5;
const C2: u64 = 0x4cf5ad432745937f;

/// MurmurHash3 函数的优化 Rust 实现
/// keys = u64.to_be_bytes();这个函数只能处理u64.
/// Performs the MurmurHash3 hash computation on a given key.
#[inline]
pub fn murmur_hash3(key: u64) -> u64 {
    let k1 = key.wrapping_mul(C1).rotate_left(31).wrapping_mul(C2);

    let mut h1 = k1;
    h1 = h1.rotate_left(27);
    h1 = h1.wrapping_mul(5).wrapping_add(0x52dce729);

    // 最终处理
    h1 ^= 8u64;
    h1 = fmix64(h1);

    h1
}

#[inline]
pub fn sea_hash(key: u64) -> u64 {
    let mut hasher = SeaHasher::default();
    hasher.write_u64(key);
    hasher.finish()
}

/// KHasher is designed to hash only u64 values.
///
/// # Safety
/// This hasher assumes that the input slice to `write` is exactly 8 bytes long,
/// representing a `u64`. Using this hasher with input that is not 8 bytes,
/// or not properly representing a `u64`, may lead to undefined behavior including
/// but not limited to memory safety violations
#[derive(Default)]
pub struct KHasher {
    hash: u64,
}

impl Hasher for KHasher {
    // fn write(&mut self, _: &[u8]) {
    // #[cfg(debug_assertions)]
    // assert_eq!(bytes.len(), 8, "KHasher input must be exactly 8 bytes.");

    // // 使用 unsafe 从字节切片直接读取 u64
    // let key: u64 = unsafe { std::ptr::read_unaligned(bytes.as_ptr() as *const u64) };

    // self.hash = murmur_hash3(key);
    // }

    // #[inline]
    // fn finish(&self) -> u64 {
    //     self.0
    // }

    #[inline]
    fn write(&mut self, _: &[u8]) {}

    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.hash = i;
    }

    #[inline]
    fn finish(&self) -> u64 {
        self.hash
    }
}

#[derive(Default)]
pub struct KBuildHasher;

impl BuildHasher for KBuildHasher {
    type Hasher = KHasher;

    fn build_hasher(&self) -> Self::Hasher {
        KHasher { hash: 0 }
    }
}

#[derive(Default)]
pub struct SBuildHasher;

impl BuildHasher for SBuildHasher {
    type Hasher = SeaHasher;

    fn build_hasher(&self) -> Self::Hasher {
        SeaHasher::default()
    }
}
