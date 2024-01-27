use criterion::{black_box, criterion_group, criterion_main, Criterion};
use kr2r::{fmix64, murmur_hash3, sea_hash};
use std::hash::Hasher;
use twox_hash::xxh3;
extern crate farmhash;

#[inline]
pub fn xx_hash(key: u64) -> u64 {
    let mut xhash = xxh3::Hash64::default();
    xhash.write_u64(key);
    xhash.finish()
}

#[inline]
pub fn farm(key: u64) -> u64 {
    // let bytes = key.to_be_bytes();
    // let byte_slce: &[u8] = &bytes;
    farmhash::hash64(&key.to_be_bytes())
}

fn criterion_benchmark(c: &mut Criterion) {
    let key = 0x12345678abcdef01u64;

    c.bench_function("fmix64", |b| b.iter(|| fmix64(black_box(key))));
    c.bench_function("murmur_hash3", |b| b.iter(|| murmur_hash3(black_box(key))));
    c.bench_function("sea_hash", |b| b.iter(|| sea_hash(black_box(key))));
    c.bench_function("xx_hash", |b| b.iter(|| xx_hash(black_box(key))));
    c.bench_function("farm", |b| b.iter(|| farm(black_box(key))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
