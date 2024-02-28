use hyperloglogplus::{HyperLogLog, HyperLogLogPlus};
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::hash::BuildHasher;

use crate::KBuildHasher;

type TaxId = u32;
pub const TAXID_MAX: TaxId = TaxId::MAX;
// 定义 taxon_counts_t 类型
pub type TaxonCounts = HashMap<TaxId, u64>;

// 定义一个错误类型
#[derive(Debug)]
pub struct UnionError;

impl fmt::Display for UnionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Failed to union HyperLogLogPlus instances")
    }
}

pub trait Unionable {
    type Error;
    fn union(&mut self, other: &Self) -> Result<(), Self::Error>;
    fn distinct_count(&mut self) -> usize;
    fn add_kmer(&mut self, kmer: u64);
}

impl<B> Unionable for HyperLogLogPlus<u64, B>
where
    // H: Hash + ?Sized,
    B: BuildHasher,
{
    type Error = UnionError;

    fn union(&mut self, other: &Self) -> Result<(), Self::Error> {
        self.merge(other).map_err(|_| UnionError)
    }

    fn distinct_count(&mut self) -> usize {
        self.count().round() as usize
    }

    fn add_kmer(&mut self, kmer: u64) {
        self.insert(&kmer);
    }
}

impl Unionable for HashSet<u64> {
    type Error = UnionError;

    fn union(&mut self, other: &Self) -> Result<(), Self::Error> {
        self.extend(other.iter().cloned());
        Ok(())
    }

    fn distinct_count(&mut self) -> usize {
        self.len()
    }

    fn add_kmer(&mut self, kmer: u64) {
        self.insert(kmer);
    }
}

pub struct ReadCounts<T>
where
    T: Unionable,
{
    n_reads: u64,
    n_kmers: u64,
    kmers: T,
}

impl<T> ReadCounts<T>
where
    T: Unionable,
{
    // pub fn new(kmers: T) -> Self {
    //     ReadCounts {
    //         n_reads: 0,
    //         n_kmers: 0,
    //         kmers,
    //     }
    // }

    pub fn with_capacity(kmers: T, n_reads: u64, n_kmers: u64) -> Self {
        ReadCounts {
            n_reads,
            n_kmers,
            kmers, // kmers: T::with_capacity(n_kmers as usize),
        }
    }

    pub fn read_count(&self) -> u64 {
        self.n_reads
    }

    pub fn increment_read_count(&mut self) {
        self.n_reads += 1;
    }

    pub fn kmer_count(&self) -> u64 {
        self.n_kmers
    }

    pub fn distinct_kmer_count(&mut self) -> usize {
        self.kmers.distinct_count()
    }

    pub fn add_kmer(&mut self, kmer: u64) {
        self.n_kmers += 1;
        self.kmers.add_kmer(kmer);
    }

    pub fn merge(&mut self, other: &ReadCounts<T>) -> Result<(), UnionError> {
        self.n_reads += other.n_reads;
        self.n_kmers += other.n_kmers;
        self.kmers.union(&other.kmers).map_err(|_| UnionError)
    }
}

// 定义 READCOUNTER 类型
#[cfg(feature = "exact_counting")]
pub type ReadCounter = ReadCounts<HashSet<u64>>;

#[cfg(feature = "exact_counting")]
impl Default for ReadCounter {
    fn default() -> Self {
        let kmers: HashSet<u64> = HashSet::new();
        ReadCounts::new(kmers)
    }
}

#[cfg(not(feature = "exact_counting"))]
pub type ReadCounter = ReadCounts<HyperLogLogPlus<u64, KBuildHasher>>;

#[cfg(not(feature = "exact_counting"))]
impl Default for ReadCounter {
    fn default() -> Self {
        let kmers: HyperLogLogPlus<u64, KBuildHasher> =
            HyperLogLogPlus::new(16, KBuildHasher::default()).unwrap();
        ReadCounts::with_capacity(kmers, 0, 0)
    }
}

#[cfg(not(feature = "exact_counting"))]
impl ReadCounter {
    pub fn new(n_reads: u64, n_kmers: u64) -> Self {
        let kmers: HyperLogLogPlus<u64, KBuildHasher> =
            HyperLogLogPlus::new(16, KBuildHasher::default()).unwrap();
        ReadCounts::with_capacity(kmers, n_reads, n_kmers)
    }
}

pub type TaxonCounters = HashMap<u64, ReadCounter>;
