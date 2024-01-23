use clap::{error::ErrorKind, Error, Parser};
use hyperloglogplus::{HyperLogLog, HyperLogLogPlus};
use kr2r::mmscanner::{MinimizerScanner, BITS_PER_CHAR, DEFAULT_SPACED_SEED_MASK};
use kr2r::utils::{expand_spaced_seed_mask, find_library_fna_files};
use kr2r::{murmur_hash3, KBuildHasher};
use seq_io::fasta::{Reader, Record};
use seq_io::parallel::read_parallel;
use std::collections::HashSet;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

#[derive(Parser, Debug)]
#[clap(
    version,
    about = "estimate capacity",
    long_about = "Estimates the size of the Kraken 2 hash table."
)]
struct Args {
    /// 构建数据库的目录
    #[arg(long = "db", default_value = "lib")]
    database: PathBuf,

    /// Set length of k-mers, k must be positive integer, k=35, k cannot be less than l
    #[clap(short, long, value_parser = clap::value_parser!(u64).range(1..), required = true)]
    k_mer: u64,

    /// Set length of minimizers, 1 <= l <= 31
    #[clap(short, long, value_parser = clap::value_parser!(u8).range(1..=31), required = true)]
    l_mer: u8,

    /// Set maximum qualifying hash code
    #[clap(short, long, default_value = "4")]
    n: usize,

    /// Spaced seed mask
    #[clap(short = 'S', long, default_value= "0", value_parser = parse_binary)]
    spaced_seed_mask: u64,

    /// Minimizer ordering toggle mask
    #[clap(short = 'T', long, value_parser = parse_binary)]
    toggle_mask: Option<u64>,

    /// Read block size
    #[clap(short = 'B', long, default_value = "31457280")]
    block_size: usize,

    /// Number of threads
    #[clap(short = 'p', long, default_value = "4")]
    threads: usize,
}

fn parse_binary(src: &str) -> Result<u64, std::num::ParseIntError> {
    u64::from_str_radix(src, 2)
}

const RANGE_SECTIONS: u64 = 1024;
const RANGE_MASK: u64 = RANGE_SECTIONS - 1;

fn main() {
    let mut args = Args::parse();
    if args.k_mer < args.l_mer as u64 {
        let err = Error::raw(ErrorKind::ValueValidation, "k cannot be less than l");
        err.exit();
    }
    if args.spaced_seed_mask != DEFAULT_SPACED_SEED_MASK {
        args.spaced_seed_mask =
            expand_spaced_seed_mask(args.spaced_seed_mask, BITS_PER_CHAR as u64);
    }
    let fna_files = find_library_fna_files(args.database);
    let mut hllp: HyperLogLogPlus<u64, _> =
        HyperLogLogPlus::new(16, KBuildHasher::default()).unwrap();

    let counter = AtomicUsize::new(0); // 初始化原子计数器

    // let sets: HashSet<u64> = HashSet::new();
    for fna_file in fna_files {
        println!("fna_file {:?}", fna_file);
        let reader = Reader::from_path(fna_file).unwrap();
        let k_mer = args.k_mer as usize;
        let l_mer = args.l_mer as usize;

        read_parallel(
            reader,
            args.threads as u32,
            args.threads - 2 as usize,
            |record_set| {
                let mut count = 0;
                let mut scanner = MinimizerScanner::default(k_mer, l_mer);
                scanner.set_spaced_seed_mask(args.spaced_seed_mask);
                if let Some(toggle_mask) = args.toggle_mask {
                    scanner.set_toggle_mask(toggle_mask);
                }
                let mut minimizer_set = HashSet::new();
                for record in record_set.into_iter() {
                    let seq = record.seq();
                    scanner.set_seq_end(seq);
                    while let Some(minimizer) = scanner.next_minimizer(seq) {
                        count += 1;
                        let hash_v = murmur_hash3(minimizer);
                        if hash_v & RANGE_MASK < args.n as u64 {
                            minimizer_set.insert(minimizer);
                        }
                    }
                    scanner.reset();
                }
                (minimizer_set, count)
            },
            |record_sets| {
                while let Some(Ok((_, (m_set, count)))) = record_sets.next() {
                    for minimizer in m_set {
                        hllp.insert(&minimizer);
                    }
                    // sets.extend(m_set);

                    counter.fetch_add(count, Ordering::SeqCst);
                }
            },
        );
    }

    let final_count = counter.load(Ordering::SeqCst); // 读取计数器的最终值

    let hllp_count = hllp.count();
    // println!("sets {:?}", sets.len() * 1024 / args.n);
    println!("Final count: {:?}", final_count);
    println!("HLLP count: {:?}", hllp_count * 1024f64 / args.n as f64);
}
