use clap::{error::ErrorKind, Error, Parser};
use hyperloglogplus::{HyperLogLog, HyperLogLogPlus};
use kr2r::mmscanner::{MinimizerScanner, BITS_PER_CHAR, DEFAULT_SPACED_SEED_MASK};
use kr2r::utils::{expand_spaced_seed_mask, find_library_fna_files};
use kr2r::KBuildHasher;
use seq_io::fasta::{Reader, Record};
use seq_io::parallel::read_parallel;
use std::path::PathBuf;
use std::sync::{
    atomic::{AtomicUsize, Ordering},
    Arc, Mutex,
};

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
    let hllp: HyperLogLogPlus<u64, _> = HyperLogLogPlus::new(16, KBuildHasher::default()).unwrap();

    let hllp = Arc::new(Mutex::new(hllp));
    let counter = Arc::new(AtomicUsize::new(0)); // 初始化原子计数器

    for fna_file in fna_files {
        println!("fna_file {:?}", fna_file);
        let reader = Reader::from_path(fna_file).unwrap();
        let counter_clone = counter.clone();
        read_parallel(
            reader,
            4,
            2,
            |record_set| {
                for record in record_set.into_iter() {
                    let k_mer = args.k_mer as usize;
                    let l_mer = args.l_mer as usize;
                    let mut scranner =
                        MinimizerScanner::default(record.seq().to_vec(), k_mer, l_mer);
                    scranner.set_spaced_seed_mask(args.spaced_seed_mask);
                    if let Some(toggle_mask) = args.toggle_mask {
                        scranner.set_toggle_mask(toggle_mask);
                    }
                    while let Some(minimizer) = scranner.next_minimizer() {
                        let mut hllp_clone = hllp.lock().unwrap();
                        hllp_clone.insert(&minimizer);
                        counter_clone.fetch_add(1, Ordering::SeqCst); // 递增计数器
                    }
                }
            },
            |_| {},
        );
    }

    let final_count = counter.load(Ordering::SeqCst); // 读取计数器的最终值

    let mut hllp_clone = hllp.lock().unwrap();
    let hllp_count = hllp_clone.count();
    println!("Final count: {:?}", final_count);
    println!("HLLP count: {:?}", hllp_count);
}
