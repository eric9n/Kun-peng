use clap::{error::ErrorKind, Error, Parser};
use hyperloglogplus::{HyperLogLog, HyperLogLogPlus};
use kr2r::mmscanner::MinimizerScanner;
use kr2r::utils::{expand_spaced_seed_mask, find_library_fna_files};
use kr2r::{fmix64 as murmur_hash3, KBuildHasher};
use kr2r::{Meros, BITS_PER_CHAR, DEFAULT_SPACED_SEED_MASK};
use seq_io::fasta::{Reader, Record};
use seq_io::parallel::read_parallel;
use serde_json;
use std::collections::HashSet;
use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

#[derive(Parser, Debug, Clone)]
#[clap(
    version,
    about = "estimate capacity",
    long_about = "Estimates the size of the Kraken 2 hash table."
)]
struct Args {
    /// build database directory or file
    #[arg(long, default_value = "lib")]
    source: PathBuf,

    /// estimate capacity from cache if exists
    #[arg(long, default_value = "false")]
    cache: bool,

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

fn build_output_path(input_path: &str, extension: &str) -> String {
    let path = Path::new(input_path);
    let parent_dir = path.parent().unwrap_or_else(|| Path::new(""));
    let stem = path.file_stem().unwrap_or_else(|| path.as_os_str());

    let mut output_path = parent_dir.join(stem);
    output_path.set_extension(extension);

    output_path.to_str().unwrap().to_owned()
}

fn process_sequence(
    fna_file: &str,
    // hllp: &mut HyperLogLogPlus<u64, KBuildHasher>,
    args: Args,
) -> HyperLogLogPlus<u64, KBuildHasher> {
    // 构建预期的 JSON 文件路径
    let json_path = build_output_path(fna_file, "hllp.json");
    // 检查是否存在 JSON 文件
    if args.cache && Path::new(&json_path).exists() {
        // 如果存在，从文件读取并反序列化
        let mut file = File::open(json_path).unwrap();
        let mut serialized_hllp = String::new();
        file.read_to_string(&mut serialized_hllp).unwrap();
        let hllp: HyperLogLogPlus<u64, KBuildHasher> =
            serde_json::from_str(&serialized_hllp).unwrap();

        return hllp;
    }

    let k_mer = args.k_mer as usize;
    let l_mer = args.l_mer as usize;
    let mut hllp: HyperLogLogPlus<u64, _> =
        HyperLogLogPlus::new(16, KBuildHasher::default()).unwrap();

    let reader = Reader::from_path(fna_file).unwrap();
    read_parallel(
        reader,
        args.threads as u32,
        args.threads - 2 as usize,
        |record_set| {
            let meros = Meros::new(
                k_mer,
                l_mer,
                Some(args.spaced_seed_mask),
                args.toggle_mask,
                None,
            );

            let mut scanner = MinimizerScanner::new(meros);

            let mut minimizer_set = HashSet::new();
            for record in record_set.into_iter() {
                let seq = record.seq();
                scanner.set_seq_end(seq);
                while let Some(minimizer) = scanner.next_minimizer(seq) {
                    let hash_v = murmur_hash3(minimizer);
                    if hash_v & RANGE_MASK < args.n as u64 {
                        minimizer_set.insert(hash_v);
                    }
                }
                scanner.reset();
            }
            minimizer_set
        },
        |record_sets| {
            while let Some(Ok((_, m_set))) = record_sets.next() {
                for minimizer in m_set {
                    hllp.insert(&minimizer);
                }
            }
        },
    );

    // 序列化 hllp 对象并将其写入文件
    let serialized_hllp = serde_json::to_string(&hllp).unwrap();
    let mut file = File::create(&json_path).unwrap();
    file.write_all(serialized_hllp.as_bytes()).unwrap();

    hllp
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

    let mut hllp: HyperLogLogPlus<u64, KBuildHasher> =
        HyperLogLogPlus::new(16, KBuildHasher::default()).unwrap();

    let source: PathBuf = args.source.clone();
    let fna_files = if source.is_file() {
        vec![source.to_string_lossy().to_string()]
    } else {
        find_library_fna_files(args.source)
    };

    for fna_file in fna_files {
        println!("fna_file {:?}", fna_file);
        let args_clone = Args {
            source: source.clone(),
            ..args
        };
        let local_hllp = process_sequence(&fna_file, args_clone);
        if let Err(e) = hllp.merge(&local_hllp) {
            println!("hllp merge err {:?}", e);
        }
    }

    // let final_count = counter.load(Ordering::SeqCst); // 读取计数器的最终值

    let hllp_count = (hllp.count() * RANGE_SECTIONS as f64 / args.n as f64).round() as u64;
    // println!("Final count: {:?}", final_count);
    println!("estimate count: {:?}", hllp_count);
}
