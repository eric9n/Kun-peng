use clap::{error::ErrorKind, Error, Parser};
use hyperloglogplus::{HyperLogLog, HyperLogLogPlus};
use kr2r::args::KLMTArgs;
use kr2r::mmscanner::MinimizerScanner;
use kr2r::utils::{find_library_fna_files, format_bytes};
use kr2r::KBuildHasher;
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
pub struct Args {
    /// build database directory or file
    #[arg(long, default_value = "lib")]
    pub database: PathBuf,

    /// 包含原始配置
    #[clap(flatten)]
    pub klmt: KLMTArgs,

    /// estimate capacity from cache if exists
    #[arg(long, default_value_t = true)]
    pub cache: bool,

    /// Set maximum qualifying hash code
    #[clap(short, long, default_value = "4")]
    pub n: usize,

    /// Proportion of the hash table to be populated
    /// (build task only; def: 0.7, must be
    ///    between 0 and 1).
    #[clap(long, long, default_value_t = 0.7)]
    pub load_factor: f64,

    /// Number of threads
    #[clap(short = 'p', long, default_value_t = 10)]
    pub threads: usize,
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

    let meros = args.klmt.as_meros();

    let mut hllp: HyperLogLogPlus<u64, _> =
        HyperLogLogPlus::new(16, KBuildHasher::default()).unwrap();

    let reader = Reader::from_path(fna_file).unwrap();
    let range_n = args.n as u64;
    read_parallel(
        reader,
        args.threads as u32,
        args.threads - 2 as usize,
        |record_set| {
            let mut minimizer_set = HashSet::new();
            for record in record_set.into_iter() {
                let seq = record.seq();
                let kmer_iter = MinimizerScanner::new(&seq, meros)
                    .into_iter()
                    .filter(|hash_key| hash_key & RANGE_MASK < range_n)
                    .collect::<HashSet<u64>>();

                minimizer_set.extend(kmer_iter);
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

pub fn run(args: Args) -> usize {
    let meros = args.klmt.as_meros();

    if meros.k_mer < meros.l_mer {
        let err = Error::raw(ErrorKind::ValueValidation, "k cannot be less than l");
        err.exit();
    }

    let mut hllp: HyperLogLogPlus<u64, KBuildHasher> =
        HyperLogLogPlus::new(16, KBuildHasher::default()).unwrap();

    let source: PathBuf = args.database.clone();
    let fna_files = if source.is_file() {
        vec![source.to_string_lossy().to_string()]
    } else {
        find_library_fna_files(args.database)
    };

    for fna_file in fna_files {
        let args_clone = Args {
            database: source.clone(),
            ..args
        };
        let local_hllp = process_sequence(&fna_file, args_clone);
        if let Err(e) = hllp.merge(&local_hllp) {
            println!("hllp merge err {:?}", e);
        }
    }

    let hllp_count = (hllp.count() * RANGE_SECTIONS as f64 / args.n as f64).round() as u64;
    let required_capacity = (hllp_count + 8192) as f64 / args.load_factor;
    println!(
        "estimate count: {:?}, required capacity: {:?}, Estimated hash table requirement: {:}",
        hllp_count,
        required_capacity.ceil(),
        format_bytes(required_capacity * 4f64)
    );
    required_capacity.ceil() as usize
}

#[allow(dead_code)]
fn main() {
    let args = Args::parse();
    run(args);
}
