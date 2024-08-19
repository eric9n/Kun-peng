use clap::{error::ErrorKind, Error, Parser};
use hyperloglogplus::{HyperLogLog, HyperLogLogPlus};
use kr2r::args::KLMTArgs;
use kr2r::utils::{find_files, format_bytes, open_file};
use kr2r::KBuildHasher;

use seqkmer::{read_parallel, BufferFastaReader};
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

fn build_output_path<P: AsRef<Path>>(input_path: &P, extension: &str) -> String {
    let path = input_path.as_ref();
    let parent_dir = path.parent().unwrap_or_else(|| Path::new(""));
    let stem = path.file_stem().unwrap_or_else(|| path.as_os_str());

    let mut output_path = parent_dir.join(stem);
    output_path.set_extension(extension);

    output_path.to_str().unwrap().to_owned()
}

fn process_sequence<P: AsRef<Path>>(
    fna_file: &P,
    // hllp: &mut HyperLogLogPlus<u64, KBuildHasher>,
    args: Args,
) -> HyperLogLogPlus<u64, KBuildHasher> {
    // 构建预期的 JSON 文件路径
    let json_path = build_output_path(fna_file, &format!("hllp_{}.json", args.n));
    // 检查是否存在 JSON 文件
    if args.cache && Path::new(&json_path).exists() {
        // 如果存在，从文件读取并反序列化
        let mut file = open_file(json_path).unwrap();
        let mut serialized_hllp = String::new();
        file.read_to_string(&mut serialized_hllp).unwrap();
        let hllp: HyperLogLogPlus<u64, KBuildHasher> =
            serde_json::from_str(&serialized_hllp).unwrap();

        return hllp;
    }

    let meros = args.klmt.as_meros();

    let mut hllp: HyperLogLogPlus<u64, _> =
        HyperLogLogPlus::new(16, KBuildHasher::default()).unwrap();

    let mut reader = BufferFastaReader::from_path(fna_file, 1)
        .expect("Failed to open the FASTA file with FastaReader");
    let range_n = args.n as u64;
    read_parallel(
        &mut reader,
        args.threads,
        &meros,
        |record_set| {
            let mut minimizer_set = HashSet::new();

            for record in record_set {
                record.body.apply_mut(|m_iter| {
                    let kmer_iter: HashSet<u64> = m_iter
                        .filter(|(_, hash_key)| *hash_key & RANGE_MASK < range_n)
                        .map(|(_, hash_key)| hash_key)
                        .collect();

                    minimizer_set.extend(kmer_iter);
                });
            }
            minimizer_set
        },
        |record_sets| {
            while let Some(data) = record_sets.next() {
                let m_set = data.unwrap();
                for minimizer in m_set {
                    hllp.insert(&minimizer);
                }
            }
        },
    )
    .expect("read parallel error");

    // 序列化 hllp 对象并将其写入文件
    let serialized_hllp = serde_json::to_string(&hllp).unwrap();

    if let Ok(mut file) = File::create(&json_path) {
        // 尝试写入数据
        if let Err(e) = file.write_all(serialized_hllp.as_bytes()) {
            eprintln!("Failed to write to file: {}", e);
        }
    } else {
        eprintln!("Failed to create file: {}", json_path);
    }

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
        vec![source.clone()]
    } else {
        let library_dir = &args.database.join("library");
        find_files(library_dir, "library", ".fna")
    };

    if fna_files.is_empty() {
        panic!("Error: No library.fna files found in the specified directory. Please ensure that the directory contains at least one library.fna file and try again.");
    }

    println!("estimate start... ");

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
