// 使用时需要引用模块路径
use crate::utils::expand_spaced_seed_mask;
use crate::{construct_seed_template, parse_binary};
use clap::Parser;
use seqkmer::Meros;
use seqkmer::{
    BITS_PER_CHAR, DEFAULT_KMER_LENGTH, DEFAULT_MINIMIZER_LENGTH, DEFAULT_MINIMIZER_SPACES,
    DEFAULT_TOGGLE_MASK,
};
use std::path::PathBuf;

pub const U32MAXPLUS: u64 = u32::MAX as u64;
pub const ONEGB: u64 = 1073741824;

#[derive(Parser, Debug, Clone)]
#[clap(version, about = "build database")]
pub struct Build {
    /// ncbi library fna database directory
    #[arg(long = "db", required = true)]
    pub database: PathBuf,

    // /// Kraken 2 options filename, default = $database/opts.k2d
    // #[clap(short = 'o')]
    // pub options_filename: Option<PathBuf>,
    /// 包含原始配置
    #[clap(flatten)]
    pub klmt: KLMTArgs,

    /// Bit storage requested for taxid 0 <= r < 31
    #[clap(short, long, value_parser = clap::value_parser!(u8).range(0..31), default_value_t = 0)]
    pub requested_bits_for_taxid: u8,

    /// Number of threads
    #[clap(short = 'p', long, default_value_t = num_cpus::get())]
    pub threads: usize,
}

const BATCH_SIZE: usize = 16 * 1024 * 1024;

/// Command line arguments for the classify program.
///
/// This structure defines the command line arguments that are accepted by the classify program.
/// It uses the `clap` crate for parsing command line arguments.
/// combines the functionality of the 'splitr', 'annotate', and 'resolve' commands into a single workflow.
/// This command streamlines the process of splitting fast(q/a) files, annotating sequences, and resolving the taxonomy tree,
/// providing a comprehensive solution for sequence classification.
#[derive(Parser, Debug, Clone)]
#[clap(
    version,
    about = "Integrates 'splitr', 'annotate', and 'resolve' into a unified workflow for sequence classification. classify a set of sequences",
    long_about = "classify a set of sequences"
)]
pub struct ClassifyArgs {
    // /// database hash chunk directory and other files
    #[arg(long = "db", required = true)]
    pub database: PathBuf,

    /// chunk directory
    #[clap(long)]
    pub chunk_dir: PathBuf,

    /// File path for outputting normal Kraken output.
    #[clap(long = "output-dir", value_parser)]
    pub kraken_output_dir: Option<PathBuf>,

    /// Enable paired-end processing.
    #[clap(short = 'P', long = "paired-end-processing", action)]
    pub paired_end_processing: bool,

    /// Process pairs with mates in the same file.
    #[clap(short = 'S', long = "single-file-pairs", action)]
    pub single_file_pairs: bool,

    /// Minimum quality score for FASTQ data.
    #[clap(
        short = 'Q',
        long = "minimum-quality-score",
        value_parser,
        default_value_t = 0
    )]
    pub minimum_quality_score: i32,

    /// The number of threads to use.
    #[clap(short = 'p', long = "num-threads", value_parser, default_value_t = num_cpus::get())]
    pub num_threads: usize,

    #[clap(long, default_value_t = BATCH_SIZE)]
    pub batch_size: usize,

    /// Confidence score threshold
    #[clap(
        short = 'T',
        long = "confidence-threshold",
        value_parser,
        default_value_t = 0.0
    )]
    pub confidence_threshold: f64,

    /// The minimum number of hit groups needed for a call.
    #[clap(
        short = 'g',
        long = "minimum-hit-groups",
        value_parser,
        default_value_t = 2
    )]
    pub minimum_hit_groups: usize,

    /// Enables use of a Kraken 2 compatible shared database.
    #[clap(long, default_value_t = false)]
    pub kraken_db_type: bool,

    /// In comb. w/ -R, provide minimizer information in report
    #[clap(short = 'K', long, value_parser, default_value_t = false)]
    pub report_kmer_data: bool,

    /// In comb. w/ -R, report taxa w/ 0 count
    #[clap(short = 'z', long, value_parser, default_value_t = false)]
    pub report_zero_counts: bool,

    /// output file contains all unclassified sequence
    #[clap(long, value_parser, default_value_t = false)]
    pub full_output: bool,

    /// A list of input file paths (FASTA/FASTQ) to be processed by the classify program.
    // #[clap(short = 'F', long = "files")]
    pub input_files: Vec<String>,
}

#[derive(Parser, Debug, Clone, Copy)]
#[clap(version, about = "k-mer")]
pub struct KLMTArgs {
    /// Set length of k-mers, k must be positive integer, k=35, k cannot be less than l
    #[clap(short, long, value_parser = clap::value_parser!(u64).range(1..), default_value_t = DEFAULT_KMER_LENGTH)]
    pub k_mer: u64,

    /// Set length of minimizers, 1 <= l <= 31
    #[clap(short, long, value_parser = clap::value_parser!(u8).range(1..=31), default_value_t = DEFAULT_MINIMIZER_LENGTH)]
    pub l_mer: u8,

    // /// Spaced seed mask
    // #[clap(short = 'S', long, default_value= "0", value_parser = parse_binary)]
    // spaced_seed_mask: u64,
    /// Number of characters in minimizer that are ignored in comparisons
    #[clap(long, default_value_t = DEFAULT_MINIMIZER_SPACES)]
    pub minimizer_spaces: u8,

    /// Minimizer ordering toggle mask
    #[clap(short = 'T', long, default_value_t = DEFAULT_TOGGLE_MASK)]
    pub toggle_mask: u64,

    #[clap(long)]
    pub min_clear_hash_value: Option<u64>,
}

impl KLMTArgs {
    pub fn as_meros(&self) -> Meros {
        let seed = construct_seed_template(self.l_mer as usize, self.minimizer_spaces as usize);
        let space_seed_mask = parse_binary(&seed).unwrap();
        let space_seed_mask = expand_spaced_seed_mask(space_seed_mask, BITS_PER_CHAR as u64);

        Meros::new(
            self.k_mer as usize,
            self.l_mer as usize,
            Some(space_seed_mask),
            Some(self.toggle_mask),
            self.min_clear_hash_value,
        )
    }
}

pub fn parse_size(s: &str) -> Result<usize, String> {
    let len = s.len();
    if len < 2 {
        return Err("Size must be at least two characters".to_string());
    }

    let (num, suffix) = s.split_at(len - 1);
    let number: f64 = num.parse().map_err(|_| "Invalid number".to_string())?;
    match suffix {
        "G" | "g" => Ok((number * 1_073_741_824.0) as usize), // 2^30
        "M" | "m" => Ok((number * 1_048_576.0) as usize),     // 2^20
        "K" | "k" => Ok((number * 1_024.0) as usize),         // 2^10
        _ => Err("Invalid size suffix. Use 'G', 'M', or 'K'".to_string()),
    }
}
