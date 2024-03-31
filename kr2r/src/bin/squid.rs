use clap::{Parser, Subcommand};
mod annotate;
mod hashshard;
mod resolve;
mod splitr;

use kr2r::utils::find_and_sort_files;
use std::io::Result;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    #[clap(subcommand)]
    cmd: Commands,
}

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
    about = "Integrates 'splitr', 'annotate', and 'resolve' into a unified workflow for sequence classification., classify a set of sequences",
    long_about = "classify a set of sequences"
)]
pub struct ClassifyArgs {
    /// database hash chunk directory and other files
    #[clap(long)]
    hash_dir: PathBuf,

    // /// The file path for the Kraken 2 options.
    // #[clap(short = 'o', long = "options-filename", value_parser, required = true)]
    // options_filename: String,
    /// Enable paired-end processing.
    #[clap(short = 'P', long = "paired-end-processing", action)]
    paired_end_processing: bool,

    /// Process pairs with mates in the same file.
    #[clap(short = 'S', long = "single-file-pairs", action)]
    single_file_pairs: bool,

    /// Minimum quality score for FASTQ data, default is 0.
    #[clap(
        short = 'Q',
        long = "minimum-quality-score",
        value_parser,
        default_value_t = 0
    )]
    minimum_quality_score: i32,

    /// The number of threads to use, default is 1.
    #[clap(short = 'p', long = "num-threads", value_parser, default_value_t = 10)]
    num_threads: i32,

    /// chunk directory
    #[clap(long)]
    chunk_dir: PathBuf,

    /// 批量处理大小 default: 8MB
    #[clap(long, default_value_t = annotate::BATCH_SIZE)]
    batch_size: usize,

    /// Confidence score threshold, default is 0.0.
    #[clap(
        short = 'T',
        long = "confidence-threshold",
        value_parser,
        default_value_t = 0.0
    )]
    confidence_threshold: f64,

    /// The minimum number of hit groups needed for a call.
    #[clap(
        short = 'g',
        long = "minimum-hit-groups",
        value_parser,
        default_value_t = 2
    )]
    minimum_hit_groups: usize,

    /// File path for outputting normal Kraken output.
    #[clap(long = "output-dir", value_parser)]
    kraken_output_dir: Option<PathBuf>,

    /// A list of input file paths (FASTA/FASTQ) to be processed by the classify program.
    // #[clap(short = 'F', long = "files")]
    input_files: Vec<String>,
}

impl From<ClassifyArgs> for splitr::Args {
    fn from(item: ClassifyArgs) -> Self {
        Self {
            hash_dir: item.hash_dir,
            paired_end_processing: item.paired_end_processing,
            single_file_pairs: item.single_file_pairs,
            minimum_quality_score: item.minimum_quality_score,
            num_threads: item.num_threads,
            chunk_dir: item.chunk_dir,
            input_files: item.input_files,
        }
    }
}

impl From<ClassifyArgs> for annotate::Args {
    fn from(item: ClassifyArgs) -> Self {
        Self {
            hash_dir: item.hash_dir,
            chunk_dir: item.chunk_dir,
            batch_size: item.batch_size,
        }
    }
}

impl From<ClassifyArgs> for resolve::Args {
    fn from(item: ClassifyArgs) -> Self {
        Self {
            hash_dir: item.hash_dir,
            chunk_dir: item.chunk_dir,
            batch_size: item.batch_size,
            confidence_threshold: item.confidence_threshold,
            minimum_hit_groups: item.minimum_hit_groups,
            kraken_output_dir: item.kraken_output_dir,
        }
    }
}

#[derive(Subcommand, Debug)]
enum Commands {
    Hashshard(hashshard::Args),
    Splitr(splitr::Args),
    Annotate(annotate::Args),
    Resolve(resolve::Args),
    Classify(ClassifyArgs),
}

fn main() -> Result<()> {
    let args = Args::parse();
    // let matches = Command::new("squid")
    //     .version("2.1")
    //     .author("Eric <eric9n@gmail.com>")
    //     .about("classify a set of sequences like Kraken 2")
    //     // .subcommand(splitr::Args::command().name("splitr"))
    //     // .subcommand(annotate::Args::command().name("annotate"))
    //     // .subcommand(resolve::Args::command().name("resolve"))
    //     // .subcommand(hashshard::Args::command().name("hashshard"))
    //     .get_matches();
    match args.cmd {
        Commands::Hashshard(cmd_args) => {
            hashshard::run(cmd_args)?;
        }
        Commands::Splitr(cmd_args) => {
            splitr::run(cmd_args)?;
        }
        Commands::Annotate(cmd_args) => {
            annotate::run(cmd_args)?;
        }
        Commands::Resolve(cmd_args) => {
            resolve::run(cmd_args)?;
        }
        Commands::Classify(cmd_args) => {
            let splitr_args = splitr::Args::from(cmd_args.clone());
            let chunk_files = find_and_sort_files(&splitr_args.chunk_dir, "sample", ".k2")?;
            if !chunk_files.is_empty() {
                panic!("{} must be empty", &splitr_args.chunk_dir.display());
            }
            splitr::run(splitr_args)?;
            let annotate_args = annotate::Args::from(cmd_args.clone());
            annotate::run(annotate_args)?;
            let resolve_args = resolve::Args::from(cmd_args.clone());
            resolve::run(resolve_args)?;
        }
    }

    // let cmd1_args = splitr::Args::from_arg_matches(sub_matches).expect("parse splitr arg error");
    Ok(())
    // match args() {
    //     Some(("splitr", sub_matches)) => {
    //         let cmd1_args =
    //             splitr::Args::from_arg_matches(sub_matches).expect("parse splitr arg error");
    //         splitr::run(cmd1_args)
    //     }
    //     Some(("annotate", sub_matches)) => {
    //         let cmd1_args =
    //             annotate::Args::from_arg_matches(sub_matches).expect("parse annotate arg error");
    //         annotate::run(cmd1_args)
    //     }
    //     Some(("resolve", sub_matches)) => {
    //         let cmd1_args =
    //             resolve::Args::from_arg_matches(sub_matches).expect("parse resolve arg error");
    //         resolve::run(cmd1_args)
    //     }
    //     Some(("hashshard", sub_matches)) => {
    //         let cmd1_args =
    //             hashshard::Args::from_arg_matches(sub_matches).expect("parse resolve arg error");
    //         hashshard::run(cmd1_args)
    //     }
    //     _ => Ok(()),
    // }
}
