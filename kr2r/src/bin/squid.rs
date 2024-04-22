use clap::{Parser, Subcommand};
mod annotate;
mod build_k2_db;
mod estimate_capacity;
mod hashshard;
mod resolve;
mod seqid2taxid;
mod splitr;

use kr2r::args::ClassifyArgs;
use kr2r::args::{Build, Taxo, ONEGB, U32MAXPLUS};
use kr2r::utils::find_and_sort_files;
use std::io::Result;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser, Debug, Clone)]
#[clap(author, version, about="build database", long_about = None)]
struct BuildArgs {
    #[clap(flatten)]
    pub build: Build,

    #[clap(flatten)]
    taxo: Taxo,

    // /// Name of Kraken 2 database
    // #[arg(short, long = "db")]
    // database: PathBuf,
    /// chunk directory
    #[clap(long)]
    chunk_dir: PathBuf,

    /// chunk size 1-4(GB) [1073741824-4294967295]
    #[clap(long, value_parser = clap::value_parser!(u64).range(ONEGB..U32MAXPLUS), default_value_t = ONEGB)]
    chunk_size: u64,

    /// estimate capacity from cache if exists
    #[arg(long, default_value_t = true)]
    cache: bool,

    /// Set maximum qualifying hash code
    #[clap(long, default_value = "4")]
    max_n: usize,

    /// Proportion of the hash table to be populated
    /// (build task only; def: 0.7, must be
    ///    between 0 and 1).
    #[clap(long, long, default_value_t = 0.7)]
    load_factor: f64,
}

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    #[clap(subcommand)]
    cmd: Commands,
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
            report_kmer_data: item.report_kmer_data,
            report_zero_counts: item.report_zero_counts,
            full_output: item.full_output,
        }
    }
}

impl From<BuildArgs> for estimate_capacity::Args {
    fn from(item: BuildArgs) -> Self {
        Self {
            source: item.build.source,
            klmt: item.build.klmt,
            cache: item.cache,
            n: item.max_n,
            load_factor: item.load_factor,
            threads: item.build.threads,
        }
    }
}

impl From<BuildArgs> for build_k2_db::Args {
    fn from(item: BuildArgs) -> Self {
        Self {
            build: item.build,
            taxo: item.taxo,
            chunk_dir: item.chunk_dir,
            chunk_size: item.chunk_size,
        }
    }
}

impl From<BuildArgs> for seqid2taxid::Args {
    fn from(item: BuildArgs) -> Self {
        Self {
            source: item.build.source,
            id_to_taxon_map_filename: item.taxo.id_to_taxon_map_filename,
        }
    }
}

#[derive(Subcommand, Debug)]
enum Commands {
    Estimate(estimate_capacity::Args),
    Seqid2taxid(seqid2taxid::Args),
    Build(BuildArgs),
    Hashshard(hashshard::Args),
    Splitr(splitr::Args),
    Annotate(annotate::Args),
    Resolve(resolve::Args),
    Classify(ClassifyArgs),
}

fn main() -> Result<()> {
    let args = Args::parse();

    match args.cmd {
        Commands::Estimate(cmd_args) => {
            estimate_capacity::run(cmd_args);
        }
        Commands::Seqid2taxid(cmd_args) => {
            seqid2taxid::run(cmd_args)?;
        }
        Commands::Build(cmd_args) => {
            let seq_args = seqid2taxid::Args::from(cmd_args.clone());
            seqid2taxid::run(seq_args)?;
            let ec_args = estimate_capacity::Args::from(cmd_args.clone());
            let required_capacity = estimate_capacity::run(ec_args);

            let build_args = build_k2_db::Args::from(cmd_args.clone());
            if let Err(e) = build_k2_db::run(build_args, required_capacity) {
                panic!("Application error: {}", e);
            }
        }
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
            let start = Instant::now();

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

            let duration = start.elapsed();
            println!("Classify took: {:?}", duration);
        }
    }

    Ok(())
}
