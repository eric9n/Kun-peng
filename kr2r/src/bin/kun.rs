use clap::{Parser, Subcommand};
mod annotate;
mod build_k2_db;
mod chunk_db;
mod direct;
mod estimate_capacity;
mod hashshard;
mod merge_fna;
mod resolve;
// mod seqid2taxid;
mod splitr;

use kr2r::args::ClassifyArgs;
use kr2r::args::{parse_size, Build};
use kr2r::utils::find_files;
// use std::io::Result;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser, Debug, Clone)]
#[clap(author, version, about="build database", long_about = None)]
struct BuildArgs {
    // /// database hash chunk directory and other files
    // #[clap(long)]
    // pub k2d_dir: Option<PathBuf>,
    /// Directory to store downloaded files
    #[arg(short, long, required = true)]
    pub download_dir: PathBuf,

    // chunk_dir: PathBuf,
    #[clap(flatten)]
    pub build: Build,

    // #[arg(short = 'm')]
    // pub id_to_taxon_map_filename: Option<PathBuf>,
    // #[clap(flatten)]
    // taxo: Taxo,
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

    /// library fna temp file max size
    #[arg(long = "max-file-size", value_parser = parse_size, default_value = "2G")]
    pub max_file_size: usize,
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
            database: item.database,
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
            database: item.database,
            chunk_dir: item.chunk_dir,
            batch_size: item.batch_size,
            buffer_size: item.buffer_size,
            num_threads: item.num_threads,
        }
    }
}

impl From<ClassifyArgs> for resolve::Args {
    fn from(item: ClassifyArgs) -> Self {
        Self {
            database: item.database,
            chunk_dir: item.chunk_dir,
            confidence_threshold: item.confidence_threshold,
            minimum_hit_groups: item.minimum_hit_groups,
            kraken_output_dir: item.kraken_output_dir,
            report_kmer_data: item.report_kmer_data,
            report_zero_counts: item.report_zero_counts,
            // full_output: item.full_output,
            num_threads: item.num_threads,
        }
    }
}

impl From<BuildArgs> for estimate_capacity::Args {
    fn from(item: BuildArgs) -> Self {
        Self {
            database: item.build.database,
            klmt: item.build.klmt,
            cache: item.cache,
            n: item.max_n,
            load_factor: item.load_factor,
            threads: item.build.threads,
        }
    }
}

impl From<BuildArgs> for chunk_db::Args {
    fn from(item: BuildArgs) -> Self {
        Self {
            build: item.build,
            hash_capacity: parse_size("1G").unwrap(),
        }
    }
}

impl From<BuildArgs> for merge_fna::Args {
    fn from(item: BuildArgs) -> Self {
        Self {
            download_dir: item.download_dir,
            database: item.build.database,
            max_file_size: item.max_file_size,
        }
    }
}

#[derive(Subcommand, Debug)]
enum Commands {
    Estimate(estimate_capacity::Args),
    // Seqid2taxid(seqid2taxid::Args),
    Build(BuildArgs),
    Hashshard(hashshard::Args),
    Splitr(splitr::Args),
    Annotate(annotate::Args),
    Resolve(resolve::Args),
    Classify(ClassifyArgs),
    Direct(direct::Args),
    MergeFna(merge_fna::Args),
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    match args.cmd {
        Commands::MergeFna(cmd_args) => {
            merge_fna::run(cmd_args)?;
        }
        Commands::Estimate(cmd_args) => {
            estimate_capacity::run(cmd_args);
        }
        Commands::Build(cmd_args) => {
            let fna_args = merge_fna::Args::from(cmd_args.clone());
            merge_fna::run(fna_args)?;
            let ec_args = estimate_capacity::Args::from(cmd_args.clone());
            let required_capacity = estimate_capacity::run(ec_args);

            let build_args = chunk_db::Args::from(cmd_args.clone());
            let database = &build_args.build.database.clone();
            chunk_db::run(build_args, required_capacity)?;
            build_k2_db::run(database)?;
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
            let chunk_files = find_files(&splitr_args.chunk_dir, "sample", ".k2");
            let sample_files = find_files(&splitr_args.chunk_dir, "sample_id", ".map");
            let bin_files = find_files(&splitr_args.chunk_dir, "sample", ".bin");
            if !chunk_files.is_empty() || !sample_files.is_empty() || !bin_files.is_empty() {
                return Err(Box::new(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!(
                        "The directory '{}' must not contain files with extensions '.k2', '.map', or '.bin' for 'sample' and 'sample_id'",
                        &splitr_args.chunk_dir.display()
                    ),
                )));
            }
            splitr::run(splitr_args)?;
            let annotate_args = annotate::Args::from(cmd_args.clone());
            annotate::run(annotate_args)?;
            let resolve_args = resolve::Args::from(cmd_args.clone());
            resolve::run(resolve_args)?;

            let duration = start.elapsed();
            println!("Classify took: {:?}", duration);
        }
        Commands::Direct(cmd_args) => {
            direct::run(cmd_args)?;
        }
    }

    Ok(())
}
