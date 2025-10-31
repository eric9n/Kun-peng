use clap::{Parser, Subcommand};
mod annotate;
mod build_db;
mod chunk_db;
mod direct;
mod estimate_capacity;
mod hashshard;
mod merge_fna;
mod resolve;
mod splitr;
mod add_library;

use kun_peng::args::ClassifyArgs;
use kun_peng::args::{parse_size, Build};
use kun_peng::utils::find_files;
use std::path::PathBuf;
use std::time::Instant;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: jemallocator::Jemalloc = jemallocator::Jemalloc;

#[derive(Parser, Debug, Clone)]
#[clap(author, version, about="Run the complete database build process", long_about = "Run the complete database build process.
This is an all-in-one command that automatically executes all steps for 'merge_fna' (merge downloaded library files) and 'build-db' (estimate, chunk, build hash tables).
If you already have a 'library/' directory and only want to build the hash tables, use the 'build-db' command instead."
)]
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
    #[clap(long, default_value_t = 0.7)]
    load_factor: f64,

    /// library fna temp file max size
    #[arg(long = "max-file-size", value_parser = parse_size, default_value = "2G")]
    pub max_file_size: usize,
}

#[derive(Parser, Debug, Clone)]
#[clap(author, version, about="Run the final database construction steps (estimate, chunk, build)", long_about = "Build-only: Run estimate, chunk, and build steps on an existing library dir.")]
struct BuildDBArgs {
    #[clap(flatten)]
    pub build: Build,

    /// Manually set the precise hash table capacity (number of slots).
    #[arg(
        short = 'c',
        long,
        value_name = "EXACT_SLOT_COUNT",
        long_help = "Manually set the precise hash table capacity (number of slots).
If this value is set, the time-consuming 'estimate_capacity' step will be skipped.

WARNING: This is a critical performance parameter.
- A value that is too low will cause the build to fail or result in a high load factor, which severely degrades classification speed.
- A value that is too high will build successfully but will waste significant memory and disk space.

It is highly recommended to run the 'estimate_capacity' command first to determine a safe and optimal value."
    )]
    pub required_capacity: Option<usize>,

    /// estimate capacity from cache if exists
    #[arg(long, default_value_t = true)]
    cache: bool,

    /// Set maximum qualifying hash code
    #[clap(long, default_value = "4")]
    max_n: usize,

    /// Proportion of the hash table to be populated
    #[clap(long, default_value_t = 0.7)]
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
            num_threads: item.num_threads,
            confidence_threshold: item.confidence_threshold,
            minimum_hit_groups: item.minimum_hit_groups,
            output_dir: item.output_dir,
            report_kmer_data: item.report_kmer_data,
            report_zero_counts: item.report_zero_counts,
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

impl From<BuildDBArgs> for estimate_capacity::Args {
    fn from(item: BuildDBArgs) -> Self {
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

impl From<BuildDBArgs> for chunk_db::Args {
    fn from(item: BuildDBArgs) -> Self {
        Self {
            build: item.build,
            hash_capacity: parse_size("1G").unwrap(),
        }
    }
}


#[derive(Subcommand, Debug)]
enum Commands {
    Estimate(estimate_capacity::Args),
    Build(BuildArgs),
    BuildDB(BuildDBArgs),
    Hashshard(hashshard::Args),
    Splitr(splitr::Args),
    Annotate(annotate::Args),
    Resolve(resolve::Args),
    Classify(ClassifyArgs),
    Direct(direct::Args),
    MergeFna(merge_fna::Args),
    AddLibrary(add_library::Args),
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
        Commands::AddLibrary(cmd_args) => {
            add_library::run(cmd_args)?;
        }
        Commands::Build(cmd_args) => {
            let fna_args = merge_fna::Args::from(cmd_args.clone());
            merge_fna::run(fna_args)?;
            let ec_args = estimate_capacity::Args::from(cmd_args.clone());
            let required_capacity = estimate_capacity::run(ec_args);

            let build_args = chunk_db::Args::from(cmd_args.clone());
            let database = &build_args.build.database.clone();
            chunk_db::run(build_args, required_capacity)?;
            build_db::run(database)?;
        }
        Commands::BuildDB(cmd_args) => {
            println!("Running: BuildDB (Building from existing library)");
            let required_capacity = match cmd_args.required_capacity {
                Some(cap) => {
                    println!("Using user-provided capacity: {}", cap);
                    cap
                }
                None => {
                    println!("Estimating capacity...");
                    let ec_args = estimate_capacity::Args::from(cmd_args.clone());
                    estimate_capacity::run(ec_args)
                }
            };

            let build_args = chunk_db::Args::from(cmd_args.clone());
            let database = &build_args.build.database.clone();
            chunk_db::run(build_args, required_capacity)?;
            build_db::run(database)?;
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
