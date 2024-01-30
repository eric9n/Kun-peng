use clap::Parser;
use kr2r::readcounts::TaxonCounters;
use kr2r::utils::IndexOptions;
use std::fs::File;
use std::io::{self, ErrorKind, Read};
use std::path::Path;
use std::sync::Mutex;
use std::time::Duration;

fn read_index_options<P: AsRef<Path>>(file_path: P) -> io::Result<IndexOptions> {
    let mut file = File::open(file_path)?;
    let mut buffer = vec![0; std::mem::size_of::<IndexOptions>()];
    file.read_exact(&mut buffer)?;

    let idx_opts: IndexOptions = unsafe {
        // 确保这种转换是安全的，这依赖于数据的确切布局和来源
        std::ptr::read(buffer.as_ptr() as *const _)
    };

    Ok(idx_opts)
}

#[derive(Parser, Debug, Clone)]
#[clap(
    version,
    about = "classify",
    long_about = "classify a set of sequences"
)]
struct Args {
    #[clap(short = 'H', long = "index-filename", value_parser, required = true)]
    index_filename: String,

    #[clap(short = 't', long = "taxonomy-filename", value_parser, required = true)]
    taxonomy_filename: String,

    #[clap(short = 'o', long = "options-filename", value_parser, required = true)]
    options_filename: String,

    #[clap(
        short = 'T',
        long = "confidence-threshold",
        value_parser,
        default_value_t = 0.0
    )]
    confidence_threshold: f64,

    #[clap(short = 'q', long = "quick-mode", action)]
    quick_mode: bool,

    #[clap(short = 'p', long = "num-threads", value_parser, default_value_t = 1)]
    num_threads: i32,

    #[clap(short = 'g', long = "minimum-hit-groups", value_parser)]
    minimum_hit_groups: Option<i32>,

    #[clap(short = 'P', long = "paired-end-processing", action)]
    paired_end_processing: bool,

    #[clap(short = 'S', long = "single-file-pairs", action)]
    single_file_pairs: bool,

    #[clap(short = 'm', long = "mpa-style-report", action)]
    mpa_style_report: bool,

    #[clap(short = 'K', long = "report-kmer-data", action)]
    report_kmer_data: bool,

    #[clap(short = 'R', long = "report-filename", value_parser)]
    report_filename: Option<String>,

    #[clap(short = 'z', long = "report-zero-counts", action)]
    report_zero_counts: bool,

    #[clap(short = 'C', long = "classified-output-filename", value_parser)]
    classified_output_filename: Option<String>,

    #[clap(short = 'U', long = "unclassified-output-filename", value_parser)]
    unclassified_output_filename: Option<String>,

    #[clap(short = 'O', long = "kraken-output-filename", value_parser)]
    kraken_output_filename: Option<String>,

    #[clap(short = 'n', long = "print-scientific-name", action)]
    print_scientific_name: bool,

    #[clap(
        short = 'Q',
        long = "minimum-quality-score",
        value_parser,
        default_value_t = 0
    )]
    minimum_quality_score: i32,

    #[clap(short = 'M', long = "use-memory-mapping", action)]
    use_memory_mapping: bool,
}

struct ClassificationStats {
    total_sequences: u64,
    total_bases: u64,
    total_classified: u64,
}

struct OutputData {
    block_id: u64,
    kraken_str: String,
    classified_out1_str: String,
    classified_out2_str: String,
    unclassified_out1_str: String,
    unclassified_out2_str: String,
}

struct OutputStreamData {
    initialized: Mutex<bool>,
    printing_sequences: bool,
    classified_output1: Option<File>,
    classified_output2: Option<File>,
    unclassified_output1: Option<File>,
    unclassified_output2: Option<File>,
    kraken_output: Option<File>,
}

fn report_stats(elapsed: Duration, stats: &ClassificationStats) {
    let seconds = elapsed.as_secs() as f64 + elapsed.subsec_micros() as f64 * 1e-6;
    let total_unclassified = stats.total_sequences - stats.total_classified;

    eprintln!(
        "{} sequences ({:.2} Mbp) processed in {:.3}s ({:.1} Kseq/m, {:.2} Mbp/m).",
        stats.total_sequences,
        stats.total_bases as f64 / 1e6,
        seconds,
        stats.total_sequences as f64 / 1e3 / (seconds / 60.0),
        stats.total_bases as f64 / 1e6 / (seconds / 60.0)
    );

    eprintln!(
        "  {} sequences classified ({:.2}%)",
        stats.total_classified,
        stats.total_classified as f64 * 100.0 / stats.total_sequences as f64
    );

    eprintln!(
        "  {} sequences unclassified ({:.2}%)",
        total_unclassified,
        total_unclassified as f64 * 100.0 / stats.total_sequences as f64
    );
}

fn initialize_outputs(args: Args, outputs: &mut OutputStreamData) -> io::Result<()> {
    let mut initialized = outputs.initialized.lock().unwrap();
    if !*initialized {
        if let Some(filename) = args.classified_output_filename {
            outputs.classified_output1 = Some(open_output_file(
                &filename,
                args.paired_end_processing,
                true,
            )?);
            if args.paired_end_processing {
                outputs.classified_output2 = Some(open_output_file(
                    &filename,
                    args.paired_end_processing,
                    false,
                )?);
            }
            outputs.printing_sequences = true;
        }
        if let Some(filename) = args.unclassified_output_filename {
            outputs.unclassified_output1 = Some(open_output_file(
                &filename,
                args.paired_end_processing,
                true,
            )?);
            if args.paired_end_processing {
                outputs.unclassified_output2 = Some(open_output_file(
                    &filename,
                    args.paired_end_processing,
                    false,
                )?);
            }
        }
        if let Some(filename) = args.kraken_output_filename {
            outputs.kraken_output = Some(File::create(filename)?);
        }

        *initialized = true;
    }
    Ok(())
}

fn open_output_file(filename: &str, is_paired: bool, is_first: bool) -> io::Result<File> {
    if is_paired {
        let fields: Vec<&str> = filename.split('#').collect();
        if fields.len() < 2 {
            return Err(io::Error::new(
                ErrorKind::InvalidInput,
                "Paired filename format missing # character",
            ));
        } else if fields.len() > 2 {
            return Err(io::Error::new(
                ErrorKind::InvalidInput,
                "Paired filename format has >1 # character",
            ));
        }
        let modified_filename = format!(
            "{}_{}{}",
            fields[0],
            if is_first { "1" } else { "2" },
            fields[1]
        );
        File::create(&modified_filename)
    } else {
        File::create(filename)
    }
}

fn main() -> io::Result<()> {
    let args = Args::parse();
    let idx_opts = read_index_options(args.options_filename)?;

    Ok(())
}
