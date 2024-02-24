use clap::Parser;
use kr2r::compact_hash::CompactHashTable;
use kr2r::iclassify::classify_seq;
use kr2r::mmscanner::MinimizerScanner;
// use kr2r::readcounts::TaxonCounters;
use kr2r::pair;
use kr2r::taxonomy::Taxonomy;
use kr2r::IndexOptions;
use seq_io::fastq::{Reader as FqReader, Record, RefRecord};
use seq_io::parallel::read_parallel;
// use std::fs::File;
use std::io::{Error, ErrorKind, Result};
// use std::sync::Mutex;
// use std::time::Duration;

/// Command line arguments for the classify program.
///
/// This structure defines the command line arguments that are accepted by the classify program.
/// It uses the `clap` crate for parsing command line arguments.
#[derive(Parser, Debug, Clone)]
#[clap(
    version,
    about = "classify",
    long_about = "classify a set of sequences"
)]
struct Args {
    /// The file path for the Kraken 2 index.
    #[clap(short = 'H', long = "index-filename", value_parser, required = true)]
    index_filename: String,

    /// The file path for the Kraken 2 taxonomy.
    #[clap(short = 't', long = "taxonomy-filename", value_parser, required = true)]
    taxonomy_filename: String,

    /// The file path for the Kraken 2 options.
    #[clap(short = 'o', long = "options-filename", value_parser, required = true)]
    options_filename: String,

    /// Confidence score threshold, default is 0.0.
    #[clap(
        short = 'T',
        long = "confidence-threshold",
        value_parser,
        default_value_t = 0.0
    )]
    confidence_threshold: f64,

    /// Enable quick mode for faster processing.
    #[clap(short = 'q', long = "quick-mode", action)]
    quick_mode: bool,

    /// The number of threads to use, default is 1.
    #[clap(short = 'p', long = "num-threads", value_parser, default_value_t = 1)]
    num_threads: i32,

    /// The minimum number of hit groups needed for a call.
    #[clap(
        short = 'g',
        long = "minimum-hit-groups",
        value_parser,
        default_value_t = 2
    )]
    minimum_hit_groups: i32,

    /// Enable paired-end processing.
    #[clap(short = 'P', long = "paired-end-processing", action)]
    paired_end_processing: bool,

    /// Process pairs with mates in the same file.
    #[clap(short = 'S', long = "single-file-pairs", action)]
    single_file_pairs: bool,

    /// Use mpa-style report format.
    #[clap(short = 'm', long = "mpa-style-report", action)]
    mpa_style_report: bool,

    /// Report k-mer data in the output.
    #[clap(short = 'K', long = "report-kmer-data", action)]
    report_kmer_data: bool,

    /// File path for outputting the report.
    #[clap(short = 'R', long = "report-filename", value_parser)]
    report_filename: Option<String>,

    /// Report taxa with zero count.
    #[clap(short = 'z', long = "report-zero-counts", action)]
    report_zero_counts: bool,

    /// File path for outputting classified sequences.
    #[clap(short = 'C', long = "classified-output-filename", value_parser)]
    classified_output_filename: Option<String>,

    /// File path for outputting unclassified sequences.
    #[clap(short = 'U', long = "unclassified-output-filename", value_parser)]
    unclassified_output_filename: Option<String>,

    /// File path for outputting normal Kraken output.
    #[clap(short = 'O', long = "kraken-output-filename", value_parser)]
    kraken_output_filename: Option<String>,

    /// Print scientific name instead of taxid in Kraken output.
    #[clap(short = 'n', long = "print-scientific-name", action)]
    print_scientific_name: bool,

    /// Minimum quality score for FASTQ data, default is 0.
    #[clap(
        short = 'Q',
        long = "minimum-quality-score",
        value_parser,
        default_value_t = 0
    )]
    minimum_quality_score: i32,

    /// Use memory mapping to access hash and taxonomy data.
    #[clap(short = 'M', long = "use-memory-mapping", action)]
    use_memory_mapping: bool,

    /// Input files for processing.
    ///
    /// A list of input file paths (FASTA/FASTQ) to be processed by the classify program.
    // #[clap(short = 'F', long = "files")]
    input_files: Vec<String>,
}

// struct ClassificationStats {
//     total_sequences: u64,
//     total_bases: u64,
//     total_classified: u64,
// }

// struct OutputData {
//     block_id: u64,
//     kraken_str: String,
//     classified_out1_str: String,
//     classified_out2_str: String,
//     unclassified_out1_str: String,
//     unclassified_out2_str: String,
// }

// struct OutputStreamData {
//     initialized: Mutex<bool>,
//     printing_sequences: bool,
//     classified_output1: Option<File>,
//     classified_output2: Option<File>,
//     unclassified_output1: Option<File>,
//     unclassified_output2: Option<File>,
//     kraken_output: Option<File>,
// }

// /// 感觉作用不大。可以删掉
// fn report_stats(elapsed: Duration, stats: &ClassificationStats) {
//     let seconds = elapsed.as_secs() as f64 + elapsed.subsec_micros() as f64 * 1e-6;
//     let total_unclassified = stats.total_sequences - stats.total_classified;

//     eprintln!(
//         "{} sequences ({:.2} Mbp) processed in {:.3}s ({:.1} Kseq/m, {:.2} Mbp/m).",
//         stats.total_sequences,
//         stats.total_bases as f64 / 1e6,
//         seconds,
//         stats.total_sequences as f64 / 1e3 / (seconds / 60.0),
//         stats.total_bases as f64 / 1e6 / (seconds / 60.0)
//     );

//     eprintln!(
//         "  {} sequences classified ({:.2}%)",
//         stats.total_classified,
//         stats.total_classified as f64 * 100.0 / stats.total_sequences as f64
//     );

//     eprintln!(
//         "  {} sequences unclassified ({:.2}%)",
//         total_unclassified,
//         total_unclassified as f64 * 100.0 / stats.total_sequences as f64
//     );
// }

// fn initialize_outputs(args: Args, outputs: &mut OutputStreamData) -> Result<()> {
//     let mut initialized = outputs.initialized.lock().unwrap();
//     if !*initialized {
//         if let Some(filename) = args.classified_output_filename {
//             outputs.classified_output1 = Some(open_output_file(
//                 &filename,
//                 args.paired_end_processing,
//                 true,
//             )?);
//             if args.paired_end_processing {
//                 outputs.classified_output2 = Some(open_output_file(
//                     &filename,
//                     args.paired_end_processing,
//                     false,
//                 )?);
//             }
//             outputs.printing_sequences = true;
//         }
//         if let Some(filename) = args.unclassified_output_filename {
//             outputs.unclassified_output1 = Some(open_output_file(
//                 &filename,
//                 args.paired_end_processing,
//                 true,
//             )?);
//             if args.paired_end_processing {
//                 outputs.unclassified_output2 = Some(open_output_file(
//                     &filename,
//                     args.paired_end_processing,
//                     false,
//                 )?);
//             }
//         }
//         if let Some(filename) = args.kraken_output_filename {
//             outputs.kraken_output = Some(File::create(filename)?);
//         }

//         *initialized = true;
//     }
//     Ok(())
// }

// fn open_output_file(filename: &str, is_paired: bool, is_first: bool) -> Result<File> {
//     if is_paired {
//         let fields: Vec<&str> = filename.split('#').collect();
//         if fields.len() < 2 {
//             return Err(Error::new(
//                 ErrorKind::InvalidInput,
//                 "Paired filename format missing # character",
//             ));
//         } else if fields.len() > 2 {
//             return Err(Error::new(
//                 ErrorKind::InvalidInput,
//                 "Paired filename format has >1 # character",
//             ));
//         }
//         let modified_filename = format!(
//             "{}_{}{}",
//             fields[0],
//             if is_first { "1" } else { "2" },
//             fields[1]
//         );
//         File::create(&modified_filename)
//     } else {
//         File::create(filename)
//     }
// }

fn check_feature(dna_db: bool) -> Result<()> {
    #[cfg(feature = "dna")]
    if !dna_db {
        return Err(Error::new(
            ErrorKind::Other,
            "Feature 'dna' is enabled but 'dna_db' is false",
        ));
    }

    #[cfg(feature = "protein")]
    if dna_db {
        return Err(Error::new(
            ErrorKind::Other,
            "Feature 'protein' is enabled but 'dna_db' is true",
        ));
    }

    Ok(())
}

fn get_record_id(ref_record: &RefRecord) -> String {
    std::str::from_utf8(ref_record.head().split(|b| *b == b' ').next().unwrap())
        .unwrap_or_default()
        .into()
}

/// 处理fastq文件
fn process_files(args: Args, idx_opts: IndexOptions, cht: &CompactHashTable, taxonomy: &Taxonomy) {
    let queue_len = if args.num_threads > 2 {
        args.num_threads as usize - 2
    } else {
        1
    };
    let meros = idx_opts.as_meros();

    if args.paired_end_processing && !args.single_file_pairs {
        // 处理成对的文件
        for file_pair in args.input_files.chunks(2) {
            let file1 = &file_pair[0];
            let file2 = &file_pair[1];
            // 对 file1 和 file2 执行分类处理
            let pair_reader = pair::PairReader::from_path(file1, file2).unwrap();
            read_parallel(
                pair_reader,
                args.num_threads as u32,
                queue_len,
                |record_set| {
                    let mut scanner = MinimizerScanner::new(idx_opts.as_meros());

                    for records in record_set.into_iter() {
                        let dna_id = get_record_id(&records.0);
                        let record_list = vec![records.0, records.1];
                        classify_seq(
                            &taxonomy,
                            &cht,
                            &mut scanner,
                            &record_list,
                            args.minimum_quality_score,
                            meros,
                            args.confidence_threshold,
                            args.minimum_hit_groups,
                            dna_id.into(),
                        );
                    }
                },
                |record_sets| {
                    while let Some(Ok((_, _))) = record_sets.next() {
                        // counter.fetch_add(count, Ordering::SeqCst);
                    }
                },
            )
        }
    } else {
        for file in args.input_files {
            // 对 file 执行分类处理
            let reader = FqReader::from_path(file).unwrap();
            read_parallel(
                reader,
                args.num_threads as u32,
                queue_len,
                |record_set| {
                    let mut scanner = MinimizerScanner::new(idx_opts.as_meros());
                    for record1 in record_set.into_iter() {
                        let dna_id = get_record_id(&record1);
                        let record_list = vec![record1];
                        classify_seq(
                            &taxonomy,
                            &cht,
                            &mut scanner,
                            &record_list,
                            args.minimum_quality_score,
                            meros,
                            args.confidence_threshold,
                            args.minimum_hit_groups,
                            dna_id.into(),
                        );
                    }
                },
                |record_sets| {
                    while let Some(Ok((_, _))) = record_sets.next() {}
                },
            )
        }
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    let idx_opts = IndexOptions::read_index_options(args.options_filename.clone())?;
    check_feature(idx_opts.dna_db)?;
    let taxo = Taxonomy::from_file(&args.taxonomy_filename)?;
    let cht = CompactHashTable::from(args.index_filename.clone())?;

    if args.paired_end_processing && !args.single_file_pairs && args.input_files.len() % 2 != 0 {
        // 验证文件列表是否为偶数个
        return Err(Error::new(
            ErrorKind::InvalidInput,
            "Paired-end processing requires an even number of input files.",
        ));
    }
    process_files(args, idx_opts, &cht, &taxo);
    Ok(())
}
