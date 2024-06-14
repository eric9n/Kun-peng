use clap::Parser;
use kr2r::classify::{add_hitlist_string, count_values, resolve_tree, trim_pair_info};
use kr2r::compact_hash::{CHTable, Compact, HashConfig, Row};
use kr2r::mmscanner::MinimizerScanner;
use kr2r::readcounts::{TaxonCounters, TaxonCountersDash};
use kr2r::report::report_kraken_style;
use kr2r::seq::{self, open_fasta_reader, SeqX};
use kr2r::taxonomy::Taxonomy;
use kr2r::utils::{
    create_sample_file, detect_file_format, find_and_sort_files, get_lastest_file_index, FileFormat,
};
use kr2r::{IndexOptions, Meros};
use seq_io::fasta::Record;
use seq_io::fastq::Record as FqRecord;
use seq_io::parallel::read_parallel;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::io::{Error, ErrorKind, Result};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use seqkmer::fastq::FastqReader;
use seqkmer::parallel::read_parallel as s_parallel;
use seqkmer::reader::Reader;
use seqkmer::Meros as SMeros;

#[derive(Parser, Debug, Clone)]
#[clap(
    version,
    about = "Directly load all hash tables for classification annotation",
    long_about = "Directly load all hash tables for classification annotation"
)]
pub struct Args {
    /// database hash chunk directory and other files
    #[clap(long)]
    pub k2d_dir: PathBuf,

    // /// The file path for the Kraken 2 options.
    // #[clap(short = 'o', long = "options-filename", value_parser, required = true)]
    // options_filename: String,
    /// Enable paired-end processing.
    #[clap(short = 'P', long = "paired-end-processing", action)]
    pub paired_end_processing: bool,

    /// Process pairs with mates in the same file.
    #[clap(short = 'S', long = "single-file-pairs", action)]
    pub single_file_pairs: bool,

    /// Minimum quality score for FASTQ data, default is 0.
    #[clap(
        short = 'Q',
        long = "minimum-quality-score",
        value_parser,
        default_value_t = 0
    )]
    pub minimum_quality_score: i32,

    /// Confidence score threshold, default is 0.0.
    #[clap(
        short = 'T',
        long = "confidence-threshold",
        value_parser,
        default_value_t = 0.0
    )]
    pub confidence_threshold: f64,

    /// In comb. w/ -R, provide minimizer information in report
    #[clap(short = 'K', long, value_parser, default_value_t = false)]
    pub report_kmer_data: bool,

    /// In comb. w/ -R, report taxa w/ 0 count
    #[clap(short = 'z', long, value_parser, default_value_t = false)]
    pub report_zero_counts: bool,

    /// The minimum number of hit groups needed for a call.
    #[clap(
        short = 'g',
        long = "minimum-hit-groups",
        value_parser,
        default_value_t = 2
    )]
    pub minimum_hit_groups: usize,

    /// The number of threads to use, default is 10.
    #[clap(short = 'p', long = "num-threads", value_parser, default_value_t = 16)]
    pub num_threads: i32,

    /// File path for outputting normal Kraken output.
    #[clap(long = "output-dir", value_parser)]
    pub kraken_output_dir: Option<PathBuf>,

    /// A list of input file paths (FASTA/FASTQ) to be processed by the classify program.
    // #[clap(short = 'F', long = "files")]
    pub input_files: Vec<String>,
}

fn process_seq(
    miner: MinimizerScanner,
    hash_config: &HashConfig,
    chtable: &CHTable,
    offset: u32,
) -> (u32, Vec<Row>) {
    let chunk_size = hash_config.hash_capacity;
    let value_bits = hash_config.value_bits;

    let mut rows = Vec::new();
    let mut kmer_count = 0;
    for (sort, hash_key) in miner.into_iter().enumerate() {
        let idx = hash_config.index(hash_key);
        let partition_index = idx / chunk_size;
        let index = idx % chunk_size;
        let taxid = chtable.get_from_page(index, hash_key, partition_index + 1);
        if taxid > 0 {
            let compacted_key = hash_key.left(value_bits) as u32;
            let high = u32::combined(compacted_key, taxid, value_bits);
            let row = Row::new(high, 0, sort as u32 + 1 + offset);
            rows.push(row);
        }
        kmer_count += 1;
    }
    (kmer_count, rows)
}

fn process_record(
    dna_id: String,
    seq1: Vec<u8>,
    seq2: Option<Vec<u8>>,
    args: &Args,
    taxonomy: &Taxonomy,
    meros: Meros,
    chtable: &CHTable,
    hash_config: &HashConfig,
    cur_taxon_counts: &TaxonCountersDash,
    classify_counter: &AtomicUsize,
) -> String {
    let value_mask = hash_config.value_mask;
    let mut seq_len_str = String::new();
    let seq1_len = seq1.len();
    seq_len_str.push_str(&seq1_len.to_string());

    let scan1 = MinimizerScanner::new(&seq1, meros);
    let (kmer_count1, mut rows) = process_seq(scan1, &hash_config, chtable, 0);
    let kmer_count2 = if let Some(seq) = seq2 {
        let scan2 = MinimizerScanner::new(&seq, meros);
        let (kmer_count2, rows2) = process_seq(scan2, &hash_config, chtable, kmer_count1);
        rows.extend_from_slice(&rows2);
        seq_len_str.push_str(format!("|{}", seq.len()).as_str());
        Some(kmer_count2)
    } else {
        None
    };
    let total_kmers: usize = (kmer_count1 + kmer_count2.unwrap_or(0)) as usize;
    let (counts, cur_counts, hit_groups) = count_values(&rows, value_mask, kmer_count1);
    let hit_string = add_hitlist_string(&rows, value_mask, kmer_count1, kmer_count2, taxonomy);
    let mut call = resolve_tree(&counts, taxonomy, total_kmers, args.confidence_threshold);
    if call > 0 && hit_groups < args.minimum_hit_groups {
        call = 0;
    };

    cur_counts.iter().for_each(|entry| {
        cur_taxon_counts
            .entry(*entry.key())
            .or_default()
            .merge(entry.value())
            .unwrap();
    });

    let ext_call = taxonomy.nodes[call as usize].external_id;
    let clasify = if call > 0 {
        classify_counter.fetch_add(1, Ordering::SeqCst);
        cur_taxon_counts
            .entry(call as u64)
            .or_default()
            .increment_read_count();

        "C"
    } else {
        "U"
    };
    // 使用锁来同步写入
    let output_line = format!(
        "{}\t{}\t{}\t{}\t{}\n",
        clasify, dna_id, ext_call, seq_len_str, hit_string
    );
    output_line
}

fn process_seq1(
    miner: Vec<u64>,
    hash_config: &HashConfig,
    chtable: &CHTable,
    offset: u32,
) -> (u32, Vec<Row>) {
    let chunk_size = hash_config.hash_capacity;
    let value_bits = hash_config.value_bits;

    let mut rows = Vec::new();
    let mut kmer_count = 0;
    for (sort, hash_key) in miner.into_iter().enumerate() {
        let idx = hash_config.index(hash_key);
        let partition_index = idx / chunk_size;
        let index = idx % chunk_size;
        let taxid = chtable.get_from_page(index, hash_key, partition_index + 1);
        if taxid > 0 {
            let compacted_key = hash_key.left(value_bits) as u32;
            let high = u32::combined(compacted_key, taxid, value_bits);
            let row = Row::new(high, 0, sort as u32 + 1 + offset);
            rows.push(row);
        }
        kmer_count += 1;
    }
    (kmer_count, rows)
}

fn process_record1(
    dna_id: String,
    seq1: Vec<u64>,
    seq2: Option<Vec<u8>>,
    args: &Args,
    taxonomy: &Taxonomy,
    meros: Meros,
    chtable: &CHTable,
    hash_config: &HashConfig,
    cur_taxon_counts: &TaxonCountersDash,
    classify_counter: &AtomicUsize,
) -> String {
    let value_mask = hash_config.value_mask;
    let mut seq_len_str = String::new();
    let seq1_len = seq1.len();
    seq_len_str.push_str(&seq1_len.to_string());

    let (kmer_count1, mut rows) = process_seq1(seq1, &hash_config, chtable, 0);
    let kmer_count2 = if let Some(seq) = seq2 {
        let scan2 = MinimizerScanner::new(&seq, meros);
        let (kmer_count2, rows2) = process_seq(scan2, &hash_config, chtable, kmer_count1);
        rows.extend_from_slice(&rows2);
        seq_len_str.push_str(format!("|{}", seq.len()).as_str());
        Some(kmer_count2)
    } else {
        None
    };
    let total_kmers: usize = (kmer_count1 + kmer_count2.unwrap_or(0)) as usize;
    let (counts, cur_counts, hit_groups) = count_values(&rows, value_mask, kmer_count1);
    let hit_string = add_hitlist_string(&rows, value_mask, kmer_count1, kmer_count2, taxonomy);
    let mut call = resolve_tree(&counts, taxonomy, total_kmers, args.confidence_threshold);
    if call > 0 && hit_groups < args.minimum_hit_groups {
        call = 0;
    };

    cur_counts.iter().for_each(|entry| {
        cur_taxon_counts
            .entry(*entry.key())
            .or_default()
            .merge(entry.value())
            .unwrap();
    });

    let ext_call = taxonomy.nodes[call as usize].external_id;
    let clasify = if call > 0 {
        classify_counter.fetch_add(1, Ordering::SeqCst);
        cur_taxon_counts
            .entry(call as u64)
            .or_default()
            .increment_read_count();

        "C"
    } else {
        "U"
    };
    // 使用锁来同步写入
    let output_line = format!(
        "{}\t{}\t{}\t{}\t{}\n",
        clasify, dna_id, ext_call, seq_len_str, hit_string
    );
    output_line
}

fn process_fasta_file(
    args: &Args,
    meros: Meros,
    hash_config: HashConfig,
    file_index: usize,
    files: &[String],
    chtable: &CHTable,
    taxonomy: &Taxonomy,
    total_taxon_counts: &mut TaxonCounters,
) -> io::Result<(usize, usize)> {
    let score = args.minimum_quality_score;
    let mut files_iter = files.iter();
    let file1 = files_iter.next().cloned().unwrap();

    let mut writer: Box<dyn Write + Send> = match &args.kraken_output_dir {
        Some(ref file_path) => {
            let filename = file_path.join(format!("output_{}.txt", file_index));
            let file = File::create(filename)?;
            Box::new(BufWriter::new(file)) as Box<dyn Write + Send>
        }
        None => Box::new(io::stdout()) as Box<dyn Write + Send>,
    };

    let cur_taxon_counts = TaxonCountersDash::new();
    let sequence_count = AtomicUsize::new(0);
    let classify_counter = AtomicUsize::new(0);

    let reader = open_fasta_reader(&file1).expect("Unable to create fasta reader from path");
    read_parallel(
        reader,
        args.num_threads as u32,
        args.num_threads as usize,
        |record_set| {
            let mut buffer = String::new();

            for records in record_set.into_iter() {
                let dna_id = trim_pair_info(records.id().unwrap_or_default());
                sequence_count.fetch_add(1, Ordering::SeqCst);

                let seq1: Vec<u8> = records.seq_x(score);
                let seq2 = None;
                let output_line = process_record(
                    dna_id,
                    seq1,
                    seq2,
                    args,
                    taxonomy,
                    meros,
                    chtable,
                    &hash_config,
                    &cur_taxon_counts,
                    &classify_counter,
                );

                buffer.push_str(&output_line);
            }
            buffer
        },
        |record_sets| {
            while let Some(Ok((_, buffer))) = record_sets.next() {
                writer
                    .write_all(buffer.as_bytes())
                    .expect("write data error");
            }
        },
    );

    let mut sample_taxon_counts: HashMap<
        u64,
        kr2r::readcounts::ReadCounts<hyperloglogplus::HyperLogLogPlus<u64, kr2r::KBuildHasher>>,
    > = HashMap::new();
    cur_taxon_counts.iter().for_each(|entry| {
        total_taxon_counts
            .entry(*entry.key())
            .or_default()
            .merge(&entry.value())
            .unwrap();
        sample_taxon_counts
            .entry(*entry.key())
            .or_default()
            .merge(&entry.value())
            .unwrap();
    });

    let thread_sequences = sequence_count.load(Ordering::SeqCst);
    let thread_classified = classify_counter.load(Ordering::SeqCst);
    if let Some(output) = &args.kraken_output_dir {
        let filename = output.join(format!("output_{}.kreport2", file_index));
        report_kraken_style(
            filename,
            args.report_zero_counts,
            args.report_kmer_data,
            &taxonomy,
            &sample_taxon_counts,
            thread_sequences as u64,
            (thread_sequences - thread_classified) as u64,
        )?;
    }

    Ok((thread_sequences, thread_sequences - thread_classified))
}

/// fastq
fn process_fastq_file(
    args: &Args,
    meros: Meros,
    hash_config: HashConfig,
    file_index: usize,
    files: &[String],
    chtable: &CHTable,
    taxonomy: &Taxonomy,
    total_taxon_counts: &mut TaxonCounters,
) -> io::Result<(usize, usize)> {
    let score = args.minimum_quality_score;
    let mut files_iter = files.iter();
    let file1 = files_iter.next().cloned().unwrap();
    let file2 = files_iter.next().cloned();

    let mut writer: Box<dyn Write + Send> = match &args.kraken_output_dir {
        Some(ref file_path) => {
            let filename = file_path.join(format!("output_{}.txt", file_index));
            let file = File::create(filename)?;
            Box::new(BufWriter::new(file)) as Box<dyn Write + Send>
        }
        None => Box::new(io::stdout()) as Box<dyn Write + Send>,
    };

    let cur_taxon_counts = TaxonCountersDash::new();

    let sequence_count = AtomicUsize::new(0);
    let classify_counter = AtomicUsize::new(0);

    let mut reader1 = FastqReader::from_path(&file1, 1, 0)?;
    let _ = s_parallel(
        &mut reader1,
        13,
        15,
        None,
        SMeros::default(),
        |seq1, seq| {
            let dna_id = trim_pair_info(&seq.id);
            sequence_count.fetch_add(1, Ordering::SeqCst);

            let seq2 = None;
            let output_line = process_record1(
                dna_id,
                seq1,
                seq2,
                args,
                taxonomy,
                meros,
                chtable,
                &hash_config,
                &cur_taxon_counts,
                &classify_counter,
            );
            None
        },
    );
    // let reader = seq::PairFastqReader::from_path(&file1, file2.as_ref())
    //     .expect("Unable to create pair reader from paths");
    // read_parallel(
    //     reader,
    //     args.num_threads as u32,
    //     args.num_threads as usize,
    //     |record_set| {
    //         let mut buffer = String::new();

    //         for records in record_set.into_iter() {
    //             let dna_id = trim_pair_info(records.0.id().unwrap_or_default());
    //             sequence_count.fetch_add(1, Ordering::SeqCst);
    //             let seq1: Vec<u8> = records.0.seq_x(score);
    //             let seq2 = records.1.map(|seq| seq.seq_x(score));
    //             let output_line = process_record(
    //                 dna_id,
    //                 seq1,
    //                 seq2,
    //                 args,
    //                 taxonomy,
    //                 meros,
    //                 chtable,
    //                 &hash_config,
    //                 &cur_taxon_counts,
    //                 &classify_counter,
    //             );

    //             buffer.push_str(&output_line);
    //         }
    //         buffer
    //     },
    //     |record_sets| {
    //         while let Some(Ok((_, buffer))) = record_sets.next() {
    //             writer
    //                 .write_all(buffer.as_bytes())
    //                 .expect("write data error");
    //         }
    //     },
    // );

    let mut sample_taxon_counts: HashMap<
        u64,
        kr2r::readcounts::ReadCounts<hyperloglogplus::HyperLogLogPlus<u64, kr2r::KBuildHasher>>,
    > = HashMap::new();
    cur_taxon_counts.iter().for_each(|entry| {
        total_taxon_counts
            .entry(*entry.key())
            .or_default()
            .merge(&entry.value())
            .unwrap();
        sample_taxon_counts
            .entry(*entry.key())
            .or_default()
            .merge(&entry.value())
            .unwrap();
    });

    let thread_sequences = sequence_count.load(Ordering::SeqCst);
    let thread_classified = classify_counter.load(Ordering::SeqCst);
    if let Some(output) = &args.kraken_output_dir {
        let filename = output.join(format!("output_{}.kreport2", file_index));
        report_kraken_style(
            filename,
            args.report_zero_counts,
            args.report_kmer_data,
            &taxonomy,
            &sample_taxon_counts,
            thread_sequences as u64,
            (thread_sequences - thread_classified) as u64,
        )?;
    }

    Ok((thread_sequences, thread_sequences - thread_classified))
}

fn process_files(
    args: Args,
    meros: Meros,
    hash_config: HashConfig,
    chtable: &CHTable,
    taxonomy: &Taxonomy,
) -> Result<()> {
    let (mut file_index, mut file_writer) = if let Some(out_dir) = &args.kraken_output_dir {
        let file_path = out_dir.join("sample_file.map");
        let file_writer = create_sample_file(&file_path);
        let file_index = get_lastest_file_index(&file_path)?;
        (file_index, Box::new(file_writer) as Box<dyn Write + Send>)
    } else {
        (0, Box::new(io::stdout()) as Box<dyn Write + Send>)
    };

    let mut process_funcs = |files: Vec<&[String]>| -> Result<()> {
        let file_bits = (((files.len() + file_index) as f64).log2().ceil() as usize).max(1);
        if file_bits > hash_config.value_bits {
            panic!("The number of files is too large to process.");
        }

        let mut total_taxon_counts = TaxonCounters::new();
        let mut total_seqs: usize = 0;
        let mut total_unclassified: usize = 0;
        for file_pair in files {
            file_index += 1;

            writeln!(file_writer, "{}\t{}", file_index, file_pair.join(","))?;
            file_writer.flush().unwrap();

            match detect_file_format(&file_pair[0])? {
                FileFormat::Fastq => {
                    let (thread_sequences, thread_unclassified) = process_fastq_file(
                        &args,
                        meros,
                        hash_config,
                        file_index,
                        file_pair,
                        chtable,
                        taxonomy,
                        &mut total_taxon_counts,
                    )?;
                    total_seqs += thread_sequences;
                    total_unclassified += thread_unclassified;
                }
                FileFormat::Fasta => {
                    let (thread_sequences, thread_unclassified) = process_fasta_file(
                        &args,
                        meros,
                        hash_config,
                        file_index,
                        file_pair,
                        chtable,
                        taxonomy,
                        &mut total_taxon_counts,
                    )?;
                    total_seqs += thread_sequences;
                    total_unclassified += thread_unclassified;
                }
            }
        }
        if let Some(output) = &args.kraken_output_dir {
            let filename = output.join("output.kreport2");
            report_kraken_style(
                filename,
                args.report_zero_counts,
                args.report_kmer_data,
                &taxonomy,
                &total_taxon_counts,
                total_seqs as u64,
                total_unclassified as u64,
            )?;
        }

        Ok(())
    };

    if args.paired_end_processing && !args.single_file_pairs {
        // 处理成对的文件
        let files = args.input_files.chunks(2).collect();
        process_funcs(files)?;
    } else {
        let files = args.input_files.chunks(1).collect();
        process_funcs(files)?;
    }

    Ok(())
}

pub fn run(args: Args) -> Result<()> {
    let options_filename = &args.k2d_dir.join("opts.k2d");
    let idx_opts = IndexOptions::read_index_options(options_filename)?;

    if args.paired_end_processing && !args.single_file_pairs && args.input_files.len() % 2 != 0 {
        // 验证文件列表是否为偶数个
        return Err(Error::new(
            ErrorKind::InvalidInput,
            "Paired-end processing requires an even number of input files.",
        ));
    }

    let taxonomy_filename = args.k2d_dir.join("taxo.k2d");
    let taxo = Taxonomy::from_file(taxonomy_filename)?;

    let hash_config = HashConfig::from_hash_header(&args.k2d_dir.join("hash_config.k2d"))?;

    println!("{:?}", hash_config);
    if hash_config.hash_capacity == 0 {
        panic!("`hash_capacity` can't be zero!");
    }
    println!("start...");
    let start = Instant::now();
    let meros = idx_opts.as_meros();
    let hash_files = find_and_sort_files(&args.k2d_dir, "hash", ".k2d")?;
    let chtable = CHTable::from_hash_files(hash_config, hash_files)?;

    process_files(args, meros, hash_config, &chtable, &taxo)?;
    let duration = start.elapsed();
    println!("classify took: {:?}", duration);
    Ok(())
}

#[allow(dead_code)]
fn main() {
    let args = Args::parse();
    if let Err(e) = run(args) {
        eprintln!("Application error: {}", e);
    }
}
