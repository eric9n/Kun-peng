use kr2r::compact_hash::{HashConfig, Slot};
use kr2r::mmscanner::MinimizerScanner;
use kr2r::seq::{self, open_fasta_reader, SeqX};
use kr2r::utils::{
    create_partition_files, create_partition_writers, create_sample_file, detect_file_format,
    get_file_limit, get_lastest_file_index, FileFormat,
};
use kr2r::{IndexOptions, Meros};
use seq_io::fasta::Record;
use seq_io::fastq::Record as FqRecord;
use seq_io::parallel::read_parallel;
use std::fs;
use std::io::{BufWriter, Write};
use std::io::{Error, ErrorKind, Result};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use clap::Parser;
/// Command line arguments for the splitr program.
///
/// This structure defines the command line arguments that are accepted by the splitr program.
/// It uses the `clap` crate for parsing command line arguments.
#[derive(Parser, Debug, Clone)]
#[clap(
    version,
    about = "Split fast(q/a) file into ranges",
    long_about = "Split fast(q/a) file into ranges"
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

    /// The number of threads to use, default is 10.
    #[clap(short = 'p', long = "num-threads", value_parser, default_value_t = 10)]
    pub num_threads: i32,

    /// chunk directory
    #[clap(long)]
    pub chunk_dir: PathBuf,

    /// A list of input file paths (FASTA/FASTQ) to be processed by the classify program.
    // #[clap(short = 'F', long = "files")]
    pub input_files: Vec<String>,
}

fn init_chunk_writers(
    args: &Args,
    partition: usize,
    chunk_size: usize,
) -> Vec<BufWriter<fs::File>> {
    let chunk_files = create_partition_files(partition, &args.chunk_dir, "sample");

    let mut writers = create_partition_writers(&chunk_files);

    writers.iter_mut().enumerate().for_each(|(index, writer)| {
        // 获取对应的文件大小
        let file_size = writer
            .get_ref()
            .metadata()
            .expect("Failed to get file metadata")
            .len();

        if file_size == 0 {
            writer
                .write_all(&index.to_le_bytes())
                .expect("Failed to write partition");

            let chunk_size_bytes = chunk_size.to_le_bytes();
            writer
                .write_all(&chunk_size_bytes)
                .expect("Failed to write chunk size");

            writer.flush().expect("Failed to flush writer");
        }
    });

    writers
}

/// 处理record
fn process_record<I>(
    iter: I,
    hash_config: &HashConfig,
    seq_id: u64,
    chunk_size: usize,
    idx_bits: usize,
    seq_index: &AtomicUsize,
) -> (usize, Vec<(usize, Slot<u64>)>)
where
    I: Iterator<Item = u64>,
{
    let mut k2_slot_list = Vec::new();
    let mut kmer_count = 0;

    for hash_key in iter.into_iter() {
        let mut slot = hash_config.slot_u64(hash_key, seq_id);
        let seq_sort = seq_index.fetch_add(1, Ordering::SeqCst);
        let partition_index = slot.idx / chunk_size;

        slot.idx = seq_sort << idx_bits | (slot.idx % chunk_size);

        k2_slot_list.push((partition_index, slot));
        kmer_count += 1;
    }
    (kmer_count, k2_slot_list)
}

fn write_data_to_file(
    k2_map: String,
    k2_slot_list: Vec<(usize, Slot<u64>)>,
    writers: &mut Vec<BufWriter<fs::File>>,
    slot_size: usize,
    sample_writer: &mut BufWriter<fs::File>,
) {
    for slot in k2_slot_list {
        let partition_index = slot.0;
        if let Some(writer) = writers.get_mut(partition_index) {
            writer.write_all(slot.1.as_slice(slot_size)).unwrap();
        }
    }

    sample_writer.write_all(k2_map.as_bytes()).unwrap();
}

fn process_fastq_file(
    args: &Args,
    meros: Meros,
    hash_config: HashConfig,
    file_index: usize,
    files: &[String],
    writers: &mut Vec<BufWriter<fs::File>>,
    sample_writer: &mut BufWriter<fs::File>,
) {
    let chunk_size = hash_config.hash_capacity;
    let idx_bits = ((chunk_size as f64).log2().ceil() as usize).max(1);
    let slot_size = std::mem::size_of::<Slot<u64>>();
    let score = args.minimum_quality_score;

    let mut files_iter = files.iter();
    let file1 = files_iter.next().cloned().unwrap();
    let file2 = files_iter.next().cloned();

    let line_index = AtomicUsize::new(0);

    let reader = seq::PairFastqReader::from_path(&file1, file2.as_ref())
        .expect("Unable to create pair reader from paths");
    read_parallel(
        reader,
        args.num_threads as u32,
        args.num_threads as usize,
        |record_set| {
            let mut k2_slot_list = Vec::new();

            let mut buffer = String::new();

            for records in record_set.into_iter() {
                let dna_id = records.0.id().unwrap_or_default().to_string();
                // 拼接seq_id
                let index = line_index.fetch_add(1, Ordering::SeqCst);
                let seq_id = (file_index << 32 | index) as u64;
                let seq_index = AtomicUsize::new(0);

                let seq1 = records.0.seq_x(score);
                let scan1 = MinimizerScanner::new(&seq1, meros);

                let (kmer_count1, slot_list1) = process_record(
                    scan1,
                    &hash_config,
                    seq_id,
                    chunk_size,
                    idx_bits,
                    &seq_index,
                );

                k2_slot_list.extend(slot_list1);
                let (kmer_count, seq_size) = if let Some(record3) = records.1 {
                    let seq2 = record3.seq_x(score);
                    let scan2 = MinimizerScanner::new(&seq2, meros);
                    let (kmer_count2, slot_list2) = process_record(
                        scan2,
                        &hash_config,
                        seq_id,
                        chunk_size,
                        idx_bits,
                        &seq_index,
                    );
                    k2_slot_list.extend(slot_list2);
                    (
                        format!("{}|{}", kmer_count1, kmer_count2),
                        format!("{}|{}", seq1.len(), seq2.len()),
                    )
                } else {
                    (kmer_count1.to_string(), format!("{}", seq1.len()))
                };

                buffer.push_str(
                    format!("{}\t{}\t{}\t{}\n", index, dna_id, seq_size, kmer_count).as_str(),
                );
            }
            (buffer, k2_slot_list)
        },
        |record_sets| {
            while let Some(Ok((_, (k2_map, k2_slot_list)))) = record_sets.next() {
                write_data_to_file(k2_map, k2_slot_list, writers, slot_size, sample_writer);
            }
        },
    )
}

fn process_fasta_file(
    args: &Args,
    meros: Meros,
    hash_config: HashConfig,
    file_index: usize,
    files: &[String],
    writers: &mut Vec<BufWriter<fs::File>>,
    sample_writer: &mut BufWriter<fs::File>,
) {
    let chunk_size = hash_config.hash_capacity;
    let idx_bits = ((chunk_size as f64).log2().ceil() as usize).max(1);
    let slot_size = std::mem::size_of::<Slot<u64>>();
    let score = args.minimum_quality_score;

    let mut files_iter = files.iter();
    let file1 = files_iter.next().cloned().unwrap();

    let line_index = AtomicUsize::new(0);

    let reader = open_fasta_reader(&file1).expect("Unable to create fasta reader from path");
    read_parallel(
        reader,
        args.num_threads as u32,
        args.num_threads as usize,
        |record_set| {
            let mut k2_slot_list = Vec::new();

            let mut buffer = String::new();

            for records in record_set.into_iter() {
                let dna_id = records.id().unwrap_or_default().to_string();
                // 拼接seq_id
                let index = line_index.fetch_add(1, Ordering::SeqCst);
                let seq_id = (file_index << 32 | index) as u64;
                let seq_index = AtomicUsize::new(0);

                let seq1 = records.seq_x(score);
                let scan1 = MinimizerScanner::new(&seq1, meros);

                let (kmer_count1, slot_list) = process_record(
                    scan1,
                    &hash_config,
                    seq_id,
                    chunk_size,
                    idx_bits,
                    &seq_index,
                );

                k2_slot_list.extend(slot_list);
                let (kmer_count, seq_size) = (kmer_count1.to_string(), format!("{}", seq1.len()));
                buffer.push_str(
                    format!("{}\t{}\t{}\t{}\n", index, dna_id, seq_size, kmer_count).as_str(),
                );
            }
            (buffer, k2_slot_list)
        },
        |record_sets| {
            while let Some(Ok((_, (k2_map, k2_slot_list)))) = record_sets.next() {
                write_data_to_file(k2_map, k2_slot_list, writers, slot_size, sample_writer);
            }
        },
    )
}

fn convert(args: Args, meros: Meros, hash_config: HashConfig) -> Result<()> {
    let partition = hash_config.partition;
    let mut writers: Vec<BufWriter<fs::File>> =
        init_chunk_writers(&args, partition, hash_config.hash_capacity);

    let file_path = args.chunk_dir.join("sample_file.map");
    let mut file_writer = create_sample_file(&file_path);
    // 如果文件内容为空，则默认最大值为0
    let mut file_index = get_lastest_file_index(&file_path)?;

    let mut process_files = |files: Vec<&[String]>| -> Result<()> {
        let file_bits = (((files.len() + file_index) as f64).log2().ceil() as usize).max(1);
        if file_bits > hash_config.value_bits {
            panic!("The number of files is too large to process.");
        }

        for file_pair in files {
            file_index += 1;

            writeln!(file_writer, "{}\t{}", file_index, file_pair.join(","))?;
            file_writer.flush().unwrap();

            create_sample_file(
                args.chunk_dir
                    .join(format!("sample_file_{}.bin", file_index)),
            );
            let mut sample_writer =
                create_sample_file(args.chunk_dir.join(format!("sample_id_{}.map", file_index)));

            match detect_file_format(&file_pair[0])? {
                FileFormat::Fastq => {
                    process_fastq_file(
                        &args,
                        meros,
                        hash_config,
                        file_index,
                        file_pair,
                        &mut writers,
                        &mut sample_writer,
                    );
                }
                FileFormat::Fasta => {
                    process_fasta_file(
                        &args,
                        meros,
                        hash_config,
                        file_index,
                        file_pair,
                        &mut writers,
                        &mut sample_writer,
                    );
                }
            }
        }
        Ok(())
    };

    if args.paired_end_processing && !args.single_file_pairs {
        // 处理成对的文件
        let files = args.input_files.chunks(2).collect();
        process_files(files)?;
    } else {
        let files = args.input_files.chunks(1).collect();
        process_files(files)?;
    }

    Ok(())
}

pub fn run(args: Args) -> Result<()> {
    // let args = Args::parse();
    let options_filename = &args.k2d_dir.join("opts.k2d");
    let idx_opts = IndexOptions::read_index_options(options_filename)?;

    if args.paired_end_processing && !args.single_file_pairs && args.input_files.len() % 2 != 0 {
        // 验证文件列表是否为偶数个
        return Err(Error::new(
            ErrorKind::InvalidInput,
            "Paired-end processing requires an even number of input files.",
        ));
    }
    let hash_config = HashConfig::from_hash_header(&args.k2d_dir.join("hash_config.k2d"))?;

    println!("hash_config {:?}", hash_config);
    if hash_config.hash_capacity == 0 {
        panic!("`hash_capacity` can't be zero!");
    }
    println!("start...");
    let file_num_limit = get_file_limit();
    if hash_config.partition >= file_num_limit {
        panic!("Exceeds File Number Limit");
    }

    let meros = idx_opts.as_meros();
    let start = Instant::now();
    convert(args, meros, hash_config)?;
    let duration = start.elapsed();
    println!("splitr took: {:?}", duration);

    Ok(())
}

#[allow(dead_code)]
fn main() {
    let args = Args::parse();
    if let Err(e) = run(args) {
        eprintln!("Application error: {}", e);
    }
}
