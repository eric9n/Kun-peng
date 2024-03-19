use clap::Parser;
use kr2r::compact_hash::{HashConfig, Slot};
use kr2r::mmscanner::MinimizerScanner;
use kr2r::seq::{self, SeqX};
use kr2r::utils::{
    create_partition_files, create_partition_writers, create_sample_map, detect_file_format,
    get_file_limit, FileFormat,
};
use kr2r::{IndexOptions, Meros};
use seq_io::fastq::{Reader, Record};
// use std::collections::HashMap;
use seq_io::parallel::read_parallel;
use std::fs;
use std::io::{self, BufWriter, Write};
use std::io::{Error, ErrorKind, Result};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

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
struct Args {
    /// The file path for the Kraken 2 index.
    #[clap(short = 'H', long = "index-filename", value_parser, required = true)]
    index_filename: String,

    /// The file path for the Kraken 2 options.
    #[clap(short = 'o', long = "options-filename", value_parser, required = true)]
    options_filename: String,

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

    // default: 1073741824(1G)
    #[clap(long, default_value_t = 1073741824)]
    chunk_size: usize,

    #[clap(long, default_value = "sample")]
    chunk_prefix: String,

    #[clap(long, default_value = "sample_file.map")]
    file_map: String,

    #[clap(long, default_value = "sample_id")]
    id_map_prefix: String,

    /// Input files for processing.
    ///
    /// A list of input file paths (FASTA/FASTQ) to be processed by the classify program.
    // #[clap(short = 'F', long = "files")]
    input_files: Vec<String>,
}

fn init_chunk_writers(args: &Args, partition: usize) -> Vec<BufWriter<fs::File>> {
    let chunk_files = create_partition_files(partition, &args.chunk_dir, &args.chunk_prefix);

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

            let chunk_size_bytes = args.chunk_size.to_le_bytes();
            writer
                .write_all(&chunk_size_bytes)
                .expect("Failed to write chunk size");

            writer.flush().expect("Failed to flush writer");
        }
    });

    writers
}

/// 获取最新的文件序号
fn get_lastest_file_index(file_path: &PathBuf) -> Result<usize> {
    let file_content = fs::read_to_string(&file_path)?;
    // 如果文件内容为空，则默认最大值为0
    let index = if file_content.is_empty() {
        0
    } else {
        file_content
            .lines() // 将内容按行分割
            .filter_map(|line| line.split('\t').next()) // 获取每行的第一列
            .filter_map(|num_str| num_str.parse::<usize>().ok()) // 尝试将第一列的字符串转换为整型
            .max() // 找到最大值
            .unwrap_or(1)
    };
    Ok(index)
}

/// 处理record
fn process_record<I>(
    iter: I,
    hash_config: &HashConfig<u32>,
    seq_id: u64,
    chunk_size: usize,
) -> (usize, Vec<(usize, Slot<u64>)>)
where
    I: Iterator<Item = u64>,
{
    let mut k2_slot_list = Vec::new();
    let mut kmer_count = 0;

    for hash_key in iter.into_iter() {
        let mut slot = hash_config.slot_u64(hash_key, seq_id);
        let partition_index = slot.idx / chunk_size;
        slot.idx = slot.idx % chunk_size;
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

fn convert(args: Args, meros: Meros, hash_config: HashConfig<u32>, partition: usize) -> Result<()> {
    let chunk_size = args.chunk_size;
    let slot_size = std::mem::size_of::<Slot<u64>>();
    let score = args.minimum_quality_score;
    let mut writers: Vec<BufWriter<fs::File>> = init_chunk_writers(&args, partition);

    let file_path = args.chunk_dir.join(args.file_map);
    let mut file_writer = create_sample_map(&file_path);
    // 如果文件内容为空，则默认最大值为0
    let mut file_index = get_lastest_file_index(&file_path)?;

    if args.paired_end_processing && !args.single_file_pairs {
        // 处理成对的文件
        for file_pair in args.input_files.chunks(2) {
            let file1 = &file_pair[0];
            let file2 = &file_pair[1];
            file_index += 1;

            writeln!(file_writer, "{}\t{},{}", file_index, file1, file2)?;
            file_writer.flush().unwrap();

            let line_index = AtomicUsize::new(0);

            let file_bits = ((file_index as f64).log2().ceil() as usize).max(1);
            if file_bits > 32 - hash_config.value_bits {
                panic!("The number of files is too large to process.");
            }

            let mut sample_writer = create_sample_map(
                args.chunk_dir
                    .join(format!("{}_{}.map", args.id_map_prefix, file_index)),
            );

            match detect_file_format(&file1)? {
                FileFormat::Fastq => {
                    let reader = seq::PairFastqReader::from_path(file1, file2)
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
                                let seq1 = records.0.seq_x(score);
                                let seq2 = records.1.seq_x(score);

                                line_index.fetch_add(1, Ordering::SeqCst);
                                let index = line_index.load(Ordering::SeqCst);

                                // 拼接seq_id
                                let seq_id = (file_index << 32 | index) as u64;

                                let scan1 = MinimizerScanner::new(&seq1, meros);
                                let scan2 = MinimizerScanner::new(&seq2, meros);

                                let (kmer_count, slot_list) = process_record(
                                    scan1.chain(scan2).into_iter(),
                                    &hash_config,
                                    seq_id,
                                    chunk_size,
                                );

                                k2_slot_list.extend(slot_list);
                                buffer.push_str(
                                    format!("{}\t{}\t{}\n", index, dna_id, kmer_count).as_str(),
                                );
                            }
                            (buffer, k2_slot_list)
                        },
                        |record_sets| {
                            while let Some(Ok((_, (k2_map, k2_slot_list)))) = record_sets.next() {
                                write_data_to_file(
                                    k2_map,
                                    k2_slot_list,
                                    &mut writers,
                                    slot_size,
                                    &mut sample_writer,
                                );
                            }
                        },
                    )
                }
                FileFormat::Fasta => {
                    return Err(io::Error::new(
                        io::ErrorKind::Other,
                        "Unrecognized file format(fasta paired)",
                    ));
                }
            }
        }
    } else {
        for file in args.input_files {
            // 对 file 执行分类处理
            file_index += 1;

            writeln!(file_writer, "{}\t{}", file_index, file)?;
            file_writer.flush().unwrap();

            let line_index = AtomicUsize::new(0);

            let file_bits = ((file_index as f64).log2().ceil() as usize).max(1);
            if file_bits > 32 - hash_config.value_bits {
                panic!("The number of files is too large to process.");
            }

            let mut sample_writer = create_sample_map(
                args.chunk_dir
                    .join(format!("{}_{}.map", args.id_map_prefix, file_index)),
            );
            match detect_file_format(&file)? {
                FileFormat::Fastq => {
                    let reader =
                        Reader::from_path(file).expect("Unable to create pair reader from paths");
                    read_parallel(
                        reader,
                        args.num_threads as u32,
                        args.num_threads as usize,
                        |record_set| {
                            let mut k2_slot_list = Vec::new();

                            let mut buffer = String::new();

                            for records in record_set.into_iter() {
                                let dna_id = records.id().unwrap_or_default().to_string();
                                let seq1 = records.seq_x(score);

                                line_index.fetch_add(1, Ordering::SeqCst);
                                let index = line_index.load(Ordering::SeqCst);

                                // 拼接seq_id
                                let seq_id = (file_index << 32 | index) as u64;

                                let scan1 = MinimizerScanner::new(&seq1, meros);

                                let (kmer_count, slot_list) = process_record(
                                    scan1.into_iter(),
                                    &hash_config,
                                    seq_id,
                                    chunk_size,
                                );

                                k2_slot_list.extend(slot_list);
                                buffer.push_str(
                                    format!("{}\t{}\t{}\n", index, dna_id, kmer_count).as_str(),
                                );
                            }
                            (buffer, k2_slot_list)
                        },
                        |record_sets| {
                            while let Some(Ok((_, (k2_map, k2_slot_list)))) = record_sets.next() {
                                write_data_to_file(
                                    k2_map,
                                    k2_slot_list,
                                    &mut writers,
                                    slot_size,
                                    &mut sample_writer,
                                );
                            }
                        },
                    )
                }
                FileFormat::Fasta => {}
            };
        }
    }

    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    let idx_opts = IndexOptions::read_index_options(args.options_filename.clone())?;

    if args.paired_end_processing && !args.single_file_pairs && args.input_files.len() % 2 != 0 {
        // 验证文件列表是否为偶数个
        return Err(Error::new(
            ErrorKind::InvalidInput,
            "Paired-end processing requires an even number of input files.",
        ));
    }
    let hash_config = HashConfig::<u32>::from(args.index_filename.clone())?;

    let partition = (hash_config.capacity + args.chunk_size - 1) / args.chunk_size;
    println!("start...");

    let file_num_limit = get_file_limit();

    if partition >= file_num_limit {
        panic!("Exceeds File Number Limit");
    }

    let meros = idx_opts.as_meros();
    // 开始计时
    let start = Instant::now();

    convert(args, meros, hash_config, partition)?;
    // 计算持续时间
    let duration = start.elapsed();

    // 打印运行时间
    println!("splitr took: {:?}", duration);

    Ok(())
}
