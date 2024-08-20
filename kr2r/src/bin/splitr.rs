use clap::Parser;
use kr2r::compact_hash::{HashConfig, Slot};
use kr2r::utils::{
    create_partition_files, create_partition_writers, create_sample_file, get_file_limit,
    get_lastest_file_index, set_fd_limit,
};
use kr2r::IndexOptions;
use seqkmer::{read_parallel, FastxReader, Meros, MinimizerIterator, OptionPair, Reader};
use std::fs;
use std::io::{BufWriter, Write};
use std::io::{Error, ErrorKind, Result};
use std::path::PathBuf;
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
pub struct Args {
    /// database hash chunk directory and other files
    #[arg(long = "db", required = true)]
    pub database: PathBuf,

    // /// The file path for the Kraken 2 options.
    // #[clap(short = 'o', long = "options-filename", value_parser, required = true)]
    // options_filename: String,
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

    /// chunk directory
    #[clap(long)]
    pub chunk_dir: PathBuf,

    /// A list of input file paths (FASTA/FASTQ) to be processed by the classify program.
    /// Supports fasta or fastq format files (e.g., .fasta, .fastq) and gzip compressed files (e.g., .fasta.gz, .fastq.gz).
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
fn process_record(
    k2_slot_list: &mut Vec<(usize, Slot<u64>)>,
    marker: &mut MinimizerIterator,
    hash_config: &HashConfig,
    chunk_size: usize,
    seq_id: u64,
    idx_bits: usize,
) {
    let offset = k2_slot_list.len();
    for (sort, hash_key) in marker {
        let mut slot = hash_config.slot_u64(hash_key, seq_id);
        let seq_sort = sort + offset;
        let partition_index = slot.idx / chunk_size;

        slot.idx = seq_sort << idx_bits | (slot.idx % chunk_size);
        k2_slot_list.push((partition_index, slot));
    }
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

fn process_fastx_file<R>(
    args: &Args,
    meros: Meros,
    hash_config: HashConfig,
    file_index: usize,
    reader: &mut R,
    writers: &mut Vec<BufWriter<fs::File>>,
    sample_writer: &mut BufWriter<fs::File>,
) -> Result<()>
where
    R: Reader,
{
    let chunk_size = hash_config.hash_capacity;
    let idx_bits = ((chunk_size as f64).log2().ceil() as usize).max(1);
    let slot_size = std::mem::size_of::<Slot<u64>>();

    read_parallel(
        reader,
        args.num_threads as usize,
        &meros,
        |seqs| {
            let mut buffer = String::new();
            let mut k2_slot_list = Vec::new();
            for seq in seqs {
                let mut init: Vec<(usize, Slot<u64>)> = Vec::new();
                let header = &seq.header;
                let index = header.reads_index;
                let dna_id = header.id.trim();
                let seq_id = (file_index << 32 | index) as u64;

                seq.body.apply_mut(|m_iter| {
                    process_record(
                        &mut init,
                        m_iter,
                        &hash_config,
                        chunk_size,
                        seq_id,
                        idx_bits,
                    );
                });
                k2_slot_list.extend_from_slice(&init);

                let size_str = seq.fmt_size();
                let seq_size_str = seq.fmt_seq_size();
                buffer.push_str(
                    format!("{}\t{}\t{}\t{}\n", index, dna_id, seq_size_str, size_str).as_str(),
                );
            }
            (buffer, k2_slot_list)
        },
        |dataset| {
            while let Some(data) = dataset.next() {
                let (buffer, k2_slot_list) = data.unwrap();
                write_data_to_file(buffer, k2_slot_list, writers, slot_size, sample_writer);
            }
        },
    )
    .expect("failed");

    Ok(())
}

/// 处理样本文件
fn process_files<F>(args: &Args, hash_config: HashConfig, mut action: F) -> Result<()>
where
    F: FnMut(usize, OptionPair<String>) -> Result<()>,
{
    let file_path = args.chunk_dir.join("sample_file.map");
    let mut file_writer = create_sample_file(&file_path);
    let mut file_index = get_lastest_file_index(&file_path)?;

    let chunk_size = if args.paired_end_processing && !args.single_file_pairs {
        2
    } else {
        1
    };
    let files = args.input_files.chunks(chunk_size).collect::<Vec<_>>();

    let file_bits = (((files.len() + file_index) as f64).log2().ceil() as usize).max(1);
    if file_bits > hash_config.value_bits {
        panic!("The number of files is too large to process.");
    }

    for file_pair in files {
        file_index += 1;
        let path_pair = OptionPair::from_slice(file_pair);
        writeln!(
            file_writer,
            "{}\t{}",
            file_index,
            path_pair.reduce_str(",", |a| a.to_string())
        )?;
        file_writer.flush().unwrap();

        action(file_index, path_pair)?;
    }

    Ok(())
}

pub fn run(args: Args) -> Result<()> {
    // let args = Args::parse();
    let options_filename = &args.database.join("opts.k2d");
    let idx_opts = IndexOptions::read_index_options(options_filename)?;

    if args.paired_end_processing && !args.single_file_pairs && args.input_files.len() % 2 != 0 {
        // 验证文件列表是否为偶数个
        return Err(Error::new(
            ErrorKind::InvalidInput,
            "Paired-end processing requires an even number of input files.",
        ));
    }
    let hash_config = HashConfig::from_hash_header(&args.database.join("hash_config.k2d"))?;

    println!("{:?}", hash_config);
    if hash_config.hash_capacity == 0 {
        panic!("`hash_capacity` can't be zero!");
    }
    println!("splitr start...");
    let file_num_limit = get_file_limit();
    if hash_config.partition >= file_num_limit {
        eprintln!(
            "file num limit {:?}, need: {:?}",
            file_num_limit, hash_config.partition
        );
        set_fd_limit(hash_config.partition as u64 + 1)
            .expect("Failed to set file descriptor limit, please run this operation with administrative/root privileges.");
        // panic!("Exceeds File Number Limit");
    }

    let meros = idx_opts.as_meros();
    let start = Instant::now();
    let partition = hash_config.partition;
    let mut writers: Vec<BufWriter<fs::File>> =
        init_chunk_writers(&args, partition, hash_config.hash_capacity);

    process_files(&args, hash_config, |file_index, path_pair| {
        let mut sample_writer =
            create_sample_file(args.chunk_dir.join(format!("sample_id_{}.map", file_index)));

        let score = args.minimum_quality_score;
        let mut reader = FastxReader::from_paths(path_pair, file_index, score)?;
        process_fastx_file(
            &args,
            meros,
            hash_config,
            file_index,
            &mut reader,
            &mut writers,
            &mut sample_writer,
        )
        .expect("process fastx file error");
        Ok(())
    })?;
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
