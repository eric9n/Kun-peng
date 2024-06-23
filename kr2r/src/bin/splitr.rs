use clap::Parser;
use kr2r::compact_hash::{HashConfig, Slot};
use kr2r::utils::{
    create_partition_files, create_partition_writers, create_sample_file, get_file_limit,
    get_lastest_file_index,
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
                let dna_id = header.id.clone();
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
            Some((buffer, k2_slot_list))
        },
        |dataset| {
            while let Some(Some((buffer, k2_slot_list))) = dataset.next() {
                write_data_to_file(buffer, k2_slot_list, writers, slot_size, sample_writer);
            }
        },
    )
    .expect("failed");

    Ok(())
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

            let score = args.minimum_quality_score;
            let paths = OptionPair::from_slice(file_pair);
            let mut reader = FastxReader::from_paths(paths, file_index, score)?;
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

    println!("hash_config {:?}", hash_config);
    if hash_config.hash_capacity == 0 {
        panic!("`hash_capacity` can't be zero!");
    }
    println!("splitr start...");
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
