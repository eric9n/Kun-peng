use clap::Parser;
use kr2r::compact_hash::{Compact, HashConfig, Slot};
use kr2r::mmscanner::MinimizerScanner;
use kr2r::seq::{self, SeqX};
use kr2r::utils::{
    create_partition_files, create_partition_writers, create_sample_map, detect_file_format,
    find_and_sort_files, get_file_limit, FileFormat,
};
use kr2r::{IndexOptions, Meros};
use seq_io::fastq::Record;
// use std::collections::HashMap;
use std::fs;
use std::fs::File;
use std::io::{self, BufReader, Error, ErrorKind, Read, Result, Seek, Write};
use std::path::Path;
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

    /// chunk directory
    #[clap(long)]
    chunk_dir: PathBuf,

    #[clap(long, default_value = "sample")]
    chunk_prefix: String,
}

// 定义每批次处理的 Slot 数量
const BATCH_SIZE: usize = 81920;

fn read_chunk_header<P: AsRef<Path>>(file_path: P) -> io::Result<(u64, u64)> {
    let file = File::open(file_path)?;
    let mut reader = BufReader::new(file);
    let mut buffer = [0u8; 16]; // u64 + u64 = 8 bytes + 8 bytes

    reader.read_exact(&mut buffer)?;

    let index = u64::from_le_bytes(
        buffer[0..8]
            .try_into()
            .expect("Failed to convert bytes to u64 for index"),
    );
    let chunk_size = u64::from_le_bytes(
        buffer[8..16]
            .try_into()
            .expect("Failed to convert bytes to u64 for chunk size"),
    );

    Ok((index, chunk_size))
}

fn process_chunk_file<P: AsRef<Path>>(chunk_file: P) -> Result<()> {
    let file = File::open(chunk_file)?;
    let mut reader = BufReader::new(file);

    let slot_size = std::mem::size_of::<Slot<u64>>();
    let batch_buffer_size = slot_size * BATCH_SIZE;
    let mut batch_buffer = vec![0u8; batch_buffer_size];

    let mut buffer = [0u8; 16]; // u64 + u64 = 8 bytes + 8 bytes
    reader.read_exact(&mut buffer)?;

    while let Ok(bytes_read) = reader.read(&mut batch_buffer) {
        if bytes_read == 0 {
            break;
        } // 文件末尾

        // 处理读取的数据批次
        let slots_in_batch = bytes_read / slot_size;

        let slots = unsafe {
            std::slice::from_raw_parts(batch_buffer.as_ptr() as *const Slot<u64>, slots_in_batch)
        };

        slots.into_iter().for_each(|slot| {
            println!("slot {:?}", slot.idx);
        });
    }
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    let idx_opts = IndexOptions::read_index_options(args.options_filename.clone())?;

    let hash_config = HashConfig::<u32>::from(args.index_filename.clone())?;

    // let partition = (hash_config.capacity + args.chunk_size - 1) / args.chunk_size;

    let chunk_files = find_and_sort_files(&args.chunk_dir, &args.chunk_prefix, ".k2")?;
    println!("chunk files {:?}", chunk_files);
    for chunk_file in chunk_files {
        let (index, chunk_size) = read_chunk_header(&chunk_file)?;
        println!("index {:?}", index);
        println!("chunk_size {:?}", chunk_size);
        process_chunk_file(chunk_file)?;
    }
    Ok(())
}
