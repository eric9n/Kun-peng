use clap::Parser;
use kr2r::compact_hash::HashConfig;
// use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{Result, Write};
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;

use memmap2::MmapOptions;

fn mmap_read_write<P: AsRef<Path>, Q: AsRef<Path>>(
    source_path: P,
    dest_path: Q,
    offset: u64,
    length: usize,
    b_size: usize,
    partition_index: Option<usize>,
) -> Result<()> {
    let file = OpenOptions::new().read(true).open(&source_path)?;

    let mmap = unsafe { MmapOptions::new().offset(offset).len(length).map(&file)? };

    // 打开目标文件，准备写入数据
    let mut dest_file = File::create(dest_path)?;

    if let Some(index) = partition_index {
        dest_file
            .write_all(&index.to_le_bytes())
            .expect("Failed to write partition");

        let size = length / b_size;
        let chunk_size_bytes = size.to_le_bytes();

        dest_file
            .write_all(&chunk_size_bytes)
            .expect("Failed to write chunk size");
    }

    // 将内存映射的数据写入目标文件
    dest_file.write_all(&mmap)?;

    Ok(())
}

/// Command line arguments for the splitr program.
///
/// This structure defines the command line arguments that are accepted by the splitr program.
/// It uses the `clap` crate for parsing command line arguments.
#[derive(Parser, Debug, Clone)]
#[clap(version, about = "HashShard", long_about = "divide hash file")]
struct Args {
    /// The file path for the Kraken 2 index.
    #[clap(short = 'H', long = "index-filename", value_parser, required = true)]
    index_filename: String,

    /// The file path for the Kraken 2 options.
    // #[clap(short = 'o', long = "options-filename", value_parser, required = true)]
    // options_filename: String,

    /// chunk directory
    #[clap(long)]
    hash_dir: PathBuf,

    // default: 1073741824(1G)
    #[clap(long, default_value_t = 1073741824)]
    hash_size: usize,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let hash_config = HashConfig::<u32>::from(args.index_filename.clone())?;

    let partition = (hash_config.capacity + args.hash_size - 1) / args.hash_size;
    println!("start...");
    // 开始计时
    let start = Instant::now();

    let file_len = hash_config.capacity * 4 + 32;
    let b_size = std::mem::size_of::<u32>();

    let config_file = args.hash_dir.join("hash_config.k2d");
    mmap_read_write(&args.index_filename, config_file, 0, 32, b_size, None)?;

    for i in 0..partition {
        let chunk_file = args.hash_dir.join(format!("hash_{}.k2d", i));
        let offset = (32 + args.hash_size * i * b_size) as u64;
        let mut length = args.hash_size * b_size;
        if (offset as usize + length) > file_len {
            length = file_len - offset as usize;
        }
        mmap_read_write(
            &args.index_filename,
            chunk_file,
            offset,
            length,
            b_size,
            Some(i),
        )?
    }

    // 计算持续时间
    let duration = start.elapsed();

    // 打印运行时间
    println!("hashshard took: {:?}", duration);

    Ok(())
}
