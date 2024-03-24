use clap::Parser;
use kr2r::compact_hash::HashConfig;
use std::fs::{self, File, OpenOptions};
use std::io::BufWriter;
use std::io::{Result, Write};
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;

use memmap2::MmapOptions;

fn mmap_read_write<P: AsRef<Path>, Q: AsRef<Path>>(
    source_path: P,
    dest_path: Q,
    partition: usize,
    cap: usize,
    offset: u64,
    length: usize,
) -> Result<()> {
    // 打开目标文件，准备写入数据
    let mut dest_file = BufWriter::new(File::create(dest_path)?);
    dest_file
        .write_all(&partition.to_le_bytes())
        .expect("Failed to write capacity");
    dest_file
        .write_all(&cap.to_le_bytes())
        .expect("Failed to write capacity");

    let file = OpenOptions::new().read(true).open(&source_path)?;
    let mmap = unsafe { MmapOptions::new().offset(offset).len(length).map(&file)? };

    // 将内存映射的数据写入目标文件
    dest_file.write_all(&mmap)?;

    Ok(())
}

#[derive(Parser, Debug, Clone)]
#[clap(version, about = "split hash file", long_about = "split hash file")]
pub struct Args {
    /// The database path for the Kraken 2 index.
    #[clap(long, value_parser, required = true)]
    db: PathBuf,

    /// database hash chunk directory and other files
    #[clap(long)]
    hash_dir: Option<PathBuf>,

    // default: 1073741824(1G)
    #[clap(long, default_value_t = 1073741824)]
    hash_size: usize,
}

pub fn run(args: Args) -> Result<()> {
    let index_filename = &args.db.join("hash.k2d");
    let hash_config = HashConfig::<u32>::from(index_filename)?;

    let partition = (hash_config.capacity + args.hash_size - 1) / args.hash_size;
    println!("start...");
    // 开始计时
    let start = Instant::now();

    let file_len = hash_config.capacity * 4 + 32;
    let b_size = std::mem::size_of::<u32>();

    let hash_dir = args.hash_dir.unwrap_or(args.db.clone());
    let config_file = hash_dir.join("hash_config.k2d");
    if config_file.exists() {
        panic!("hash config is exists!!!");
    }

    mmap_read_write(
        &index_filename,
        config_file,
        partition,
        args.hash_size,
        0,
        32,
    )?;

    for i in 1..=partition {
        let chunk_file = hash_dir.join(format!("hash_{}.k2d", i));
        let offset = (32 + args.hash_size * (i - 1) * b_size) as u64;
        let mut length = args.hash_size * b_size;
        if (offset as usize + length) > file_len {
            length = file_len - offset as usize;
        }
        let cap = length / b_size;
        mmap_read_write(&index_filename, chunk_file, i, cap, offset, length)?
    }

    // 计算持续时间
    let duration = start.elapsed();

    // 打印运行时间
    println!("hashshard took: {:?}", duration);

    let source_taxo_file = &args.db.join("taxo.k2d");
    let dst_tax_file = hash_dir.join("taxo.k2d");
    if !dst_tax_file.exists() {
        fs::copy(source_taxo_file, dst_tax_file)?;
    }

    let source_opts_file = &args.db.join("opts.k2d");
    let dst_opts_file = hash_dir.join("opts.k2d");
    if !dst_opts_file.exists() {
        fs::copy(source_opts_file, dst_opts_file)?;
    }

    Ok(())
}

#[allow(dead_code)]
fn main() {
    let args = Args::parse();
    if let Err(e) = run(args) {
        eprintln!("Application error: {}", e);
    }
}
