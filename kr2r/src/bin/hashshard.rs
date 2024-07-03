use clap::Parser;
use kr2r::args::parse_size;
use kr2r::compact_hash::HashConfig;
// use memmap2::MmapOptions;
use std::fs::{self, create_dir_all, File, OpenOptions};
use std::io::BufWriter;
use std::io::{self, BufReader, Read, Result as IOResult, Seek, Write};
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;

fn mmap_read_write<P: AsRef<Path>, Q: AsRef<Path>>(
    source_path: P,
    dest_path: Q,
    partition: usize,
    cap: usize,
    offset: u64,
    length: usize,
) -> IOResult<()> {
    // 打开目标文件，准备写入数据
    let mut dest_file = BufWriter::new(File::create(dest_path)?);
    dest_file
        .write_all(&partition.to_le_bytes())
        .expect("Failed to write capacity");
    dest_file
        .write_all(&cap.to_le_bytes())
        .expect("Failed to write capacity");

    // 打开源文件，并创建一个缓冲读取器
    let mut file = OpenOptions::new().read(true).open(&source_path)?;
    file.seek(io::SeekFrom::Start(offset))?;
    let mut reader = BufReader::new(file);

    // 创建一个缓冲区，用于存储读取的数据
    let mut buffer = vec![0; length];
    reader.read_exact(&mut buffer)?;

    // 将读取的数据写入目标文件
    dest_file.write_all(&buffer)?;
    // let file = OpenOptions::new().read(true).open(&source_path)?;
    // let mmap = unsafe { MmapOptions::new().offset(offset).len(length).map(&file)? };

    // // 将内存映射的数据写入目标文件
    // dest_file.write_all(&mmap)?;

    Ok(())
}

#[derive(Parser, Debug, Clone)]
#[clap(
    version,
    about = "Convert Kraken2 database files to Kun-peng database format for efficient processing and analysis."
)]
pub struct Args {
    /// The database directory for the Kraken 2 index. contains index files(hash.k2d opts.k2d taxo.k2d)
    #[clap(long = "db", value_parser, required = true)]
    database: PathBuf,

    // /// database hash chunk directory and other files
    // #[clap(long)]
    // k2d_dir: Option<PathBuf>,
    /// Specifies the hash file capacity. Acceptable formats include numeric values followed by 'K', 'M', or 'G' (e.g., '1.5G', '250M', '1024K').
    /// Note: The specified capacity affects the index size, with a factor of 4 applied. For example, specifying '1G' results in an index size of '4G'.
    /// Default: 1G (capacity 1G = file size 4G)
    #[clap(long = "hash-capacity", value_parser = parse_size, default_value = "1G", help = "Specifies the hash file capacity.\nAcceptable formats include numeric values followed by 'K', 'M', or 'G' (e.g., '1.5G', '250M', '1024K').\nNote: The specified capacity affects the index size, with a factor of 4 applied.\nFor example, specifying '1G' results in an index size of '4G'.\nDefault: 1G (capacity 1G = file size 4G)")]
    hash_capacity: usize,
}

pub fn run(args: Args) -> IOResult<()> {
    let index_filename = &args.database.join("hash.k2d");
    let hash_config = HashConfig::from_hash_header(index_filename)?;

    let partition = (hash_config.capacity + args.hash_capacity - 1) / args.hash_capacity;
    println!("hashshard start...");
    // 开始计时
    let start = Instant::now();

    let file_len = hash_config.capacity * 4 + 32;
    let b_size = std::mem::size_of::<u32>();

    let k2d_dir = args.database.clone();

    create_dir_all(&k2d_dir).expect(&format!("create hash dir error {:?}", k2d_dir));

    let config_file = k2d_dir.join("hash_config.k2d");
    if config_file.exists() {
        panic!("hash config is exists!!!");
    }

    mmap_read_write(
        &index_filename,
        config_file,
        partition,
        args.hash_capacity,
        0,
        32,
    )?;

    for i in 1..=partition {
        let chunk_file = k2d_dir.join(format!("hash_{}.k2d", i));
        let offset = (32 + args.hash_capacity * (i - 1) * b_size) as u64;
        let mut length = args.hash_capacity * b_size;
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

    let source_taxo_file = &args.database.join("taxo.k2d");
    let dst_tax_file = k2d_dir.join("taxo.k2d");
    if !dst_tax_file.exists() {
        fs::copy(source_taxo_file, dst_tax_file)?;
    }

    let source_opts_file = &args.database.join("opts.k2d");
    let dst_opts_file = k2d_dir.join("opts.k2d");
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
