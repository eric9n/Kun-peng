use clap::Parser;
use kr2r::compact_hash::{CHTable, Compact, HashConfig, Row, Slot};
use kr2r::utils::{find_and_sort_files, open_file};
// use std::collections::HashMap;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{self, BufReader, BufWriter, Read, Result, Write};
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;

// 定义每批次处理的 Slot 数量
pub const BATCH_SIZE: usize = 8 * 1024 * 1024;

/// Command line arguments for the splitr program.
///
/// This structure defines the command line arguments that are accepted by the splitr program.
/// It uses the `clap` crate for parsing command line arguments.
#[derive(Parser, Debug, Clone)]
#[clap(
    version,
    about = "annotate a set of sequences",
    long_about = "annotate a set of sequences"
)]
pub struct Args {
    /// database hash chunk directory and other files
    #[clap(long)]
    pub k2d_dir: PathBuf,

    /// Enables use of a Kraken 2 compatible shared database. Default is false.
    #[clap(long, default_value_t = false)]
    pub kraken_db_type: bool,

    /// chunk directory
    #[clap(long)]
    pub chunk_dir: PathBuf,

    /// 批量处理大小 default: 8MB
    #[clap(long, default_value_t = BATCH_SIZE)]
    pub batch_size: usize,
}

fn read_chunk_header<R: Read>(reader: &mut R) -> io::Result<(usize, usize)> {
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

    Ok((index as usize, chunk_size as usize))
}

fn write_to_file(
    file_index: u64,
    bytes: &[u8],
    last_file_index: &mut Option<u64>,
    writer: &mut Option<BufWriter<File>>,
    chunk_dir: &PathBuf,
) -> io::Result<()> {
    if last_file_index.is_none() || last_file_index.unwrap() != file_index {
        if let Some(mut w) = writer.take() {
            w.flush()?;
        }

        let file_name = format!("sample_file_{}.bin", file_index);
        let file_path = chunk_dir.join(file_name);
        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&file_path)?;
        *writer = Some(BufWriter::new(file));

        *last_file_index = Some(file_index);
    }

    if let Some(w) = writer.as_mut() {
        w.write_all(bytes)?;
    }

    Ok(())
}

fn process_batch<R>(
    reader: &mut R,
    hash_config: &HashConfig,
    chtm: &CHTable,
    chunk_dir: PathBuf,
    batch_size: usize,
    page_index: usize,
) -> std::io::Result<()>
where
    R: Read + Send,
{
    let slot_size = std::mem::size_of::<Slot<u64>>();
    let row_size = std::mem::size_of::<Row>();
    let mut batch_buffer = vec![0u8; slot_size * batch_size];
    let mut last_file_index: Option<u64> = None;
    let mut writer: Option<BufWriter<File>> = None;

    let value_mask = hash_config.get_value_mask();
    let value_bits = hash_config.get_value_bits();
    let idx_mask = hash_config.get_idx_mask();
    let idx_bits = hash_config.get_idx_bits();

    while let Ok(bytes_read) = reader.read(&mut batch_buffer) {
        if bytes_read == 0 {
            break;
        } // 文件末尾

        // 处理读取的数据批次
        let slots_in_batch = bytes_read / slot_size;

        let slots = unsafe {
            std::slice::from_raw_parts(batch_buffer.as_ptr() as *const Slot<u64>, slots_in_batch)
        };

        let result: HashMap<u64, Vec<u8>> = slots
            .into_par_iter()
            .filter_map(|slot| {
                let indx = slot.idx & idx_mask;
                let compacted = slot.value.left(value_bits) as u32;
                let taxid = chtm.get_from_page(indx, compacted, page_index);

                if taxid > 0 {
                    let kmer_id = slot.idx >> idx_bits;
                    let file_index = slot.value.right(value_mask) >> 32;
                    let seq_id = slot.get_seq_id() as u32;
                    let left = slot.value.left(value_bits) as u32;
                    let high = u32::combined(left, taxid, value_bits);
                    let row = Row::new(high, seq_id, kmer_id as u32);
                    // let value = slot.to_b(high);
                    // let value_bytes = value.to_le_bytes(); // 将u64转换为[u8; 8]
                    let value_bytes = row.as_slice(row_size);
                    Some((file_index, value_bytes.to_vec()))
                } else {
                    None
                }
            })
            .fold(
                || HashMap::new(),
                |mut acc: HashMap<u64, Vec<u8>>, (file_index, value_bytes)| {
                    acc.entry(file_index)
                        .or_insert_with(Vec::new)
                        .extend(value_bytes);
                    acc
                },
            )
            .reduce(
                || HashMap::new(),
                |mut acc, h| {
                    for (k, mut v) in h {
                        acc.entry(k).or_insert_with(Vec::new).append(&mut v);
                    }
                    acc
                },
            );

        let mut file_indices: Vec<_> = result.keys().cloned().collect();
        file_indices.sort_unstable(); // 对file_index进行排序

        for file_index in file_indices {
            if let Some(bytes) = result.get(&file_index) {
                write_to_file(
                    file_index,
                    bytes,
                    &mut last_file_index,
                    &mut writer,
                    &chunk_dir,
                )?;
            }
        }
    }

    if let Some(w) = writer.as_mut() {
        w.flush()?;
    }

    Ok(())
}

fn process_chunk_file<P: AsRef<Path>>(
    args: &Args,
    chunk_file: P,
    hash_files: &Vec<PathBuf>,
) -> Result<()> {
    let file = open_file(chunk_file)?;
    let mut reader = BufReader::new(file);

    let (page_index, _) = read_chunk_header(&mut reader)?;

    let start = Instant::now();

    let config = HashConfig::from_hash_header(&args.k2d_dir.join("hash_config.k2d"))?;
    let parition = hash_files.len();
    let chtm = if args.kraken_db_type {
        CHTable::from_pair(
            config,
            &hash_files[page_index],
            &hash_files[(page_index + 1) % parition],
        )?
    } else {
        CHTable::from(config, &hash_files[page_index])?
    };

    // 计算持续时间
    let duration = start.elapsed();
    // 打印运行时间
    println!("load table took: {:?}", duration);
    process_batch(
        &mut reader,
        &config,
        &chtm,
        args.chunk_dir.clone(),
        args.batch_size,
        page_index,
    )?;

    Ok(())
}

pub fn run(args: Args) -> Result<()> {
    let chunk_files = find_and_sort_files(&args.chunk_dir, "sample", ".k2")?;

    let hash_files = find_and_sort_files(&args.k2d_dir, "hash", ".k2d")?;

    // 开始计时
    let start = Instant::now();
    println!("start...");
    for chunk_file in chunk_files {
        println!("chunk_file {:?}", chunk_file);
        process_chunk_file(&args, chunk_file, &hash_files)?;
    }
    // 计算持续时间
    let duration = start.elapsed();
    // 打印运行时间
    println!("annotate took: {:?}", duration);

    Ok(())
}

#[allow(dead_code)]
fn main() {
    let args = Args::parse();
    if let Err(e) = run(args) {
        eprintln!("Application error: {}", e);
    }
}
