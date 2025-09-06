use clap::Parser;
use kun_peng::compact_hash::{read_next_page, Compact, HashConfig, Page, Row, Slot};
use kun_peng::utils::{find_and_sort_files, open_file};
use seqkmer::buffer_read_parallel;
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{self, BufReader, BufWriter, Read, Result, Write};
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;

// 定义每批次处理的 Slot 数量
pub const BUFFER_SIZE: usize = 48 * 1024 * 1024;

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
    #[arg(long = "db", required = true)]
    pub database: PathBuf,

    /// chunk directory
    #[clap(long)]
    pub chunk_dir: PathBuf,

    #[clap(long, default_value_t = BUFFER_SIZE)]
    pub buffer_size: usize,

    /// The size of each batch for processing taxid match results, used to control memory usage
    #[clap(long, value_parser = clap::value_parser!(u32).range(1..=32), default_value_t = 4)]
    pub batch_size: u32,

    /// The number of threads to use.
    #[clap(short = 'p', long = "num-threads", value_parser, default_value_t = num_cpus::get())]
    pub num_threads: usize,
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

fn _write_to_file(
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

fn write_to_file(
    file_index: u64,
    seq_id_mod: u32,
    bytes: &[u8],
    writers: &mut HashMap<(u64, u32), BufWriter<File>>,
    chunk_dir: &PathBuf,
) -> io::Result<()> {
    // 检查是否已经有该文件的 writer，没有则创建一个新的
    let writer = writers.entry((file_index, seq_id_mod)).or_insert_with(|| {
        let file_name = format!("sample_file_{}_{}.bin", file_index, seq_id_mod);
        let file_path = chunk_dir.join(file_name);
        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&file_path)
            .expect("failed to open file");
        BufWriter::new(file)
    });

    writer.write_all(bytes)?;

    Ok(())
}

fn clean_up_writers(
    writers: &mut HashMap<(u64, u32), BufWriter<File>>,
    current_file_index: u64,
) -> io::Result<()> {
    let keys_to_remove: Vec<(u64, u32)> = writers
        .keys()
        .cloned()
        .filter(|(idx, _)| *idx != current_file_index)
        .collect();

    for key in keys_to_remove {
        if let Some(mut writer) = writers.remove(&key) {
            writer.flush()?; // 刷新并清理
        }
    }

    Ok(())
}

fn process_batch<R>(
    reader: &mut R,
    hash_config: &HashConfig,
    page: &Page,
    chunk_dir: PathBuf,
    buffer_size: usize,
    bin_threads: u32,
    // page_index: usize,
    num_threads: usize,
) -> std::io::Result<()>
where
    R: Read + Send,
{
    let row_size = std::mem::size_of::<Row>();
    let mut writers: HashMap<(u64, u32), BufWriter<File>> = HashMap::new();
    let mut current_file_index: Option<u64> = None;

    let value_mask = hash_config.get_value_mask();
    let value_bits = hash_config.get_value_bits();
    let idx_mask = hash_config.get_idx_mask();
    let idx_bits = hash_config.get_idx_bits();

    buffer_read_parallel(
        reader,
        num_threads,
        buffer_size,
        |dataset: Vec<Slot<u64>>| {
            let mut results: HashMap<(u64, u32), Vec<u8>> = HashMap::new();
            for slot in dataset {
                let indx = slot.idx & idx_mask;
                let compacted = slot.value.left(value_bits) as u32;
                // let taxid = chtm.get_from_page(indx, compacted, page_index);
                let taxid = page.find_index(indx, compacted, value_bits, value_mask);

                if taxid > 0 {
                    let kmer_id = slot.idx >> idx_bits;
                    let file_index = slot.value.right(value_mask) >> 32;
                    let seq_id = slot.get_seq_id() as u32;
                    let left = slot.value.left(value_bits) as u32;
                    let high = u32::combined(left, taxid, value_bits);
                    let row = Row::new(high, seq_id, kmer_id as u32);
                    let value_bytes = row.as_slice(row_size);
                    let seq_id_mod = seq_id % bin_threads;

                    results
                        .entry((file_index, seq_id_mod))
                        .or_insert_with(Vec::new)
                        .extend(value_bytes);
                }
            }
            results
        },
        |result| {
            while let Some(data) = result.next() {
                let res = data.unwrap();
                let mut file_keys: Vec<_> = res.keys().cloned().collect();
                file_keys.sort_unstable(); // 对 (file_index, seq_id_mod) 进行排序

                for (file_index, seq_id_mod) in file_keys {
                    if let Some(bytes) = res.get(&(file_index, seq_id_mod)) {
                        // 如果当前处理的 file_index 改变了，清理非当前的 writers
                        if current_file_index != Some(file_index) {
                            clean_up_writers(&mut writers, file_index).expect("clean writer");
                            current_file_index = Some(file_index);
                        }

                        write_to_file(file_index, seq_id_mod, bytes, &mut writers, &chunk_dir)
                            .expect("write to file error");
                    }
                }
            }
        },
    )
    .expect("failed");

    // 最终批次处理完成后，刷新所有的 writer
    for writer in writers.values_mut() {
        writer.flush()?;
    }

    Ok(())
}

fn process_chunk_file<P: AsRef<Path>>(
    args: &Args,
    chunk_file: P,
    hash_files: &Vec<PathBuf>,
    large_page: &mut Page,
) -> Result<()> {
    let file = open_file(chunk_file)?;
    let mut reader = BufReader::new(file);

    let (page_index, _) = read_chunk_header(&mut reader)?;

    let start = Instant::now();

    println!("start load table...");
    let config = HashConfig::from_hash_header(&args.database.join("hash_config.k2d"))?;

    read_next_page(large_page, hash_files, page_index, config)?;
    // 计算持续时间
    let duration = start.elapsed();
    // 打印运行时间
    println!("load table took: {:?}", duration);
    process_batch(
        &mut reader,
        &config,
        &large_page,
        args.chunk_dir.clone(),
        args.buffer_size,
        args.batch_size,
        // page_index,
        args.num_threads,
    )?;

    Ok(())
}

pub fn run(args: Args) -> Result<()> {
    let chunk_files = find_and_sort_files(&args.chunk_dir, "sample", ".k2", true)?;
    let hash_files = find_and_sort_files(
        &args.database, "hash", ".k2d", true,
    )
    .expect("Invalid or incomplete database: missing hash files.");

    // 开始计时
    let start = Instant::now();
    println!("annotate start...");
    let config = HashConfig::from_hash_header(&args.database.join("hash_config.k2d"))
        .expect("Invalid or incomplete database: missing hash_config.k2d.");
    let mut large_page = Page::with_capacity(0, config.hash_capacity);
    for chunk_file in &chunk_files {
        process_chunk_file(&args, chunk_file, &hash_files, &mut large_page)?;
        let _ = std::fs::remove_file(chunk_file);
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
