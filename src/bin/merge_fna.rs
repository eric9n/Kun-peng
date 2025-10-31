use clap::Parser;
use flate2::read::GzDecoder;
use kun_peng::args::parse_size;
use kun_peng::utils::{find_files, open_file};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::{create_dir_all, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Read, Result, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

#[derive(Parser, Debug, Clone)]
#[clap(version, about = "A tool for processing genomic files")]
pub struct Args {
    /// Directory to store downloaded files
    #[arg(short, long, required = true)]
    pub download_dir: PathBuf,

    /// ncbi library fna database directory
    #[arg(long = "db", required = true)]
    pub database: PathBuf,

    /// library fna temp file max size
    #[arg(long = "max-file-size", value_parser = parse_size, default_value = "2G")]
    pub max_file_size: usize,
}

struct SizedWriter {
    writer: BufWriter<File>,
    bytes_written: u64,
    thread_index: usize,
    file_suffix: AtomicUsize,
    library_dir: PathBuf,
    max_file_size: u64,
}

impl SizedWriter {
    fn new(library_dir: &PathBuf, thread_index: usize, max_file_size: u64) -> Result<Self> {
        let file_suffix = AtomicUsize::new(0);
        let path = Self::get_file_path(library_dir, thread_index, 0);
        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .write(true)
            .open(&path)?;
        let writer = BufWriter::new(file);

        Ok(Self {
            writer,
            bytes_written: 0,
            thread_index,
            file_suffix,
            library_dir: library_dir.to_path_buf(),
            max_file_size,
        })
    }

    fn get_file_path(library_dir: &PathBuf, thread_index: usize, suffix: usize) -> PathBuf {
        library_dir.join(format!("library_{}_{}.fna", thread_index, suffix))
    }

    fn is_kraken_taxid_start(buf: &[u8]) -> bool {
        buf.starts_with(b">taxid")
    }

    fn write(&mut self, buf: &[u8]) -> Result<usize> {
        if Self::is_kraken_taxid_start(buf)
            && self.bytes_written + buf.len() as u64 > self.max_file_size
        {
            self.writer.flush()?;
            let new_suffix = self.file_suffix.fetch_add(1, Ordering::SeqCst) + 1;
            let new_path = Self::get_file_path(&self.library_dir, self.thread_index, new_suffix);
            let new_file = OpenOptions::new()
                .create(true)
                .append(true)
                .write(true)
                .open(&new_path)?;
            self.writer = BufWriter::new(new_file);
            self.bytes_written = 0;
        }

        let bytes_written = self.writer.write(buf)?;
        self.bytes_written += bytes_written as u64;
        Ok(bytes_written)
    }

    fn flush(&mut self) -> Result<()> {
        self.writer.flush()
    }
}

fn parse_assembly_fna(assembly_file: &PathBuf, site: &str) -> Result<Vec<(String, String)>> {
    let mut gz_files = Vec::new();
    let file = open_file(&assembly_file)?;
    let reader = BufReader::new(file);
    let lines = reader.lines();

    let parent_path = assembly_file
        .parent()
        .expect("Can't find assembly file parent directory");
    for line in lines {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() > 19 {
            let (taxid, _, ftp_path) = (fields[5], fields[11], fields[19]);

            if ftp_path == "na" {
                continue;
            }

            // let levels = vec!["Complete Genome", "Chromosome"];
            // if !levels.contains(&asm_level) {
            //     continue;
            // }

            let fna_file_name = format!(
                "{}/{}/{}_genomic.fna.gz",
                parent_path.to_string_lossy(),
                site,
                ftp_path.split('/').last().unwrap_or_default()
            );
            gz_files.push((fna_file_name, taxid.into()));
        }
    }
    Ok(gz_files)
}

fn process_gz_file(
    gz_file: &PathBuf,
    map_writer: &mut BufWriter<File>,
    fna_writer: &mut SizedWriter,
    fna_start: &regex::Regex,
    taxid: &str,
) -> Result<()> {
    let file = open_file(gz_file)?;
    let decompressor = GzDecoder::new(BufReader::new(file));
    let mut reader = BufReader::new(decompressor);

    let mut line = String::new();
    
    // fna_buffer 用于存储一个 *完整* 的 FASTA 记录
    let mut fna_buffer = String::new(); 

    while reader.read_line(&mut line)? != 0 {
        if let Some(caps) = fna_start.captures(&line) {
            // 找到了一个新的 FASTA 记录头 (>)

            // 1. 写入 *上一个* 完整的 FASTA 记录 (如果它存在)
            // SizedWriter 现在会接收一个以 ">taxid" 开头的完整记录，
            // 它的文件分割逻辑可以正确运行了。
            if !fna_buffer.is_empty() {
                fna_writer.write(fna_buffer.as_bytes())?;
                fna_buffer.clear();
            }

            // 2. 为 *新* 记录处理 map 和 fna 头
            let seqid = &caps[1];
            
            // 直接写入 map_writer (它本身有缓冲)
            map_writer.write_all(
                format!("taxid|{}|{}\t{}\n", taxid, seqid, taxid).as_bytes()
            )?;

            // 开始在 fna_buffer 中累积 *新* 的 FASTA 记录
            fna_buffer.push_str(&format!(">taxid|{}|{}", taxid, &line[1..]));
        } else {
            // This is a sequence line, append it to the current fna_buffer
            fna_buffer.push_str(&line);
        }

        line.clear();
    }

    // After the loop ends, don't forget to write out the last accumulated FASTA record
    if !fna_buffer.is_empty() {
        fna_writer.write(fna_buffer.as_bytes())?;
    }

    // Flush once at the end of the function
    fna_writer.flush()?;
    map_writer.flush()?;

    Ok(())
}

const PREFIX: &'static str = "assembly_summary";
const SUFFIX: &'static str = "txt";

fn merge_fna_parallel(
    assembly_files: &Vec<PathBuf>,
    database: &PathBuf,
    library_dir: &PathBuf,
    max_file_size: u64,
) -> Result<()> {
    let pattern = format!(r"{}_(\S+)\.{}", PREFIX, SUFFIX);
    let file_site = regex::Regex::new(&pattern).unwrap();

    let fna_start: regex::Regex = regex::Regex::new(r"^>(\S+)").unwrap();
    let is_empty = AtomicBool::new(true);
    let writers: Arc<Mutex<HashMap<usize, SizedWriter>>> = Arc::new(Mutex::new(HashMap::new()));

    for assembly_file in assembly_files {
        if let Some(caps) = file_site.captures(assembly_file.to_string_lossy().as_ref()) {
            if let Some(matched) = caps.get(1) {
                let gz_files = parse_assembly_fna(assembly_file, matched.as_str())?;

                gz_files.par_iter().for_each(|(gz_path, taxid)| {
                    let gz_file = PathBuf::from(&gz_path);
                    if !gz_file.exists() {
                        // eprintln!("{} does not exist", gz_file.to_string_lossy());
                        return;
                    }

                    let thread_index = rayon::current_thread_index().unwrap_or(0);
                    let mut writers = writers.lock().unwrap();
                    let mut fna_writer = writers.entry(thread_index).or_insert_with(|| {
                        SizedWriter::new(&library_dir, thread_index, max_file_size).unwrap()
                    });
                    let seqid2taxid_path =
                        database.join(format!("seqid2taxid_{}.map", thread_index));
                    let mut map_writer = BufWriter::new(
                        OpenOptions::new()
                            .create(true)
                            .write(true)
                            .append(true)
                            .open(&seqid2taxid_path)
                            .unwrap(),
                    );

                    if let Err(e) = process_gz_file(
                        &gz_file,
                        &mut map_writer,
                        &mut fna_writer,
                        &fna_start,
                        &taxid,
                    ) {
                        eprintln!("process_gz_file error: {}", e);
                    } else {
                        fna_writer.flush().unwrap();
                        map_writer.flush().unwrap();
                        is_empty.fetch_and(false, Ordering::Relaxed);
                    }
                });
            }
        }
    }

    let seqid_files = find_files(database, "seqid2taxid_", "map");
    let seqid2taxid_path = database.join("seqid2taxid.map");
    merge_files(&seqid_files, &seqid2taxid_path)?;
    if is_empty.load(Ordering::Relaxed) {
        panic!("genimics fna files is empty! please check download dir");
    }
    Ok(())
}

fn merge_files(paths: &Vec<PathBuf>, output_path: &PathBuf) -> Result<()> {
    let output = Arc::new(Mutex::new(BufWriter::new(File::create(output_path)?)));
    let buffer_size = 64 * 1024 * 1024; // Increased buffer size for better performance

    paths.par_iter().try_for_each(|path| -> Result<()> {
        let input = File::open(path)?;
        let mut reader = BufReader::new(input);
        let mut buffer = vec![0; buffer_size];

        loop {
            let bytes_read = reader.read(&mut buffer)?;
            if bytes_read == 0 {
                break; // 文件读取完毕
            }

            let mut output = output.lock().unwrap();
            output.write_all(&buffer[..bytes_read])?;
        }

        std::fs::remove_file(path)?;
        Ok(())
    })?;

    // Ensure all buffered data is flushed to the output file
    let mut output = output.lock().unwrap();
    output.flush()?;
    Ok(())
}
pub fn run(args: Args) -> Result<()> {
    // 开始计时
    let start = Instant::now();
    println!("merge fna start...");
    let download_dir = args.download_dir;
    let database = &args.database;
    let max_file_size = &args.max_file_size;

    let dst_tax_dir = database.join("taxonomy");
    create_dir_all(&dst_tax_dir)?;

    let library_dir = database.join("library");
    // create_dir_all(&library_dir)?;

    let source_names_file = &download_dir.join("taxonomy").join("names.dmp");
    assert!(source_names_file.exists());
    let dst_name_file = &dst_tax_dir.join("names.dmp");
    if !dst_name_file.exists() {
        std::fs::copy(source_names_file, dst_name_file)?;
    }

    let source_nodes_file = &download_dir.join("taxonomy").join("nodes.dmp");
    assert!(source_nodes_file.exists());
    let dst_nodes_file = &dst_tax_dir.join("nodes.dmp");
    if !dst_nodes_file.exists() {
        std::fs::copy(source_nodes_file, dst_nodes_file)?;
    }

    let seqid2taxid_path = database.join("seqid2taxid.map");
    if seqid2taxid_path.exists() {
        if let Ok(mut entries) = std::fs::read_dir(&library_dir) {
            if entries.next().is_some() {
                // 如果 library 目录至少有一个文件，我们就认为构建已完成
                println!("Build appears complete (seqid2taxid.map, taxo.k2d, and library files exist). Skipping.");
                return Ok(());
            }
        }
    }

    if library_dir.exists() {
        std::fs::remove_dir_all(&library_dir)?;
    }

    create_dir_all(&library_dir)?; // 重新创建干净的目录
    if seqid2taxid_path.exists() {
        std::fs::remove_file(seqid2taxid_path)?;
    }
    let assembly_files = find_files(&download_dir, &PREFIX, &SUFFIX);

    merge_fna_parallel(
        &assembly_files,
        &args.database,
        &library_dir,
        *max_file_size as u64,
    )?;

    // 计算持续时间
    let duration = start.elapsed();
    println!("merge fna took: {:?}", duration);
    Ok(())
}

#[allow(dead_code)]
fn main() {
    let args = Args::parse();
    if let Err(e) = run(args) {
        eprintln!("Application error: {}", e);
    }
}
