use clap::Parser;
use flate2::bufread::MultiGzDecoder; // 支持 .gz 和 .fna
use kun_peng::args::parse_size;
use kun_peng::db::generate_taxonomy;
use kun_peng::utils::{find_files, read_id_to_taxon_map};
use rayon::prelude::*;
use regex::Regex;
use std::collections::{HashMap, HashSet}; 
use std::fs::{create_dir_all, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Read, Result, Write}; 
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Instant, SystemTime, UNIX_EPOCH}; 
use walkdir::WalkDir;

use md5::Context; 

#[derive(Parser, Debug, Clone)]
#[clap(version, about = "Add new FASTA files to an existing Kun-Peng database library")]
pub struct Args {
    /// Main database directory (must contain existing library/ and taxonomy/ dirs)
    #[arg(long = "db", required = true)]
    pub database: PathBuf,

    /// Input files or directories (containing .fa, .fna, .fasta, .fsa, *.gz files)
    #[arg(short = 'i', long, required = true, num_args = 1..)]
    pub input_library: Vec<PathBuf>,

    /// library fna temp file max size
    #[arg(long = "max-file-size", value_parser = parse_size, default_value = "2G")]
    pub max_file_size: usize,
}

// ... (SizedWriter 结构体保持不变) ...
struct SizedWriter {
    writer: BufWriter<File>,
    bytes_written: u64,
    thread_index: usize,
    file_suffix: AtomicUsize,
    library_dir: PathBuf,
    max_file_size: u64,
    file_prefix: String,
}

impl SizedWriter {
    fn new(
        library_dir: &PathBuf,
        file_prefix: String,
        thread_index: usize,
        max_file_size: u64,
    ) -> Result<Self> {
        let file_suffix = AtomicUsize::new(0);
        let path = Self::get_file_path(library_dir, &file_prefix, thread_index, 0);
        let file = OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(true) // 这是安全的，因为 file_prefix 是唯一的
            .open(&path)?;
        let writer = BufWriter::new(file);

        Ok(Self {
            writer,
            bytes_written: 0,
            thread_index,
            file_suffix,
            library_dir: library_dir.to_path_buf(),
            max_file_size,
            file_prefix,
        })
    }

    fn get_file_path(
        library_dir: &PathBuf,
        prefix: &str,
        thread_index: usize,
        suffix: usize,
    ) -> PathBuf {
        library_dir.join(format!("{}_{}_{}.fna", prefix, thread_index, suffix))
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
            let new_path = Self::get_file_path(
                &self.library_dir,
                &self.file_prefix,
                self.thread_index,
                new_suffix,
            );
            let new_file = OpenOptions::new()
                .create(true)
                .write(true)
                .truncate(true) 
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


// ... (find_fasta_files 保持不变) ...
/// 查找所有 FASTA 文件 (包括 .gz)，递归地
fn find_fasta_files(paths: &[PathBuf]) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    let extensions = ["fa", "fna", "fasta", "fsa"];

    for path in paths {
        if path.is_dir() {
            for entry in WalkDir::new(path).into_iter().filter_map(|e| e.ok()) {
                if entry.file_type().is_file() {
                    let p = entry.path();
                    let ext = p.extension().and_then(|s| s.to_str()).unwrap_or("");
                    if extensions.contains(&ext) {
                        files.push(p.to_path_buf());
                    } else if ext == "gz" {
                        if let Some(stem) = p.file_stem() {
                            if let Some(inner_ext) =
                                Path::new(stem).extension().and_then(|s| s.to_str())
                            {
                                if extensions.contains(&inner_ext) {
                                    files.push(p.to_path_buf());
                                }
                            }
                        }
                    }
                }
            }
        } else if path.is_file() {
            files.push(path.clone());
        }
    }
    Ok(files)
}

// ... (parse_header_to_map_entry 保持不变, 它正确地返回 None) ...
/// 解析 FASTA 标题以提取 seqid 和 taxid
/// 稳健地查找 "taxid|TAXID" 或 "kraken:taxid|TAXID" 这样的模式
/// 输出: "full_seq_id\tTAXID"
fn parse_header_to_map_entry(header: &str) -> Option<String> {
    // 1. 去除 '>' 并获取第一个 "单词" (ID 部分)
    let id_part = header
        .strip_prefix('>')
        .unwrap_or(header)
        .split_whitespace()
        .next()
        .unwrap_or("");

    if id_part.is_empty() {
        return None;
    }

    // 2. 按 '|' 分割
    let parts: Vec<&str> = id_part.split('|').collect();

    // 3. 迭代查找 "taxid" 键值对
    for i in 0..parts.len().saturating_sub(1) {
        let key = parts[i];
        let value = parts[i + 1];

        // 检查键是否包含 "taxid"
        if key.contains("taxid") {
            // 检查值是否为非空且纯数字
            if !value.is_empty() && value.chars().all(char::is_numeric) {
                // 找到了！返回 "完整的ID部分\t分类ID"
                return Some(format!("{}\t{}", id_part, value));
            }
        }
    }
    None
}

// --- 已修改 ---
/// 处理单个 FASTA 文件 (gz 或 plain)
fn process_fasta_file(
    fasta_file: &PathBuf,
    map_writer: &mut BufWriter<File>,
    fna_writer: &mut SizedWriter,
    fna_start: &Regex,
) -> Result<()> { // <-- 这个 Result 可以是 Box<dyn Error>
    let file = File::open(fasta_file)?;
    let is_gzipped = fasta_file.extension().and_then(|s| s.to_str()) == Some("gz");
    
    let reader: Box<dyn BufRead> = if is_gzipped {
        Box::new(BufReader::new(MultiGzDecoder::new(BufReader::new(file))))
    } else {
        Box::new(BufReader::new(file))
    };
    let mut reader = BufReader::new(reader);

    let mut line = String::new();
    let mut fna_buffer = String::new();

    while reader.read_line(&mut line)? != 0 {
        if fna_start.is_match(&line) {
            // 找到了一个新的 FASTA 记录头 (>)
            if !fna_buffer.is_empty() {
                fna_writer.write(fna_buffer.as_bytes())?;
                fna_buffer.clear();
            }

            // --- 这是新的错误处理逻辑 ---
            if let Some(map_entry) = parse_header_to_map_entry(&line) {
                // 成功: 写入 map, 准备 fna_buffer
                map_writer.write_all(map_entry.as_bytes())?;
                map_writer.write_all(b"\n")?;
                fna_buffer.push_str(&line);
            } else {
                // 失败: 构造错误消息并返回 Err
                let error_message = format!(
                    "Error in file '{}': Could not parse a valid taxid from FASTA header. \
                     \nPlease ensure the header contains a 'taxid|123' format.\
                     \nProblematic header: \"{}\"",
                    fasta_file.display(),
                    line.trim()
                );
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    error_message,
                ));
            }
            // --- 结束新逻辑 ---
        } else {
            // 这是序列行
            if !fna_buffer.is_empty() {
                fna_buffer.push_str(&line);
            }
        }
        line.clear();
    }

    if !fna_buffer.is_empty() {
        fna_writer.write(fna_buffer.as_bytes())?;
    }

    fna_writer.flush()?;
    map_writer.flush()?;
    Ok(())
}

// ... (merge_files 保持不变) ...
/// 复制自 merge_fna.rs
fn merge_files(paths: &Vec<PathBuf>, output_path: &PathBuf) -> Result<()> {
    let output = Arc::new(Mutex::new(BufWriter::new(File::create(output_path)?)));
    let buffer_size = 64 * 1024 * 1024;

    paths.par_iter().try_for_each(|path| -> Result<()> {
        let input = File::open(path)?;
        let mut reader = BufReader::new(input);
        let mut buffer = vec![0; buffer_size];

        loop {
            let bytes_read = reader.read(&mut buffer)?;
            if bytes_read == 0 {
                break;
            }
            let mut output = output.lock().unwrap();
            output.write_all(&buffer[..bytes_read])?;
        }
        std::fs::remove_file(path)?;
        Ok(())
    })?;

    let mut output = output.lock().unwrap();
    output.flush()?;
    Ok(())
}


// ... (load_processed_log 和 hash_file_content 保持不变) ...
/// 从 added.md5 日志文件加载已处理文件的哈希值
fn load_processed_log(log_path: &Path) -> Result<HashSet<String>> {
    if !log_path.exists() {
        return Ok(HashSet::new());
    }
    let file = File::open(log_path)?;
    let reader = BufReader::new(file);
    let mut hashes = HashSet::new();
    for line in reader.lines() {
        if let Ok(line_content) = line {
            if let Some(hash) = line_content.split('\t').nth(1) {
                hashes.insert(hash.to_string());
            }
        }
    }
    Ok(hashes)
}

/// 计算文件的 md5 哈希值
fn hash_file_content(file_path: &Path) -> Result<String> {
    let mut file = File::open(file_path)?;
    
    // 使用 'Context::new()'
    let mut hasher = Context::new(); 
    std::io::copy(&mut file, &mut hasher)?;
    let hash_bytes = hasher.finalize(); 
    Ok(format!("{:x}", hash_bytes)) 
}


// --- 已修改 ---
/// 并行添加 FASTA 文件，基于 merge_fna.rs::merge_fna_parallel
fn add_fna_parallel(
    fasta_files: &Vec<PathBuf>, // <-- 列表现在是预先过滤过的
    database: &PathBuf,
    library_dir: &PathBuf,
    max_file_size: u64,
    run_prefix: String, 
) -> Result<()> { // <-- 这个 Result 会从 try_for_each 传播上来
    let fna_start: Regex = Regex::new(r"^>").unwrap(); 
    let writers: Arc<Mutex<HashMap<usize, SizedWriter>>> = Arc::new(Mutex::new(HashMap::new()));

    // --- 改为 .try_for_each() 以便能中途退出 ---
    let result = fasta_files.par_iter().try_for_each(|fasta_file| -> Result<()> {
        let thread_index = rayon::current_thread_index().unwrap_or(0);
        let mut writers = writers.lock().unwrap();
        
        let run_prefix_clone = run_prefix.clone();

        let fna_writer = writers.entry(thread_index).or_insert_with(|| {
            SizedWriter::new(
                &library_dir,
                run_prefix_clone, 
                thread_index,
                max_file_size
            ).unwrap()
        });

        let seqid2taxid_path =
            database.join(format!("add_seqid2taxid_{}.map", thread_index));
        let mut map_writer = BufWriter::new(
            OpenOptions::new()
                .create(true)
                .write(true)
                .append(true) 
                .open(&seqid2taxid_path)
                .unwrap(),
        );

        // --- '?' 将在出错时立即传播 Err, 停止 .try_for_each ---
        process_fasta_file(&fasta_file, &mut map_writer, fna_writer, &fna_start)?;

        Ok(()) // <-- 此文件成功
    });

    // 检查 .try_for_each 的结果
    if let Err(e) = result {
        return Err(e);
    }

    // 确保所有 SizedWriter 都已 flush
    let mut writers = writers.lock().unwrap();
    for (_, writer) in writers.iter_mut() {
        writer.flush()?;
    }
    Ok(())
}

pub fn run(args: Args) -> Result<()> {
    let start = Instant::now();
    println!("Adding files to library...");
    let database = &args.database;
    let max_file_size = &args.max_file_size;

    // 1. 准备目录
    let dst_tax_dir = database.join("taxonomy");
    create_dir_all(&dst_tax_dir)?;
    let library_dir = database.join("library");
    create_dir_all(&library_dir)?; 

    let names_file = dst_tax_dir.join("names.dmp");
    let nodes_file = dst_tax_dir.join("nodes.dmp");
    assert!(names_file.exists(), "names.dmp not found in taxonomy directory");
    assert!(nodes_file.exists(), "nodes.dmp not found in taxonomy directory");

    // 3. 哈希校验和文件过滤
    let log_path = library_dir.join("added.md5");
    println!("Loading processed file log from: {}", log_path.display());

    let processed_hashes_arc = Arc::new(load_processed_log(&log_path)?);
    println!("Loaded {} previously processed file hashes.", processed_hashes_arc.len());

    let new_hashes_for_this_run_arc = Arc::new(Mutex::new(HashSet::new()));

    let all_fasta_files = find_fasta_files(&args.input_library)?;
    println!("Found {} total FASTA files. Checking for duplicates...", all_fasta_files.len());

    // (文件路径, 哈希值) 的 Vec
    let files_to_process_with_hash: Vec<(PathBuf, String)> = all_fasta_files
        .into_par_iter() 
        .filter_map(|file_path| {
            let hash = match hash_file_content(&file_path) {
                Ok(h) => h,
                Err(e) => {
                    eprintln!("Error hashing file {}: {}. Skipping.", file_path.display(), e);
                    return None;
                }
            };

            // 检查是否 *之前* 处理过
            if processed_hashes_arc.contains(&hash) {
                println!("Skipping (already processed): {}", file_path.display());
                return None;
            }
            
            // 检查是否在 *本轮* 运行中已经见过 (处理重复输入)
            let mut new_hashes = new_hashes_for_this_run_arc.lock().unwrap();
            if new_hashes.contains(&hash) {
                println!("Skipping (duplicate in this run): {}", file_path.display());
                return None;
            }

            // 这是个新文件
            new_hashes.insert(hash.clone()); // 标记为本轮已见
            Some((file_path, hash)) // <-- 保留路径和哈希
        })
        .collect();
    
    // 从元组 Vec 中提取文件路径
    let files_to_process: Vec<PathBuf> = files_to_process_with_hash
        .iter()
        .map(|(path, _hash)| path.clone())
        .collect();
    
    // 4. 检查是否有新文件需要处理
    if files_to_process.is_empty() {
        println!("No new files to add. Database is up-to-date.");
        return Ok(());
    }
    println!("Processing {} new files...", files_to_process.len());


    // 5. 生成唯一的运行前缀
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("Time went backwards")
        .as_secs();
    let run_prefix = format!("library_add_{}", timestamp);
    println!("Using unique run prefix: {}", run_prefix);

    // 6. 传递 *过滤后* 的列表到并行处理器
    // --- '?' 将捕获来自 'add_fna_parallel' 的任何错误并停止 'run' ---
    add_fna_parallel(
        &files_to_process, // <-- 使用过滤后的列表
        &args.database,
        &library_dir,
        *max_file_size as u64,
        run_prefix,
    )?; 

    // 7. 合并并追加 map 文件
    let add_map_files = find_files(database, "add_seqid2taxid_", "map");
    if !add_map_files.is_empty() {
        let temp_map_merge_path = database.join("add_seqid2taxid.map.tmp");
        
        // a. 将所有 add_*.map 文件合并成一个临时文件
        merge_files(&add_map_files, &temp_map_merge_path)?;
        
        // b. 加载现有的主 map 文件
        let main_map_path = database.join("seqid2taxid.map");
        let mut existing_map = if main_map_path.exists() {
            // 使用我们之前讨论过的 read_id_to_taxon_map 函数
            read_id_to_taxon_map(&main_map_path)?
        } else {
            HashMap::new()
        };
        
        println!("De-duplicating and appending new seqid entries...");

        // c. 打开主 map 文件用于追加
        let mut main_map_writer = BufWriter::new(OpenOptions::new()
            .create(true) 
            .append(true)
            .write(true)
            .open(&main_map_path)?);
        
        // d. 逐行读取临时文件并进行校验
        let temp_map_file = File::open(&temp_map_merge_path)?;
        let reader = BufReader::new(temp_map_file);

        for (line_number, line_result) in reader.lines().enumerate() {
            let line = line_result?;
            let parts: Vec<&str> = line.trim().split_whitespace().collect();
            if parts.len() < 2 { continue; }
            let seq_id = parts[0].to_string();
            if let Ok(taxid) = parts[1].parse::<u64>() {
                
                match existing_map.get(&seq_id) {
                    Some(old_taxid) => {
                        // --- 情况 A (冲突) ---
                        if *old_taxid != taxid {
                            let error_message = format!(
                                "Error: Inconsistent mapping for sequence ID '{}'. \
                                 The main database maps it to '{}', but a new file maps it to '{}'. \
                                 (Source: '{}', Line: {})",
                                seq_id, old_taxid, taxid, temp_map_merge_path.display(), line_number + 1
                            );
                            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, error_message));
                        }
                    }
                    None => {
                        // 'insert' 返回 None 意味着它是真的新条目 (既不在旧 map 也不在本轮)
                        if existing_map.insert(seq_id.clone(), taxid).is_none() {
                             // 写入主 map 文件
                            writeln!(main_map_writer, "{}\t{}", seq_id, taxid)?;
                        }
                    }
                }
            }
        }
        main_map_writer.flush()?;
        
        // e. 清理临时文件
        std::fs::remove_file(&temp_map_merge_path)?;
        println!("Appended new entries to seqid2taxid.map");
    }

    // 8. 重新生成分类法
    println!("Regenerating taxonomy (taxo.k2d)...");
    let id_to_taxon_map_filename = args.database.join("seqid2taxid.map");
    let id_to_taxon_map = read_id_to_taxon_map(&id_to_taxon_map_filename)?;
    
    let taxonomy_filename = args.database.join("taxo.k2d");
    let ncbi_taxonomy_directory = &dst_tax_dir;

    let _ = generate_taxonomy(
        &ncbi_taxonomy_directory,
        &taxonomy_filename,
        &id_to_taxon_map,
    )?;

    // 9. 成功后，将新哈希值写入日志 (TSV 格式)
    // --- 只有在前面所有步骤都成功后 (没有 Err) 才会执行到这里 ---
    println!("Updating processed file log: {}", log_path.display());
    let log_file = OpenOptions::new()
        .append(true)
        .create(true)
        .write(true)
        .open(&log_path)?;
    let mut log_writer = BufWriter::new(log_file);

    // 迭代 (路径, 哈希) 元组
    for (path, hash) in &files_to_process_with_hash {
        writeln!(log_writer, "{}\t{}", path.display(), hash)?;
    }
    log_writer.flush()?;

    let duration = start.elapsed();
    println!("Finished adding {} new files in: {:?}", files_to_process.len(), duration);
    Ok(())
}

#[allow(dead_code)]
fn main() {
    let args = Args::parse();
    if let Err(e) = run(args) {
        eprintln!("Application error: {}", e);
    }
}