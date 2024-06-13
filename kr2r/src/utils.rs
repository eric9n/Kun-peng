use std::collections::HashMap;
use std::fs::{self, create_dir_all, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Result, Seek, Write};
use std::path::{Path, PathBuf};
use walkdir::WalkDir;

/// 读取 seqid2taxid.map 文件。为了裁剪 ncbi 的 taxonomy 树
pub fn read_id_to_taxon_map<P: AsRef<Path>>(filename: P) -> Result<HashMap<String, u64>> {
    let file = open_file(filename)?;
    let reader = BufReader::new(file);
    let mut id_map = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.trim().split_whitespace().collect();
        if parts.len() < 2 {
            continue;
        }
        let seq_id = parts[0].to_string();
        if let Ok(taxid) = parts[1].parse::<u64>() {
            id_map.insert(seq_id, taxid);
        }
    }

    Ok(id_map)
}

/// Expands a spaced seed mask based on the given bit expansion factor.
///
/// # Examples
///
/// Basic usage:
///
/// ```
/// # use kr2r::utils::expand_spaced_seed_mask; // Replace with the appropriate crate name
/// // Expanding 0b1010 (binary for 10) with a factor of 2
/// assert_eq!(expand_spaced_seed_mask(0b1010, 2), 204);
///
/// // Expanding 0b0101 (binary for 5) with a factor of 1
/// assert_eq!(expand_spaced_seed_mask(0b0101, 1), 5);
/// ```
///
/// When the bit expansion factor is zero or greater than 64:
///
/// ```
/// # use kr2r::utils::expand_spaced_seed_mask;
/// // No expansion, factor is 0
/// assert_eq!(expand_spaced_seed_mask(0b1010, 0), 0b1010);
///
/// // No expansion, factor is greater than 64
/// assert_eq!(expand_spaced_seed_mask(0b1010, 65), 0b1010);
/// ```
pub fn expand_spaced_seed_mask(spaced_seed_mask: u64, bit_expansion_factor: u64) -> u64 {
    // 检查 bit_expansion_factor 是否在有效范围内
    if bit_expansion_factor == 0 || bit_expansion_factor > 64 {
        return spaced_seed_mask;
    }

    let mut new_mask = 0;
    let bits = (1 << bit_expansion_factor) - 1;

    for i in (0..64 / bit_expansion_factor).rev() {
        new_mask <<= bit_expansion_factor;
        if (spaced_seed_mask >> i) & 1 == 1 {
            new_mask |= bits;
        }
    }

    new_mask
}

pub fn find_library_fna_files<P: AsRef<Path>>(path: P) -> Vec<String> {
    WalkDir::new(path)
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().file_name() == Some("library.fna".as_ref()))
        .map(|e| e.path().to_string_lossy().into_owned())
        .collect()
}

pub fn summary_prelim_map_files<P: AsRef<Path>>(data_dir: P) -> Result<PathBuf> {
    let lib_path = data_dir.as_ref().join("library");

    let prelim_map_files: Vec<PathBuf> = WalkDir::new(&lib_path)
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().file_name() == Some("prelim_map.txt".as_ref()))
        .map(|e| e.path().to_path_buf())
        .collect();

    let summary_file_path = data_dir.as_ref().join("prelim_map.txt");
    if summary_file_path.exists() && summary_file_path.is_file() {
        // 删除文件
        fs::remove_file(&summary_file_path)?;
    }

    let mut summary_file = File::create(&summary_file_path)?;
    // 遍历找到的文件并汇总内容
    for path in prelim_map_files {
        let file = open_file(path)?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            writeln!(summary_file, "{}", line)?;
        }
    }

    Ok(summary_file_path)
}

pub fn create_seqid2taxid_file<P: AsRef<Path>>(prelim_map_file: P, output_file: P) -> Result<()> {
    let file = open_file(prelim_map_file)?;
    let reader = BufReader::new(file);
    let mut output = File::create(output_file).unwrap();

    for line in reader.lines() {
        let line = line?;
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() >= 3 {
            // 写入第二列和第三列到输出文件
            writeln!(output, "{}\t{}", columns[1], columns[2])?;
        }
    }

    Ok(())
}

pub fn format_bytes(size: f64) -> String {
    let suffixes = ["B", "KB", "MB", "GB", "TB", "PB", "EB"];
    let mut size = size;
    let mut current_suffix = &suffixes[0];

    for suffix in &suffixes[1..] {
        if size >= 1024.0 {
            current_suffix = suffix;
            size /= 1024.0;
        } else {
            break;
        }
    }

    format!("{:.2}{}", size, current_suffix)
}

#[derive(Debug)]
pub enum FileFormat {
    Fasta,
    Fastq,
}

use flate2::read::GzDecoder;
use std::io::{self, Read};

pub fn is_gzipped(file: &mut File) -> io::Result<bool> {
    let mut buffer = [0; 2];
    file.read_exact(&mut buffer)?;
    file.rewind()?; // 重置文件指针到开头
    Ok(buffer == [0x1F, 0x8B])
}

pub fn detect_file_format<P: AsRef<Path>>(path: P) -> io::Result<FileFormat> {
    let mut file = open_file(path)?;
    let read1: Box<dyn io::Read + Send> = if is_gzipped(&mut file)? {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };

    let reader = BufReader::new(read1);
    let mut lines = reader.lines();

    if let Some(first_line) = lines.next() {
        let line = first_line?;

        if line.starts_with('>') {
            return Ok(FileFormat::Fasta);
        } else if line.starts_with('@') {
            let _ = lines.next();
            if let Some(third_line) = lines.next() {
                let line: String = third_line?;
                if line.starts_with('+') {
                    return Ok(FileFormat::Fastq);
                }
            }
        } else {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "Unrecognized fasta(fastq) file format",
            ));
        }
    }

    Err(io::Error::new(
        io::ErrorKind::Other,
        "Unrecognized fasta(fastq) file format",
    ))
    // let mut buffer = [0; 1]; // 仅分配一个字节的缓冲区

    // // 读取文件的第一个字节
    // let bytes_read = reader.read(&mut buffer)?;

    // if bytes_read == 0 {
    //     return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "Empty file"));
    // }

    // match buffer[0] {
    //     b'>' => Ok(FileFormat::Fasta),
    //     b'@' => Ok(FileFormat::Fastq),
    //     _ => Err(io::Error::new(
    //         io::ErrorKind::Other,
    //         "Unrecognized file format",
    //     )),
    // }
}

#[cfg(unix)]
extern crate libc;

#[cfg(unix)]
use libc::{getrlimit, rlimit, RLIMIT_NOFILE};

#[cfg(unix)]
pub fn get_file_limit() -> usize {
    let mut limits = rlimit {
        rlim_cur: 0, // 当前（软）限制
        rlim_max: 0, // 最大（硬）限制
    };

    // 使用unsafe块调用getrlimit，因为这是一个外部C函数
    let result = unsafe { getrlimit(RLIMIT_NOFILE, &mut limits) };

    if result == 0 {
        // 如果成功，返回当前软限制转换为usize
        limits.rlim_cur as usize
    } else {
        // 如果失败，输出错误并可能返回一个默认值或panic
        eprintln!("Failed to get file limit");
        0
    }
}

#[cfg(windows)]
pub fn get_file_limit() -> usize {
    8192
}

pub fn create_partition_files(partition: usize, base_path: &PathBuf, prefix: &str) -> Vec<PathBuf> {
    create_dir_all(&base_path).expect(&format!("create dir error {:?}", base_path));
    let file_path = base_path.clone();
    (1..=partition)
        .into_iter()
        .map(|item| file_path.join(format!("{}_{}.k2", prefix, item)))
        .collect()
}

pub fn create_partition_writers(partition_files: &Vec<PathBuf>) -> Vec<BufWriter<File>> {
    partition_files
        .into_iter()
        .map(|item| {
            // 尝试创建文件，如果失败则直接返回错误
            let file = OpenOptions::new()
                .write(true)
                // .append(true) // 确保以追加模式打开文件
                .create(true) // 如果文件不存在，则创建
                .open(item)
                .unwrap();
            BufWriter::new(file)
        })
        .collect()
}

pub fn create_sample_file<P: AsRef<Path>>(filename: P) -> BufWriter<File> {
    let file = OpenOptions::new()
        .write(true)
        .append(true) // 确保以追加模式打开文件
        .create(true) // 如果文件不存在，则创建
        .open(filename)
        .unwrap();
    BufWriter::new(file)
}

use regex::Regex;

// 函数定义
pub fn find_and_sort_files(
    directory: &Path,
    prefix: &str,
    suffix: &str,
) -> io::Result<Vec<PathBuf>> {
    // 构建正则表达式以匹配文件名中的数字
    let pattern = format!(r"{}_(\d+){}", prefix, suffix);
    let re = Regex::new(&pattern).unwrap();

    // 读取指定目录下的所有条目
    let entries = fs::read_dir(directory)?
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| {
            path.is_file()
                && path
                    .file_name()
                    .unwrap()
                    .to_str()
                    .map_or(false, |s| s.starts_with(prefix) && s.ends_with(suffix))
        })
        .collect::<Vec<PathBuf>>();

    // 使用正则表达式提取数字并排序
    let mut sorted_entries = entries
        .into_iter()
        .filter_map(|path| {
            re.captures(path.file_name()?.to_str()?)
                .and_then(|caps| caps.get(1).map(|m| m.as_str().parse::<i32>().ok()))
                .flatten()
                .map(|num| (path, num))
        })
        .collect::<Vec<(PathBuf, i32)>>();

    sorted_entries.sort_by_key(|k| k.1);

    // 检查数字是否从0开始连续
    for (i, (_, num)) in sorted_entries.iter().enumerate() {
        let a_idx = i + 1;
        if a_idx as i32 != *num {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                "File numbers are not continuous starting from 0.",
            ));
        }
    }

    // 返回排序后的文件路径
    Ok(sorted_entries.into_iter().map(|(path, _)| path).collect())
}

pub fn open_file<P: AsRef<Path>>(path: P) -> io::Result<File> {
    File::open(&path).map_err(|e| {
        if e.kind() == io::ErrorKind::NotFound {
            io::Error::new(e.kind(), format!("File not found: {:?}", path.as_ref()))
        } else {
            e
        }
    })
}

/// 获取最新的文件序号
pub fn get_lastest_file_index(file_path: &PathBuf) -> Result<usize> {
    let file_content = fs::read_to_string(&file_path)?;
    // 如果文件内容为空，则默认最大值为0
    let index = if file_content.is_empty() {
        0
    } else {
        file_content
            .lines() // 将内容按行分割
            .filter_map(|line| line.split('\t').next()) // 获取每行的第一列
            .filter_map(|num_str| num_str.parse::<usize>().ok()) // 尝试将第一列的字符串转换为整型
            .max() // 找到最大值
            .unwrap_or(1)
    };
    Ok(index)
}
