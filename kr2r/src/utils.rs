use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Result, Write};
use std::path::{Path, PathBuf};
use walkdir::WalkDir;

/// 读取 seqid2taxid.map 文件。为了裁剪 ncbi 的 taxonomy 树
pub fn read_id_to_taxon_map<P: AsRef<Path>>(filename: P) -> Result<HashMap<String, u64>> {
    let file = File::open(filename)?;
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
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            writeln!(summary_file, "{}", line)?;
        }
    }

    Ok(summary_file_path)
}

pub fn create_seqid2taxid_file<P: AsRef<Path>>(prelim_map_file: P, output_file: P) -> Result<()> {
    let file = File::open(prelim_map_file)?;
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
