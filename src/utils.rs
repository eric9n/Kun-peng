use std::collections::{BTreeMap as Map, HashMap};
use std::fs::{self, create_dir_all, File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Result};
use std::path::{Path, PathBuf};
use walkdir::WalkDir;

/// Reads the seqid2taxid.map file to create a mapping for trimming the NCBI taxonomy tree.
///
/// This function reads a file containing sequence IDs and their corresponding taxon IDs,
/// and creates a HashMap mapping sequence IDs to taxon IDs.
///
/// # Arguments
///
/// * `filename` - A path-like object representing the file to be read.
///
/// # Returns
///
/// Returns a `Result` containing a `HashMap<String, u64>` where the keys are sequence IDs
/// and the values are taxon IDs, or an error if the file cannot be read or parsed.
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
/// # use kun_peng::utils::expand_spaced_seed_mask; // Replace with the appropriate crate name
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
/// # use kun_peng::utils::expand_spaced_seed_mask;
/// // No expansion, factor is 0
/// assert_eq!(expand_spaced_seed_mask(0b1010, 0), 0b1010);
///
/// // No expansion, factor is greater than 64
/// assert_eq!(expand_spaced_seed_mask(0b1010, 65), 0b1010);
/// ```
pub fn expand_spaced_seed_mask(spaced_seed_mask: u64, bit_expansion_factor: u64) -> u64 {
    // Check if bit_expansion_factor is within the valid range
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

pub fn find_files<P: AsRef<Path>>(path: P, prefix: &str, suffix: &str) -> Vec<PathBuf> {
    let mut files: Vec<PathBuf> = WalkDir::new(path)
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|e| {
            e.path()
                .file_name()
                .and_then(|name| name.to_str())
                .map(|name| name.starts_with(prefix) && name.ends_with(suffix))
                .unwrap_or(false)
        })
        .map(|e| e.path().to_path_buf())
        .collect();
    files.sort_unstable();
    files
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

#[cfg(unix)]
extern crate libc;

#[cfg(unix)]
use libc::{getrlimit, rlimit, setrlimit, RLIMIT_NOFILE};

/// Get the current file descriptor limit for the process.
///
/// This function retrieves the soft limit on the number of open file descriptors
/// allowed for the current process on Unix-like systems.
///
/// # Returns
///
/// Returns the current soft limit as a `usize`, or 0 if the limit couldn't be retrieved.
///
/// # Platform-specific behavior
///
/// This function is only available on Unix-like systems. On other platforms,
/// it may return a default value or not be compiled.
#[cfg(unix)]
pub fn get_file_limit() -> usize {
    let mut limits = rlimit {
        rlim_cur: 0, // Current (soft) limit
        rlim_max: 0, // Maximum (hard) limit
    };

    // Use an unsafe block to call getrlimit, as it's an external C function
    let result = unsafe { getrlimit(RLIMIT_NOFILE, &mut limits) };

    if result == 0 {
        // If successful, return the current soft limit converted to usize
        limits.rlim_cur as usize
    } else {
        // If failed, print an error and return 0
        eprintln!("Failed to get file limit");
        0
    }
}

#[cfg(unix)]
pub fn set_fd_limit(new_limit: u64) -> io::Result<()> {
    let rlim = rlimit {
        rlim_cur: new_limit,
        rlim_max: new_limit,
    };

    let ret = unsafe { setrlimit(RLIMIT_NOFILE, &rlim) };
    if ret != 0 {
        return Err(io::Error::last_os_error());
    }
    Ok(())
}

#[cfg(windows)]
pub fn get_file_limit() -> usize {
    8192
}

#[cfg(windows)]
pub fn set_fd_limit(_new_limit: u64) -> io::Result<()> {
    Ok(())
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
            // Try to create the file, return an error if failed
            let file = OpenOptions::new()
                .write(true)
                .append(true) // Ensure opening the file in append mode
                .create(true) // Create the file if it doesn't exist
                .open(item)
                .unwrap();
            BufWriter::new(file)
        })
        .collect()
}

pub fn create_sample_file<P: AsRef<Path>>(filename: P) -> BufWriter<File> {
    let file = OpenOptions::new()
        .write(true)
        .append(true) // Ensure opening the file in append mode
        .create(true) // Create the file if it doesn't exist
        .open(filename)
        .unwrap();
    BufWriter::new(file)
}

use regex::Regex;

pub fn find_and_trans_bin_files(
    directory: &Path,
    prefix: &str,
    suffix: &str,
    check: bool,
) -> io::Result<Map<usize, Vec<PathBuf>>> {
    // Aggregate file paths with the same number
    // Build a regular expression to match the first number in the filename
    let pattern = format!(r"{}_(\d+)_\d+{}", prefix, suffix);
    let re = Regex::new(&pattern).expect("Invalid regex pattern");

    // Read all entries in the specified directory
    let mut map_entries = Map::new();
    for entry in fs::read_dir(directory)? {
        let path = entry?.path();

        if path.is_file() {
            if let Some(file_name) = path.file_name().and_then(|name| name.to_str()) {
                // Use the regular expression to match the filename and extract the first number part
                if let Some(cap) = re.captures(file_name) {
                    if let Some(m) = cap.get(1) {
                        if let Ok(num) = m.as_str().parse::<usize>() {
                            map_entries.entry(num).or_insert_with(Vec::new).push(path);
                        }
                    }
                }
            }
        }
    }

    if check {
        // Check if the numbers are continuous starting from 1
        let mut keys: Vec<_> = map_entries.keys().cloned().collect();
        keys.sort_unstable();
        for (i, &key) in keys.iter().enumerate() {
            if i + 1 != key {
                return Err(io::Error::new(
                    io::ErrorKind::NotFound,
                    "File numbers are not continuous starting from 1.",
                ));
            }
        }
    }

    // Return the aggregated file paths
    Ok(map_entries)
}

pub fn find_and_trans_files(
    directory: &Path,
    prefix: &str,
    suffix: &str,
    check: bool,
) -> io::Result<Map<usize, PathBuf>> {
    // Build a regular expression to match the number in the filename
    let pattern = format!(r"{}_(\d+){}", prefix, suffix);
    let re = Regex::new(&pattern).unwrap();

    // Read all entries in the specified directory
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

    // Extract numbers using the regular expression and store them in a BTreeMap
    let mut map_entries = Map::new();
    for path in entries {
        if let Some(fname) = path.file_name().and_then(|name| name.to_str()) {
            if let Some(cap) = re.captures(fname) {
                if let Some(m) = cap.get(1) {
                    if let Ok(num) = m.as_str().parse::<usize>() {
                        map_entries.insert(num, path);
                    }
                }
            }
        }
    }

    if check {
        // Check if the numbers are continuous starting from 0
        let mut keys: Vec<_> = map_entries.keys().cloned().collect();
        keys.sort_unstable();
        for (i, key) in keys.iter().enumerate() {
            if i + 1 != *key {
                return Err(io::Error::new(
                    io::ErrorKind::NotFound,
                    "File numbers are not continuous starting from 1.",
                ));
            }
        }
    }

    // Return the sorted file paths
    Ok(map_entries)
}

// Function definition
pub fn find_and_sort_files(
    directory: &Path,
    prefix: &str,
    suffix: &str,
    check: bool,
) -> io::Result<Vec<PathBuf>> {
    // Build a regular expression to match the number in the filename
    let pattern = format!(r"{}_(\d+){}", prefix, suffix);
    let re = Regex::new(&pattern).unwrap();

    // Read all entries in the specified directory
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

    // Extract numbers using the regular expression and sort
    let mut sorted_entries = entries
        .into_iter()
        .filter_map(|path| {
            re.captures(path.file_name()?.to_str()?)
                .and_then(|caps| caps.get(1).map(|m| m.as_str().parse::<usize>().ok()))
                .flatten()
                .map(|num| (path, num))
        })
        .collect::<Vec<(PathBuf, usize)>>();

    sorted_entries.sort_by_key(|k| k.1);

    if check {
        // Check if the numbers are continuous starting from 0
        for (i, (_, num)) in sorted_entries.iter().enumerate() {
            let a_idx = i + 1;
            if a_idx != *num {
                return Err(io::Error::new(
                    io::ErrorKind::NotFound,
                    "File numbers are not continuous starting from 1.",
                ));
            }
        }
    }

    // Return the sorted file paths
    Ok(sorted_entries
        .iter()
        .map(|(path, _)| path.clone())
        .collect())
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

/// Get the latest file index
pub fn get_lastest_file_index(file_path: &PathBuf) -> Result<usize> {
    let file_content = fs::read_to_string(&file_path)?;
    // If the file content is empty, default the maximum value to 0
    let index = if file_content.is_empty() {
        0
    } else {
        file_content
            .lines() // Split the content by lines
            .filter_map(|line| line.split('\t').next()) // Get the first column of each line
            .filter_map(|num_str| num_str.parse::<usize>().ok()) // Try to convert the first column string to integer
            .max() // Find the maximum value
            .unwrap_or(1)
    };
    Ok(index)
}
