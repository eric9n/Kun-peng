use std::path::Path;
use walkdir::WalkDir;

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
