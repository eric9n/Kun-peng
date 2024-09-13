use crate::compact_hash::Row;
use crate::utils::open_file;
use seqkmer::Meros;
use seqkmer::OptionPair;
use seqkmer::CURRENT_REVCOM_VERSION;
use std::fs::File;
use std::io::{Read, Result as IoResult, Write};
use std::mem;
use std::path::Path;

/// Parses a binary string into a u64
///
/// # Arguments
///
/// * `src` - The binary string to parse
///
/// # Returns
///
/// A Result containing the parsed u64 value or a parsing error
pub fn parse_binary(src: &str) -> Result<u64, std::num::ParseIntError> {
    u64::from_str_radix(src, 2)
}

/// Constructs a seed template based on minimizer length and spaces
///
/// # Arguments
///
/// * `minimizer_len` - The length of the minimizer
/// * `minimizer_spaces` - The number of spaces in the minimizer
///
/// # Returns
///
/// A String representing the constructed seed template
pub fn construct_seed_template(minimizer_len: usize, minimizer_spaces: usize) -> String {
    if minimizer_len / 4 < minimizer_spaces {
        panic!(
            "number of minimizer spaces ({}) exceeds max for minimizer len ({}); max: {}",
            minimizer_spaces,
            minimizer_len,
            minimizer_len / 4
        );
    }
    let core = "1".repeat(minimizer_len - 2 * minimizer_spaces);
    let spaces = "01".repeat(minimizer_spaces);
    format!("{}{}", core, spaces)
}

/// Converts a u64 value to Option<u64>, filtering out zero values
///
/// # Arguments
///
/// * `value` - The u64 value to convert
///
/// # Returns
///
/// An Option<u64> containing the value if it's non-zero, or None otherwise
pub fn u64_to_option(value: u64) -> Option<u64> {
    Option::from(value).filter(|&x| x != 0)
}

/// Represents a group of hits with associated rows and range
pub struct HitGroup {
    pub rows: Vec<Row>,
    /// Range example: (0..10], left-open right-closed
    pub range: OptionPair<(usize, usize)>,
}

impl HitGroup {
    /// Creates a new HitGroup
    pub fn new(rows: Vec<Row>, range: OptionPair<(usize, usize)>) -> Self {
        Self { rows, range }
    }

    /// Calculates the capacity of the HitGroup
    pub fn capacity(&self) -> usize {
        self.range.reduce(0, |acc, range| acc + range.1 - range.0)
    }

    /// Calculates the required score based on a confidence threshold
    pub fn required_score(&self, confidence_threshold: f64) -> u64 {
        (confidence_threshold * self.capacity() as f64).ceil() as u64
    }
}

/// Represents options for indexing
#[repr(C)]
#[derive(Debug)]
pub struct IndexOptions {
    pub k: usize,
    pub l: usize,
    pub spaced_seed_mask: u64,
    pub toggle_mask: u64,
    pub dna_db: bool,
    pub minimum_acceptable_hash_value: u64,
    pub revcom_version: i32, // Throws an error if equal to 0
    pub db_version: i32,     // Reserved for future database structure changes
    pub db_type: i32,        // Reserved for future use of other data structures
}

impl IndexOptions {
    /// Creates a new IndexOptions instance
    pub fn new(
        k: usize,
        l: usize,
        spaced_seed_mask: u64,
        toggle_mask: u64,
        dna_db: bool,
        minimum_acceptable_hash_value: u64,
    ) -> Self {
        Self {
            k,
            l,
            spaced_seed_mask,
            toggle_mask,
            dna_db,
            minimum_acceptable_hash_value,
            revcom_version: CURRENT_REVCOM_VERSION as i32,
            db_version: 0,
            db_type: 0,
        }
    }

    /// Reads IndexOptions from a file
    ///
    /// # Arguments
    ///
    /// * `file_path` - The path to the file containing IndexOptions
    ///
    /// # Returns
    ///
    /// An IoResult containing the read IndexOptions
    pub fn read_index_options<P: AsRef<Path>>(file_path: P) -> IoResult<Self> {
        let mut file = open_file(file_path)?;
        let mut buffer = vec![0; std::mem::size_of::<Self>()];
        file.read_exact(&mut buffer)?;

        let idx_opts = unsafe {
            // Ensure this conversion is safe, depending on the exact layout and source of the data
            std::ptr::read(buffer.as_ptr() as *const Self)
        };
        if idx_opts.revcom_version != CURRENT_REVCOM_VERSION as i32 {
            // Trigger a panic if the version is 0
            panic!("Unsupported version (revcom_version == 0)");
        }

        Ok(idx_opts)
    }

    /// Writes IndexOptions to a file
    ///
    /// # Arguments
    ///
    /// * `file_path` - The path to the file where IndexOptions will be written
    ///
    /// # Returns
    ///
    /// An IoResult indicating success or failure of the write operation
    pub fn write_to_file<P: AsRef<Path>>(&self, file_path: P) -> IoResult<()> {
        let mut file = File::create(file_path)?;

        // Convert the struct to a byte slice. This is an unsafe operation as we're
        // forcing memory content to be interpreted as bytes. This requires IndexOptions
        // to be #[repr(C)], and all fields must be safely copyable in their raw memory representation.
        let bytes: &[u8] = unsafe {
            std::slice::from_raw_parts(
                (self as *const IndexOptions) as *const u8,
                mem::size_of::<IndexOptions>(),
            )
        };

        file.write_all(bytes)?;
        Ok(())
    }

    /// Creates IndexOptions from a Meros instance
    pub fn from_meros(meros: Meros) -> Self {
        Self::new(
            meros.k_mer,
            meros.l_mer,
            meros.spaced_seed_mask,
            meros.toggle_mask,
            true,
            meros.min_clear_hash_value.unwrap_or_default(),
        )
    }

    /// Converts IndexOptions to a Meros instance
    pub fn as_meros(&self) -> Meros {
        Meros::new(
            self.k,
            self.l,
            u64_to_option(self.spaced_seed_mask),
            u64_to_option(self.toggle_mask),
            u64_to_option(self.minimum_acceptable_hash_value),
        )
    }
}
