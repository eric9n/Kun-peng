#[cfg(feature = "dna")]
pub mod constants {
    pub const DEFAULT_KMER_LENGTH: u64 = 35;
    pub const DEFAULT_MINIMIZER_LENGTH: u8 = 31;
    pub const DEFAULT_MINIMIZER_SPACES: u8 = 7;

    pub const BITS_PER_CHAR: usize = 2;
}

#[cfg(feature = "protein")]
pub mod constants {
    pub const DEFAULT_KMER_LENGTH: u64 = 15;
    pub const DEFAULT_MINIMIZER_LENGTH: u8 = 12;
    pub const DEFAULT_MINIMIZER_SPACES: u8 = 0;

    pub const BITS_PER_CHAR: usize = 4;
}

pub const DEFAULT_TOGGLE_MASK: u64 = 0xe37e28c4271b5a2d;
pub const DEFAULT_SPACED_SEED_MASK: u64 = 0;
pub const CURRENT_REVCOM_VERSION: u8 = 1;
