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

#[cfg(feature = "dna")]
#[inline]
pub fn char_to_value(c: u8) -> Option<u64> {
    match c {
        b'A' | b'a' => Some(0x00),
        b'C' | b'c' => Some(0x01),
        b'G' | b'g' => Some(0x02),
        b'T' | b't' => Some(0x03),
        _ => None,
    }
}

#[cfg(feature = "protein")]
#[inline]
pub fn char_to_value(c: u8) -> Option<64> {
    match c {
        // stop codons/rare amino acids
        b'*' | b'U' | b'u' | b'O' | b'o' => Some(0x00),
        // alanine
        b'A' | b'a' => Some(0x01),
        // asparagine, glutamine, serine
        b'N' | b'n' | b'Q' | b'q' | b'S' | b's' => Some(0x02),
        // cysteine
        b'C' | b'c' => Some(0x03),
        // aspartic acid, glutamic acid
        b'D' | b'd' | b'E' | b'e' => Some(0x04),
        // phenylalanine
        b'F' | b'f' => Some(0x05),
        // glycine
        b'G' | b'g' => Some(0x06),
        // histidine
        b'H' | b'h' => Some(0x07),
        // isoleucine, leucine
        b'I' | b'i' | b'L' | b'l' => Some(0x08),
        // lysine
        b'K' | b'k' => Some(0x09),
        // proline
        b'P' | b'p' => Some(0x0a),
        // arginine
        b'R' | b'r' => Some(0x0b),
        // methionine, valine
        b'M' | b'm' | b'V' | b'v' => Some(0x0c),
        // threonine
        b'T' | b't' => Some(0x0d),
        // tryptophan
        b'W' | b'w' => Some(0x0e),
        // tyrosine
        b'Y' | b'y' => Some(0x0f),
        _ => None,
    }
}

#[inline]
fn reverse_complement(mut kmer: u64, n: usize) -> u64 {
    // Reverse bits while leaving bit pairs (nucleotides) intact.

    // Swap consecutive pairs of bits
    kmer = (kmer >> 2 & 0x3333333333333333) | (kmer << 2 & 0xCCCCCCCCCCCCCCCC);

    // Swap consecutive nibbles (4-bit groups)
    kmer = (kmer >> 4 & 0x0F0F0F0F0F0F0F0F) | (kmer << 4 & 0xF0F0F0F0F0F0F0F0);

    // Swap consecutive bytes
    kmer = (kmer >> 8 & 0x00FF00FF00FF00FF) | (kmer << 8 & 0xFF00FF00FF00FF00);

    // Swap consecutive pairs of bytes
    kmer = (kmer >> 16 & 0x0000FFFF0000FFFF) | (kmer << 16 & 0xFFFF0000FFFF0000);

    // Swap the two halves of the 64-bit word
    kmer = (kmer >> 32) | (kmer << 32);

    // Complement the bits, shift to the right length, and mask to get the desired length
    (!kmer >> (64 - n * 2)) & ((1u64 << (n * 2)) - 1)

    // if revcom_version == 0 {
    //     // Complement the bits and mask to get the desired length
    //     !kmer & ((1u64 << (n * 2)) - 1)
    // } else {
    //     // Complement the bits, shift to the right length, and mask to get the desired length
    //     (!kmer >> (64 - n * 2)) & ((1u64 << (n * 2)) - 1)
    // }
}

#[cfg(feature = "dna")]
#[inline]
pub fn canonical_representation(kmer: u64, n: usize) -> u64 {
    let revcom = reverse_complement(kmer, n);
    if kmer < revcom {
        kmer
    } else {
        revcom
    }
}

#[cfg(feature = "protein")]
#[inline]
pub fn canonical_representation(kmer: u64, n: usize, revcom_version: u8) -> u64 {
    kmer
}

pub const DEFAULT_TOGGLE_MASK: u64 = 0xe37e28c4271b5a2d;
pub const DEFAULT_SPACED_SEED_MASK: u64 = 0;
pub const CURRENT_REVCOM_VERSION: u8 = 1;
