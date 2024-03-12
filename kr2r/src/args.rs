// 使用时需要引用模块路径
use crate::utils::expand_spaced_seed_mask;
use crate::{construct_seed_template, parse_binary, Meros, BITS_PER_CHAR};
use crate::{
    DEFAULT_KMER_LENGTH, DEFAULT_MINIMIZER_LENGTH, DEFAULT_MINIMIZER_SPACES, DEFAULT_TOGGLE_MASK,
};
use clap::Parser;

use std::path::PathBuf;

#[derive(Parser, Debug, Clone)]
#[clap(version, about = "build database")]
pub struct Build {
    /// build database directory or file
    #[arg(long, required = true)]
    pub source: PathBuf,

    /// Kraken 2 hash table filename
    #[clap(short = 'H', required = true)]
    pub hashtable_filename: PathBuf,

    /// Kraken 2 options filename
    #[clap(short = 'o', required = true)]
    pub options_filename: PathBuf,

    /// Set length of k-mers, k must be positive integer, k=35, k cannot be less than l
    #[clap(short, long, value_parser = clap::value_parser!(u64).range(1..), default_value_t = DEFAULT_KMER_LENGTH)]
    pub k_mer: u64,

    /// Set length of minimizers, 1 <= l <= 31
    #[clap(short, long, value_parser = clap::value_parser!(u8).range(1..=31), default_value_t = DEFAULT_MINIMIZER_LENGTH)]
    pub l_mer: u8,

    /// Bit storage requested for taxid 0 <= r < 31
    #[clap(short, long, value_parser = clap::value_parser!(u8).range(0..31), default_value_t = 0)]
    pub requested_bits_for_taxid: u8,

    /// Minimizer ordering toggle mask
    #[clap(short = 'T', long, default_value_t = DEFAULT_TOGGLE_MASK)]
    pub toggle_mask: u64,

    /// Number of characters in minimizer that are ignored in comparisons
    #[clap(long, default_value_t = DEFAULT_MINIMIZER_SPACES)]
    pub minimizer_spaces: u8,

    // /// Name of Kraken 2 database
    // #[arg(short, long = "db")]
    // database: PathBuf,
    #[arg(short = 'c', long, required = true)]
    pub required_capacity: u64,

    /// Number of threads
    #[clap(short = 'p', long, default_value_t = 4)]
    pub threads: usize,
}

impl Build {
    pub fn as_meros(&self) -> Meros {
        let seed = construct_seed_template(self.l_mer as usize, self.minimizer_spaces as usize);
        let space_seed_mask = parse_binary(&seed).unwrap();
        let space_seed_mask = expand_spaced_seed_mask(space_seed_mask, BITS_PER_CHAR as u64);

        Meros::new(
            self.k_mer as usize,
            self.l_mer as usize,
            Some(space_seed_mask),
            Some(self.toggle_mask),
            Some(0),
        )
    }
}

#[derive(Parser, Debug, Clone)]
#[clap(version, about = "taxonomy")]
pub struct Taxo {
    /// Kraken 2 taxonomy filename
    #[clap(short = 't', required = true)]
    pub taxonomy_filename: PathBuf,

    /// Sequence ID to taxon map filename
    #[clap(short = 'm', required = true)]
    pub id_to_taxon_map_filename: PathBuf,

    /// NCBI taxonomy directory name
    #[clap(short, long, required = true)]
    pub ncbi_taxonomy_directory: PathBuf,
}
