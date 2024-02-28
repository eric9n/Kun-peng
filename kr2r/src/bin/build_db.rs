// 使用时需要引用模块路径
use clap::Parser;
use kr2r::compact_hash::CompactHashTable;
use kr2r::utils::{expand_spaced_seed_mask, find_library_fna_files, read_id_to_taxon_map};
use kr2r::{construct_seed_template, parse_binary, Meros};
use kr2r::{
    IndexOptions, BITS_PER_CHAR, DEFAULT_KMER_LENGTH, DEFAULT_MINIMIZER_LENGTH,
    DEFAULT_MINIMIZER_SPACES, DEFAULT_TOGGLE_MASK,
};
// use std::collections::HashMap;
use kr2r::db::{generate_taxonomy, get_bits_for_taxid, process_sequence};
use std::path::PathBuf;
use std::sync::atomic::AtomicU32;

#[derive(Parser, Debug, Clone)]
#[clap(version, about = "build database")]
struct Build {
    /// build database directory or file
    #[arg(long, required = true)]
    source: PathBuf,

    /// Kraken 2 hash table filename
    #[clap(short = 'H', required = true)]
    hashtable_filename: PathBuf,

    /// Kraken 2 taxonomy filename
    #[clap(short = 't', required = true)]
    taxonomy_filename: PathBuf,

    /// Sequence ID to taxon map filename
    #[clap(short = 'm', required = true)]
    id_to_taxon_map_filename: PathBuf,

    /// Kraken 2 options filename
    #[clap(short = 'o', required = true)]
    options_filename: PathBuf,

    /// NCBI taxonomy directory name
    #[clap(short, long, required = true)]
    ncbi_taxonomy_directory: PathBuf,

    /// Set length of k-mers, k must be positive integer, k=35, k cannot be less than l
    #[clap(short, long, value_parser = clap::value_parser!(u64).range(1..), default_value_t = DEFAULT_KMER_LENGTH)]
    k_mer: u64,

    /// Set length of minimizers, 1 <= l <= 31
    #[clap(short, long, value_parser = clap::value_parser!(u8).range(1..=31), default_value_t = DEFAULT_MINIMIZER_LENGTH)]
    l_mer: u8,

    /// Bit storage requested for taxid 0 <= r < 31
    #[clap(short, long, value_parser = clap::value_parser!(u8).range(0..31), default_value_t = 0)]
    requested_bits_for_taxid: u8,

    /// Minimizer ordering toggle mask
    #[clap(short = 'T', long, default_value_t = DEFAULT_TOGGLE_MASK)]
    toggle_mask: u64,

    /// Number of characters in minimizer that are ignored in comparisons
    #[clap(long, default_value_t = DEFAULT_MINIMIZER_SPACES)]
    minimizer_spaces: u8,

    // /// Name of Kraken 2 database
    // #[arg(short, long = "db")]
    // database: PathBuf,
    #[arg(short = 'c', long, required = true)]
    required_capacity: u64,

    /// Number of threads
    #[clap(short = 'p', long, default_value_t = 4)]
    threads: usize,
}

impl Build {
    pub fn as_meros(&self, spaced_seed_mask: u64) -> Meros {
        Meros::new(
            self.k_mer as usize,
            self.l_mer as usize,
            Some(spaced_seed_mask),
            Some(self.toggle_mask),
            Some(0),
        )
    }
}

// fn set_minimizer_lca(hash: &mut CompactHashTableMut, minimizer: u64, taxid: u64, tax: &Taxonomy) {
//     let mut old_value: u32 = 0;
//     let mut new_value: u64 = taxid;
//     // while !hash.compare_and_set(minimizer, new_value as u32, &mut old_value) {
//     //     new_value = tax.lowest_common_ancestor(old_value as u64, taxid);
//     // }
// }

// // This function exists to deal with NCBI's use of \x01 characters to denote
// // the start of a new FASTA header in the same line (for non-redundant DBs).
// // We return all sequence IDs in a header line, not just the first.
// fn extract_ncbi_sequence_ids(header: &[u8]) -> Vec<String> {
//     let mut list = Vec::new();
//     let mut current_str = Vec::new();

//     // Skip the first byte if it's '>'.
//     for &c in header.iter().skip(1) {
//         match c {
//             b'\x01' | b' ' | b'\t' | b'\n' | b'\r' => {
//                 if !current_str.is_empty() {
//                     // Convert current_str bytes to String and push to list
//                     if let Ok(s) = String::from_utf8(current_str.clone()) {
//                         list.push(s);
//                     }
//                     current_str.clear();
//                 }
//             }
//             _ => {
//                 current_str.push(c);
//             }
//         }
//     }

//     // Check if there's a remaining ID to be added
//     if !current_str.is_empty() {
//         if let Ok(s) = String::from_utf8(current_str) {
//             list.push(s);
//         }
//     }

//     list
// }

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Build::parse();
    let seed = construct_seed_template(
        args.l_mer.clone() as usize,
        args.minimizer_spaces.clone() as usize,
    );
    let space_seed_mask = parse_binary(&seed).unwrap();
    let space_seed_mask = expand_spaced_seed_mask(space_seed_mask, BITS_PER_CHAR as u64);

    let meros = args.as_meros(space_seed_mask);

    let id_to_taxon_map = read_id_to_taxon_map(&args.id_to_taxon_map_filename)?;

    let taxonomy = generate_taxonomy(
        &args.ncbi_taxonomy_directory,
        &args.taxonomy_filename,
        &id_to_taxon_map,
    )?;

    let bits_for_taxid = get_bits_for_taxid(
        args.requested_bits_for_taxid as usize,
        taxonomy.node_count() as f64,
    )
    .expect("more bits required for storing taxid");

    // 1211893248
    let capacity = args.required_capacity as usize;

    let chtm =
        CompactHashTable::<AtomicU32>::new(args.hashtable_filename, capacity, bits_for_taxid)?;

    let source: PathBuf = args.source.clone();
    let fna_files = if source.is_file() {
        vec![source.to_string_lossy().to_string()]
    } else {
        find_library_fna_files(args.source)
    };

    println!("start...");
    for fna_file in fna_files {
        process_sequence(
            fna_file,
            meros,
            &chtm,
            &taxonomy,
            &id_to_taxon_map,
            args.threads as u32,
        )
    }
    println!("success...");

    let idx_opts = IndexOptions::new(
        args.k_mer as usize,
        args.l_mer as usize,
        space_seed_mask,
        args.toggle_mask,
        true,
        0,
    );
    idx_opts.write_to_file(args.options_filename)?;

    Ok(())
}
