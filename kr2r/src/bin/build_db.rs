// 使用时需要引用模块路径
use kr2r::compact_hash::CompactHashTableMut;
use kr2r::fmix64 as murmur_hash3;
use kr2r::mmscanner::MinimizerScanner;
use kr2r::taxonomy::{NCBITaxonomy, Taxonomy};
use kr2r::utils::{expand_spaced_seed_mask, find_library_fna_files, read_id_to_taxon_map};
use kr2r::{construct_seed_template, parse_binary, Meros};
use kr2r::{
    IndexOptions, BITS_PER_CHAR, DEFAULT_KMER_LENGTH, DEFAULT_MINIMIZER_LENGTH,
    DEFAULT_MINIMIZER_SPACES, DEFAULT_TOGGLE_MASK,
};
use seq_io::fasta::{Reader, Record};
use seq_io::parallel::read_parallel;
use std::collections::{BTreeSet, HashMap};
use std::path::{Path, PathBuf};

use clap::{error::ErrorKind, Error, Parser};

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

fn generate_taxonomy(
    args: &Build,
    id_map: &HashMap<String, u64>,
) -> Result<Taxonomy, Box<dyn std::error::Error>> {
    let nodes_filename = args.ncbi_taxonomy_directory.join("nodes.dmp");
    let names_filename = args.ncbi_taxonomy_directory.join("names.dmp");
    let mut ncbi = NCBITaxonomy::from_ncbi(nodes_filename, names_filename)?;

    for (_, id) in id_map.into_iter() {
        ncbi.mark_node(*id);
    }
    let taxo = ncbi.convert_to_kraken_taxonomy();
    taxo.write_to_disk(&args.taxonomy_filename)?;

    Ok(taxo)
}

fn get_bits_for_taxid(requested_bits_for_taxid: usize, node_count: f64) -> usize {
    // 计算存储taxonomy节点数量所需的最小位数
    let bits_needed_for_value = (node_count.log2().ceil() as usize).max(1);

    // 检查是否需要更多位来存储taxid
    if requested_bits_for_taxid > 0 && bits_needed_for_value > requested_bits_for_taxid as usize {
        Error::raw(
            ErrorKind::ValueValidation,
            "more bits required for storing taxid",
        )
        .exit();
    }

    bits_needed_for_value.max(requested_bits_for_taxid)
}

fn set_minimizer_lca(hash: &mut CompactHashTableMut, minimizer: u64, taxid: u64, tax: &Taxonomy) {
    let mut old_value: u32 = 0;
    let mut new_value: u64 = taxid;
    while !hash.compare_and_set(minimizer, new_value as u32, &mut old_value) {
        new_value = tax.lowest_common_ancestor(old_value as u64, taxid);
    }
}

// This function exists to deal with NCBI's use of \x01 characters to denote
// the start of a new FASTA header in the same line (for non-redundant DBs).
// We return all sequence IDs in a header line, not just the first.
fn extract_ncbi_sequence_ids(header: &[u8]) -> Vec<String> {
    let mut list = Vec::new();
    let mut current_str = Vec::new();

    // Skip the first byte if it's '>'.
    for &c in header.iter().skip(1) {
        match c {
            b'\x01' | b' ' | b'\t' | b'\n' | b'\r' => {
                if !current_str.is_empty() {
                    // Convert current_str bytes to String and push to list
                    if let Ok(s) = String::from_utf8(current_str.clone()) {
                        list.push(s);
                    }
                    current_str.clear();
                }
            }
            _ => {
                current_str.push(c);
            }
        }
    }

    // Check if there's a remaining ID to be added
    if !current_str.is_empty() {
        if let Ok(s) = String::from_utf8(current_str) {
            list.push(s);
        }
    }

    list
}

fn process_sequence<P: AsRef<Path>>(
    fna_file: P,
    meros: Meros,
    chtm: &mut CompactHashTableMut,
    taxonomy: &Taxonomy,
    id_to_taxon_map: &HashMap<String, u64>,
    threads: u32,
) {
    let reader = Reader::from_path(fna_file).unwrap();
    let queue_len = (threads - 2) as usize;
    read_parallel(
        reader,
        threads,
        queue_len,
        |record_set| {
            let mut scanner = MinimizerScanner::new(meros);

            let mut hash_taxid = HashMap::<u64, BTreeSet<u64>>::new();

            for record in record_set.into_iter() {
                let seq = record.seq();
                let all_sequence_ids = extract_ncbi_sequence_ids(record.head());
                let mut taxid = 0;
                let mut hash_set = BTreeSet::<u64>::new();
                for seq_id in all_sequence_ids {
                    if id_to_taxon_map.contains_key(&seq_id) {
                        continue;
                    }
                    let ext_taxid = id_to_taxon_map.get(&seq_id).unwrap();
                    taxid = taxonomy
                        .lowest_common_ancestor(taxid, taxonomy.get_internal_id(*ext_taxid));
                }
                scanner.set_seq_end(seq);
                while let Some(minimizer) = scanner.next_minimizer(seq) {
                    let hash_key = murmur_hash3(minimizer);
                    hash_set.insert(hash_key);
                }
                hash_taxid.insert(taxid, hash_set);
                scanner.reset();
            }
            hash_taxid
        },
        |record_sets| {
            while let Some(Ok((_, m_set))) = record_sets.next() {
                for (taxid, hash_set) in m_set.into_iter() {
                    for hash_key in hash_set {
                        set_minimizer_lca(chtm, hash_key, taxid, taxonomy)
                    }
                }
            }
        },
    );
}

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

    let taxonomy = generate_taxonomy(&args, &id_to_taxon_map)?;

    let bits_for_taxid = get_bits_for_taxid(
        args.requested_bits_for_taxid as usize,
        taxonomy.node_count() as f64,
    );

    // 1211893248
    let capacity = args.required_capacity as usize;

    let mut chtm = CompactHashTableMut::new(args.hashtable_filename, capacity, bits_for_taxid)?;

    let source: PathBuf = args.source.clone();
    let fna_files = if source.is_file() {
        vec![source.to_string_lossy().to_string()]
    } else {
        find_library_fna_files(args.source)
    };

    for fna_file in fna_files {
        process_sequence(
            fna_file,
            meros,
            &mut chtm,
            &taxonomy,
            &id_to_taxon_map,
            args.threads as u32,
        )
    }

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
