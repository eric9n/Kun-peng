// 使用时需要引用模块路径
use crate::compact_hash::{CellIndex, Compact, CompactHashTable};
use crate::fmix64 as murmur_hash3;
use crate::mmscanner::MinimizerScanner;
use crate::taxonomy::{NCBITaxonomy, Taxonomy};
use crate::Meros;

use rayon::prelude::*;
use seq_io::fasta::{Reader, Record};
use seq_io::parallel::read_parallel;
use std::collections::{BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU32, AtomicUsize, Ordering};

/// 处理fasta文件,构建数据库
pub fn process_sequence<P: AsRef<Path>>(
    fna_file: P,
    meros: Meros,
    chtm: &CompactHashTable<AtomicU32>,
    taxonomy: &Taxonomy,
    id_to_taxon_map: &HashMap<String, u64>,
    threads: u32,
) {
    let reader = Reader::from_path(fna_file).unwrap();
    let queue_len = (threads - 2) as usize;

    let total_counter = AtomicUsize::new(0);
    let size_counter = AtomicUsize::new(0);
    let seq_counter = AtomicUsize::new(0);
    let update_counter = AtomicUsize::new(0);

    let capacity = chtm.capacity;
    let value_bits = chtm.value_bits;

    read_parallel(
        reader,
        threads,
        queue_len,
        |record_set| {
            let mut scanner = MinimizerScanner::new(meros);
            let mut hash_set = BTreeSet::<CellIndex>::new();

            for record in record_set.into_iter() {
                seq_counter.fetch_add(1, Ordering::SeqCst);
                let seq = record.seq();
                if let Ok(seq_id) = record.id() {
                    if let Some(ext_taxid) = id_to_taxon_map.get(seq_id) {
                        let taxid = taxonomy.get_internal_id(*ext_taxid);

                        scanner.set_seq_end(seq);
                        while let Some(minimizer) = scanner.next_minimizer(seq) {
                            let hash_key = murmur_hash3(minimizer);

                            let ci = CellIndex::new(
                                hash_key.index(capacity),
                                hash_key.compacted(value_bits),
                                taxid as u32,
                            );

                            hash_set.insert(ci);

                            total_counter.fetch_add(1, Ordering::SeqCst);
                        }

                        scanner.reset();
                    };
                }
            }
            hash_set
        },
        |record_sets| {
            while let Some(Ok((_, hash_set))) = record_sets.next() {
                hash_set.into_par_iter().for_each(|item| {
                    if let Some(mut ci) = &chtm.set_cell(item) {
                        let new_taxid = taxonomy.lca(item.cell.taxid, ci.cell.taxid);
                        if ci.cell.taxid != new_taxid {
                            ci.cell.taxid = new_taxid;
                            chtm.update_cell(ci);
                            update_counter.fetch_add(1, Ordering::SeqCst);
                        }
                    } else {
                        size_counter.fetch_add(1, Ordering::SeqCst);
                    }
                });
            }
        },
    );
    let size_count = size_counter.load(Ordering::SeqCst);
    let seq_count = seq_counter.load(Ordering::SeqCst);
    let total_count = total_counter.load(Ordering::SeqCst);
    let update_count = update_counter.load(Ordering::SeqCst);
    println!("seq_count {:?}", seq_count);
    println!("size_count {:?}", size_count);
    println!("total_count {:?}", total_count);
    println!("update_count {:?}", update_count);
    chtm.update_size(size_count);
}

pub fn generate_taxonomy(
    ncbi_taxonomy_directory: &PathBuf,
    taxonomy_filename: &PathBuf,
    id_map: &HashMap<String, u64>,
) -> Result<Taxonomy, Box<dyn std::error::Error>> {
    let nodes_filename = ncbi_taxonomy_directory.join("nodes.dmp");
    let names_filename = ncbi_taxonomy_directory.join("names.dmp");
    let mut ncbi = NCBITaxonomy::from_ncbi(nodes_filename, names_filename)?;

    for (_, id) in id_map.into_iter() {
        ncbi.mark_node(*id);
    }
    let mut taxo = ncbi.convert_to_kraken_taxonomy();
    taxo.generate_external_to_internal_id_map();
    taxo.build_path_cache();
    taxo.write_to_disk(&taxonomy_filename)?;

    Ok(taxo)
}

pub fn get_bits_for_taxid(
    requested_bits_for_taxid: usize,
    node_count: f64,
) -> Result<usize, String> {
    // 计算存储taxonomy节点数量所需的最小位数
    let bits_needed_for_value = (node_count.log2().ceil() as usize).max(1);

    // 检查是否需要更多位来存储taxid
    if requested_bits_for_taxid > 0 && bits_needed_for_value > requested_bits_for_taxid as usize {
        return Err("more bits required for storing taxid".to_string());
    }

    Ok(bits_needed_for_value.max(requested_bits_for_taxid))
}
