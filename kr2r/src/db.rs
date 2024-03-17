// 使用时需要引用模块路径
use crate::compact_hash::{CHTableMut, Cell, Compact, HashConfig, Slot};
use crate::mmscanner::MinimizerScanner;
use crate::taxonomy::{NCBITaxonomy, Taxonomy};
use crate::Meros;

use seq_io::fasta::{Reader, Record};
use seq_io::parallel::read_parallel;
use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Result as IOResult, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

// 定义每批次处理的 Cell 数量
const BATCH_SIZE: usize = 81920;

/// 处理k2格式的临时文件,构建数据库
pub fn process_k2file<P: AsRef<Path>, B: Compact>(
    chunk_file: P,
    partition_index: usize,
    chtm: &mut CHTableMut<B>,
    taxonomy: &Taxonomy,
) -> IOResult<()> {
    let total_counter = AtomicUsize::new(0);
    let size_counter = AtomicUsize::new(0);

    let value_mask = chtm.config.value_mask;
    let value_bits = chtm.config.value_bits;

    let file = File::open(chunk_file)?;
    let mut reader = BufReader::new(file);

    let cell_size = std::mem::size_of::<Cell<B>>();
    let batch_buffer_size = cell_size * BATCH_SIZE;
    let mut batch_buffer = vec![0u8; batch_buffer_size];

    while let Ok(bytes_read) = reader.read(&mut batch_buffer) {
        if bytes_read == 0 {
            break;
        } // 文件末尾

        // 处理读取的数据批次
        let cells_in_batch = bytes_read / cell_size;

        let cells = unsafe {
            std::slice::from_raw_parts(batch_buffer.as_ptr() as *const Cell<B>, cells_in_batch)
        };

        cells.into_iter().for_each(|cell| {
            let item = cell.as_slot();
            let item_taxid: u32 = item.value.right(value_mask).to_u32();
            if let Some((flag, mut slot)) = &chtm.set_page_cell(item, partition_index) {
                let slot_taxid = slot.value.right(value_mask).to_u32();
                let new_taxid = taxonomy.lca(item_taxid, slot_taxid);
                if slot_taxid != new_taxid {
                    slot.update_right(B::from_u32(new_taxid), value_bits);
                    chtm.update_cell(flag, slot);
                }
            } else {
                size_counter.fetch_add(1, Ordering::SeqCst);
            }
        });
        total_counter.fetch_add(cells.len(), Ordering::SeqCst);
    }

    chtm.copy_from_page();

    let size_count = size_counter.load(Ordering::SeqCst);
    let total_count = total_counter.load(Ordering::SeqCst);
    println!("size_count {:?}", size_count);
    println!("total_count {:?}", total_count);
    chtm.add_size(size_count);

    Ok(())
}

/// 直接处理fna文件构建数据库
pub fn process_fna<P: AsRef<Path>>(
    fna_file: P,
    meros: Meros,
    chtm: &mut CHTableMut<u32>,
    taxonomy: &Taxonomy,
    id_to_taxon_map: &HashMap<String, u64>,
    threads: u32,
) {
    let reader = Reader::from_path(fna_file).unwrap();
    let queue_len = (threads - 2) as usize;

    let total_counter = AtomicUsize::new(0);
    let size_counter = AtomicUsize::new(0);
    let seq_counter = AtomicUsize::new(0);
    // let update_counter = AtomicUsize::new(0);

    // let capacity = chtm.config.capacity;
    let value_bits = chtm.config.value_bits;
    let hash_config = chtm.config;
    let value_mask = hash_config.value_mask;

    read_parallel(
        reader,
        threads,
        queue_len,
        |record_set| {
            let mut hash_set = BTreeSet::<Slot<u32>>::new();

            for record in record_set.into_iter() {
                seq_counter.fetch_add(1, Ordering::SeqCst);
                if let Ok(seq_id) = record.id() {
                    if let Some(ext_taxid) = id_to_taxon_map.get(seq_id) {
                        let taxid = taxonomy.get_internal_id(*ext_taxid);
                        let kmer_iter = MinimizerScanner::new(record.seq(), meros)
                            .into_iter()
                            .map(|hash_key| hash_config.slot(hash_key, taxid as u32))
                            .collect::<BTreeSet<Slot<u32>>>();

                        total_counter.fetch_add(kmer_iter.len(), Ordering::SeqCst);
                        hash_set.extend(kmer_iter);
                    };
                }
            }
            hash_set
        },
        |record_sets| {
            while let Some(Ok((_, hash_set))) = record_sets.next() {
                hash_set.into_iter().for_each(|item| {
                    let item_taxid = item.value.right(value_mask);
                    if let Some(mut slot) = &chtm.set_table_cell(item.idx, item.value) {
                        let slot_taxid = slot.value.right(value_mask);
                        let new_taxid = taxonomy.lca(item_taxid, slot_taxid);
                        if slot_taxid != new_taxid {
                            slot.update_right(new_taxid, value_bits);
                            chtm.update_cell(&1, slot);
                            // update_counter.fetch_add(1, Ordering::SeqCst);
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
    // let update_count = update_counter.load(Ordering::SeqCst);
    println!("seq_count {:?}", seq_count);
    println!("size_count {:?}", size_count);
    println!("total_count {:?}", total_count);
    chtm.add_size(size_count);
}

/// 生成taxonomy树文件
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

/// 获取需要存储最大内部taxid的bit数量
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

/// 将fna文件转换成k2格式的临时文件
pub fn convert_fna_to_k2_format<P: AsRef<Path>, B: Compact>(
    fna_file: P,
    meros: Meros,
    taxonomy: &Taxonomy,
    id_to_taxon_map: &HashMap<String, u64>,
    hash_config: HashConfig<B>,
    writers: &mut Vec<BufWriter<File>>,
    chunk_size: usize,
    threads: u32,
) {
    let reader = Reader::from_path(fna_file).unwrap();
    let queue_len = (threads - 2) as usize;
    let value_bits = hash_config.value_bits;

    read_parallel(
        reader,
        threads,
        queue_len,
        |record_set| {
            let mut k2_cell_list = Vec::new();
            for record in record_set.into_iter() {
                if let Ok(seq_id) = record.id() {
                    if let Some(ext_taxid) = id_to_taxon_map.get(seq_id) {
                        let taxid = taxonomy.get_internal_id(*ext_taxid);
                        for hash_key in MinimizerScanner::new(record.seq(), meros).into_iter() {
                            let index: usize = hash_config.index(hash_key);
                            let idx = index % chunk_size;
                            let partition_index = index / chunk_size;
                            let cell = Cell::new(
                                idx as u32,
                                B::hash_value(hash_key, value_bits, B::from_u32(taxid)),
                            );
                            k2_cell_list.push((partition_index, cell));
                        }
                    };
                }
            }
            k2_cell_list
        },
        |record_sets| {
            while let Some(Ok((_, k2_cell_list))) = record_sets.next() {
                let cell_size = std::mem::size_of::<Cell<B>>();
                for cell in k2_cell_list {
                    let partition_index = cell.0;
                    if let Some(writer) = writers.get_mut(partition_index) {
                        writer.write_all(cell.1.as_slice(cell_size)).unwrap();
                    }
                }

                // for writer in writers.into_iter() {
                //     writer.flush().expect("io error");
                // }
            }
        },
    );
}

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
