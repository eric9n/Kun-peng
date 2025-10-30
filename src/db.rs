use crate::compact_hash::{Compact, HashConfig, Slot};
// use crate::mmscanner::MinimizerScanner;
use crate::taxonomy::{NCBITaxonomy, Taxonomy};
use seqkmer::{read_parallel, BufferFastaReader, Meros};

use crate::utils::open_file;
use byteorder::{LittleEndian, WriteBytesExt};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Result as IOResult, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU32, AtomicUsize, Ordering};

// Define the number of Cells processed per batch
const BATCH_SIZE: usize = 81920;

/// Sets a cell in the page with the given item, handling collisions and LCA calculations
///
/// # Arguments
///
/// * `taxonomy` - The taxonomy used for LCA calculations
/// * `page` - The page of AtomicU32 cells
/// * `item` - The Slot item to be set
/// * `page_size` - The size of the page
/// * `value_bits` - The number of bits used for the value
/// * `value_mask` - The mask used to extract the value
fn set_page_cell(
    taxonomy: &Taxonomy,
    page: &[AtomicU32],
    item: &Slot<u32>,
    page_size: usize,
    value_bits: usize,
    value_mask: usize,
) {
    let mut idx = item.idx % page_size;
    let item_taxid: u32 = item.value.right(value_mask).to_u32();
    let compact_key = item.value.left(value_bits);
    let first_idx = idx;

    loop {
        let result = page[idx].fetch_update(Ordering::SeqCst, Ordering::Relaxed, |current| {
            let current_taxid = current.right(value_mask).to_u32();
            let current_key = current.left(value_bits);

            if current == 0 || current_taxid == u32::default() {
                Some(item.value)
            } else if current_key == compact_key {
                let new_taxid = taxonomy.lca(item_taxid, current_taxid);
                Some(u32::combined(compact_key, new_taxid, value_bits))
            } else {
                None // Current value doesn't match, try the next index
            }
        });

        match result {
            Ok(_) => {
                // `fetch_update` 成功 (我们写入了数据或LCA)
                // 任务完成，退出循环
                break;
            }
            Err(_) => {
                // `fetch_update` 失败 (返回 None)，意味着 slot 被
                // 另一个 *不同的 key* 占用了。我们必须进行线性探测。
                idx = (idx + 1) % page_size;
                if idx == first_idx {
                    // 我们已经绕了完整的一圈，没有找到空位
                    // TODO: 在这里添加哈希页已满的日志或错误处理
                    break; // 放弃
                }
                // 继续下一次 `loop` 迭代，尝试新的 `idx`
            }
        }
    }
}

/// Writes the hash table to a file
///
/// # Arguments
///
/// * `page` - The page of AtomicU32 cells to write
/// * `file_path` - The path to the output file
/// * `page_index` - The index of the current page
/// * `capacity` - The capacity of the page
///
/// # Returns
///
/// The number of non-zero items written to the file
fn write_hashtable_to_file(
    page: &Vec<AtomicU32>,
    file_path: &PathBuf,
    page_index: u64,
    capacity: u64,
) -> IOResult<usize> {
    // Open the file for writing
    let file = File::create(file_path)?;
    let mut writer = BufWriter::new(file);
    let mut count = 0;
    writer.write_u64::<LittleEndian>(page_index)?;
    writer.write_u64::<LittleEndian>(capacity)?;

    for item in page {
        let value = item.load(Ordering::Relaxed);
        if value != 0 {
            count += 1;
        }

        writer.write_u32::<LittleEndian>(value)?;
    }

    writer.flush()?; // Ensure all content is written to the file
    Ok(count)
}

/// Processes a k2 file and updates the hash table
///
/// # Arguments
///
/// * `config` - The HashConfig for the process
/// * `database` - The path to the database
/// * `chunk_file` - The path to the chunk file
/// * `taxonomy` - The taxonomy used for processing
/// * `page_size` - The size of each page
/// * `page_index` - The index of the current page
///
/// # Returns
///
/// The number of items processed
pub fn process_k2file(
    config: HashConfig,
    database: &PathBuf,
    chunk_file: &PathBuf,
    taxonomy: &Taxonomy,
    page_size: usize,
    page_index: usize,
) -> IOResult<usize> {
    let total_counter = AtomicUsize::new(0);

    let value_mask = config.value_mask;
    let value_bits = config.value_bits;

    let start_index = (page_index - 1) * page_size;
    let end_index = std::cmp::min(page_index * page_size, config.capacity);

    let capacity = end_index - start_index;
    let page_file = database.join(format!("hash_{}.k2d", page_index));

    let page: Vec<AtomicU32> = (0..capacity).map(|_| AtomicU32::new(0)).collect();

    let file = open_file(&chunk_file)?;
    let mut reader = BufReader::new(file);

    let cell_size = std::mem::size_of::<Slot<u32>>();
    let batch_buffer_size = cell_size * BATCH_SIZE;
    let mut batch_buffer = vec![0u8; batch_buffer_size];

    while let Ok(bytes_read) = reader.read(&mut batch_buffer) {
        if bytes_read == 0 {
            break;
        } // End of file

        // Process the read data batch
        let cells_in_batch = bytes_read / cell_size;

        let cells = unsafe {
            std::slice::from_raw_parts(batch_buffer.as_ptr() as *const Slot<u32>, cells_in_batch)
        };
        cells.par_iter().for_each(|item| {
            set_page_cell(taxonomy, &page, item, capacity, value_bits, value_mask);
        });
        total_counter.fetch_add(cells.len(), Ordering::SeqCst);
    }

    let size_count =
        write_hashtable_to_file(&page, &page_file, page_index as u64, capacity as u64)?;
    Ok(size_count)
}

/// Generates a taxonomy tree file
///
/// # Arguments
///
/// * `ncbi_taxonomy_directory` - The directory containing NCBI taxonomy files
/// * `taxonomy_filename` - The output filename for the generated taxonomy
/// * `id_map` - A map of string IDs to u64 IDs
///
/// # Returns
///
/// The generated Taxonomy
pub fn generate_taxonomy(
    ncbi_taxonomy_directory: &PathBuf,
    taxonomy_filename: &PathBuf,
    id_map: &HashMap<String, u64>,
) -> IOResult<Taxonomy> {
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

/// Calculates the number of bits required to store the maximum internal taxid
///
/// # Arguments
///
/// * `requested_bits_for_taxid` - The requested number of bits for taxid storage
/// * `node_count` - The number of nodes in the taxonomy
///
/// # Returns
///
/// The number of bits required for taxid storage
pub fn get_bits_for_taxid(
    requested_bits_for_taxid: usize,
    node_count: f64,
) -> Result<usize, String> {
    // Calculate the minimum number of bits needed to store the taxonomy node count
    let bits_needed_for_value = (node_count.log2().ceil() as usize).max(1);

    // Check if more bits are required for storing taxid
    if requested_bits_for_taxid > 0 && bits_needed_for_value > requested_bits_for_taxid as usize {
        return Err("more bits required for storing taxid".to_string());
    }

    Ok(bits_needed_for_value.max(requested_bits_for_taxid))
}

/// Converts an FNA file to the k2 format temporary file
///
/// # Arguments
///
/// * `fna_file` - The input FNA file path
/// * `meros` - The Meros instance for k-mer processing
/// * `taxonomy` - The taxonomy used for processing
/// * `id_to_taxon_map` - A map of string IDs to taxon IDs
/// * `hash_config` - The HashConfig for the process
/// * `writers` - A vector of BufWriters for output
/// * `chunk_size` - The size of each chunk
/// * `threads` - The number of threads to use for processing
pub fn convert_fna_to_k2_format<P: AsRef<Path>>(
    fna_file: P,
    meros: Meros,
    taxonomy: &Taxonomy,
    id_to_taxon_map: &HashMap<String, u64>,
    hash_config: HashConfig,
    writers: &mut Vec<BufWriter<File>>,
    chunk_size: usize,
    threads: usize,
) {
    let mut reader = BufferFastaReader::from_path(fna_file, 1).unwrap();
    let value_bits = hash_config.value_bits;
    let cell_size = std::mem::size_of::<Slot<u32>>();

    read_parallel(
        &mut reader,
        threads,
        &meros,
        |seqs| {
            let mut k2_cell_list = Vec::new();

            for record in seqs {
                let header = &record.header;
                record.body.apply_mut(|m_iter| {
                    if let Some(ext_taxid) = id_to_taxon_map.get(&header.id) {
                        let taxid = taxonomy.get_internal_id(*ext_taxid);
                        let k2_cell: Vec<(usize, Slot<u32>)> = m_iter
                            .map(|(_, hash_key)| {
                                let index: usize = hash_config.index(hash_key);
                                let idx = index % chunk_size;
                                let partition_index = index / chunk_size;
                                let cell =
                                    Slot::new(idx, u32::hash_value(hash_key, value_bits, taxid));
                                (partition_index, cell)
                            })
                            .collect();

                        k2_cell_list.extend_from_slice(&k2_cell);
                    }
                });
            }

            k2_cell_list
        },
        |record_sets| {
            while let Some(data) = record_sets.next() {
                let k2_cell_map = data.unwrap();
                for cell in k2_cell_map {
                    let partition_index = cell.0;
                    if let Some(writer) = writers.get_mut(partition_index) {
                        writer.write_all(&cell.1.as_slice(cell_size)).unwrap();
                    }
                }
            }
        },
    )
    .expect("failed");
}
