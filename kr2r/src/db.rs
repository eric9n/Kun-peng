// 使用时需要引用模块路径
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
// 定义每批次处理的 Cell 数量
const BATCH_SIZE: usize = 81920;

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
                None // 当前值不匹配，尝试下一个索引
            }
        });

        if result.is_ok() || idx == first_idx {
            break;
        }

        idx = (idx + 1) % page_size;
        if idx == first_idx {
            break;
        }
    }
}

fn write_hashtable_to_file(
    page: &Vec<AtomicU32>,
    file_path: &PathBuf,
    page_index: u64,
    capacity: u64,
) -> IOResult<usize> {
    // 打开文件用于写入
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

    writer.flush()?; // 确保所有内容都被写入文件
    Ok(count)
}

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
        } // 文件末尾

        // 处理读取的数据批次
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

/// 生成taxonomy树文件
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

            Some(k2_cell_list)
        },
        |record_sets| {
            while let Some(Some(k2_cell_map)) = record_sets.next() {
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
