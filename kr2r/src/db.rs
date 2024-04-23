// 使用时需要引用模块路径
use crate::compact_hash::{CHTableMut, Cell, Compact, HashConfig, Slot};
use crate::mmscanner::MinimizerScanner;
use crate::taxonomy::{NCBITaxonomy, Taxonomy};
use crate::Meros;

use byteorder::{LittleEndian, WriteBytesExt};
use rayon::prelude::*;
use seq_io::fasta::{Reader, Record};
use seq_io::parallel::read_parallel;
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
    item: Slot<u32>,
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
            // 成功更新或遍历一圈后回到起点
            break;
        }

        idx = (idx + 1) % page_size; // 移动到下一个索引
        if idx == first_idx {
            // 如果遍历完一整圈还没有找到插入点，可能需要处理溢出或者重新哈希等策略
            break;
        }
    }
}

fn write_u32_to_file(page: &Vec<AtomicU32>, file_path: &PathBuf) -> IOResult<()> {
    // 打开文件用于写入
    let file = File::create(file_path)?;
    let mut writer = BufWriter::new(file);
    let mut count = 0;
    // 遍历 Vec 并写入每个 u32，使用小端字节序
    for item in page {
        let value = item.load(Ordering::Relaxed);
        if value != 0 {
            count += 1;
        }
        writer.write_u32::<LittleEndian>(value)?;
    }

    println!("count {:?}", count);
    writer.flush()?; // 确保所有内容都被写入文件
    Ok(())
}

use memmap2::{Mmap, MmapMut, MmapOptions};
use std::fs::OpenOptions;

pub fn process_k2file1(
    config: HashConfig<u32>,
    chunk_file: &PathBuf,
    taxonomy: &Taxonomy,
    page_size: usize,
    page_index: usize,
) -> IOResult<()> {
    let total_counter = AtomicUsize::new(0);
    // let size_counter = AtomicUsize::new(0);

    let value_mask = config.value_mask;
    let value_bits = config.value_bits;

    let start_index = page_index * page_size;
    let end_index = std::cmp::min((page_index + 1) * page_size, config.capacity);

    let capacity = end_index - start_index;
    let page_file = &chunk_file
        .parent()
        .unwrap()
        .join(format!("page_{}", page_index));
    let hash_file = &chunk_file
        .parent()
        .unwrap()
        .join(format!("hash_{}", page_index));

    let file_len = std::mem::size_of::<u32>() * config.capacity;
    let file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .open(&hash_file)?;
    file.set_len(file_len as u64)?;

    let mut mut_mmap = unsafe { MmapOptions::new().len(file_len).map_mut(&file)? };
    let table = unsafe {
        std::slice::from_raw_parts_mut(mut_mmap.as_mut_ptr() as *mut u32, config.capacity)
    };
    let page: Vec<AtomicU32> = (0..capacity).map(|_| AtomicU32::new(0)).collect();
    // 使用 table 中的值更新 page 直到遇到第一个 0
    for (i, &value) in table.iter().enumerate() {
        if value == 0 {
            break;
        }
        page[i].store(value, Ordering::Relaxed);
    }

    let file = File::open(&chunk_file)?;
    let mut reader = BufReader::new(file);

    let cell_size = std::mem::size_of::<Cell>();
    let batch_buffer_size = cell_size * BATCH_SIZE;
    let mut batch_buffer = vec![0u8; batch_buffer_size];

    while let Ok(bytes_read) = reader.read(&mut batch_buffer) {
        if bytes_read == 0 {
            break;
        } // 文件末尾

        // 处理读取的数据批次
        let cells_in_batch = bytes_read / cell_size;

        let cells = unsafe {
            std::slice::from_raw_parts(batch_buffer.as_ptr() as *const Cell, cells_in_batch)
        };
        cells.par_iter().for_each(|cell| {
            let item = cell.as_slot();
            set_page_cell(taxonomy, &page, item, capacity, value_bits, value_mask);
        });
        total_counter.fetch_add(cells.len(), Ordering::SeqCst);
    }

    write_u32_to_file(&page, &page_file)?;
    println!("total_counter {:?}", total_counter.load(Ordering::SeqCst));
    Ok(())
}

/// 处理k2格式的临时文件,构建数据库
pub fn process_k2file<P: AsRef<Path>>(
    chunk_file: P,
    chtm: &mut CHTableMut<u32>,
    taxonomy: &Taxonomy,
) -> IOResult<()> {
    let total_counter = AtomicUsize::new(0);
    let size_counter = AtomicUsize::new(0);

    let value_mask = chtm.config.value_mask;
    let value_bits = chtm.config.value_bits;

    let file = File::open(chunk_file)?;
    let mut reader = BufReader::new(file);

    let cell_size = std::mem::size_of::<Cell>();
    let batch_buffer_size = cell_size * BATCH_SIZE;
    let mut batch_buffer = vec![0u8; batch_buffer_size];

    while let Ok(bytes_read) = reader.read(&mut batch_buffer) {
        if bytes_read == 0 {
            break;
        } // 文件末尾

        // 处理读取的数据批次
        let cells_in_batch = bytes_read / cell_size;

        let cells = unsafe {
            std::slice::from_raw_parts(batch_buffer.as_ptr() as *const Cell, cells_in_batch)
        };

        cells.iter().for_each(|cell| {
            let item = cell.as_slot();
            let item_taxid: u32 = item.value.right(value_mask).to_u32();

            if let Some((flag, mut slot)) = &chtm.set_page_cell(item) {
                let slot_taxid = slot.value.right(value_mask).to_u32();
                let new_taxid = taxonomy.lca(item_taxid, slot_taxid);
                if slot_taxid != new_taxid {
                    slot.update_right(u32::from_u32(new_taxid), value_bits);
                    chtm.update_cell(flag, slot);
                }
            } else {
                size_counter.fetch_add(1, Ordering::SeqCst);
            }
        });
        total_counter.fetch_add(cells.len(), Ordering::SeqCst);
    }

    println!("total_counter {:?}", total_counter.load(Ordering::SeqCst));
    chtm.copy_from_page();

    let size_count = size_counter.load(Ordering::SeqCst);
    let total_count = total_counter.load(Ordering::SeqCst);
    println!("size_count {:?}", size_count);
    println!("total_count {:?}", total_count);
    chtm.add_size(size_count);

    Ok(())
}

// /// 直接处理fna文件构建数据库
// pub fn process_fna<P: AsRef<Path>>(
//     fna_file: P,
//     meros: Meros,
//     chtm: &mut CHTableMut<u32>,
//     taxonomy: &Taxonomy,
//     id_to_taxon_map: &HashMap<String, u64>,
//     threads: u32,
// ) {
//     let reader = Reader::from_path(fna_file).unwrap();
//     let queue_len = (threads - 2) as usize;

//     let total_counter = AtomicUsize::new(0);
//     let size_counter = AtomicUsize::new(0);
//     let seq_counter = AtomicUsize::new(0);
//     // let update_counter = AtomicUsize::new(0);

//     // let capacity = chtm.config.capacity;
//     let value_bits = chtm.config.value_bits;
//     let hash_config = chtm.config;
//     let value_mask = hash_config.value_mask;

//     read_parallel(
//         reader,
//         threads,
//         queue_len,
//         |record_set| {
//             let mut hash_set = BTreeSet::<Slot<u32>>::new();

//             for record in record_set.into_iter() {
//                 seq_counter.fetch_add(1, Ordering::SeqCst);
//                 if let Ok(seq_id) = record.id() {
//                     if let Some(ext_taxid) = id_to_taxon_map.get(seq_id) {
//                         let taxid = taxonomy.get_internal_id(*ext_taxid);
//                         let kmer_iter = MinimizerScanner::new(record.seq(), meros)
//                             .into_iter()
//                             .map(|hash_key| hash_config.slot(hash_key, taxid as u32))
//                             .collect::<BTreeSet<Slot<u32>>>();

//                         total_counter.fetch_add(kmer_iter.len(), Ordering::SeqCst);
//                         hash_set.extend(kmer_iter);
//                     };
//                 }
//             }
//             hash_set
//         },
//         |record_sets| {
//             while let Some(Ok((_, hash_set))) = record_sets.next() {
//                 hash_set.into_iter().for_each(|item| {
//                     let item_taxid = item.value.right(value_mask);
//                     if let Some(mut slot) = &chtm.set_table_cell(item.idx, item.value) {
//                         let slot_taxid = slot.value.right(value_mask);
//                         let new_taxid = taxonomy.lca(item_taxid, slot_taxid);
//                         if slot_taxid != new_taxid {
//                             slot.update_right(new_taxid, value_bits);
//                             chtm.update_cell(&1, slot);
//                             // update_counter.fetch_add(1, Ordering::SeqCst);
//                         }
//                     } else {
//                         size_counter.fetch_add(1, Ordering::SeqCst);
//                     }
//                 });
//             }
//         },
//     );
//     let size_count = size_counter.load(Ordering::SeqCst);
//     let seq_count = seq_counter.load(Ordering::SeqCst);
//     let total_count = total_counter.load(Ordering::SeqCst);
//     // let update_count = update_counter.load(Ordering::SeqCst);
//     println!("seq_count {:?}", seq_count);
//     println!("size_count {:?}", size_count);
//     println!("total_count {:?}", total_count);
//     chtm.add_size(size_count);
// }

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
    let cell_size = std::mem::size_of::<Cell>();

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
                            let cell =
                                Cell::new(idx as u32, u32::hash_value(hash_key, value_bits, taxid));
                            k2_cell_list.push((partition_index, cell));
                        }
                    };
                }
            }
            k2_cell_list
        },
        |record_sets| {
            while let Some(Ok((_, k2_cell_map))) = record_sets.next() {
                for cell in k2_cell_map {
                    let partition_index = cell.0;
                    if let Some(writer) = writers.get_mut(partition_index) {
                        writer.write_all(&cell.1.as_slice(cell_size)).unwrap();
                    }
                }
            }
        },
    );
}
