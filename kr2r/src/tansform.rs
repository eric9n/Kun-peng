// 使用时需要引用模块路径
use crate::compact_hash::{CellIndex, Compact, K2Cell};
use crate::mmscanner::MinimizerScanner;
use crate::table::CompactHashTable;
use crate::taxonomy::Taxonomy;
use crate::Meros;
use rayon::prelude::*;
use seq_io::fasta::{Reader, Record};
use seq_io::parallel::read_parallel;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Result, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

pub fn find_partition(index: usize, chunk_size: usize) -> usize {
    // 将index转换为usize（假设在这个上下文中这样做是安全的）
    // 注意：在实际应用中，确保这种转换不会造成数据丢失

    // 计算index落在哪个分区内
    let partition_index = index / chunk_size;

    // 计算index在其所属分区内的位置
    // let position_in_partition = index % chunk_size;

    partition_index
}

fn to_u8_vec(cell: &K2Cell) -> Vec<u8> {
    let cell_ptr = cell as *const K2Cell as *const u8;
    let cell_size = std::mem::size_of::<K2Cell>();
    let cell_bytes = unsafe {
        // 将K2Cell实例的内存表示转换为字节切片
        std::slice::from_raw_parts(cell_ptr, cell_size)
    };
    cell_bytes.to_vec()
}

pub fn convert_fna_to_k2_format<P: AsRef<Path>>(
    fna_file: P,
    meros: Meros,
    taxonomy: &Taxonomy,
    id_to_taxon_map: &HashMap<String, u64>,
    threads: u32,
    capacity: usize,
    value_bits: usize,
    writers: &mut Vec<BufWriter<File>>,
    chunk_size: usize,
) {
    let reader = Reader::from_path(fna_file).unwrap();
    let queue_len = (threads - 2) as usize;

    println!("size {:?}", std::mem::size_of::<K2Cell>());

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
                            let index: usize = hash_key.index(capacity);
                            let idx = index % chunk_size;
                            let partition_index = index / chunk_size;
                            let cell = K2Cell::new(
                                idx as u32,
                                hash_key.combined_value(value_bits, taxid as u32),
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
                for cell in k2_cell_list {
                    let partition_index = cell.0;
                    if let Some(writer) = writers.get_mut(partition_index) {
                        writer.write_all(to_u8_vec(&cell.1).as_slice()).unwrap();
                    }
                }

                // for writer in writers.iter_mut() {
                //     writer.flush().expect("Failed to flush writer");
                // }
            }
        },
    );
}

const BATCH_SIZE: usize = 81920; // 定义每批次处理的 K2Cell 数量

pub fn process_sequence<P: AsRef<Path>>(
    chunk_file: P,
    chtm: &mut CompactHashTable,
    taxonomy: &Taxonomy,
    value_bits: usize,
) -> Result<()> {
    let total_counter = AtomicUsize::new(0);
    let size_counter = AtomicUsize::new(0);
    let seq_counter = AtomicUsize::new(0);

    let value_mask = (1 << value_bits) - 1;

    let file = File::open(chunk_file)?;
    let mut reader = BufReader::new(file);

    let cell_size = std::mem::size_of::<K2Cell>();
    let batch_buffer_size = cell_size * BATCH_SIZE;
    let mut batch_buffer = vec![0u8; batch_buffer_size];

    while let Ok(bytes_read) = reader.read(&mut batch_buffer) {
        if bytes_read == 0 {
            break;
        } // 文件末尾

        // 处理读取的数据批次
        let cells_in_batch = bytes_read / cell_size;

        let cells: Vec<K2Cell> = (0..cells_in_batch)
            .into_par_iter()
            .map(|i| {
                let start = i * cell_size;
                let cell =
                    unsafe { std::ptr::read(batch_buffer[start..].as_ptr() as *const K2Cell) };

                cell
            })
            .collect();

        cells.into_iter().for_each(|cell| {
            let item = cell.to_cellindex(value_bits, value_mask);
            if let Some(mut ci) = &chtm.set_page_cell(item) {
                let new_taxid = taxonomy.lca(item.cell.taxid, ci.cell.taxid);
                if ci.cell.taxid != new_taxid {
                    ci.cell.taxid = new_taxid;
                    chtm.update_cell(ci);
                }
            } else {
                size_counter.fetch_add(1, Ordering::SeqCst);
            }
        });
        // for i in 0..cells_in_batch {
        //     let start = i * cell_size;
        //     let cell = unsafe { std::ptr::read(batch_buffer[start..].as_ptr() as *const K2Cell) };
        //     let item = cell.to_cellindex(value_bits, value_mask);

        //     if let Some(mut ci) = &chtm.set_page_cell(item) {
        //         let new_taxid = taxonomy.lca(item.cell.taxid, ci.cell.taxid);
        //         if ci.cell.taxid != new_taxid {
        //             ci.cell.taxid = new_taxid;
        //             chtm.update_cell(ci);
        //         }
        //     } else {
        //         size_counter.fetch_add(1, Ordering::SeqCst);
        //     }
        // }
    }

    chtm.copy_from_page();

    let size_count = size_counter.load(Ordering::SeqCst);
    let seq_count = seq_counter.load(Ordering::SeqCst);
    let total_count = total_counter.load(Ordering::SeqCst);
    // let update_count = update_counter.load(Ordering::SeqCst);
    println!("seq_count {:?}", seq_count);
    println!("size_count {:?}", size_count);
    println!("total_count {:?}", total_count);
    chtm.add_size(size_count);

    Ok(())
}
