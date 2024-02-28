use kr2r::compact_hash::{Cell, CellIndex, Compact, CompactHash, CompactHashTable};
use memmap2::MmapOptions;
use rayon::prelude::*;
use std::fs::OpenOptions;

fn main() -> std::io::Result<()> {
    let capacity = 2usize;
    let file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .open("make_mut.md")?;
    file.set_len(4 * capacity as u64)?;

    // 分配足够的空间来存储两个 CompactHashCell
    let mut mut_mmap = unsafe { MmapOptions::new().len(4 * capacity).map_mut(&file)? };

    // 使用 unsafe 获取 CompactHashCell 的切片引用
    let cells: &mut [u32] =
        unsafe { std::slice::from_raw_parts_mut(mut_mmap.as_mut_ptr() as *mut u32, capacity) };

    // for i in 0..capacity {
    //     println!("cell {:?}", cells[i]);
    // }

    // 直接在内存映射上操作 CompactHashCell 实例
    // cells[0] = CompactHashCell(4);
    // cells[1] = CompactHashCell(5);

    // let chtm = CompactHashTable::from("test.k2d")?;
    // println!("chtm {:?}", chtm);
    // let cell = chtm.table[0];

    // println!("cell {:?}", cell.compacted_key(10));
    // println!("cell taxid {:?}", cell.taxid(chtm.value_mask));

    let chtm1 = CompactHashTable::from("lib/hash1.k2d")?;
    println!("chtm1 {:?}", chtm1);
    println!("value::: {:?}", chtm1.get(2361267427824489423));
    let key: u64 = 2361267427824489423;

    println!("compacted key {:?}", key.compacted(10));
    println!("index: {:?}", key.index(560818468));

    let table = chtm1.table;
    let cell = &table[24550723];
    println!("cell {:?}", cell.compacted_key(10));
    println!("cell taxid {:?}", cell.taxid(chtm1.value_mask));

    let chtm1 = CompactHashTable::from("lib/hash.k2d")?;
    println!("chtm1 {:?}", chtm1);
    println!("value::: {:?}", chtm1.get(2361267427824489423));

    // let chtm = CompactHashTable::new("test.k2d", 1, 10)?;
    // chtm.set_cell(CellIndex::new(0, 1, 32));
    // let table = chtm.table;
    // let cell = &table[24550723];
    // println!("cell {:?}", cell.compacted_key(10));
    // println!("cell taxid {:?}", cell.taxid(chtm.value_mask));

    // let chtm1 = CompactHashTableMut::new("hash_file.k2d", 1731287771, 10)?;

    // let mut_table = chtm1.table;

    // let table = chtm.table;
    // mut_table.copy_from_slice(table);
    // let mut size = 0;

    // println!("chtm11 {:?}", chtm1);

    // chtm1.save_size(393344595)?;
    // let hash_file = File::open("lib/hash.k2d")?;
    // let hash_mmp = unsafe { Mmap::map(&file)? };

    // mut_mmap.flush();
    // 直接修改第 2-5 字节
    // 注意：这里的索引范围是 2..6，因为范围是左闭右开的
    // mut_mmap[2..6].copy_from_slice(b"test");

    Ok(())
}
