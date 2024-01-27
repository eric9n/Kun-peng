use kr2r::compact_hash::CompactHashTable;
use kr2r::compact_hash::CompactHashTableMut;
use memmap2::{Mmap, MmapOptions};
use std::fs::File;
use std::fs::OpenOptions;
use std::io::Write;
use std::ops::DerefMut;

use kr2r::compact_hash::CompactHashCell;

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
    let cells: &mut [CompactHashCell] = unsafe {
        std::slice::from_raw_parts_mut(mut_mmap.as_mut_ptr() as *mut CompactHashCell, capacity)
    };

    for i in 0..capacity {
        println!("cell {:?}", cells[i]);
    }

    // 直接在内存映射上操作 CompactHashCell 实例
    cells[0] = CompactHashCell(4);
    cells[1] = CompactHashCell(5);

    let chtm = CompactHashTable::from("lib/hash.k2d")?;
    println!("chtm {:?}", chtm);

    let chtm1 = CompactHashTableMut::new("hash_file.k2d", 1731287771, 10)?;

    let mut_table = chtm1.table;
    let table = chtm.table;
    mut_table.copy_from_slice(table);

    // chtm1.save_size(393344595)?;
    // let hash_file = File::open("lib/hash.k2d")?;
    // let hash_mmp = unsafe { Mmap::map(&file)? };

    // mut_mmap.flush();
    // 直接修改第 2-5 字节
    // 注意：这里的索引范围是 2..6，因为范围是左闭右开的
    // mut_mmap[2..6].copy_from_slice(b"test");

    Ok(())
}
