#[repr(C)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Cell {
    pub idx: u64,
    pub reads_id: u32,
    pub data: u32,
}
use kr2r::utils::format_bytes;

fn main() {
    let cell_size = std::mem::size_of::<Cell>();

    println!("cell size {:?}", cell_size);

    println!("format {:?}", format_bytes(557655920.0));
}
