// 使用时需要引用模块路径
use clap::Parser;
use kr2r::table::{CHTableMut, HashConfig};
use kr2r::utils::{find_library_fna_files, read_id_to_taxon_map};
use kr2r::IndexOptions;
// use std::collections::HashMap;
use kr2r::args::Build;
use kr2r::db::{generate_taxonomy, get_bits_for_taxid, process_fna};
use std::path::PathBuf;
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Build::parse();
    let meros = args.as_meros();

    let id_to_taxon_map = read_id_to_taxon_map(&args.id_to_taxon_map_filename)?;

    let taxonomy = generate_taxonomy(
        &args.ncbi_taxonomy_directory,
        &args.taxonomy_filename,
        &id_to_taxon_map,
    )?;

    let value_bits = get_bits_for_taxid(
        args.requested_bits_for_taxid as usize,
        taxonomy.node_count() as f64,
    )
    .expect("more bits required for storing taxid");

    // 1211893248
    let capacity = args.required_capacity as usize;

    let hash_config = HashConfig::<u32>::new(capacity, value_bits, 0);
    let mut chtm = CHTableMut::<u32>::new(args.hashtable_filename, hash_config, 0, capacity)?;

    let source: PathBuf = args.source.clone();
    let fna_files = if source.is_file() {
        vec![source.to_string_lossy().to_string()]
    } else {
        find_library_fna_files(args.source)
    };

    println!("start...");
    // 开始计时
    let start = Instant::now();
    for fna_file in fna_files {
        process_fna(
            fna_file,
            meros,
            &mut chtm,
            &taxonomy,
            &id_to_taxon_map,
            args.threads as u32,
        )
    }
    // 计算持续时间
    let duration = start.elapsed();
    // 打印运行时间
    println!("build db took: {:?}", duration);

    let idx_opts = IndexOptions::new(
        args.k_mer as usize,
        args.l_mer as usize,
        meros.spaced_seed_mask,
        args.toggle_mask,
        true,
        0,
    );
    idx_opts.write_to_file(args.options_filename)?;

    Ok(())
}
