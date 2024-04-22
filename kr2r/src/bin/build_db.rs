// 使用时需要引用模块路径
use clap::Parser;
use kr2r::compact_hash::{CHTableMut, HashConfig};
use kr2r::utils::{find_library_fna_files, read_id_to_taxon_map};
use kr2r::IndexOptions;
// use std::collections::HashMap;
use kr2r::args::{Build, Taxo};
use kr2r::db::{generate_taxonomy, get_bits_for_taxid, process_fna};
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// 包含原始配置
    #[clap(flatten)]
    build: Build,

    #[clap(flatten)]
    taxo: Taxo,

    // /// Name of Kraken 2 database
    // #[arg(short, long = "db")]
    // database: PathBuf,
    #[arg(short = 'c', long, required = true)]
    pub required_capacity: u64,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let meros = args.build.klmt.as_meros();

    let id_to_taxon_map_filename = args
        .taxo
        .id_to_taxon_map_filename
        .unwrap_or(args.build.source.join("seqid2taxid.map"));

    let id_to_taxon_map = read_id_to_taxon_map(&id_to_taxon_map_filename)?;

    let taxonomy_filename = args
        .taxo
        .taxonomy_filename
        .unwrap_or(args.build.source.join("taxo.k2d"));

    let ncbi_taxonomy_directory = args
        .taxo
        .ncbi_taxonomy_directory
        .unwrap_or(args.build.source.join("taxonomy"));

    let taxonomy = generate_taxonomy(
        &ncbi_taxonomy_directory,
        &taxonomy_filename,
        &id_to_taxon_map,
    )?;

    let value_bits = get_bits_for_taxid(
        args.build.requested_bits_for_taxid as usize,
        taxonomy.node_count() as f64,
    )
    .expect("more bits required for storing taxid");

    let capacity = args.required_capacity as usize;

    let hashtable_filename = args
        .build
        .hashtable_filename
        .unwrap_or(args.build.source.clone().join("hash.k2d"));
    let hash_config = HashConfig::<u32>::new(capacity, value_bits, 0, 0, 0);
    let mut chtm = CHTableMut::new(hashtable_filename, hash_config, 0, capacity)?;

    let source: PathBuf = args.build.source.clone();
    let fna_files = if source.is_file() {
        vec![source.to_string_lossy().to_string()]
    } else {
        find_library_fna_files(args.build.source)
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
            args.build.threads as u32,
        )
    }
    // 计算持续时间
    let duration = start.elapsed();
    // 打印运行时间
    println!("build db took: {:?}", duration);

    let idx_opts = IndexOptions::from_meros(meros);
    let options_filename = args
        .build
        .options_filename
        .unwrap_or(source.join("opts.k2d"));
    idx_opts.write_to_file(options_filename)?;

    Ok(())
}
