// 使用时需要引用模块路径
use clap::Parser;
use kun_peng::args::{parse_size, Build};
use kun_peng::compact_hash::HashConfig;
use kun_peng::db::{convert_fna_to_k2_format, get_bits_for_taxid};
use kun_peng::taxonomy::Taxonomy;
use kun_peng::utils::{
    create_partition_files, create_partition_writers, find_files, get_file_limit,
    read_id_to_taxon_map, set_fd_limit,
};
use kun_peng::IndexOptions;
use std::time::Instant;

#[derive(Parser, Debug, Clone)]
#[clap(author, version, about="prebuild database", long_about = None)]
pub struct Args {
    #[clap(long, value_parser = parse_size, default_value = "1G", help = "Specifies the hash file capacity.\nAcceptable formats include numeric values followed by 'K', 'M', or 'G' (e.g., '1.5G', '250M', '1024K').\nNote: The specified capacity affects the index size, with a factor of 4 applied.\nFor example, specifying '1G' results in an index size of '4G'.\nDefault: 1G (capacity 1G = file size 4G)")]
    pub hash_capacity: usize,

    /// 包含原始配置
    #[clap(flatten)]
    pub build: Build,
}

pub fn run(args: Args, required_capacity: usize) -> Result<(), Box<dyn std::error::Error>> {
    let file_num_limit = get_file_limit();
    let meros = args.build.klmt.as_meros();
    let k2d_dir = &args.build.database;

    let id_to_taxon_map_filename = args.build.database.join("seqid2taxid.map");
    let id_to_taxon_map = read_id_to_taxon_map(&id_to_taxon_map_filename)?;

    let taxonomy_filename = k2d_dir.join("taxo.k2d");
    let taxonomy = Taxonomy::from_file(taxonomy_filename)?;

    let value_bits = get_bits_for_taxid(
        args.build.requested_bits_for_taxid as usize,
        taxonomy.node_count() as f64,
    )
    .expect("more bits required for storing taxid");

    let capacity = required_capacity;
    let partition = (capacity + args.hash_capacity - 1) / args.hash_capacity;
    let hash_config = HashConfig::new(1, capacity, value_bits, 0, partition, args.hash_capacity);

    // 开始计时
    let start = Instant::now();

    let chunk_size = args.hash_capacity as usize;

    if partition >= file_num_limit {
        set_fd_limit(partition as u64 + 1).expect("Failed to set file descriptor limit");
        // panic!("Exceeds File Number Limit");
    }

    let chunk_files = create_partition_files(partition, &k2d_dir, "chunk");
    let mut writers = create_partition_writers(&chunk_files);

    let library_dir = &args.build.database.join("library");
    let fna_files = find_files(&library_dir, "library", ".fna");

    for fna_file in fna_files {
        println!("convert fna file {:?}", fna_file);
        convert_fna_to_k2_format(
            fna_file,
            meros,
            &taxonomy,
            &id_to_taxon_map,
            hash_config,
            &mut writers,
            chunk_size,
            args.build.threads,
        );
    }

    let hash_filename = k2d_dir.join("hash_config.k2d");
    hash_config.write_to_file(&hash_filename)?;
    // 计算持续时间
    let duration = start.elapsed();
    // 打印运行时间
    println!("chunk db took: {:?}", duration);

    let options_filename = k2d_dir.join("opts.k2d");
    let idx_opts = IndexOptions::from_meros(meros);
    idx_opts.write_to_file(options_filename)?;

    Ok(())
}

#[derive(Parser, Debug, Clone)]
#[clap(author, version, about, long_about = None)]
struct ChunkArgs {
    #[clap(flatten)]
    db_args: Args,

    #[arg(short = 'c', long, required = true)]
    pub required_capacity: usize,
}

#[allow(dead_code)]
fn main() {
    let args = ChunkArgs::parse();
    if let Err(e) = run(args.db_args, args.required_capacity) {
        eprintln!("Application error: {}", e);
    }
}
