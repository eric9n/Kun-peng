// 使用时需要引用模块路径
use clap::Parser;
use kr2r::args::{parse_size, Build};
use kr2r::compact_hash::HashConfig;
use kr2r::db::{
    convert_fna_to_k2_format, generate_taxonomy, get_bits_for_taxid, process_k2file,
    write_config_to_file,
};
use kr2r::utils::{
    create_partition_files, create_partition_writers, find_library_fna_files, get_file_limit,
    read_id_to_taxon_map,
};
use kr2r::IndexOptions;
use std::fs::remove_file;
use std::time::Instant;

#[derive(Parser, Debug, Clone)]
#[clap(author, version, about="build database", long_about = None)]
pub struct Args {
    // /// database hash chunk directory and other files
    // #[clap(long)]
    // pub k2d_dir: Option<PathBuf>,
    #[clap(long, value_parser = parse_size, default_value = "1G", help = "Specifies the hash file capacity.\nAcceptable formats include numeric values followed by 'K', 'M', or 'G' (e.g., '1.5G', '250M', '1024K').\nNote: The specified capacity affects the index size, with a factor of 4 applied.\nFor example, specifying '1G' results in an index size of '4G'.\nDefault: 1G (capacity 1G = file size 4G)")]
    pub hash_capacity: usize,

    // chunk temp directory
    // #[clap(long)]
    // pub chunk_dir: PathBuf,
    /// 包含原始配置
    #[clap(flatten)]
    pub build: Build,
    // #[arg(short = 'm')]
    // pub id_to_taxon_map_filename: Option<PathBuf>,
}

pub fn run(args: Args, required_capacity: usize) -> Result<(), Box<dyn std::error::Error>> {
    let file_num_limit = get_file_limit();
    let meros = args.build.klmt.as_meros();

    let id_to_taxon_map_filename = args.build.database.join("seqid2taxid.map");

    let id_to_taxon_map = read_id_to_taxon_map(&id_to_taxon_map_filename)?;

    let k2d_dir = &args.build.database;

    let taxonomy_filename = k2d_dir.join("taxo.k2d");

    let ncbi_taxonomy_directory = &args.build.database.join("taxonomy");

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

    let capacity = required_capacity;
    let partition = (capacity + args.hash_capacity - 1) / args.hash_capacity;
    let hash_config = HashConfig::new(capacity, value_bits, 0, partition, args.hash_capacity);

    // 开始计时
    let start = Instant::now();

    let chunk_size = args.hash_capacity as usize;

    if partition >= file_num_limit {
        panic!("Exceeds File Number Limit");
    }

    let chunk_files = create_partition_files(partition, &k2d_dir, "chunk");
    let mut writers = create_partition_writers(&chunk_files);

    let fna_files = find_library_fna_files(&args.build.database);

    for fna_file in &fna_files {
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
    let partition = chunk_files.len();
    let mut size: u64 = 0;

    for i in 1..=partition {
        // 计算持续时间
        let count = process_k2file(
            hash_config,
            &k2d_dir,
            &chunk_files[i - 1],
            &taxonomy,
            chunk_size,
            i,
        )?;
        size += count as u64;
        let duration = start.elapsed();
        println!(
            "process chunk file {:?}/{:}: duration: {:?}",
            i, partition, duration
        );
    }

    write_config_to_file(
        &hash_filename,
        partition as u64,
        args.hash_capacity as u64,
        capacity as u64,
        size,
        32 - hash_config.value_bits as u64,
        hash_config.value_bits as u64,
    )?;

    // 计算持续时间
    let duration = start.elapsed();
    // 打印运行时间
    println!("build k2 db took: {:?}", duration);

    let options_filename = k2d_dir.join("opts.k2d");
    let idx_opts = IndexOptions::from_meros(meros);
    idx_opts.write_to_file(options_filename)?;

    for chunk_file in chunk_files {
        remove_file(chunk_file)?;
    }
    Ok(())
}

#[derive(Parser, Debug, Clone)]
#[clap(author, version, about, long_about = None)]
struct BuildArgs {
    #[clap(flatten)]
    db_args: Args,

    #[arg(short = 'c', long, required = true)]
    pub required_capacity: usize,
}

#[allow(dead_code)]
fn main() {
    let args = BuildArgs::parse();
    if let Err(e) = run(args.db_args, args.required_capacity) {
        eprintln!("Application error: {}", e);
    }
}
