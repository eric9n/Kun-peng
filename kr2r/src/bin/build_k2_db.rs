// 使用时需要引用模块路径
use clap::Parser;
use kr2r::args::{Build, Taxo, ONEGB, U32MAXPLUS};
use kr2r::compact_hash::{CHTableMut, HashConfig};
use kr2r::db::{convert_fna_to_k2_format, generate_taxonomy, get_bits_for_taxid, process_k2file};
use kr2r::utils::{
    create_partition_files, create_partition_writers, find_library_fna_files, format_bytes,
    get_file_limit, read_id_to_taxon_map,
};
use kr2r::IndexOptions;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser, Debug, Clone)]
#[clap(author, version, about="build database", long_about = None)]
pub struct Args {
    /// 包含原始配置
    #[clap(flatten)]
    pub build: Build,

    #[clap(flatten)]
    pub taxo: Taxo,

    // // /// Name of Kraken 2 database
    // // #[arg(short, long = "db")]
    // // database: PathBuf,
    // #[arg(short = 'c', long, required = true)]
    // pub required_capacity: u64,
    /// chunk directory
    #[clap(long)]
    pub chunk_dir: PathBuf,

    /// chunk size 1-4(GB) [1073741824-4294967295] default: 1GB
    #[clap(long, value_parser = clap::value_parser!(u64).range(ONEGB..U32MAXPLUS), default_value_t = ONEGB)]
    pub chunk_size: u64,
}

pub fn run(args: Args, required_capacity: usize) -> Result<(), Box<dyn std::error::Error>> {
    let file_num_limit = get_file_limit();
    let meros = args.build.klmt.as_meros();

    let id_to_taxon_map = read_id_to_taxon_map(&args.taxo.id_to_taxon_map_filename)?;

    let taxonomy_filename = args
        .taxo
        .taxonomy_filename
        .unwrap_or(args.build.source.join("taxo.k2d"));

    let taxonomy = generate_taxonomy(
        &args.taxo.ncbi_taxonomy_directory,
        &taxonomy_filename,
        &id_to_taxon_map,
    )?;

    let value_bits = get_bits_for_taxid(
        args.build.requested_bits_for_taxid as usize,
        taxonomy.node_count() as f64,
    )
    .expect("more bits required for storing taxid");

    let capacity = required_capacity;
    let hash_config = HashConfig::<u32>::new(capacity, value_bits, 0, 0, 0);

    // 开始计时
    let start = Instant::now();

    let chunk_size = args.chunk_size as usize;

    let partition = (capacity + chunk_size - 1) / chunk_size;
    println!("start...");

    if partition >= file_num_limit {
        panic!("Exceeds File Number Limit");
    }

    let chunk_files = create_partition_files(partition, &args.chunk_dir, "chunk");
    let mut writers = create_partition_writers(&chunk_files);

    println!("chunk_size {:?}", format_bytes(chunk_size as f64));

    let source: PathBuf = args.build.source.clone();
    let fna_files = if source.is_file() {
        vec![source.to_string_lossy().to_string()]
    } else {
        find_library_fna_files(args.build.source)
    };

    for fna_file in &fna_files {
        convert_fna_to_k2_format(
            fna_file,
            meros,
            &taxonomy,
            &id_to_taxon_map,
            hash_config,
            &mut writers,
            chunk_size,
            args.build.threads as u32,
        );
    }

    let hash_filename = args
        .build
        .hashtable_filename
        .unwrap_or(source.join("hash.k2d"))
        .clone();
    let partition = chunk_files.len();
    for i in 0..partition {
        let mut chtm = CHTableMut::new(&hash_filename, hash_config, i, chunk_size)?;
        process_k2file(&chunk_files[i], &mut chtm, &taxonomy)?;
    }
    // 计算持续时间
    let duration = start.elapsed();
    // 打印运行时间
    println!("build k2 db took: {:?}", duration);

    let options_filename = args
        .build
        .options_filename
        .unwrap_or(source.clone().join("opts.k2d"));
    let idx_opts = IndexOptions::from_meros(meros);
    idx_opts.write_to_file(options_filename)?;

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
