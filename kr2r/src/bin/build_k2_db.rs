// 使用时需要引用模块路径
use clap::Parser;
use kr2r::args::Build;
use kr2r::compact_hash::{CHTableMut, HashConfig};
use kr2r::db::{
    convert_fna_to_k2_format, create_partition_files, generate_taxonomy, get_bits_for_taxid,
    process_k2file,
};
use kr2r::db::{create_partition_writers, get_file_limit};
use kr2r::utils::{find_library_fna_files, read_id_to_taxon_map};
use kr2r::IndexOptions;
use std::path::PathBuf;
use std::time::Instant;

fn format_bytes(size: f64) -> String {
    let suffixes = ["B", "kB", "MB", "GB", "TB", "PB", "EB"];
    let mut size = size;
    let mut current_suffix = &suffixes[0];

    for suffix in &suffixes[1..] {
        if size >= 1024.0 {
            current_suffix = suffix;
            size /= 1024.0;
        } else {
            break;
        }
    }

    format!("{:.2}{}", size, current_suffix)
}

pub const U32MAX: u64 = u32::MAX as u64;
pub const ONEGB: u64 = 1073741824;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// 包含原始配置
    #[clap(flatten)]
    build: Build,

    /// 新增的参数
    #[clap(long)]
    chunk_dir: PathBuf,

    #[clap(long, value_parser = clap::value_parser!(u64).range(ONEGB..U32MAX), default_value_t = ONEGB as usize)]
    chunk_size: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let file_num_limit = get_file_limit();

    let args = Args::parse();
    let meros = args.build.as_meros();

    let id_to_taxon_map = read_id_to_taxon_map(&args.build.id_to_taxon_map_filename)?;

    let taxonomy = generate_taxonomy(
        &args.build.ncbi_taxonomy_directory,
        &args.build.taxonomy_filename,
        &id_to_taxon_map,
    )?;

    let value_bits = get_bits_for_taxid(
        args.build.requested_bits_for_taxid as usize,
        taxonomy.node_count() as f64,
    )
    .expect("more bits required for storing taxid");

    let capacity = args.build.required_capacity as usize;

    let chunk_size = args.chunk_size;

    // 使用整数数学向上取整
    let partition = (capacity + chunk_size - 1) / chunk_size;
    if partition >= file_num_limit {
        panic!("Exceeds File Number Limit");
    }

    let chunk_files = create_partition_files(partition, &args.chunk_dir);
    let mut writers = create_partition_writers(&chunk_files);

    println!("chunk_size {:?}", format_bytes(chunk_size as f64));

    let source: PathBuf = args.build.source.clone();
    let fna_files = if source.is_file() {
        vec![source.to_string_lossy().to_string()]
    } else {
        find_library_fna_files(args.build.source)
    };

    println!("start...");
    // 开始计时
    let start = Instant::now();
    let hash_config = HashConfig::new(capacity, value_bits, 0);
    for fna_file in fna_files {
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

    let hash_filename = args.build.hashtable_filename.clone();
    for i in 0..2 {
        let mut chtm = CHTableMut::new(&hash_filename, hash_config, i, chunk_size)?;
        process_k2file(&chunk_files[i], &mut chtm, &taxonomy)?;
    }
    // 计算持续时间
    let duration = start.elapsed();
    // 打印运行时间
    println!("build k2 db took: {:?}", duration);

    let idx_opts = IndexOptions::from_meros(meros);
    idx_opts.write_to_file(args.build.options_filename)?;

    Ok(())
}
