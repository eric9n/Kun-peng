// 使用时需要引用模块路径
use clap::Parser;
use kr2r::args::{Build, Taxo};
use kr2r::compact_hash::{CHTableMut, HashConfig};
use kr2r::db::{convert_fna_to_k2_format, generate_taxonomy, get_bits_for_taxid, process_k2file};
use kr2r::utils::{
    create_partition_files, create_partition_writers, find_and_sort_files, find_library_fna_files,
    format_bytes, get_file_limit, read_id_to_taxon_map,
};
use kr2r::IndexOptions;
use std::path::PathBuf;
use std::time::Instant;

pub const U32MAXPLUS: u64 = u32::MAX as u64 + 2;
pub const ONEGB: u64 = 1073741824;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// 包含原始配置
    #[clap(flatten)]
    build: Build,

    #[clap(flatten)]
    taxo: Taxo,

    /// chunk directory
    #[clap(long)]
    chunk_dir: PathBuf,

    /// chunk size 1-4(GB)
    #[clap(long, value_parser = clap::value_parser!(u64).range(ONEGB..U32MAXPLUS), default_value_t = ONEGB)]
    chunk_size: u64,

    #[clap(long, default_value = "chunk")]
    chunk_prefix: String,

    /// process k2 file only
    #[clap(long)]
    only_k2: bool,

    /// process convert fna file to k2 format only
    #[clap(long)]
    only_convert: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let file_num_limit = get_file_limit();

    let args = Args::parse();
    let meros = args.build.as_meros();

    let id_to_taxon_map = read_id_to_taxon_map(&args.taxo.id_to_taxon_map_filename)?;

    let taxonomy = generate_taxonomy(
        &args.taxo.ncbi_taxonomy_directory,
        &args.taxo.taxonomy_filename,
        &id_to_taxon_map,
    )?;

    let value_bits = get_bits_for_taxid(
        args.build.requested_bits_for_taxid as usize,
        taxonomy.node_count() as f64,
    )
    .expect("more bits required for storing taxid");

    let capacity = args.build.required_capacity as usize;
    let hash_config = HashConfig::<u32>::new(capacity, value_bits, 0);

    // 开始计时
    let start = Instant::now();

    let chunk_size = args.chunk_size as usize;

    let chunk_files = if !args.only_k2 {
        let partition = (capacity + chunk_size - 1) / chunk_size;
        println!("start...");

        if partition >= file_num_limit {
            panic!("Exceeds File Number Limit");
        }

        let chunk_files = create_partition_files(partition, &args.chunk_dir, &args.chunk_prefix);
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
        println!("convert finished {:?}", &fna_files);

        chunk_files
    } else {
        find_and_sort_files(&args.chunk_dir, &args.chunk_prefix, ".k2", 0)?
    };
    println!("chunk_files {:?}", chunk_files);

    if !args.only_convert {
        let hash_filename = args.build.hashtable_filename.clone();
        let partition = chunk_files.len();
        for i in 0..partition {
            println!("chunk file {:?}", &chunk_files[i]);
            let mut chtm = CHTableMut::new(&hash_filename, hash_config, i, chunk_size)?;
            process_k2file(&chunk_files[i], &mut chtm, &taxonomy)?;
        }
        // 计算持续时间
        let duration = start.elapsed();
        // 打印运行时间
        println!("build k2 db took: {:?}", duration);

        let idx_opts = IndexOptions::from_meros(meros);
        idx_opts.write_to_file(args.build.options_filename)?;
    }

    Ok(())
}
