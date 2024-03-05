// 使用时需要引用模块路径
use clap::Parser;
use kr2r::table::CompactHashTable;
use kr2r::utils::{expand_spaced_seed_mask, find_library_fna_files, read_id_to_taxon_map};
use kr2r::{construct_seed_template, parse_binary, Meros};
use kr2r::{
    IndexOptions, BITS_PER_CHAR, DEFAULT_KMER_LENGTH, DEFAULT_MINIMIZER_LENGTH,
    DEFAULT_MINIMIZER_SPACES, DEFAULT_TOGGLE_MASK,
};
// use std::collections::HashMap;
use kr2r::db::{generate_taxonomy, get_bits_for_taxid};
use kr2r::tansform::{convert_fna_to_k2_format, process_sequence};
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;
extern crate libc;

use libc::{getrlimit, rlimit, RLIMIT_NOFILE};

#[derive(Parser, Debug, Clone)]
#[clap(version, about = "build database")]
struct Build {
    /// build database directory or file
    #[arg(long, required = true)]
    source: PathBuf,

    /// Kraken 2 hash table filename
    #[clap(short = 'H', required = true)]
    hashtable_filename: PathBuf,

    /// Kraken 2 taxonomy filename
    #[clap(short = 't', required = true)]
    taxonomy_filename: PathBuf,

    /// Sequence ID to taxon map filename
    #[clap(short = 'm', required = true)]
    id_to_taxon_map_filename: PathBuf,

    /// Kraken 2 options filename
    #[clap(short = 'o', required = true)]
    options_filename: PathBuf,

    /// NCBI taxonomy directory name
    #[clap(short, long, required = true)]
    ncbi_taxonomy_directory: PathBuf,

    /// Set length of k-mers, k must be positive integer, k=35, k cannot be less than l
    #[clap(short, long, value_parser = clap::value_parser!(u64).range(1..), default_value_t = DEFAULT_KMER_LENGTH)]
    k_mer: u64,

    /// Set length of minimizers, 1 <= l <= 31
    #[clap(short, long, value_parser = clap::value_parser!(u8).range(1..=31), default_value_t = DEFAULT_MINIMIZER_LENGTH)]
    l_mer: u8,

    /// Bit storage requested for taxid 0 <= r < 31
    #[clap(short, long, value_parser = clap::value_parser!(u8).range(0..31), default_value_t = 0)]
    requested_bits_for_taxid: u8,

    /// Minimizer ordering toggle mask
    #[clap(short = 'T', long, default_value_t = DEFAULT_TOGGLE_MASK)]
    toggle_mask: u64,

    /// Number of characters in minimizer that are ignored in comparisons
    #[clap(long, default_value_t = DEFAULT_MINIMIZER_SPACES)]
    minimizer_spaces: u8,

    // /// Name of Kraken 2 database
    // #[arg(short, long = "db")]
    // database: PathBuf,
    #[arg(short = 'c', long, required = true)]
    required_capacity: u64,

    /// Number of threads
    #[clap(short = 'p', long, default_value_t = 4)]
    threads: usize,
}

impl Build {
    pub fn as_meros(&self, spaced_seed_mask: u64) -> Meros {
        Meros::new(
            self.k_mer as usize,
            self.l_mer as usize,
            Some(spaced_seed_mask),
            Some(self.toggle_mask),
            Some(0),
        )
    }
}

fn get_limit() -> usize {
    let mut limits = rlimit {
        rlim_cur: 0, // 当前（软）限制
        rlim_max: 0, // 最大（硬）限制
    };

    // 使用unsafe块调用getrlimit，因为这是一个外部C函数
    let result = unsafe { getrlimit(RLIMIT_NOFILE, &mut limits) };

    if result == 0 {
        // 如果成功，返回当前软限制转换为usize
        limits.rlim_cur as usize
    } else {
        // 如果失败，输出错误并可能返回一个默认值或panic
        eprintln!("Failed to get file limit");
        // 这里你可以根据需要返回一个默认值，或者调用panic!来终止程序
        // 这里返回0作为示例，但你可能想根据你的应用逻辑来选择更合适的处理方式
        0
    }
}

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

fn create_partition_writers(partition: usize, base_path: &str) -> io::Result<Vec<BufWriter<File>>> {
    let mut writers = Vec::with_capacity(partition); // 预先分配足够的空间

    for i in 0..partition {
        // 构建每个分区文件的完整路径和名称
        let mut file_path = PathBuf::from(base_path);
        file_path.push(format!("chunk_{}.k2", i));

        // 尝试创建文件，如果失败则直接返回错误
        let file = OpenOptions::new()
            .write(true)
            .append(true) // 确保以追加模式打开文件
            .create(true) // 如果文件不存在，则创建
            .open(file_path)
            .unwrap();

        let writer = BufWriter::new(file);
        writers.push(writer);

        // 可选：打印日志或信息，确认文件已创建
        // println!("Created partition file: {:?}", file_path);
    }

    Ok(writers)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let file_num_limit = get_limit();

    let args = Build::parse();
    let seed = construct_seed_template(
        args.l_mer.clone() as usize,
        args.minimizer_spaces.clone() as usize,
    );
    let space_seed_mask = parse_binary(&seed).unwrap();
    let space_seed_mask = expand_spaced_seed_mask(space_seed_mask, BITS_PER_CHAR as u64);

    let meros = args.as_meros(space_seed_mask);

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
    // 1073741824
    // 1732968825
    let capacity = args.required_capacity as usize;

    let chunk_size = 1073741824usize;

    // 使用整数数学向上取整
    let partition = (capacity + chunk_size - 1) / chunk_size;
    if partition >= file_num_limit {
        panic!("Exceeds File Number Limit");
    }

    // let mut writers = create_partition_writers(partition, ".")?;

    println!("chunk_size {:?}", format_bytes(chunk_size as f64));

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
        // convert_fna_to_k2_format(
        //     fna_file,
        //     meros,
        //     &taxonomy,
        //     &id_to_taxon_map,
        //     args.threads as u32,
        //     capacity,
        //     value_bits,
        //     &mut writers,
        //     chunk_size,
        // );
    }

    let chunk_files = vec!["chunk_0.k2", "chunk_1.k2"];
    let hash_filename = args.hashtable_filename.clone();
    for i in 0..2 {
        println!("chunk_files[i] {:?}", chunk_files[i]);
        let mut chtm = CompactHashTable::new(&hash_filename, capacity, value_bits, i, chunk_size)?;
        process_sequence(chunk_files[i], &mut chtm, &taxonomy, value_bits)?;
    }
    // 计算持续时间
    let duration = start.elapsed();
    // 打印运行时间
    println!("build db took: {:?}", duration);

    let idx_opts = IndexOptions::new(
        args.k_mer as usize,
        args.l_mer as usize,
        space_seed_mask,
        args.toggle_mask,
        true,
        0,
    );
    idx_opts.write_to_file(args.options_filename)?;

    Ok(())
}
