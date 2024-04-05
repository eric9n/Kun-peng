use clap::Parser;
use dashmap::DashMap;
use kr2r::compact_hash::{Compact, HashConfig};
use kr2r::iclassify::{resolve_tree, trim_pair_info};
use kr2r::readcounts::{TaxonCounters, TaxonCountersDash};
use kr2r::report::report_kraken_style;
use kr2r::taxonomy::Taxonomy;
use kr2r::utils::find_and_sort_files;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Result, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

const BATCH_SIZE: usize = 8 * 1024 * 1024;

pub fn read_id_to_seq_map<P: AsRef<Path>>(filename: P) -> Result<DashMap<u32, (String, usize)>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let id_map = DashMap::new();

    reader.lines().par_bridge().for_each(|line| {
        let line = line.expect("Could not read line");
        let parts: Vec<&str> = line.trim().split_whitespace().collect();
        if parts.len() >= 3 {
            // 解析序号为u32类型的键
            if let Ok(id) = parts[0].parse::<u32>() {
                // 第二列是序列标识符，直接作为字符串
                let seq_id = parts[1].to_string();
                if let Ok(count) = parts[2].parse::<usize>() {
                    // 插入到DashMap中
                    id_map.insert(id, (seq_id, count));
                }
            }
        }
    });

    Ok(id_map)
}

pub fn count_values(
    vec: Vec<u32>,
    value_mask: usize,
) -> (HashMap<u32, u64>, TaxonCountersDash, usize) {
    let mut counts = HashMap::new();

    let mut unique_elements = HashSet::new();

    let cur_taxon_counts = TaxonCountersDash::new();

    for value in vec {
        // 使用entry API处理计数
        // entry返回的是一个Entry枚举，它代表了可能存在也可能不存在的值
        // or_insert方法在键不存在时插入默认值（在这里是0）
        // 然后无论哪种情况，我们都对计数器加1
        let key = value.right(value_mask);
        *counts.entry(key).or_insert(0) += 1;
        if !unique_elements.contains(&value) {
            cur_taxon_counts
                .entry(key as u64)
                .or_default()
                .add_kmer(value as u64);
        }
        unique_elements.insert(value);
    }

    (counts, cur_taxon_counts, unique_elements.len())
}

#[derive(Parser, Debug, Clone)]
#[clap(
    version,
    about = "resolve taxonomy tree",
    long_about = "resolve taxonomy tree"
)]
pub struct Args {
    /// database hash chunk directory and other files
    #[clap(long)]
    pub hash_dir: PathBuf,

    // chunk directory
    #[clap(long, value_parser, required = true)]
    pub chunk_dir: PathBuf,

    // /// The file path for the Kraken 2 index.
    // #[clap(short = 'H', long = "index-filename", value_parser, required = true)]
    // index_filename: PathBuf,

    // /// The file path for the Kraken 2 taxonomy.
    // #[clap(short = 't', long = "taxonomy-filename", value_parser, required = true)]
    // taxonomy_filename: String,
    /// Confidence score threshold, default is 0.0.
    #[clap(
        short = 'T',
        long = "confidence-threshold",
        value_parser,
        default_value_t = 0.0
    )]
    pub confidence_threshold: f64,

    /// In comb. w/ -R, provide minimizer information in report
    #[clap(short = 'K', long, value_parser, default_value_t = false)]
    pub report_kmer_data: bool,

    /// In comb. w/ -R, report taxa w/ 0 count
    #[clap(short = 'z', long, value_parser, default_value_t = false)]
    pub report_zero_counts: bool,

    /// The minimum number of hit groups needed for a call.
    #[clap(
        short = 'g',
        long = "minimum-hit-groups",
        value_parser,
        default_value_t = 2
    )]
    pub minimum_hit_groups: usize,

    /// 批量处理大小 default: 8MB
    #[clap(long, default_value_t = BATCH_SIZE)]
    pub batch_size: usize,

    /// File path for outputting normal Kraken output.
    #[clap(long = "output-dir", value_parser)]
    pub kraken_output_dir: Option<PathBuf>,
}

fn process_batch<P: AsRef<Path>, B: Compact>(
    sample_file: P,
    args: &Args,
    taxonomy: &Taxonomy,
    id_map: DashMap<u32, (String, usize)>,
    writer: Box<dyn Write + Send>,
    value_mask: usize,
) -> Result<(TaxonCountersDash, usize)> {
    let file = File::open(sample_file)?;
    let mut reader = BufReader::new(file);
    let size = std::mem::size_of::<B>();
    let mut batch_buffer = vec![0u8; size * BATCH_SIZE];

    let hit_counts = DashMap::new();
    let confidence_threshold = args.confidence_threshold;
    let minimum_hit_groups = args.minimum_hit_groups;

    while let Ok(bytes_read) = reader.read(&mut batch_buffer) {
        if bytes_read == 0 {
            break;
        } // 文件末尾

        // 处理读取的数据批次
        let slots_in_batch = bytes_read / size;

        let slots = unsafe {
            std::slice::from_raw_parts(batch_buffer.as_ptr() as *const B, slots_in_batch)
        };

        slots.into_par_iter().for_each(|item| {
            let cell = item.left(0).to_u32();
            let seq_id = item.right(0).to_u32();
            hit_counts.entry(seq_id).or_insert_with(Vec::new).push(cell)
        });
    }

    let writer = Mutex::new(writer);
    let classify_counter = AtomicUsize::new(0);
    let cur_taxon_counts = TaxonCountersDash::new();

    hit_counts.into_par_iter().for_each(|(k, cells)| {
        if let Some(item) = id_map.get(&k) {
            let total_kmers: usize = item.1;
            let dna_id = trim_pair_info(&item.0);
            let (counts, cur_counts, hit_groups) = count_values(cells, value_mask);
            let mut call = resolve_tree(&counts, taxonomy, total_kmers, confidence_threshold);
            if call > 0 && hit_groups < minimum_hit_groups {
                call = 0;
            };
            cur_counts.iter().for_each(|entry| {
                cur_taxon_counts
                    .entry(*entry.key())
                    .or_default()
                    .merge(entry.value())
                    .unwrap();
            });

            let ext_call = taxonomy.nodes[call as usize].external_id;
            if call > 0 {
                let output_line = format!("C\t{}\t{}\n", dna_id, ext_call);
                // 使用锁来同步写入
                let mut file = writer.lock().unwrap();
                file.write_all(output_line.as_bytes()).unwrap();
                classify_counter.fetch_add(1, Ordering::SeqCst);
                cur_taxon_counts
                    .entry(call as u64)
                    .or_default()
                    .increment_read_count();
            }
        }
    });
    Ok((cur_taxon_counts, classify_counter.load(Ordering::SeqCst)))
}

pub fn run(args: Args) -> Result<()> {
    let hash_dir = &args.hash_dir;
    let taxonomy_filename = hash_dir.join("taxo.k2d");
    let taxo = Taxonomy::from_file(taxonomy_filename)?;

    let sample_files = find_and_sort_files(&args.chunk_dir, "sample_file", ".bin")?;
    let sample_id_files = find_and_sort_files(&args.chunk_dir, "sample_id", ".map")?;

    let partition = sample_files.len();
    let hash_config = HashConfig::<u32>::from_hash_header(&args.hash_dir.join("hash_config.k2d"))?;
    let value_mask = hash_config.value_mask;

    let mut total_taxon_counts = TaxonCounters::new();
    let mut total_seqs = 0;
    let mut total_unclassified = 0;

    for i in 0..partition {
        let sample_file = &sample_files[i];
        let sample_id_map = read_id_to_seq_map(&sample_id_files[i])?;
        let thread_sequences = sample_id_map.len();
        let writer: Box<dyn Write + Send> = match &args.kraken_output_dir {
            Some(ref file_path) => {
                let filename = file_path.join(format!("output_{}.txt", i + 1));
                let file = File::create(filename)?;
                Box::new(BufWriter::new(file)) as Box<dyn Write + Send>
            }
            None => Box::new(io::stdout()) as Box<dyn Write + Send>,
        };
        let (thread_taxon_counts, thread_classified) = process_batch::<&PathBuf, u64>(
            sample_file,
            &args,
            &taxo,
            sample_id_map,
            writer,
            value_mask,
        )?;
        thread_taxon_counts.iter().for_each(|entry| {
            total_taxon_counts
                .entry(*entry.key())
                .or_default()
                .merge(&entry.value())
                .unwrap();
        });
        total_seqs += thread_sequences;
        total_unclassified += thread_sequences - thread_classified;
    }

    if let Some(output) = &args.kraken_output_dir {
        let filename = output.join("output.krepot2");
        report_kraken_style(
            filename,
            args.report_zero_counts,
            args.report_kmer_data,
            &taxo,
            &total_taxon_counts,
            total_seqs as u64,
            total_unclassified as u64,
        )?;
    }

    Ok(())
}

#[allow(dead_code)]
fn main() {
    let args = Args::parse();
    if let Err(e) = run(args) {
        eprintln!("Application error: {}", e);
    }
}
