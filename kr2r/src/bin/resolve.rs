use clap::Parser;
use dashmap::{DashMap, DashSet};
use kr2r::compact_hash::{Compact, HashConfig, Row};
// use kr2r::iclassify::{resolve_tree, trim_pair_info};
use kr2r::readcounts::{TaxonCounters, TaxonCountersDash};
use kr2r::report::report_kraken_style;
use kr2r::taxonomy::Taxonomy;
use kr2r::utils::{find_and_sort_files, open_file};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Result, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

const BATCH_SIZE: usize = 8 * 1024 * 1024;

pub fn read_id_to_seq_map<P: AsRef<Path>>(
    filename: P,
) -> Result<DashMap<u32, (String, String, u32, Option<u32>)>> {
    let file = open_file(filename)?;
    let reader = BufReader::new(file);
    let id_map = DashMap::new();

    reader.lines().par_bridge().for_each(|line| {
        let line = line.expect("Could not read line");
        let parts: Vec<&str> = line.trim().split_whitespace().collect();
        if parts.len() >= 4 {
            // 解析序号为u32类型的键
            if let Ok(id) = parts[0].parse::<u32>() {
                // 第二列是序列标识符，直接作为字符串
                let seq_id = parts[1].to_string();
                let seq_size = parts[2].to_string();
                let count_parts: Vec<&str> = parts[3].split('|').collect();
                let kmer_count1 = count_parts[0].parse::<u32>().unwrap();
                let kmer_count2 = if count_parts.len() > 1 {
                    count_parts[1].parse::<u32>().map_or(None, |i| Some(i))
                } else {
                    None
                };
                id_map.insert(id, (seq_id, seq_size, kmer_count1, kmer_count2));
            }
        }
    });

    Ok(id_map)
}

fn generate_hit_string(
    count: u32,
    rows: &Vec<Row>,
    taxonomy: &Taxonomy,
    value_mask: usize,
    offset: u32,
) -> String {
    let mut result = Vec::new();
    let mut last_pos = 0;

    for row in rows {
        if row.kmer_id < offset || row.kmer_id >= offset + count {
            continue;
        }
        let adjusted_pos = row.kmer_id - offset;

        let value = row.value;
        let key = value.right(value_mask);
        let ext_code = taxonomy.nodes[key as usize].external_id;

        if last_pos == 0 && adjusted_pos > 0 {
            result.push((0, adjusted_pos)); // 在开始处添加0
        } else if adjusted_pos - last_pos > 1 {
            result.push((0, adjusted_pos - last_pos - 1)); // 在两个特定位置之间添加0
        }
        if let Some(last) = result.last_mut() {
            if last.0 == ext_code {
                last.1 += 1;
                last_pos = adjusted_pos;
                continue;
            }
        }

        // 添加当前key的计数
        result.push((ext_code, 1));
        last_pos = adjusted_pos;
    }

    // 填充尾随0
    if last_pos < count - 1 {
        if last_pos == 0 {
            result.push((0, count - last_pos));
        } else {
            result.push((0, count - last_pos - 1));
        }
    }

    result
        .iter()
        .map(|i| format!("{}:{}", i.0, i.1))
        .collect::<Vec<String>>()
        .join(" ")
}

pub fn trim_pair_info(id: &str) -> String {
    let sz = id.len();
    if sz <= 2 {
        return id.to_string();
    }
    if id.ends_with("/1") || id.ends_with("/2") {
        return id[0..sz - 2].to_string();
    }
    id.to_string()
}

// &HashMap<u32, u64>,
pub fn resolve_tree(
    hit_counts: &HashMap<u32, u64>,
    taxonomy: &Taxonomy,
    total_minimizers: usize,
    confidence_threshold: f64,
) -> u32 {
    let required_score = (confidence_threshold * total_minimizers as f64).ceil() as u64;

    let mut max_taxon = 0u32;
    let mut max_score = 0;

    for (&taxon, _) in hit_counts {
        let mut score = 0;

        for (&taxon2, &count2) in hit_counts {
            if taxonomy.is_a_ancestor_of_b(taxon2, taxon) {
                score += count2;
            }
        }

        if score > max_score {
            max_score = score;
            max_taxon = taxon;
        } else if score == max_score {
            max_taxon = taxonomy.lca(max_taxon, taxon);
        }
    }

    max_score = *hit_counts.get(&max_taxon).unwrap_or(&0);

    while max_taxon != 0 && max_score < required_score {
        max_score = hit_counts
            .iter()
            .filter(|(&taxon, _)| taxonomy.is_a_ancestor_of_b(max_taxon, taxon))
            .map(|(_, &count)| count)
            .sum();

        if max_score >= required_score {
            break;
        }
        max_taxon = taxonomy.nodes[max_taxon as usize].parent_id as u32;
    }

    max_taxon
}

pub fn add_hitlist_string(
    rows: &Vec<Row>,
    value_mask: usize,
    kmer_count1: u32,
    kmer_count2: Option<u32>,
    taxonomy: &Taxonomy,
) -> String {
    let result1 = generate_hit_string(kmer_count1, &rows, taxonomy, value_mask, 0);
    if let Some(count) = kmer_count2 {
        let result2 = generate_hit_string(count, &rows, taxonomy, value_mask, kmer_count1);
        format!("{} |:| {}", result1, result2)
    } else {
        format!("{}", result1)
    }
}

pub fn count_values(
    rows: &Vec<Row>,
    value_mask: usize,
    kmer_count1: u32,
) -> (HashMap<u32, u64>, TaxonCountersDash, usize) {
    let mut counts = HashMap::new();

    let mut hit_count: usize = 0;

    let mut last_row: Row = Row::new(0, 0, 0);
    let cur_taxon_counts = TaxonCountersDash::new();

    for row in rows {
        let value = row.value;
        let key = value.right(value_mask);
        *counts.entry(key).or_insert(0) += 1;

        // 如果切换到第2条seq,就重新计算
        if last_row.kmer_id < kmer_count1 && row.kmer_id > kmer_count1 {
            last_row = Row::new(0, 0, 0);
        }
        if !(last_row.value == value && row.kmer_id - last_row.kmer_id == 1) {
            cur_taxon_counts
                .entry(key as u64)
                .or_default()
                .add_kmer(value as u64);
            hit_count += 1;
        }

        last_row = *row;
    }

    (counts, cur_taxon_counts, hit_count)
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
    pub k2d_dir: PathBuf,

    /// chunk directory
    #[clap(long, value_parser, required = true)]
    pub chunk_dir: PathBuf,

    /// output file contains all unclassified sequence
    #[clap(long, value_parser, default_value_t = false)]
    pub full_output: bool,
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

fn process_batch<P: AsRef<Path>>(
    sample_file: P,
    args: &Args,
    taxonomy: &Taxonomy,
    id_map: &DashMap<u32, (String, String, u32, Option<u32>)>,
    writer: &Mutex<Box<dyn Write + Send>>,
    value_mask: usize,
) -> Result<(TaxonCountersDash, usize, DashSet<u32>)> {
    let file = open_file(sample_file)?;
    let mut reader = BufReader::new(file);
    let size = std::mem::size_of::<Row>();
    let mut batch_buffer = vec![0u8; size * BATCH_SIZE];

    let hit_counts = DashMap::new();
    let hit_seq_id_set = DashSet::new();
    let confidence_threshold = args.confidence_threshold;
    let minimum_hit_groups = args.minimum_hit_groups;

    while let Ok(bytes_read) = reader.read(&mut batch_buffer) {
        if bytes_read == 0 {
            break;
        } // 文件末尾

        // 处理读取的数据批次
        let slots_in_batch = bytes_read / size;
        let slots = unsafe {
            std::slice::from_raw_parts(batch_buffer.as_ptr() as *const Row, slots_in_batch)
        };

        slots.into_par_iter().for_each(|item| {
            let seq_id = item.seq_id;
            hit_seq_id_set.insert(seq_id);
            hit_counts
                .entry(seq_id)
                .or_insert_with(Vec::new)
                .push(*item)
        });
    }

    // let writer = Mutex::new(writer);
    let classify_counter = AtomicUsize::new(0);
    let cur_taxon_counts = TaxonCountersDash::new();

    hit_counts.into_par_iter().for_each(|(k, mut rows)| {
        if let Some(item) = id_map.get(&k) {
            rows.sort_unstable();
            let total_kmers: usize = item.2 as usize + item.3.unwrap_or(0) as usize;
            let dna_id = trim_pair_info(&item.0);
            let (counts, cur_counts, hit_groups) = count_values(&rows, value_mask, item.2);
            let hit_string = add_hitlist_string(&rows, value_mask, item.2, item.3, taxonomy);
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
            let clasify = if call > 0 {
                classify_counter.fetch_add(1, Ordering::SeqCst);
                cur_taxon_counts
                    .entry(call as u64)
                    .or_default()
                    .increment_read_count();

                "C"
            } else {
                "U"
            };
            // 使用锁来同步写入
            let output_line = format!(
                "{}\t{}\t{}\t{}\t{}\n",
                clasify, dna_id, ext_call, item.1, hit_string
            );
            let mut file = writer.lock().unwrap();
            file.write_all(output_line.as_bytes()).unwrap();
        }
    });
    Ok((
        cur_taxon_counts,
        classify_counter.load(Ordering::SeqCst),
        hit_seq_id_set,
    ))
}

pub fn run(args: Args) -> Result<()> {
    let k2d_dir = &args.k2d_dir;
    let taxonomy_filename = k2d_dir.join("taxo.k2d");
    let taxo = Taxonomy::from_file(taxonomy_filename)?;

    let sample_files = find_and_sort_files(&args.chunk_dir, "sample_file", ".bin")?;
    let sample_id_files = find_and_sort_files(&args.chunk_dir, "sample_id", ".map")?;

    let partition = sample_files.len();
    let hash_config = HashConfig::from_hash_header(&args.k2d_dir.join("hash_config.k2d"))?;
    let value_mask = hash_config.value_mask;

    let mut total_taxon_counts = TaxonCounters::new();
    let mut total_seqs = 0;
    let mut total_unclassified = 0;

    // 开始计时
    let start = Instant::now();
    println!("start...");

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
        let writer = Mutex::new(writer);
        let (thread_taxon_counts, thread_classified, hit_seq_set) = process_batch::<&PathBuf>(
            sample_file,
            &args,
            &taxo,
            &sample_id_map,
            &writer,
            value_mask,
        )?;

        if args.full_output {
            sample_id_map
                .iter()
                .filter(|item| !hit_seq_set.contains(item.key()))
                .for_each(|item| {
                    let dna_id = trim_pair_info(&item.0);
                    let hit_string = add_hitlist_string(&vec![], value_mask, item.2, item.3, &taxo);
                    let output_line = format!("U\t{}\t0\t{}\t{}\n", dna_id, item.1, hit_string);
                    let mut file = writer.lock().unwrap();
                    file.write_all(output_line.as_bytes()).unwrap();
                });
        }

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
        let filename = output.join("output.kreport2");
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

    // 计算持续时间
    let duration = start.elapsed();
    // 打印运行时间
    println!("resolve took: {:?}", duration);

    Ok(())
}

#[allow(dead_code)]
fn main() {
    let args = Args::parse();
    if let Err(e) = run(args) {
        eprintln!("Application error: {}", e);
    }
}
