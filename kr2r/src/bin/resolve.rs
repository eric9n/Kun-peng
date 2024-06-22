use clap::Parser;
use dashmap::{DashMap, DashSet};
use kr2r::classify::process_hitgroup;
use kr2r::compact_hash::{HashConfig, Row};
use kr2r::readcounts::{TaxonCounters, TaxonCountersDash};
use kr2r::report::report_kraken_style;
use kr2r::taxonomy::Taxonomy;
use kr2r::utils::{find_and_sort_files, open_file};
use kr2r::HitGroup;
use rayon::prelude::*;
use seqkmer::{trim_pair_info, OptionPair};
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
) -> Result<DashMap<u32, (String, String, usize, Option<usize>)>> {
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
                let kmer_count1 = count_parts[0].parse::<usize>().unwrap();
                let kmer_count2 = if count_parts.len() > 1 {
                    count_parts[1].parse::<usize>().map_or(None, |i| Some(i))
                } else {
                    None
                };
                id_map.insert(id, (seq_id, seq_size, kmer_count1, kmer_count2));
            }
        }
    });

    Ok(id_map)
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
    id_map: &DashMap<u32, (String, String, usize, Option<usize>)>,
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
            let dna_id = trim_pair_info(&item.0);
            let range = OptionPair::from(((0, item.2), item.3.map(|size| (item.2, size + item.2))));
            let hits = HitGroup::new(rows, range);

            let hit_data = process_hitgroup(
                &hits,
                taxonomy,
                &classify_counter,
                hits.required_score(confidence_threshold),
                minimum_hit_groups,
                value_mask,
            );
            // let (counts, cur_counts, hit_groups) = count_values(&rows, value_mask, item.2);
            // let hit_string = add_hitlist_string(&rows, value_mask, item.2, item.3, taxonomy);
            // let require_score = (confidence_threshold * total_kmers as f64).ceil() as u64;
            // let mut call = resolve_tree(&counts, taxonomy, require_score);
            // if call > 0 && hit_groups < minimum_hit_groups {
            //     call = 0;
            // };

            hit_data.3.iter().for_each(|(key, value)| {
                cur_taxon_counts
                    .entry(*key)
                    .or_default()
                    .merge(value)
                    .unwrap();
            });

            // let ext_call = taxonomy.nodes[call as usize].external_id;
            // let clasify = if call > 0 {
            //     classify_counter.fetch_add(1, Ordering::SeqCst);
            //     cur_taxon_counts
            //         .entry(call as u64)
            //         .or_default()
            //         .increment_read_count();

            //     "C"
            // } else {
            //     "U"
            // };
            // 使用锁来同步写入
            let output_line = format!(
                "{}\t{}\t{}\t{}\t{}\n",
                hit_data.0, dna_id, hit_data.1, item.1, hit_data.2
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
    println!("resolve start...");

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
                    let output_line = format!(
                        "U\t{}\t0\t{}\t{}\n",
                        dna_id,
                        item.1,
                        if item.3.is_none() { "" } else { " |:| " }
                    );

                    let mut file = writer.lock().unwrap();
                    file.write_all(output_line.as_bytes()).unwrap();
                });
        }

        let mut sample_taxon_counts: HashMap<
            u64,
            kr2r::readcounts::ReadCounts<hyperloglogplus::HyperLogLogPlus<u64, kr2r::KBuildHasher>>,
        > = HashMap::new();
        thread_taxon_counts.iter().for_each(|entry| {
            total_taxon_counts
                .entry(*entry.key())
                .or_default()
                .merge(&entry.value())
                .unwrap();
            sample_taxon_counts
                .entry(*entry.key())
                .or_default()
                .merge(&entry.value())
                .unwrap();
        });
        if let Some(output) = &args.kraken_output_dir {
            let filename = output.join(format!("output_{}.kreport2", i + 1));
            report_kraken_style(
                filename,
                args.report_zero_counts,
                args.report_kmer_data,
                &taxo,
                &sample_taxon_counts,
                thread_sequences as u64,
                (thread_sequences - thread_classified) as u64,
            )?;
        }

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
