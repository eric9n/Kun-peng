use clap::Parser;
use kr2r::compact_hash::Compact;
use kr2r::iclassify::{count_values, resolve_tree};
use kr2r::taxonomy::Taxonomy;
use kr2r::utils::find_and_sort_files;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Result, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

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

#[derive(Parser, Debug, Clone)]
#[clap(
    version,
    about = "resolve taxonomy tree",
    long_about = "resolve taxonomy tree"
)]
struct Args {
    // chunk directory
    #[clap(long)]
    chunk_dir: PathBuf,

    /// The file path for the Kraken 2 taxonomy.
    #[clap(short = 't', long = "taxonomy-filename", value_parser, required = true)]
    taxonomy_filename: String,

    /// Confidence score threshold, default is 0.0.
    #[clap(
        short = 'T',
        long = "confidence-threshold",
        value_parser,
        default_value_t = 0.0
    )]
    confidence_threshold: f64,

    /// The minimum number of hit groups needed for a call.
    #[clap(
        short = 'g',
        long = "minimum-hit-groups",
        value_parser,
        default_value_t = 2
    )]
    minimum_hit_groups: usize,

    /// File path for outputting normal Kraken output.
    #[clap(long = "output-dir", value_parser)]
    kraken_output_dir: Option<PathBuf>,
}

const BATCH_SIZE: usize = 80 * 1024;
use dashmap::DashMap;

fn process_batch<P: AsRef<Path>, B: Compact>(
    sample_file: P,
    args: &Args,
    taxonomy: &Taxonomy,
    id_map: DashMap<u32, (String, usize)>,
    writer: Box<dyn Write + Send>,
) -> Result<()> {
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
            let taxid = item.left(0).to_u32();
            let seq_id = item.right(0).to_u32();
            hit_counts
                .entry(seq_id)
                .or_insert_with(Vec::new)
                .push(taxid)
        });
    }

    let writer = Mutex::new(writer);

    hit_counts.into_par_iter().for_each(|(k, v)| {
        if let Some(item) = id_map.get(&k) {
            let total_kmers = item.1;
            let minimizer_hit_groups = v.len();
            let mut call = resolve_tree(
                &count_values(v),
                taxonomy,
                total_kmers,
                confidence_threshold,
            );
            if call > 0 && minimizer_hit_groups < minimum_hit_groups {
                call = 0;
            };

            let ext_call = taxonomy.nodes[call as usize].external_id;
            let classify = if call > 0 { "C" } else { "U" };
            let output_line = format!("{}\t{}\t{}\n", classify, item.0, ext_call);
            // 使用锁来同步写入
            let mut file = writer.lock().unwrap();
            file.write_all(output_line.as_bytes()).unwrap();
        }
    });
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    let taxo = Taxonomy::from_file(&args.taxonomy_filename)?;

    let sample_files = find_and_sort_files(&args.chunk_dir, "sample_file", ".bin", 1)?;
    let sample_id_files = find_and_sort_files(&args.chunk_dir, "sample_id", ".map", 1)?;

    let partition = sample_files.len();

    for i in 0..partition {
        let sample_file = &sample_files[i];
        let sample_id_map = read_id_to_seq_map(&sample_id_files[i])?;
        let writer: Box<dyn Write + Send> = match &args.kraken_output_dir {
            Some(ref file_path) => {
                let filename = file_path.join(format!("output_{}.txt", i));
                let file = File::create(filename)?;
                Box::new(BufWriter::new(file)) as Box<dyn Write + Send>
            }
            None => Box::new(io::stdout()) as Box<dyn Write + Send>,
        };
        process_batch::<&PathBuf, u64>(sample_file, &args, &taxo, sample_id_map, writer)?;
    }
    Ok(())
}
