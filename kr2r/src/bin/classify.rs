use clap::Parser;
use kr2r::classify::process_hitgroup;
use kr2r::compact_hash::{CHTable, Compact, HashConfig, Row};
use kr2r::readcounts::{TaxonCounters, TaxonCountersDash};
use kr2r::report::report_kraken_style;
use kr2r::taxonomy::Taxonomy;
use kr2r::utils::{create_sample_file, find_and_sort_files, get_lastest_file_index};
use kr2r::{HitGroup, IndexOptions};
use seqkmer::{read_parallel, Base, FastxReader, Meros, MinimizerIterator, OptionPair, Reader};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::io::{Error, ErrorKind, Result};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

#[derive(Parser, Debug, Clone)]
#[clap(
    version,
    about = "Directly load all hash tables for classification annotation",
    long_about = "Directly load all hash tables for classification annotation"
)]
pub struct Args {
    /// database hash chunk directory and other files
    #[clap(long)]
    pub k2d_dir: PathBuf,

    // /// The file path for the Kraken 2 options.
    // #[clap(short = 'o', long = "options-filename", value_parser, required = true)]
    // options_filename: String,
    /// Enable paired-end processing.
    #[clap(short = 'P', long = "paired-end-processing", action)]
    pub paired_end_processing: bool,

    /// Process pairs with mates in the same file.
    #[clap(short = 'S', long = "single-file-pairs", action)]
    pub single_file_pairs: bool,

    /// Minimum quality score for FASTQ data.
    #[clap(
        short = 'Q',
        long = "minimum-quality-score",
        value_parser,
        default_value_t = 0
    )]
    pub minimum_quality_score: i32,

    /// Confidence score threshold.
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

    /// The number of threads to use.
    #[clap(short = 'p', long = "num-threads", value_parser, default_value_t = 16)]
    pub num_threads: i32,

    /// File path for outputting normal Kraken output.
    #[clap(long = "output-dir", value_parser)]
    pub kraken_output_dir: Option<PathBuf>,

    /// A list of input file paths (FASTA/FASTQ) to be processed by the classify program.
    // #[clap(short = 'F', long = "files")]
    pub input_files: Vec<String>,
}

fn process_seq(
    rows: &mut Vec<Row>,
    m_iter: &mut MinimizerIterator,
    hash_config: &HashConfig,
    chtable: &CHTable,
    offset: usize,
) -> usize {
    let chunk_size = hash_config.hash_capacity;
    let value_bits = hash_config.value_bits;
    let data: Vec<(usize, u64)> = m_iter.collect();
    for (sort, hash_key) in data {
        let (idx, compacted) = hash_config.compact(hash_key);
        let partition_index = idx / chunk_size;
        let index = idx % chunk_size;

        let taxid = chtable.get_from_page(index, compacted, partition_index + 1);
        if taxid > 0 {
            let high = u32::combined(compacted, taxid, value_bits);
            let row = Row::new(high, 0, sort as u32 + 1 + offset as u32);
            rows.push(row);
        }
    }
    m_iter.size + offset
}

fn process_record(
    marker: &mut Base<MinimizerIterator>,
    args: &Args,
    taxonomy: &Taxonomy,
    chtable: &CHTable,
    hash_config: &HashConfig,
    cur_taxon_counts: &TaxonCountersDash,
    classify_counter: &AtomicUsize,
) -> String {
    let id = &marker.header.id.clone();
    let rows: Vec<Row> = marker
        .fold(|rows, m_iter, offset| process_seq(rows, m_iter, &hash_config, chtable, offset));

    let hits = HitGroup::new(rows, marker.range());

    let seq_len_str = marker.fmt_seq_size();

    let required_score = hits.required_score(args.confidence_threshold);
    let hit_data = process_hitgroup(
        &hits,
        taxonomy,
        classify_counter,
        required_score,
        args.minimum_hit_groups,
        hash_config.value_mask,
    );

    hit_data.3.iter().for_each(|(key, value)| {
        cur_taxon_counts
            .entry(*key)
            .or_default()
            .merge(value)
            .unwrap();
    });
    format!(
        "{}\t{}\t{}\t{}\t{}\n",
        hit_data.0, id, hit_data.1, seq_len_str, hit_data.2
    )
}

fn process_fastx_file<R>(
    args: &Args,
    meros: Meros,
    hash_config: HashConfig,
    file_index: usize,
    reader: &mut R,
    chtable: &CHTable,
    taxonomy: &Taxonomy,
    total_taxon_counts: &mut TaxonCounters,
) -> io::Result<(usize, usize)>
where
    R: Reader,
{
    let mut writer: Box<dyn Write + Send> = match &args.kraken_output_dir {
        Some(ref file_path) => {
            let filename = file_path.join(format!("output_{}.txt", file_index));
            let file = File::create(filename)?;
            Box::new(BufWriter::new(file)) as Box<dyn Write + Send>
        }
        None => Box::new(io::stdout()) as Box<dyn Write + Send>,
    };

    let cur_taxon_counts = TaxonCountersDash::new();

    let seq_counter = AtomicUsize::new(0);
    let classify_counter = AtomicUsize::new(0);

    let _ = read_parallel(
        reader,
        args.num_threads as usize,
        &meros,
        |seqs| {
            let mut buffer = String::new();
            for record in seqs {
                seq_counter.fetch_add(1, Ordering::SeqCst);
                let output_line = process_record(
                    record,
                    args,
                    taxonomy,
                    chtable,
                    &hash_config,
                    &cur_taxon_counts,
                    &classify_counter,
                );
                buffer.push_str(&output_line);
            }

            Some(buffer)
        },
        |dataset| {
            while let Some(Some(res)) = dataset.next() {
                writer
                    .write_all(res.as_bytes())
                    .expect("Failed to write date to file");
            }
        },
    );

    let mut sample_taxon_counts: HashMap<
        u64,
        kr2r::readcounts::ReadCounts<hyperloglogplus::HyperLogLogPlus<u64, kr2r::KBuildHasher>>,
    > = HashMap::new();
    cur_taxon_counts.iter().for_each(|entry| {
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

    let thread_sequences = seq_counter.load(Ordering::SeqCst);
    let thread_classified = classify_counter.load(Ordering::SeqCst);
    if let Some(output) = &args.kraken_output_dir {
        let filename = output.join(format!("output_{}.kreport2", file_index));
        report_kraken_style(
            filename,
            args.report_zero_counts,
            args.report_kmer_data,
            &taxonomy,
            &sample_taxon_counts,
            thread_sequences as u64,
            (thread_sequences - thread_classified) as u64,
        )?;
    }

    Ok((thread_sequences, thread_sequences - thread_classified))
}

fn process_files(
    args: Args,
    meros: Meros,
    hash_config: HashConfig,
    chtable: &CHTable,
    taxonomy: &Taxonomy,
) -> Result<()> {
    let (mut file_index, mut file_writer) = if let Some(out_dir) = &args.kraken_output_dir {
        let file_path = out_dir.join("sample_file.map");
        let file_writer = create_sample_file(&file_path);
        let file_index = get_lastest_file_index(&file_path)?;
        (file_index, Box::new(file_writer) as Box<dyn Write + Send>)
    } else {
        (0, Box::new(io::stdout()) as Box<dyn Write + Send>)
    };

    let mut process_funcs = |files: Vec<&[String]>| -> Result<()> {
        let file_bits = (((files.len() + file_index) as f64).log2().ceil() as usize).max(1);
        if file_bits > hash_config.value_bits {
            panic!("The number of files is too large to process.");
        }

        let mut total_taxon_counts = TaxonCounters::new();
        let mut total_seqs: usize = 0;
        let mut total_unclassified: usize = 0;
        for file_pair in files {
            file_index += 1;

            writeln!(file_writer, "{}\t{}", file_index, file_pair.join(","))?;
            file_writer.flush().unwrap();

            let score = args.minimum_quality_score;
            let paths = OptionPair::from_slice(file_pair);
            let mut reader = FastxReader::from_paths(paths, file_index, score)?;
            // let mut reader = create_reader(file_pair, file_index, score)?;
            let (thread_sequences, thread_unclassified) = process_fastx_file(
                &args,
                meros,
                hash_config,
                file_index,
                &mut reader,
                chtable,
                taxonomy,
                &mut total_taxon_counts,
            )?;
            total_seqs += thread_sequences;
            total_unclassified += thread_unclassified;
        }
        if let Some(output) = &args.kraken_output_dir {
            let filename = output.join("output.kreport2");
            report_kraken_style(
                filename,
                args.report_zero_counts,
                args.report_kmer_data,
                &taxonomy,
                &total_taxon_counts,
                total_seqs as u64,
                total_unclassified as u64,
            )?;
        }

        Ok(())
    };

    if args.paired_end_processing && !args.single_file_pairs {
        // 处理成对的文件
        let files = args.input_files.chunks(2).collect();
        process_funcs(files)?;
    } else {
        let files = args.input_files.chunks(1).collect();
        process_funcs(files)?;
    }

    Ok(())
}

pub fn run(args: Args) -> Result<()> {
    let options_filename = &args.k2d_dir.join("opts.k2d");
    let idx_opts = IndexOptions::read_index_options(options_filename)?;

    if args.paired_end_processing && !args.single_file_pairs && args.input_files.len() % 2 != 0 {
        // 验证文件列表是否为偶数个
        return Err(Error::new(
            ErrorKind::InvalidInput,
            "Paired-end processing requires an even number of input files.",
        ));
    }

    let taxonomy_filename = args.k2d_dir.join("taxo.k2d");
    let taxo = Taxonomy::from_file(taxonomy_filename)?;

    let hash_config = HashConfig::from_hash_header(&args.k2d_dir.join("hash_config.k2d"))?;

    println!("{:?}", hash_config);
    if hash_config.hash_capacity == 0 {
        panic!("`hash_capacity` can't be zero!");
    }
    println!("classify start...");
    let start = Instant::now();
    let meros = idx_opts.as_meros();
    let hash_files = find_and_sort_files(&args.k2d_dir, "hash", ".k2d")?;
    let chtable = CHTable::from_hash_files(hash_config, hash_files)?;

    process_files(args, meros, hash_config, &chtable, &taxo)?;
    let duration = start.elapsed();
    println!("classify took: {:?}", duration);
    Ok(())
}

#[allow(dead_code)]
fn main() {
    let args = Args::parse();
    if let Err(e) = run(args) {
        eprintln!("Application error: {}", e);
    }
}
