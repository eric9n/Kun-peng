use clap::Parser;

use flate2::read::GzDecoder;
use kr2r::utils::{find_files, open_file};
use std::collections::HashMap;
use std::fs::{create_dir_all, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Result, Write};
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser, Debug, Clone)]
#[clap(version, about = "A tool for processing genomic files")]
pub struct Args {
    /// Directory to store downloaded files
    #[arg(short, long, default_value = "lib")]
    pub download_dir: PathBuf,

    /// ncbi library fna database directory
    #[arg(long = "db", required = true)]
    pub database: PathBuf,
    // /// seqid2taxid.map file path, default = $database/seqid2taxid.map
    // #[arg(short = 'm', long)]
    // pub id_to_taxon_map_filename: Option<PathBuf>,
}

fn parse_assembly_fna(assembly_file: &PathBuf, site: &str) -> Result<HashMap<String, String>> {
    let mut gz_files: HashMap<String, String> = HashMap::new();
    let file = open_file(&assembly_file)?;
    let reader = BufReader::new(file);
    let lines = reader.lines();

    let parent_path = assembly_file
        .parent()
        .expect("Can't find assembly file parent directory");
    for line in lines {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() > 19 {
            let (taxid, _, ftp_path) = (fields[5], fields[11], fields[19]);

            if ftp_path == "na" {
                continue;
            }

            // let levels = vec!["Complete Genome", "Chromosome"];
            // if !levels.contains(&asm_level) {
            //     continue;
            // }

            let fna_file_name = format!(
                "{}/{}/{}_genomic.fna.gz",
                parent_path.to_string_lossy(),
                site,
                ftp_path.split('/').last().unwrap_or_default()
            );
            gz_files.insert(fna_file_name, taxid.into());
        }
    }
    Ok(gz_files)
}

fn process_gz_file(
    gz_file: &PathBuf,
    map_writer: &mut BufWriter<File>,
    fna_writer: &mut BufWriter<File>,
    fna_start: &regex::Regex,
    taxid: &str,
) -> Result<()> {
    let file = open_file(gz_file)?;
    let decompressor = GzDecoder::new(BufReader::new(file));
    let mut reader = BufReader::new(decompressor);

    let mut line = String::new();
    let mut map_buffer = String::new(); // Buffer for map writer
    let mut fna_buffer = String::new(); // Buffer for fna writer

    while reader.read_line(&mut line)? != 0 {
        if let Some(caps) = fna_start.captures(&line) {
            let seqid = &caps[1];
            map_buffer.push_str(&format!("kraken:taxid|{}|{}\t{}\n", taxid, seqid, taxid));
            fna_buffer.push_str(&format!(">kraken:taxid|{}|{}", taxid, &line[1..]));
        } else {
            fna_buffer.push_str(&line);
        }

        // Write to the writers if the buffer size exceeds a certain threshold
        if map_buffer.len() > 10000 {
            map_writer.write_all(map_buffer.as_bytes())?;
            map_buffer.clear();
        }

        if fna_buffer.len() > 10000 {
            fna_writer.write_all(fna_buffer.as_bytes())?;
            fna_buffer.clear();
        }

        line.clear();
    }

    // Write any remaining buffered content
    if !map_buffer.is_empty() {
        map_writer.write_all(map_buffer.as_bytes())?;
    }

    if !fna_buffer.is_empty() {
        fna_writer.write_all(fna_buffer.as_bytes())?;
    }

    fna_writer.flush()?;
    map_writer.flush()?;

    Ok(())
}

const PREFIX: &'static str = "assembly_summary";
const SUFFIX: &'static str = "txt";

fn merge_fna(assembly_files: &Vec<PathBuf>, database: &PathBuf) -> Result<()> {
    let pattern = format!(r"{}_(\S+)\.{}", PREFIX, SUFFIX);
    let file_site = regex::Regex::new(&pattern).unwrap();

    let library_fna_path = database.join("library.fna");
    let seqid2taxid_path = database.join("seqid2taxid.map");
    let mut fna_writer = BufWriter::new(
        OpenOptions::new()
            .create(true)
            .write(true)
            .open(&library_fna_path)?,
    );
    let mut map_writer = BufWriter::new(
        OpenOptions::new()
            .create(true)
            .write(true)
            .open(&seqid2taxid_path)?,
    );

    let fna_start: regex::Regex = regex::Regex::new(r"^>(\S+)").unwrap();
    for assembly_file in assembly_files {
        if let Some(caps) = file_site.captures(assembly_file.to_string_lossy().as_ref()) {
            if let Some(matched) = caps.get(1) {
                let gz_files = parse_assembly_fna(assembly_file, matched.as_str())?;

                for (gz_path, taxid) in gz_files {
                    let gz_file = PathBuf::from(&gz_path);
                    if !gz_file.exists() {
                        // eprintln!("{} does not exist", gz_file.to_string_lossy());
                        continue;
                    }

                    process_gz_file(
                        &gz_file,
                        &mut map_writer,
                        &mut fna_writer,
                        &fna_start,
                        &taxid,
                    )?;
                }

                fna_writer.flush()?;
                map_writer.flush()?;
            }
        }
    }

    Ok(())
}

pub fn run(args: Args) -> Result<()> {
    // 开始计时
    let start = Instant::now();
    println!("merge fna start...");
    let download_dir = args.download_dir;
    let database = &args.database;

    let dst_tax_dir = database.join("taxonomy");
    create_dir_all(&dst_tax_dir)?;

    let source_names_file = &download_dir.join("taxonomy").join("names.dmp");
    assert!(source_names_file.exists());
    let dst_name_file = &dst_tax_dir.join("names.dmp");
    if !dst_name_file.exists() {
        std::fs::copy(source_names_file, dst_name_file)?;
    }

    let source_nodes_file = &download_dir.join("taxonomy").join("nodes.dmp");
    assert!(source_nodes_file.exists());
    let dst_nodes_file = &dst_tax_dir.join("nodes.dmp");
    if !dst_nodes_file.exists() {
        std::fs::copy(source_nodes_file, dst_nodes_file)?;
    }

    let library_fna_path = database.join("library.fna");
    let seqid2taxid_path = database.join("seqid2taxid.map");
    if library_fna_path.exists() && seqid2taxid_path.exists() {
        println!("library.fna and seqid2taxid.map exists!");
        return Ok(());
    }

    if library_fna_path.exists() {
        std::fs::remove_file(library_fna_path)?;
    }
    if seqid2taxid_path.exists() {
        std::fs::remove_file(seqid2taxid_path)?;
    }
    let assembly_files = find_files(&download_dir, &PREFIX, &SUFFIX);

    merge_fna(&assembly_files, &args.database)?;

    // 计算持续时间
    let duration = start.elapsed();
    println!("merge fna took: {:?}", duration);
    Ok(())
}

#[allow(dead_code)]
fn main() {
    let args = Args::parse();
    if let Err(e) = run(args) {
        eprintln!("Application error: {}", e);
    }
}
