use async_compression::tokio::bufread::GzipDecoder;
use regex::Regex;
use std::collections::HashMap;
use tokio::fs::OpenOptions;
use tokio::{
    fs::File,
    io::{AsyncBufReadExt, AsyncWriteExt, BufReader, BufWriter},
};

use anyhow::Result;
use std::path::PathBuf;

pub async fn parse_assembly_fna(site: &str, data_dir: &PathBuf) -> Result<HashMap<String, String>> {
    let mut gz_files: HashMap<String, String> = HashMap::new();
    let file_name = format!("assembly_summary_{}.txt", site);
    let file_path = data_dir.join(file_name);
    let file = File::open(&file_path).await?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    while let Some(line) = lines.next_line().await? {
        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() > 19 {
            let (taxid, asm_level, ftp_path) = (fields[5], fields[11], fields[19]);

            if !["Complete Genome", "Chromosome"].contains(&asm_level) || ftp_path == "na" {
                continue;
            }

            let fna_file_name = format!(
                "{}_genomic.fna.gz",
                ftp_path.split('/').last().unwrap_or_default()
            );
            gz_files.insert(fna_file_name, taxid.into());
        }
    }
    Ok(gz_files)
}

pub async fn write_to_fna(site: &str, data_dir: &PathBuf) -> Result<()> {
    log::info!("write to fna...");

    let gz_files = parse_assembly_fna(site, data_dir).await?;
    let library_fna_path = data_dir.join(format!("library_{}.fna", &site));
    let prelim_map_path = data_dir.join(format!("prelim_map_{}.txt", &site));

    let mut fna_writer = BufWriter::new(
        OpenOptions::new()
            .create(true)
            .write(true)
            .open(&library_fna_path)
            .await?,
    );
    let mut map_writer = BufWriter::new(
        OpenOptions::new()
            .create(true)
            .write(true)
            .open(&prelim_map_path)
            .await?,
    );

    let re: Regex = Regex::new(r"^>(\S+)").unwrap();

    for (gz_path, taxid) in gz_files {
        let gz_file = data_dir.join(&site).join(gz_path);
        let file = File::open(gz_file).await?;
        let decompressor = GzipDecoder::new(BufReader::new(file));
        let mut reader = BufReader::new(decompressor);

        let mut line = String::new();
        while reader.read_line(&mut line).await? != 0 {
            if let Some(caps) = re.captures(&line) {
                let seqid = &caps[1];
                map_writer
                    .write_all(format!("{}\t{}\n", seqid, taxid).as_bytes())
                    .await?;
                fna_writer
                    .write_all(format!(">kraken:taxid|{}|{}", taxid, &line[1..]).as_bytes())
                    .await?;
            } else {
                fna_writer.write_all(line.as_bytes()).await?;
            }
            line.clear();
        }
    }

    fna_writer.flush().await?;
    map_writer.flush().await?;

    log::info!("write to fna finished");
    Ok(())
}
