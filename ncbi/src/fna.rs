use async_compression::tokio::bufread::GzipDecoder;
use regex::Regex;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use tokio::fs::OpenOptions;
use tokio::{
    fs::File,
    io::{AsyncBufReadExt, AsyncReadExt, AsyncWriteExt, BufReader, BufWriter},
};

use anyhow::Result;
use tar::Archive;

pub async fn decompress_and_extract_tar_gz(
    gz_path: &PathBuf,
    out_path: &PathBuf,
    files_to_extract: Vec<String>,
) -> std::io::Result<()> {
    // Open the .tar.gz file
    let file = File::open(gz_path).await?;
    let buf_reader = BufReader::new(file);

    // Use GzipDecoder for decompression
    let gzip_decoder = GzipDecoder::new(buf_reader);

    // Create an async reader
    let mut async_reader = BufReader::new(gzip_decoder);

    // Read all decompressed data into memory
    let mut decompressed_data = Vec::new();
    async_reader.read_to_end(&mut decompressed_data).await?; // Works after importing AsyncReadExt

    // Use the tar crate to decompress the TAR archive
    let mut archive = Archive::new(&decompressed_data[..]);
    // archive.unpack(out_path)?;

    // 遍历 TAR 归档中的每个条目
    for entry in archive.entries()? {
        let mut entry = entry?;
        let path = entry.path()?.to_string_lossy().to_string();

        // 检查是否为需要提取的文件
        if files_to_extract.contains(&path) {
            let out_file_path = out_path.join(&path);

            // 创建输出文件夹
            if let Some(parent) = out_file_path.parent() {
                tokio::fs::create_dir_all(parent).await?;
            }

            // 解压缩并写入文件
            entry.unpack(out_file_path)?;
        }
    }

    Ok(())
}

async fn delete_hllp_json_files(out_dir: &Path) -> Result<()> {
    let entries = fs::read_dir(out_dir)?;

    for entry in entries {
        let entry = entry?;
        let path = entry.path();
        if path.is_file() {
            if let Some(file_name) = path.file_name() {
                if let Some(file_name_str) = file_name.to_str() {
                    if file_name_str.starts_with("hllp_") && file_name_str.ends_with(".json") {
                        tokio::fs::remove_file(&path).await?;
                        println!("Deleted: {:?}", path);
                    }
                }
            }
        }
    }

    Ok(())
}

pub async fn parse_assembly_fna(
    site: &str,
    data_dir: &PathBuf,
    asm_levels: &Vec<&str>,
) -> Result<HashMap<String, String>> {
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

            if ftp_path == "na" {
                continue;
            }

            if !asm_levels.contains(&asm_level) {
                continue;
            }

            let fna_file_name = format!(
                "{}/{}_genomic.fna.gz",
                site,
                ftp_path.split('/').last().unwrap_or_default()
            );
            gz_files.insert(fna_file_name, taxid.into());
        }
    }
    Ok(gz_files)
}

pub async fn write_to_fna(
    site: &str,
    group: &str,
    asm_levels: &Vec<&str>,
    data_dir: &PathBuf,
    out_dir: &PathBuf,
) -> Result<()> {
    log::info!("{} {} write to fna...", group, site);

    let gz_files = if site == "all" {
        let mut gz_files = parse_assembly_fna("genbank", data_dir, asm_levels).await?;
        let ref_gz_files = parse_assembly_fna("refseq", data_dir, asm_levels).await?;
        gz_files.extend(ref_gz_files);
        gz_files
    } else {
        parse_assembly_fna(site, data_dir, asm_levels).await?
    };
    let library_fna_path = out_dir.join("library.fna");
    let prelim_map_path = out_dir.join("prelim_map.txt");

    let mut fna_writer = BufWriter::new(
        OpenOptions::new()
            .create(true)
            .write(true)
            .open(&library_fna_path)
            .await?,
    );
    delete_hllp_json_files(&out_dir).await?;
    let mut map_writer = BufWriter::new(
        OpenOptions::new()
            .create(true)
            .write(true)
            .open(&prelim_map_path)
            .await?,
    );

    let re: Regex = Regex::new(r"^>(\S+)").unwrap();

    for (gz_path, taxid) in gz_files {
        let gz_file = data_dir.join(gz_path);
        if !gz_file.exists() {
            continue;
        }

        let file = File::open(gz_file).await?;
        let decompressor = GzipDecoder::new(BufReader::new(file));
        let mut reader = BufReader::new(decompressor);

        let mut line = String::new();
        while reader.read_line(&mut line).await? != 0 {
            if let Some(caps) = re.captures(&line) {
                let seqid = &caps[1];
                let full_tax_id = format!("taxid|{}", taxid);
                map_writer
                    .write_all(format!("TAXID\t{}|{}\t{}\n", full_tax_id, seqid, taxid).as_bytes())
                    .await?;
                fna_writer
                    .write_all(format!(">{}|{}", full_tax_id, &line[1..]).as_bytes())
                    .await?;
            } else {
                fna_writer.write_all(line.as_bytes()).await?;
            }
            line.clear();
        }
    }

    fna_writer.flush().await?;
    map_writer.flush().await?;

    log::info!("{} {} write to fna finished", group, site);
    Ok(())
}
