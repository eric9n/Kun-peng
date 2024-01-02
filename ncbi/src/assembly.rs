use crate::download::{Task, Tasks};
use crate::site::{NCBI_GEN_URL, NCBI_SITES};
use anyhow::Result;
use md5::{Digest, Md5};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::PathBuf;

#[allow(dead_code)]
pub struct Assembly {
    group: String,
    data_dir: PathBuf,
    meta: HashMap<String, String>,
    md5sum: HashMap<String, String>,
    taxid_map: HashMap<String, String>,
}

impl Assembly {
    pub fn new(group: &str, data_dir: &PathBuf) -> Self {
        Self {
            group: group.to_string(),
            data_dir: data_dir.clone(),
            meta: HashMap::new(),
            md5sum: HashMap::new(),
            taxid_map: HashMap::new(),
        }
    }

    pub async fn download_assembly_file(&self) -> Result<()> {
        let mut tasks: Vec<Task> = vec![];
        for site in NCBI_SITES {
            let url = format!(
                "{}{}/{}/assembly_summary.txt",
                NCBI_GEN_URL, site, self.group
            );
            let output_path = self
                .data_dir
                .clone()
                .join(format!("assembly_summary_{}.txt", site));
            tasks.push(Task::new(url, output_path));
        }
        let task_handle = Tasks::new(tasks, 4);
        task_handle.run().await?;
        Ok(())
    }

    pub fn parse_assembly_file(&mut self) {
        for site in NCBI_SITES {
            let filename = self
                .data_dir
                .clone()
                .join(format!("assembly_summary_{}.txt", site));

            let file = File::open(filename).unwrap();
            let reader = BufReader::new(file);

            for line in reader.lines() {
                let line = line.unwrap();
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

                    let fna_url = format!("{}/{}", ftp_path, fna_file_name);
                    let md5_url = format!("{}/md5checksums.txt", ftp_path);
                    self.meta.insert(fna_file_name.clone(), fna_url);
                    self.md5sum.insert(fna_file_name.clone(), md5_url);
                    self.taxid_map.insert(fna_file_name.clone(), taxid.into());
                }
            }
        }
    }

    /// 根据 assembly 文件下载基因fna文件
    pub async fn download_genomic(&self, parallel: usize) -> Result<()> {
        let mut tasks = vec![];
        for (fna_file, fna_url) in self.meta.iter() {
            let output_path = self.data_dir.clone().join(fna_file);
            tasks.push(Task::new(fna_url.into(), output_path));
        }
        let task_handle = Tasks::new(tasks, parallel);
        task_handle.run().await?;
        Ok(())
    }

    /// 根据 assembly 文件下载md5checksums文件
    pub async fn download_md5_file(&self, parallel: usize) -> Result<()> {
        let mut tasks = vec![];
        for (fna_file, md5_url) in self.md5sum.iter() {
            let md5_file = format!("{}.md5checksums.txt", fna_file);
            let output_path = self.data_dir.clone().join(md5_file);
            tasks.push(Task::new(md5_url.into(), output_path));
        }
        let task_handle = Tasks::new(tasks, parallel);
        task_handle.run().await?;
        Ok(())
    }
}

fn check_md5sum_file(fna_file: PathBuf, md5_value: String) -> Result<bool> {
    // 读取 fna 文件并计算其 MD5 散列值
    let mut hasher = Md5::new();
    let mut fna = File::open(fna_file)?;
    let mut buffer = [0; 8192]; // 使用缓冲区以减少内存占用
    loop {
        let n = fna.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        hasher.update(&buffer[..n]);
    }

    let fna_md5 = hasher.finalize();

    Ok(format!("{:x}", fna_md5) == md5_value)
}

fn write_failed_to_file(file_path: PathBuf, items: Vec<String>) -> Result<()> {
    let mut file = File::create(file_path)?;

    for item in items {
        writeln!(file, "{}", item)?;
    }

    Ok(())
}

fn get_md5_for_file(md5_file: &PathBuf, target_file_name: &str) -> Result<String> {
    let file = File::open(md5_file)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue; // 跳过注释和空行
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let file_name = parts.last().unwrap();
            if file_name.ends_with(target_file_name) {
                let md5 = parts[0].trim().to_string();
                return Ok(md5);
            }
        }
    }

    Ok("".into()) // 如果没有找到对应的文件名，则返回 None
}

pub fn check_md5sum(ably: Assembly, delete: bool) -> Result<()> {
    let mut failed_files: Vec<String> = vec![];
    for (fna_file_name, _) in ably.md5sum.iter() {
        let fna_file = ably.data_dir.clone().join(fna_file_name);
        let md5_file_name = format!("{}.md5checksums.txt", fna_file_name);
        let md5_file = ably.data_dir.clone().join(md5_file_name);
        if let Ok(md5_value) = get_md5_for_file(&md5_file, &fna_file_name) {
            if let Ok(res) = check_md5sum_file(fna_file.clone(), md5_value) {
                if !res {
                    failed_files.push(fna_file_name.to_string());
                    if delete {
                        std::fs::remove_file(&fna_file)?;
                        std::fs::remove_file(&md5_file)?;
                    }
                }
            }
        }
    }
    let file_path = ably.data_dir.clone().join("md5_failed.txt");

    if failed_files.is_empty() {
        log::info!("all genomics file success...");
        if file_path.exists() {
            std::fs::remove_file(&file_path)?;
        }
    } else {
        log::error!("{} files are failed...", failed_files.len());

        write_failed_to_file(file_path.clone(), failed_files)?;
        log::info!("saved as {:?}", &file_path.to_str());
    }
    Ok(())
}
