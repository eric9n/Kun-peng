use anyhow::{anyhow, Result};
use md5::{Digest, Md5};
use std::path::PathBuf;
use tokio::fs::File;
use tokio::io::{AsyncBufReadExt, AsyncReadExt, BufReader};

async fn get_md5_for_file(target_file: &PathBuf, md5_file: &PathBuf) -> Result<String> {
    let file = File::open(md5_file).await.map_err(|e| {
        tokio::io::Error::new(
            tokio::io::ErrorKind::NotFound,
            format!("File operation failed: {:?}-{:?}", md5_file, e),
        )
    })?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let file_name = target_file
        .file_name()
        .unwrap()
        .to_string_lossy()
        .to_string();
    while let Some(line) = lines.next_line().await? {
        if line.starts_with('#') || line.trim().is_empty() {
            continue; // 跳过注释和空行
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let hash_key = parts.last().unwrap();
            if hash_key.ends_with(&file_name) {
                let md5 = parts[0].trim().to_string();
                return Ok(md5);
            }
        }
    }

    Err(anyhow!("MD5 checksum not found for file {}", file_name))
}

pub async fn check_md5sum_file(target_file: &PathBuf, md5_file: &PathBuf) -> Result<bool> {
    let md5_value = get_md5_for_file(target_file, md5_file).await?;

    // 异步读取目标文件并计算其 MD5 散列值
    let mut hasher = Md5::new();
    let mut fna = File::open(target_file).await?;
    let mut buffer = vec![0; 8192]; // 使用缓冲区以减少内存占用

    loop {
        let n = fna.read(&mut buffer).await?;
        if n == 0 {
            break;
        }
        hasher.update(&buffer[..n]);
    }

    let fna_md5 = hasher.finalize();
    Ok(format!("{:x}", fna_md5) == md5_value)
}

// fn write_failed_to_file(file_path: PathBuf, items: Vec<String>) -> Result<()> {#[derive(Clone)]

//     let mut file = File::create(file_path)?;

//     for item in items {
//         writeln!(file, "{}", item)?;
//     }

//     Ok(())
// }
