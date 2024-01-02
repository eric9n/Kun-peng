use anyhow::Result;
use std::fs::{create_dir_all, File};
use std::path::PathBuf;
use std::str::FromStr;

// 获取 url 地址的最后一层目录
pub fn get_last_segment_of_url(url: &str) -> &str {
    url.trim_end_matches('/').split('/').last().unwrap_or("")
}

// 创建文件，并检查文件所在目录是否存在，不存在则创建
pub fn create_file_in_dir(filename: &str) -> Result<File> {
    let path = PathBuf::from_str(filename).unwrap();

    if let Some(parent) = path.parent() {
        create_dir_all(parent)?;
    }

    let file: File = File::create(filename)?;

    Ok(file)
}

pub fn create_data_dir(dirname: &str) -> Result<PathBuf, anyhow::Error> {
    let path: PathBuf = PathBuf::from_str(dirname).unwrap();
    if !path.exists() {
        create_dir_all(&path)?;
    }
    Ok(path)
}

pub fn create_dir(dirname: &PathBuf) -> Result<(), anyhow::Error> {
    if !dirname.exists() {
        create_dir_all(&dirname)?;
    }
    Ok(())
}
