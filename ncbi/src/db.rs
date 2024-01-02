use anyhow::Result;
use lazy_static::lazy_static;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::str::FromStr;
use std::sync::Arc;
use tokio::sync::Mutex;

#[allow(dead_code)]
pub struct DataBase {
    /// 包含的种类, archaea, bacteria...
    groups: Vec<String>,
    /// 文件路径
    pub dir: PathBuf,
    // 元数据,文件路径: etag
    pub meta: HashMap<String, String>,
}

pub async fn save_metadata() -> Result<()> {
    let db = DATABASE.lock().await;
    let meta_file: PathBuf = db.dir.clone().join(".metadata");
    let mut file = File::create(meta_file)?;

    for (url, etag) in db.meta.iter() {
        writeln!(file, "{},{}", url, etag)?;
    }

    Ok(())
}

fn parse_metadata(filename: &PathBuf) -> HashMap<String, String> {
    let mut meta: HashMap<String, String> = HashMap::new();
    if !filename.exists() {
        File::create(filename).expect("创建 metadata 文件失败");
    }

    let file = File::open(filename).expect("读取 metadata 文件失败");
    let reader = BufReader::new(file);
    for line in reader.lines() {
        if let Ok(item) = line {
            if item.starts_with('#') {
                continue; // 跳过以 '#' 开头的行
            }

            let parts: Vec<&str> = item.split(',').collect();
            let url = parts[0].to_string();
            let etag = parts[1].to_string();
            meta.insert(url, etag);
        }
    }
    meta
}

impl DataBase {
    pub fn new() -> Self {
        Self {
            dir: PathBuf::from_str("lib").unwrap(),
            groups: vec![],
            meta: HashMap::new(),
        }
    }

    pub async fn init(&mut self, dir: PathBuf) {
        let meta_file: PathBuf = dir.clone().join(".metadata");
        self.dir = dir;
        self.meta = parse_metadata(&meta_file);
    }

    pub fn insert_or_update(&mut self, key: String, latest_etag: String) {
        match self.meta.get_mut(&key) {
            Some(etag) => *etag = latest_etag,
            None => {
                self.meta.insert(key, latest_etag);
            }
        }
    }
}

lazy_static! {
    pub static ref DATABASE: Arc<Mutex<DataBase>> = Arc::new(Mutex::new(DataBase::new()));
}
