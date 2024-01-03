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
    // 元数据,文件路径: (etag, bool)
    pub meta: HashMap<String, (String, bool)>,
}

pub async fn save_metadata() -> Result<()> {
    let db = DATABASE.lock().await;
    let meta_file: PathBuf = db.dir.clone().join(".metadata");
    let mut file = File::create(meta_file)?;

    for (url, value) in db.meta.iter() {
        writeln!(file, "{},{},{}", url, value.0, value.1)?;
    }

    Ok(())
}

fn parse_metadata(filename: &PathBuf) -> HashMap<String, (String, bool)> {
    let mut meta = HashMap::new();
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
            if parts.len() < 3 {
                continue; // 如果数据不完整，跳过该行
            }
            let url = parts[0].to_string();
            let etag = parts[1].to_string();
            let downloaded = parts[2].parse::<bool>().unwrap_or(false);
            meta.insert(url, (etag, downloaded));
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

    pub fn insert_or_update(
        &mut self,
        key: String,
        latest_etag: Option<String>,
        downloaded: Option<bool>,
    ) {
        match self.meta.get_mut(&key) {
            Some((ref mut etag, ref mut down)) => {
                if let Some(tag) = latest_etag {
                    *etag = tag;
                }
                if let Some(flag) = downloaded {
                    *down = flag;
                }
            }
            None => {
                self.meta.insert(
                    key,
                    (latest_etag.unwrap_or_default(), downloaded.unwrap_or(false)),
                );
            }
        }
    }
}

lazy_static! {
    pub static ref DATABASE: Arc<Mutex<DataBase>> = Arc::new(Mutex::new(DataBase::new()));
}

pub async fn update_db_state(key: String, latest_etag: Option<String>, downloaded: Option<bool>) {
    let mut db = DATABASE.lock().await;
    db.insert_or_update(key, latest_etag, downloaded);
}

pub async fn check_db_state(key: String, etag: String) -> bool {
    let db = DATABASE.lock().await;
    if let Some((h_etag, _)) = db.meta.get(&key) {
        if &etag != "" && &etag == h_etag {
            return true;
        }
    }
    false
}
