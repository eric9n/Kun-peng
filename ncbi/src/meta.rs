use anyhow::Result;
use lazy_static::lazy_static;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;
use tokio::fs::File;
use tokio::io::{AsyncBufReadExt, AsyncWriteExt, BufReader};
use tokio::sync::Mutex;

const META_FILE_NAME: &'static str = ".metadata";

async fn parse_metadata(filename: &PathBuf) -> Result<HashMap<String, String>> {
    let mut meta: HashMap<String, String> = HashMap::new();
    if !filename.exists() {
        File::create(filename).await?;
    }

    let file = File::open(filename).await?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    while let Some(line) = lines.next_line().await? {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 2 {
            continue;
        }
        let url = parts[0].to_string();
        let etag = parts[1].to_string();
        meta.insert(url, etag);
    }

    Ok(meta)
}

pub struct Meta {
    inner: HashMap<String, String>,
}

impl Meta {
    pub fn new() -> Self {
        Self {
            inner: HashMap::new(),
        }
    }

    pub async fn init(&mut self, db_path: &PathBuf) {
        let meta_file: PathBuf = db_path.clone().join(META_FILE_NAME);
        if let Ok(meta) = parse_metadata(&meta_file).await {
            self.inner = meta;
        }
    }

    pub fn get_etag(&self, key: &str) -> Option<String> {
        self.inner.get(key).map(|s| s.to_string())
    }

    pub fn insert_or_update(&mut self, key: String, new_etag: String) {
        self.inner.insert(key, new_etag);
    }
}

pub async fn init_meta(db_path: &PathBuf) {
    let mut meta = META.lock().await;
    meta.init(db_path).await;
}

pub async fn get_local_etag(key: &str) -> Option<String> {
    let meta = META.lock().await;
    meta.get_etag(key)
}

pub async fn insert_local_etag(key: String, new_etag: String) {
    let mut meta = META.lock().await;
    meta.insert_or_update(key, new_etag);
}

pub async fn save_meta(db_path: &PathBuf) -> Result<()> {
    let meta = META.lock().await;
    let meta_file: PathBuf = db_path.join(META_FILE_NAME);
    let mut file = File::create(meta_file).await?;

    for (url, etag) in meta.inner.iter() {
        file.write_all(format!("{},{}\n", url, etag).as_bytes())
            .await?;
    }

    Ok(())
}

lazy_static! {
    pub static ref META: Arc<Mutex<Meta>> = Arc::new(Mutex::new(Meta::new()));
}
