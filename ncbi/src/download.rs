use crate::client::retry_client;
use anyhow::{Context, Result};
use futures::future::join_all;
use futures_util::StreamExt;
use indicatif::{ProgressBar, ProgressStyle};
use std::path::PathBuf;

use crate::db::DATABASE;
use std::sync::{Arc, Mutex};
use tokio::fs::OpenOptions;
use tokio::io::AsyncWriteExt;
use tokio::sync::Semaphore;
use tokio::time::Duration;

#[derive(Clone)]
pub struct Task {
    url: String,
    output_file: PathBuf,
}

impl Task {
    pub fn new(url: String, output_file: PathBuf) -> Self {
        Self { url, output_file }
    }

    pub async fn run(&self) -> Result<()> {
        let client = retry_client();

        // 如果文件存在，并且检查文件的 etag，如果 etag 没有变化，就跳过
        if self.output_file.exists() && self.output_file.is_file() {
            let head_resp: reqwest::Response = client.head(&self.url).send().await?;
            let h_etag = head_resp
                .headers()
                .get("etag")
                .map(|tag| tag.to_str().unwrap_or(""))
                .unwrap_or("");

            {
                let db = DATABASE.lock().await;
                if let Some(etag) = db.meta.get(&self.url) {
                    if h_etag != "" && h_etag == etag {
                        log::debug!("url {:?} does not need to be updated.", &self.url);
                        return Ok(());
                    }
                }
            }
        }

        let mut file = OpenOptions::new()
            .create(true)
            .write(true)
            .open(&self.output_file)
            .await
            .context("Failed to open file")?;

        let response = client
            .get(&self.url)
            .send()
            .await
            .context("Failed to send request")?;

        let latest_etag = response
            .headers()
            .get("etag")
            .map(|tag| tag.to_str().unwrap_or(""))
            .unwrap_or("")
            .to_string();

        let mut stream = response.bytes_stream();
        let mut buffer = Vec::new();

        while let Some(chunk) = stream.next().await {
            match chunk {
                Ok(result) => {
                    buffer.extend_from_slice(&result);
                    if buffer.len() >= 1 * 1024 * 1024 {
                        // 例如，缓冲区大小为 1MB
                        file.write_all(&buffer)
                            .await
                            .context("Error while writing to file")?;
                        buffer.clear();
                    }
                }
                Err(e) => {
                    if e.is_body() {
                        continue;
                    }
                    log::error!("Error while reading chunk: {}", e);
                    return Err(e.into());
                }
            }
        }

        if !buffer.is_empty() {
            file.write_all(&buffer)
                .await
                .context("Error while writing to file")?;
        }

        {
            let mut db = DATABASE.lock().await;
            db.insert_or_update(self.url.to_string(), latest_etag);
        }

        Ok(())
    }
}

pub struct Tasks {
    tag: String,
    max_concurrent_tasks: usize,
    tasks: Vec<Task>,
}

impl Tasks {
    pub fn new(tag: String, tasks: Vec<Task>, max_concurrent_tasks: usize) -> Self {
        Self {
            tag,
            max_concurrent_tasks,
            tasks,
        }
    }

    pub async fn run(&self) -> Result<()> {
        let semaphore = Arc::new(Semaphore::new(self.max_concurrent_tasks));
        let progress = Arc::new(Mutex::new(ProgressBar::new(self.tasks.len() as u64)));

        let mut task_list = vec![];

        // 设置进度条样式
        progress.lock().unwrap().set_style(
            ProgressStyle::default_bar()
                .template(
                    "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
                )
                .unwrap()
                .progress_chars("##-"),
        );

        for task in &self.tasks {
            let permit = semaphore.clone().acquire_owned().await?;
            let task_clone = task.clone(); // 克隆任务数据
            let progress_clone = progress.clone();

            let task_future = tokio::spawn(async move {
                let mut retries = 3;
                let mut result;
                loop {
                    result = task_clone.run().await;
                    if result.is_ok() || retries == 0 {
                        break;
                    }
                    let _ = std::fs::remove_file(&task_clone.output_file);
                    retries -= 1;
                    tokio::time::sleep(Duration::from_secs(2)).await; // 重试前等待
                }
                drop(permit); // 释放信号量
                progress_clone.lock().unwrap().inc(1); // 增加进度
                result
            });
            task_list.push(task_future);
        }

        for task_future in join_all(task_list).await {
            match task_future {
                Ok(Ok(_)) => {
                    // 任务成功完成，无需执行任何操作
                    // log::info!("{} {} tasks finished", &self.tag, self.tasks.len());
                }
                Ok(Err(e)) => {
                    // 任务失败，记录或处理错误
                    log::error!("{} Task failed: {}", &self.tag, e);
                }
                Err(e) => {
                    // 任务 panic 或无法加入，记录或处理错误
                    log::error!("{} Task panicked or could not be joined: {}", &self.tag, e);
                }
            }
        }
        progress
            .lock()
            .unwrap()
            .finish_with_message(format!("{} All tasks completed", self.tag));
        Ok(())
    }
}
