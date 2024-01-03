use crate::client::retry_client;
use crate::db::{check_db_state, update_db_state};
use anyhow::Result;
use futures::future::join_all;
use futures_util::StreamExt;
use indicatif::{ProgressBar, ProgressStyle};
use reqwest::{self, header, StatusCode};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use tokio::fs::{self, OpenOptions};
use tokio::io::AsyncWriteExt;
use tokio::sync::Semaphore;
use tokio::time::{timeout, Duration};

#[derive(Clone)]
pub struct Task {
    url: String,
    output_file: PathBuf,
}

#[allow(dead_code)]
impl Task {
    pub fn new(url: String, output_file: PathBuf) -> Self {
        Self { url, output_file }
    }

    fn get_etag(&self, response: &reqwest::Response) -> String {
        response
            .headers()
            .get(header::ETAG)
            .map(|tag| tag.to_str().unwrap_or(""))
            .unwrap_or("")
            .to_string()
    }

    async fn is_latest_etag(
        &self,
        client: &reqwest_middleware::ClientWithMiddleware,
    ) -> Result<bool> {
        // 如果文件存在，并且检查文件的 etag，如果 etag 没有变化，就跳过
        if self.output_file.exists() && self.output_file.is_file() {
            let head_resp: reqwest::Response = client.head(&self.url).send().await?;
            let h_etag = self.get_etag(&head_resp);

            if check_db_state(self.url.to_string(), h_etag.clone()).await {
                return Ok(true);
            }
        }
        Ok(false)
    }

    pub async fn run_resume(&self) -> Result<()> {
        let client: &reqwest_middleware::ClientWithMiddleware = retry_client();

        // 如果文件存在，并且检查文件的 etag，如果 etag 没有变化，就跳过
        if self.is_latest_etag(client).await? {
            return Ok(());
        }
        update_db_state(self.url.to_string(), None, Some(false)).await;

        // 检查本地文件大小
        let file_size = if self.output_file.exists() {
            fs::metadata(&self.output_file).await?.len()
        } else {
            0
        };

        // 设置 Range 头部
        let mut headers = header::HeaderMap::new();
        headers.insert(header::RANGE, format!("bytes={}-", file_size).parse()?);

        // 发送带 Range 头部的 GET 请求
        let response = client.get(&self.url).headers(headers).send().await?;
        let status = response.status();
        let latest_etag = self.get_etag(&response);

        // 仅当服务器响应 206 Partial Content 时处理响应体
        if status == StatusCode::PARTIAL_CONTENT {
            // 打开文件用于追加
            let mut file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(&self.output_file)
                .await?;

            // 读取响应内容并写入文件
            let mut stream = response.bytes_stream();
            while let Some(chunk) = timeout(Duration::from_secs(30), stream.next()).await? {
                let chunk = chunk?;
                file.write_all(&chunk).await?;
            }

            update_db_state(self.url.to_string(), Some(latest_etag), Some(true)).await;
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
                let mut retries = 6;
                let mut result;
                loop {
                    result = task_clone.run_resume().await;
                    if result.is_ok() || retries == 0 {
                        break;
                    }
                    retries -= 1;
                    tokio::time::sleep(Duration::from_secs(2)).await; // 重试前等待
                }
                drop(permit); // 释放信号量
                if result.is_ok() {
                    progress_clone.lock().unwrap().inc(1); // 增加进度
                }
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
