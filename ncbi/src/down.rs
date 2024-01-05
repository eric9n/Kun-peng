use crate::client::retry_client;
use anyhow::Result;
use futures_util::StreamExt;
use reqwest::{self, header, StatusCode};
use std::path::PathBuf;
use std::str::FromStr;
use tokio::fs::{self, OpenOptions};
use tokio::io::AsyncWriteExt;
use tokio::time::{sleep, timeout, Duration};

fn get_etag(response: &reqwest::Response) -> String {
    response
        .headers()
        .get(header::ETAG)
        .map(|tag| tag.to_str().unwrap_or(""))
        .unwrap_or("")
        .to_string()
}

fn get_content_length(response: &reqwest::Response) -> u64 {
    // 如果解析不到这个字段，设置成 1，为了让下载继续
    response
        .headers()
        .get(header::CONTENT_LENGTH)
        .and_then(|value| value.to_str().ok())
        .and_then(|value_str| u64::from_str(value_str).ok())
        .unwrap_or(1)
}

/// 断点续传
async fn download_resume(url: &str, file_name: &PathBuf, etag: &str) -> Result<String> {
    let client = retry_client();

    // 检查本地文件大小
    let file_size = if file_name.exists() {
        fs::metadata(file_name).await?.len()
    } else {
        0
    };

    let head_resp: reqwest::Response = client.head(url).send().await?;
    let latest_etag = get_etag(&head_resp);
    let content_length = get_content_length(&head_resp);

    // 如果文件下载完成，并且 etag 没有更新，则不下载
    if file_size == content_length && latest_etag != "" && latest_etag == etag {
        return Ok(etag.into());
    }

    // 设置 Range 头部
    let mut headers = header::HeaderMap::new();
    headers.insert(header::RANGE, format!("bytes={}-", file_size).parse()?);

    // 发送带 Range 头部的 GET 请求
    let response = client.get(url).headers(headers).send().await?;
    let status = response.status();

    // 仅当服务器响应 206 Partial Content 时处理响应体
    if status == StatusCode::PARTIAL_CONTENT {
        // 打开文件用于追加
        let mut file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&file_name)
            .await?;

        // 读取响应内容并写入文件
        let mut stream = response.bytes_stream();
        while let Some(chunk) = timeout(Duration::from_secs(30), stream.next()).await? {
            let chunk = chunk?;
            file.write_all(&chunk).await?;
        }
    }

    Ok(latest_etag)
}

/// 下载文件
pub async fn retry_download(
    url: &str,
    file_name: &PathBuf,
    etag: &str,
    retry: i32,
) -> Result<String> {
    let mut count = retry;
    let mut err = anyhow::anyhow!("download failed");
    while count > 0 {
        let result = download_resume(url, file_name, etag).await;
        match result {
            Ok(retrieved_etag) => return Ok(retrieved_etag),
            Err(e) => err = e,
        }
        count -= 1;
        sleep(Duration::from_secs(3)).await;
    }
    Err(err)
}
