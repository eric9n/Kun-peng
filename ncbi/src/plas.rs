use std::path::PathBuf;

use crate::client::retry_client;
use crate::load::{DownTuple, NcbiFile, NCBI_GEN_URL};
use crate::meta::get_local_etag;
use anyhow::Result;
use regex::Regex;

pub async fn download_plas_files(data_dir: PathBuf, plas_type: &str) -> Result<()> {
    let plas_url = format!("{}refseq/{}/", NCBI_GEN_URL, plas_type);
    let client = retry_client();
    let response: reqwest::Response = client.get(&plas_url).send().await?;
    if response.status().is_success() {
        let contents = response.text().await?;
        // println!("contents {:?}", contents);
        // 正则表达式匹配所有 href 属性
        let re = Regex::new(r#"href="([^"]+\.genomic\.fna\.gz)""#)?;

        // 查找并打印所有匹配的 href
        for cap in re.captures_iter(&contents) {
            let filename = &cap[1];
            let url = format!("{}{}", plas_url, filename);
            println!("url {:?}", url);
            let etag = get_local_etag(&url).await;
            let output_path = data_dir.join(filename);
            let ncbi_file = NcbiFile::Summary(DownTuple::new(
                url.clone(),
                output_path,
                etag.unwrap_or("".into()),
            ));
            ncbi_file.run().await?;
        }
    } else {
        log::error!("Failed to fetch the webpage.");
    }
    Ok(())
}
