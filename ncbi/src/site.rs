use std::collections::HashMap;

use crate::client::retry_client;
use anyhow::Result;
use regex::Regex;
use std::fmt;

// #[derive(Debug)]
pub struct Item {
    name: String,
    date: String,
}

impl fmt::Debug for Item {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {})", self.name, self.date)
    }
}

pub const NCBI_SITES: &[&str] = &["refseq"];
pub const NCBI_GEN_URL: &'static str = "https://ftp.ncbi.nlm.nih.gov/genomes/";

// 从https://ftp.ncbi.nlm.nih.gov/genomes/下的 genbank或者 refreq 目录获取解析 assembly_summary的信息。
pub fn extract_info_from_html(html: &str) -> Result<Vec<Item>> {
    let re = Regex::new(r#"<a href="[^"]*">([^<]*)</a>\s*([\d-]{10} \d{2}:\d{2})\s*([-\dKM]*)"#)
        .unwrap();
    let mut items = Vec::new();

    for cap in re.captures_iter(html) {
        let item = Item {
            name: cap[1].trim_end_matches("/").to_string(),
            date: cap[2].to_string(),
        };
        if item.name.starts_with("README") || item.name.contains("assembly") {
            continue;
        } else {
            items.push(item);
        }
    }

    Ok(items)
}

pub fn site_display(map: &HashMap<String, Vec<Item>>) {
    for (key, item) in map {
        log::info!("{}: {:?}", key, item);
    }
}

pub fn check_group(group: &str, map: &HashMap<String, Vec<Item>>) -> bool {
    for (_, items) in map {
        for item in items {
            if item.name == group {
                return true;
            }
        }
    }
    false
}

pub async fn open_ncbi_site(display: bool) -> Result<HashMap<String, Vec<Item>>> {
    let mut map: HashMap<String, Vec<Item>> = HashMap::new();
    let client = retry_client();

    for site in NCBI_SITES {
        let url: String = format!("{}{}/", NCBI_GEN_URL, site);

        let response: reqwest::Response = client.get(&url).send().await?;

        if response.status().is_success() {
            let contents = response.text().await?;
            let item_vec = extract_info_from_html(&contents).unwrap();
            map.insert(site.to_string(), item_vec);
        } else {
            log::error!("Failed to fetch the webpage.");
        }
    }

    if display {
        site_display(&map);
    }
    Ok(map)
}
