use crate::down::retry_download;
use crate::md5sum::check_md5sum_file;
use crate::meta::get_local_etag;
use crate::meta::insert_local_etag;
use anyhow::Result;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use tokio::fs::File;
use tokio::io::{AsyncBufReadExt, BufReader};
use tokio::sync::mpsc;

pub const NCBI_GEN_URL: &'static str = "https://ftp.ncbi.nlm.nih.gov/genomes/";
pub const NCBI_TAXO_URL: &'static str = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy";

#[derive(Debug)]
pub struct DownTuple {
    url: String,
    etag: String,
    pub file: PathBuf,
}

impl DownTuple {
    pub fn new(url: String, file: PathBuf, etag: String) -> Self {
        Self { url, file, etag }
    }

    pub async fn new_taxo(url_path: String, taxo_dir: &PathBuf) -> Self {
        let taxdump_url = format!("{}/{}", NCBI_TAXO_URL, url_path);
        let file_name = if url_path.contains("/") {
            url_path.split("/").last().unwrap().to_string()
        } else {
            url_path
        };
        let output_path = taxo_dir.join(file_name);
        let etag = get_local_etag(&taxdump_url).await;
        DownTuple::new(taxdump_url, output_path, etag.unwrap_or_default())
    }

    pub async fn run(&self) -> Result<()> {
        let result = retry_download(&self.url, &self.file, &self.etag, 3).await?;
        if result != self.etag {
            insert_local_etag(self.url.clone(), result).await;
        }
        Ok(())
    }
}

#[derive(Debug)]
pub enum NcbiFile {
    Summary(DownTuple),
    Genomic(DownTuple, DownTuple),
    Taxonomy(DownTuple, DownTuple),
}

impl NcbiFile {
    pub async fn from_group(group: &str, data_dir: &PathBuf, site: &str) -> Self {
        let url = format!("{}{}/{}/assembly_summary.txt", NCBI_GEN_URL, site, group);
        let output_path = data_dir.join(format!("assembly_summary_{}.txt", site));

        let etag = get_local_etag(&url).await;
        NcbiFile::Summary(DownTuple::new(
            url.clone(),
            output_path,
            etag.unwrap_or("".into()),
        ))
    }

    pub async fn new_taxo(taxo_dir: &PathBuf) -> Vec<Self> {
        let mut down_tasks = vec![];
        let files = [
            "taxdump.tar.gz",
            "accession2taxid/nucl_gb.accession2taxid.gz",
            "accession2taxid/nucl_wgs.accession2taxid.gz",
        ];
        for url_path in files.iter() {
            let taxo = DownTuple::new_taxo(url_path.to_string(), taxo_dir).await;
            let md5_file = format!("{}.md5", url_path);
            let taxo_md5 = DownTuple::new_taxo(md5_file, taxo_dir).await;

            down_tasks.push(NcbiFile::Taxonomy(taxo, taxo_md5));
        }
        down_tasks
    }

    pub async fn run(&self) -> Result<()> {
        match self {
            NcbiFile::Summary(dt) => dt.run().await?,
            NcbiFile::Genomic(dt1, dt2) => {
                dt1.run().await?;
                dt2.run().await?;
            }
            NcbiFile::Taxonomy(dt1, dt2) => {
                dt1.run().await?;
                dt2.run().await?;
            }
        }
        Ok(())
    }

    pub async fn check(&self) -> Result<()> {
        match self {
            NcbiFile::Summary(_) => unreachable!(),
            NcbiFile::Genomic(dt1, dt2) => {
                let result = check_md5sum_file(&dt1.file, &dt2.file).await?;
                if !result {
                    return Err(anyhow::anyhow!("mismatch"));
                }
                Ok(())
            }
            NcbiFile::Taxonomy(dt1, dt2) => {
                let result = check_md5sum_file(&dt1.file, &dt2.file).await?;
                if !result {
                    return Err(anyhow::anyhow!("mismatch"));
                }
                Ok(())
            }
        }
    }
}

impl NcbiFile {
    pub async fn parse_assembly_file(
        &self,
        site: &str,
        data_dir: &PathBuf,
        sender: mpsc::Sender<NcbiFile>,
        counter: Arc<AtomicUsize>,
    ) -> Result<()> {
        match self {
            NcbiFile::Genomic(_, _) => {}
            NcbiFile::Taxonomy(_, _) => {}
            NcbiFile::Summary(ncbi) => {
                let file = File::open(&ncbi.file).await?;
                let reader = BufReader::new(file);
                let mut lines = reader.lines();

                while let Some(line) = lines.next_line().await? {
                    if line.starts_with('#') {
                        continue;
                    }

                    let fields: Vec<&str> = line.split('\t').collect();
                    if fields.len() > 19 {
                        let (_, asm_level, ftp_path) = (fields[5], fields[11], fields[19]);

                        if !["Complete Genome", "Chromosome"].contains(&asm_level)
                            || ftp_path == "na"
                        {
                            continue;
                        }

                        let fna_file_name = format!(
                            "{}_genomic.fna.gz",
                            ftp_path.split('/').last().unwrap_or_default()
                        );

                        let fna_url = format!("{}/{}", ftp_path, fna_file_name);
                        let fna_file = data_dir.join(site).join(&fna_file_name);
                        let fna_etag = get_local_etag(&fna_url).await.unwrap_or_default();
                        let genomic_dt = DownTuple::new(fna_url, fna_file, fna_etag);
                        let md5_url = format!("{}/md5checksums.txt", ftp_path);
                        let md5_file = data_dir
                            .join(site)
                            .join(format!("{}_md5checksums.txt", &fna_file_name));
                        let md5_etag = get_local_etag(&md5_url).await.unwrap_or_default();
                        let md5_dt = DownTuple::new(md5_url, md5_file, md5_etag);
                        let new_file = NcbiFile::Genomic(genomic_dt, md5_dt);
                        let _ = sender.send(new_file).await;
                        counter.fetch_add(1, Ordering::SeqCst);
                    }
                }
            }
        }
        Ok(())
    }
}
