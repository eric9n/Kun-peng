use crate::load::{DownTuple, NcbiFile};
use crate::meta::get_local_etag;
use anyhow::Result;
use std::path::PathBuf;
use tokio::fs::File;
use tokio::io::{AsyncBufReadExt, BufReader};

use tokio::sync::mpsc;

pub async fn parse_assembly_file(
    ncbi_file: NcbiFile,
    data_dir: &PathBuf,
    sender: mpsc::Sender<NcbiFile>,
) -> Result<()> {
    match ncbi_file {
        NcbiFile::Genomic(_, _) => {}
        NcbiFile::Summary(ncbi) => {
            let file = File::open(ncbi.file).await?;
            let reader = BufReader::new(file);
            let mut lines = reader.lines();

            while let Some(line) = lines.next_line().await? {
                if line.starts_with('#') {
                    continue;
                }

                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() > 19 {
                    let (_, asm_level, ftp_path) = (fields[5], fields[11], fields[19]);

                    if !["Complete Genome", "Chromosome"].contains(&asm_level) || ftp_path == "na" {
                        continue;
                    }

                    let fna_file_name = format!(
                        "{}_genomic.fna.gz",
                        ftp_path.split('/').last().unwrap_or_default()
                    );

                    let fna_url = format!("{}/{}", ftp_path, fna_file_name);
                    let fna_file = data_dir.join(&fna_file_name);
                    let fna_etag = get_local_etag(&fna_url).await.unwrap_or_default();
                    let genomic_dt = DownTuple::new(fna_url, fna_file, fna_etag);
                    let md5_url = format!("{}/md5checksums.txt", ftp_path);
                    let md5_file = data_dir.join(format!("{}_md5checksums.txt", &fna_file_name));
                    let md5_etag = get_local_etag(&md5_url).await.unwrap_or_default();
                    let md5_dt = DownTuple::new(md5_url, md5_file, md5_etag);
                    let new_file = NcbiFile::Genomic(genomic_dt, md5_dt);
                    let _ = sender.send(new_file).await;
                }
            }
        }
    }
    Ok(())
}
