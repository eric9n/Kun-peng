use crate::load::NcbiFile;
use anyhow::Result;
use futures::stream::StreamExt;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use tokio::sync::{mpsc, Semaphore};

async fn process_tasks(
    task_type: String,
    num_threads: usize,
    mut receiver: mpsc::Receiver<NcbiFile>,
    next_tx: Option<mpsc::Sender<NcbiFile>>,
) -> Result<usize> {
    let semaphore = Arc::new(Semaphore::new(num_threads));
    let mut futures = vec![];
    let counter: Arc<AtomicUsize> = Arc::new(AtomicUsize::new(0));

    while let Some(task) = receiver.recv().await {
        let permit = semaphore.clone().acquire_owned().await?;
        let next_tx_clone = next_tx.clone();
        let task_type_clone = task_type.clone();
        let counter_clone = counter.clone();

        let task_future = tokio::spawn(async move {
            let result = match task_type_clone.as_str() {
                "run" => task.run().await,
                "check" => task.check().await,
                _ => unreachable!(),
            };
            drop(permit);
            if result.is_ok() {
                counter_clone.fetch_add(1, Ordering::SeqCst);
                if let Some(tx) = next_tx_clone {
                    let _ = tx.send(task).await;
                }
            }
            result
        });
        futures.push(task_future);
    }

    let mut stream = futures::stream::iter(futures).buffer_unordered(num_threads);

    while let Some(result) = stream.next().await {
        match result {
            Ok(Ok(_)) => {} // 任务成功完成
            Ok(Err(e)) => log::error!("Task failed: {}", e),
            Err(e) => log::error!("Task panicked or could not be joined: {}", e),
        }
    }

    Ok(counter.load(Ordering::SeqCst))
}

/// 处理 assembly 文件
async fn process_assembly_tasks(
    site: &str,
    group: &str,
    data_dir: &PathBuf,
    asm_levels: &Vec<&str>,
    tx: mpsc::Sender<NcbiFile>,
) -> Result<usize> {
    let counter = Arc::new(AtomicUsize::new(0));
    let assembly = NcbiFile::from_group(group, data_dir, site).await;
    match assembly.run().await {
        Ok(_) => {
            let result = assembly
                .process_summary_and_apply(site, data_dir, asm_levels, |file: NcbiFile| {
                    let tx_clone = tx.clone();
                    let counter_clone = counter.clone();
                    async move {
                        let _ = tx_clone.send(file).await;
                        counter_clone.fetch_add(1, Ordering::SeqCst);
                    }
                })
                .await;
            if result.is_err() {
                log::error!("Error parsing assembly file: {:?}", result);
            }
        }
        Err(e) => {
            log::info!("{}", e);
        }
    }

    Ok(counter.load(Ordering::SeqCst))
}

pub async fn run_task(
    site: &str,
    group: &str,
    data_dir: &PathBuf,
    asm_levels: &Vec<&str>,
    num_threads: usize,
) -> Result<()> {
    log::info!("{} {} download assembly file start...", group, site);
    let (tx, rx) = mpsc::channel(4096); // 通道大小可以根据需要调整
    let (tx1, rx1) = mpsc::channel(4096); // 通道大小可以根据需要调整
    let assembly_tasks = process_assembly_tasks(site, group, data_dir, asm_levels, tx);
    let download_handle = process_tasks("run".to_string(), num_threads, rx, Some(tx1));
    let md5_handle = process_tasks("check".to_string(), num_threads, rx1, None);
    // // 等待处理任务完成
    let (ably_res, down_res, md5_res) = tokio::join!(assembly_tasks, download_handle, md5_handle);
    log::info!(
        "{} {} file total count: {}, downloaded: {}, md5match: {}",
        group,
        site,
        ably_res?,
        down_res?,
        md5_res?
    );
    log::info!("{} {} file finished...", group, site);
    Ok(())
}

pub async fn run_check(
    site: &str,
    group: &str,
    data_dir: &PathBuf,
    asm_levels: &Vec<&str>,
    num_threads: usize,
) -> Result<()> {
    log::info!("{} {} check md5 start...", group, site);
    let (tx, rx) = mpsc::channel(4096); // 通道大小可以根据需要调整
    let assembly_tasks = process_assembly_tasks(site, group, data_dir, asm_levels, tx);
    let md5_handle = process_tasks("check".to_string(), num_threads, rx, None);
    // // 等待处理任务完成
    let (ably_res, md5_res) = tokio::join!(assembly_tasks, md5_handle);
    log::info!(
        "{} {} file total count: {}, md5match: {}",
        group,
        site,
        ably_res?,
        md5_res?
    );
    Ok(())
}

pub async fn run_taxo(taxo_dir: &PathBuf) -> Result<()> {
    log::info!("download taxonomy...");
    let files = [
        "taxdump.tar.gz",
        "accession2taxid/nucl_gb.accession2taxid.gz",
        "accession2taxid/nucl_wgs.accession2taxid.gz",
    ];
    for url_path in files.iter() {
        let ncbi_file = NcbiFile::new_taxo(taxo_dir, &url_path).await;
        let result = ncbi_file.run().await;
        if result.is_ok() && url_path.to_string() == "taxdump.tar.gz" {
            let _ = ncbi_file.decompress(taxo_dir).await;
        }
    }
    log::info!("download taxonomy finished...");
    Ok(())
}

pub async fn run_download_file(site: &str, data_dir: &PathBuf, fna_url: &str) -> Result<()> {
    let ncbi_file = NcbiFile::from_file(site, data_dir, fna_url).await;
    log::info!("{} download file start...", fna_url);
    ncbi_file.clear().await;
    ncbi_file.run().await?;
    ncbi_file.check().await?;
    log::info!("{} download file end...", fna_url);
    Ok(())
}

pub async fn run_assembly(
    site: &str,
    group: &str,
    asm_levels: &Vec<&str>,
    data_dir: &PathBuf,
) -> Result<()> {
    let assembly = NcbiFile::from_group(group, data_dir, site).await;
    if !assembly.file_exists() {
        let _ = assembly.run().await;
    }
    let total_counter = Arc::new(AtomicUsize::new(0));
    let counter = Arc::new(AtomicUsize::new(0));
    let result = assembly
        .process_summary_and_apply(site, data_dir, asm_levels, |file: NcbiFile| {
            let counter_clone = counter.clone();
            let total_counter_clone = total_counter.clone();
            async move {
                total_counter_clone.fetch_add(1, Ordering::SeqCst);
                if file.file_exists() {
                    counter_clone.fetch_add(1, Ordering::SeqCst);
                }
            }
        })
        .await;
    if result.is_err() {
        log::error!("Error parsing assembly file: {:?}", result);
    }

    let total_count = total_counter.load(Ordering::SeqCst);
    let count = counter.load(Ordering::SeqCst);
    log::info!(
        "{} {} 总文件数: {}, 本地文件数: {}",
        group,
        site,
        total_count,
        count
    );
    Ok(())
}
