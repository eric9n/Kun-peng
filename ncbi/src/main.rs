use anyhow::Result;
use clap::Parser;
use ncbi::meta::{init_meta, save_meta};
use ncbi::task;
use ncbi::utils;
use std::path::PathBuf;
use tokio::runtime::Builder;

#[derive(Parser, Debug)]
#[clap(
    version,
    about = "ncbi download resource",
    long_about = "从 ncbi 网站上下载 genomes 资源"
)]
struct Args {
    /// 列出 NCBI 站点上的种类列表信息，实时拉取
    #[arg(short, long)]
    list: bool,

    /// 构建数据库的目录
    #[arg(short, long, default_value = "lib")]
    database: PathBuf,

    /// 从 NCBI 站点上下载某个种类的数据信息，必须是列表中所列名称，archaea,bacteria,fungi...
    #[arg(short, long)]
    group: Option<String>,
    /// 仅检查文件的 md5
    #[arg(short, long, default_value = "false")]
    md5: bool,

    /// 下载时的并行大小
    #[arg(short, long, default_value = "8")]
    threads: usize,
}

async fn async_run(args: Args) -> Result<()> {
    let db_path = utils::create_data_dir(&args.database).unwrap();
    init_meta(&db_path).await;

    if let Some(group_arg) = args.group {
        let groups = utils::parse_comma_separated_list(&group_arg);
        for group in groups {
            let data_dir: PathBuf = db_path.join(group.clone());
            utils::create_dir(&data_dir)?;
            if args.md5 {
                let _ = task::run_check(&group, &data_dir, args.threads).await;
            } else {
                let _ = task::run_task(&group, &data_dir, args.threads).await;
            }
        }
    }
    save_meta(&db_path).await?;
    Ok(())
}

fn main() -> Result<()> {
    env_logger::Builder::new()
        .filter_level(log::LevelFilter::Info)
        .init();

    let args = Args::parse();
    let num_thread = args.threads.clone();
    // 创建一个 Runtime 实例，并配置线程数
    let runtime = Builder::new_multi_thread()
        .enable_all()
        .thread_name("ncbi")
        // .max_blocking_threads(100)
        .worker_threads(num_thread) // 设置所需的工作线程数
        .build()
        .expect("Failed to create runtime");

    // 使用 Runtime 运行异步代码
    runtime.block_on(async_run(args))?;

    Ok(())
}
