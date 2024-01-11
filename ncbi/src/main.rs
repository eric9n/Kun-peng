use anyhow::Result;
use clap::{Parser, ValueEnum};
use ncbi::fna::write_to_fna;
use ncbi::meta::{init_meta, save_meta};
use ncbi::task;
use ncbi::utils;
use std::fmt;
use std::path::PathBuf;
use tokio::runtime::Builder;

const NCBI_LIBRARY: &'static [&str] = &[
    "archaea", "bacteria", "viral", "fungi", "plant", "human", "protozoa",
];

fn validate_group(group: &str) -> Result<String, String> {
    let groups = utils::parse_comma_separated_list(&group);
    for grp in &groups {
        if !NCBI_LIBRARY.contains(&grp.as_str()) {
            return Err(format!("group not in ncbi library"));
        }
    }
    Ok(group.to_string())
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum Site {
    /// 下载 genbank 资源
    Genbank,
    /// 下载 refseq 资源
    Refseq,
}

impl fmt::Display for Site {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Site::Genbank => "genbank",
                Site::Refseq => "refseq",
            }
        )
    }
}

#[derive(Parser, Debug)]
#[clap(
    version,
    about = "ncbi download resource",
    long_about = "从 ncbi 网站上下载 genomes 资源"
)]
struct Args {
    #[arg(value_enum)]
    site: Site,

    /// 构建数据库的目录
    #[arg(short, long, default_value = "lib")]
    database: PathBuf,

    /// 下载 taxonomy 文件
    #[arg(short, long, default_value = "false")]
    taxonomy: bool,

    /// 从 NCBI 站点上下载某个种类的数据信息，必须是列表中所列名称，archaea,bacteria,fungi...
    #[arg(short, long, value_parser = validate_group)]
    group: Option<String>,
    /// 仅检查文件的 md5
    #[arg(short, long, default_value = "false")]
    md5: bool,

    /// 下载时的并行大小
    #[arg(short, long, default_value = "8")]
    num_threads: usize,
}

async fn async_run(args: Args) -> Result<()> {
    let db_path = utils::create_data_dir(&args.database).unwrap();
    init_meta(&db_path).await;

    let site = args.site.to_string();

    if let Some(group_arg) = args.group {
        let groups = utils::parse_comma_separated_list(&group_arg);
        for group in groups {
            let data_dir: PathBuf = db_path.join("library").join(group.clone());
            utils::create_dir(&data_dir.join(&site))?;

            let trans_group = if &group == "human" {
                "vertebrate_mammalian/Homo_sapiens".to_string()
            } else {
                group.to_string()
            };

            let result = if args.md5 {
                task::run_check(&site, &trans_group, &data_dir, args.num_threads).await
            } else {
                task::run_task(&site, &trans_group, &data_dir, args.num_threads).await
            };
            if result.is_ok() {
                let _ = write_to_fna(&site, &data_dir).await;
            }
        }
    }
    if args.taxonomy {
        let data_dir: PathBuf = db_path.join("taxonomy");
        utils::create_dir(&data_dir)?;
        let _ = task::run_taxo(&data_dir).await;
    }
    save_meta(&db_path).await?;
    Ok(())
}

fn main() -> Result<()> {
    env_logger::Builder::new()
        .filter_level(log::LevelFilter::Info)
        .init();

    let args = Args::parse();
    let num_thread = args.num_threads.clone();
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
