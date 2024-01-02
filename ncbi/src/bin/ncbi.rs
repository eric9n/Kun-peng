use anyhow::anyhow;
use clap::Parser;
use ncbi::assembly::{check_md5sum, Assembly};
use ncbi::{site, utils};

use env_logger;
use ncbi::db::{save_metadata, DATABASE};

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
    database: Option<String>,

    /// 从 NCBI 站点上下载某个种类的数据信息，必须是列表中所列名称
    #[arg(short, long)]
    group: Option<String>,
    /// 检查文件的 md5
    #[arg(short, long, default_value = "true")]
    check_md5: Option<bool>,

    /// 下载时的并行大小
    #[arg(short, long, default_value = "8")]
    parallel: usize,
}

#[tokio::main]
async fn main() -> Result<(), anyhow::Error> {
    env_logger::Builder::new()
        .filter_level(log::LevelFilter::Info)
        .init();

    let args = Args::parse();
    let db_path = utils::create_data_dir(args.database.unwrap_or("lib".into()).as_str()).unwrap();
    {
        let mut db = DATABASE.lock().await;
        db.init(db_path.clone()).await;
    }
    if let Some(group) = &args.group {
        let items = site::open_ncbi_site(args.list).await.unwrap();
        if !site::check_group(&group, &items) {
            site::site_display(&items);
            return Err(anyhow!("group 参数错误, 请选择上述所列种类进行操作"));
        }
        // 创建对应的数据库子目录
        let data_dir = db_path.join(group);
        utils::create_dir(&data_dir)?;
        let mut ably = Assembly::new(group, &data_dir);
        ably.download_assembly_file().await?;
        ably.parse_assembly_file();
        ably.download_genomic(args.parallel).await?;
        ably.download_md5_file(args.parallel).await?;
        check_md5sum(ably, true)?;
    }

    save_metadata().await?;
    Ok(())
}
