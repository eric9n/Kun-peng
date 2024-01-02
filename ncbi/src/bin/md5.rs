use anyhow::Result;
use clap::Parser;
use env_logger;
use ncbi::assembly::{check_md5sum, Assembly};
use std::path::PathBuf;
use std::str::FromStr;

#[derive(Parser, Debug)]
#[clap(
    version,
    about = "ncbi check genomics file md5sum",
    long_about = "检查 ncbi genomics 文件的 md5 值"
)]
struct Args {
    /// 数据库的路径
    #[arg(short, long, default_value = "lib")]
    database: Option<String>,

    /// 从 NCBI 站点上下载某个种类的数据信息，必须是列表中所列名称
    #[arg(short, long)]
    group: Option<String>,

    /// 删除校验错误的文件
    #[arg(long, default_value = "true")]
    delete: bool,
}

#[tokio::main]
async fn main() -> Result<()> {
    env_logger::Builder::new()
        .filter_level(log::LevelFilter::Info)
        .init();

    let args = Args::parse();
    let db_path =
        PathBuf::from_str(args.database.unwrap_or("lib".into()).as_str()).expect("找不到数据库");

    if let Some(group) = &args.group {
        let data_dir = db_path.join(group);
        let mut ably = Assembly::new(group, &data_dir);
        ably.parse_assembly_file();
        check_md5sum(ably, args.delete)?;
    }
    Ok(())
}
