use anyhow::Result;
use clap::{Parser, Subcommand, ValueEnum};
use lazy_static::lazy_static;
use ncbi::fna::write_to_fna;
use ncbi::meta::{init_meta, save_meta};
use ncbi::task;
use ncbi::utils;
use std::collections::HashMap;
use std::fmt;
use std::path::PathBuf;
use tokio::runtime::Builder;

const NCBI_LIBRARY: &'static [&str] = &[
    "archaea",
    "bacteria",
    "viral",
    "fungi",
    "plant",
    "human",
    "protozoa",
    "vertebrate_mammalian",
    "vertebrate_other",
    "invertebrate",
];

lazy_static! {
    static ref NCBI_ASM_LEVELS: HashMap<String, Vec<&'static str>> = {
        let mut m = HashMap::new();
        m.insert("complete_genome".to_string(), vec!["Complete Genome"]);
        m.insert("chromosome".to_string(), vec!["Chromosome"]);
        m.insert("scaffold".to_string(), vec!["Scaffold"]);
        m.insert("contig".into(), vec!["Contig"]);
        m.insert("basic".into(), vec!["Complete Genome", "Chromosome"]);
        m.insert("uncomplete".into(), vec!["Scaffold", "Contig"]);
        m.insert(
            "all".into(),
            vec!["Complete Genome", "Chromosome", "Scaffold", "Contig"],
        );
        m
    };
}

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
    /// genbank and refseq
    All,
}

impl fmt::Display for Site {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Site::Genbank => "genbank",
                Site::Refseq => "refseq",
                Site::All => "all",
            }
        )
    }
}

#[derive(Subcommand, Debug)]
enum Mode {
    /// 仅检查文件的 md5
    Md5,
    /// 解析 genomic 文件，并且生成 library fna 文件
    Fna {
        /// library fna 文件存储目录，为了不和原始文件混淆
        #[clap(value_parser)]
        out_dir: PathBuf,
    },
    /// 仅下载和解析 assembly 文件
    Assembly,
    /// 单独下载 genomic 文件，指定 url 地址
    Url {
        #[clap(value_parser)]
        url: String,
    },
}

#[derive(Parser, Debug)]
#[clap(
    version,
    about = "ncbi download resource",
    long_about = "从 ncbi 网站上下载 genomes 资源"
)]
struct Args {
    /// 构建数据库的目录
    #[arg(short, long = "db", default_value = "lib")]
    database: PathBuf,

    /// 下载时的并行大小
    #[arg(short, long, default_value = "8")]
    num_threads: usize,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// 从 NCBI 下载 taxonomy 文件 (alias: tax)
    #[command(alias = "tax")]
    Taxonomy,

    /// 从 NCBI 下载 genomes 数据 (alias: gen)
    #[command(alias = "gen")]
    Genomes {
        /// 从 NCBI 哪个站点目录下载（RefSeq或GenBank）
        #[arg(long, value_enum, default_value_t = Site::Refseq)]
        site: Site,

        /// Assembly level: the highest level of assembly for any object in the genome
        /// all, complete_genome, chromosome, scaffold, contig. basic: [complete_genome, chromosome]
        #[arg(long, default_value = "basic")]
        asm_level: String,

        /// 从 NCBI 站点上下载某个种类的数据信息，可以是逗号分隔的多个, archaea,bacteria,viral,fungi,plant,human,protozoa,vertebrate_mammalian,vertebrate_other,invertebrate
        #[arg(short, long, value_parser = validate_group)]
        group: String,

        /// 子命令，使用 md5 校验和生成 fna 文件
        #[command(subcommand)]
        mode: Option<Mode>,
    },
}

async fn async_run(args: Args) -> Result<()> {
    let db_path = utils::create_data_dir(&args.database).unwrap();
    init_meta(&db_path).await;

    match args.command {
        Commands::Taxonomy => {
            let data_dir: PathBuf = db_path.join("taxonomy");
            utils::create_dir(&data_dir)?;
            let _ = task::run_taxo(&data_dir).await;
        }
        Commands::Genomes {
            site,
            group,
            asm_level,
            mode,
        } => {
            // let site_str = site.to_string();
            let groups = utils::parse_comma_separated_list(&group);
            for grp in groups {
                let data_dir: PathBuf = db_path.join("library").join(grp.clone());
                match site {
                    Site::All => {
                        for s in [Site::Genbank, Site::Refseq].iter() {
                            utils::create_dir(&data_dir.join(&s.to_string()))?;
                        }
                    }
                    _ => {
                        utils::create_dir(&data_dir.join(&site.to_string()))?;
                    }
                }

                let trans_group = if &grp == "human" {
                    "vertebrate_mammalian/Homo_sapiens".to_string()
                } else {
                    grp.to_string()
                };

                let levels = NCBI_ASM_LEVELS.get(&asm_level).unwrap();

                match &mode {
                    Some(Mode::Md5) => match site {
                        Site::All => {
                            for site in [Site::Genbank, Site::Refseq].iter() {
                                let _ = task::run_check(
                                    &site.to_string(),
                                    &trans_group,
                                    &data_dir,
                                    &levels,
                                    args.num_threads,
                                )
                                .await;
                            }
                        }
                        _ => {
                            let _ = task::run_check(
                                &site.to_string(),
                                &trans_group,
                                &data_dir,
                                &levels,
                                args.num_threads,
                            )
                            .await;
                        }
                    },
                    Some(Mode::Fna { out_dir }) => {
                        let fna_out_dir = out_dir.join("library").join(grp.clone());
                        utils::create_dir(&fna_out_dir)?;
                        let _ = write_to_fna(
                            &site.to_string(),
                            &trans_group,
                            &levels,
                            &data_dir,
                            &fna_out_dir,
                        )
                        .await;
                    }
                    Some(Mode::Assembly) => match site {
                        Site::All => {
                            for s in [Site::Genbank, Site::Refseq].iter() {
                                let _ = task::run_assembly(
                                    &s.to_string(),
                                    &trans_group,
                                    &levels,
                                    &data_dir,
                                )
                                .await;
                            }
                        }
                        _ => {
                            let _ = task::run_assembly(
                                &site.to_string(),
                                &trans_group,
                                &levels,
                                &data_dir,
                            )
                            .await;
                        }
                    },
                    Some(Mode::Url { url }) => {
                        if site == Site::All {
                            log::error!("必须指定合适的site");
                        } else {
                            let result =
                                task::run_download_file(&site.to_string(), &data_dir, &url).await;
                            if result.is_err() {
                                log::error!("下载文件失败... {:?}", result);
                            }
                        }
                    }
                    None => match site {
                        Site::All => {
                            for s in [Site::Genbank, Site::Refseq].iter() {
                                let _ = task::run_task(
                                    &s.to_string(),
                                    &trans_group,
                                    &data_dir,
                                    &&levels,
                                    args.num_threads,
                                )
                                .await;
                            }
                        }
                        _ => {
                            let _ = task::run_task(
                                &site.to_string(),
                                &trans_group,
                                &data_dir,
                                &&levels,
                                args.num_threads,
                            )
                            .await;
                        }
                    },
                }
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
