use clap::Parser;

use kr2r::utils::{create_seqid2taxid_file, summary_prelim_map_files};
use std::io::Result;
use std::path::PathBuf;

#[derive(Parser, Debug, Clone)]
#[clap(version, about = "seqid to taxid map file")]
pub struct Args {
    /// the database directory
    #[arg(long, required = true)]
    pub source: PathBuf,

    /// seqid2taxid.map file path, default = $source/seqid2taxid.map
    #[arg(short = 'm', long)]
    pub id_to_taxon_map_filename: Option<PathBuf>,
}

pub fn run(args: Args) -> Result<()> {
    let prelim_file = summary_prelim_map_files(&args.source)?;
    let map_file = args
        .id_to_taxon_map_filename
        .unwrap_or(args.source.join("seqid2taxid.map"));
    create_seqid2taxid_file(prelim_file, map_file.clone())?;

    println!("finished {:?}", &map_file);
    Ok(())
}

#[allow(dead_code)]
fn main() {
    let args = Args::parse();
    if let Err(e) = run(args) {
        eprintln!("Application error: {}", e);
    }
}
