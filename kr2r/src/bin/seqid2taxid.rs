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
    #[arg(short = 'f', long)]
    map_file: Option<PathBuf>,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let prelim_file = summary_prelim_map_files(&args.source)?;

    let map_file = if let Some(seqid_file) = args.map_file {
        seqid_file
    } else {
        args.source.join("seqid2taxid.map")
    };
    create_seqid2taxid_file(prelim_file, map_file.clone())?;

    println!("finished {:?}", &map_file);
    Ok(())
}
