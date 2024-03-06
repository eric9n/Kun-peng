use clap::Parser;
use kr2r::table::HashConfig;
use kr2r::taxonomy::Taxonomy;
use kr2r::IndexOptions;
use std::io::Result;

/// Command line arguments for the classify program.
///
/// This structure defines the command line arguments that are accepted by the classify program.
/// It uses the `clap` crate for parsing command line arguments.
#[derive(Parser, Debug, Clone)]
#[clap(version, about = "inspect")]
struct Args {
    /// The file path for the Kraken 2 index.
    #[clap(short = 'H', long = "index-filename", value_parser, required = true)]
    index_filename: String,

    /// The file path for the Kraken 2 taxonomy.
    #[clap(short = 't', long = "taxonomy-filename", value_parser, required = true)]
    taxonomy_filename: String,

    /// The file path for the Kraken 2 options.
    #[clap(short = 'o', long = "options-filename", value_parser, required = true)]
    options_filename: String,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let idx_opts = IndexOptions::read_index_options(args.options_filename.clone())?;
    let taxo = Taxonomy::from_file(&args.taxonomy_filename)?;
    let config = HashConfig::<u32>::from(args.index_filename.clone())?;
    println!("index option {:?}", idx_opts);
    println!("taxonomy node count {:?}", taxo.node_count());
    println!("compact hash table {:?}", config);

    Ok(())
}
