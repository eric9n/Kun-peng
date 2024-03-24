use clap::{Parser, Subcommand};
mod annotate;
mod hashshard;
mod resolve;
mod splitr;

use std::io::Result;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    #[clap(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    Hashshard(hashshard::Args),
    Splitr(splitr::Args),
    Annotate(annotate::Args),
    Resolve(resolve::Args),
}

fn main() -> Result<()> {
    let args = Args::parse();
    // let matches = Command::new("squid")
    //     .version("2.1")
    //     .author("Eric <eric9n@gmail.com>")
    //     .about("classify a set of sequences like Kraken 2")
    //     // .subcommand(splitr::Args::command().name("splitr"))
    //     // .subcommand(annotate::Args::command().name("annotate"))
    //     // .subcommand(resolve::Args::command().name("resolve"))
    //     // .subcommand(hashshard::Args::command().name("hashshard"))
    //     .get_matches();
    match args.command {
        Commands::Hashshard(cmd_args) => {
            hashshard::run(cmd_args)?;
        }
        Commands::Splitr(cmd_args) => {
            splitr::run(cmd_args)?;
        }
        Commands::Annotate(cmd_args) => {
            annotate::run(cmd_args)?;
        }
        Commands::Resolve(cmd_args) => {
            resolve::run(cmd_args)?;
        }
    }

    // let cmd1_args = splitr::Args::from_arg_matches(sub_matches).expect("parse splitr arg error");
    Ok(())
    // match args() {
    //     Some(("splitr", sub_matches)) => {
    //         let cmd1_args =
    //             splitr::Args::from_arg_matches(sub_matches).expect("parse splitr arg error");
    //         splitr::run(cmd1_args)
    //     }
    //     Some(("annotate", sub_matches)) => {
    //         let cmd1_args =
    //             annotate::Args::from_arg_matches(sub_matches).expect("parse annotate arg error");
    //         annotate::run(cmd1_args)
    //     }
    //     Some(("resolve", sub_matches)) => {
    //         let cmd1_args =
    //             resolve::Args::from_arg_matches(sub_matches).expect("parse resolve arg error");
    //         resolve::run(cmd1_args)
    //     }
    //     Some(("hashshard", sub_matches)) => {
    //         let cmd1_args =
    //             hashshard::Args::from_arg_matches(sub_matches).expect("parse resolve arg error");
    //         hashshard::run(cmd1_args)
    //     }
    //     _ => Ok(()),
    // }
}
