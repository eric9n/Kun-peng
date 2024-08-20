## ncbi tool

The ncbi binary is used to download resources from the NCBI website. Here is the help manual for the ncbi binary:

``` sh
./target/release/ncbi -h
ncbi download resource

Usage: ncbi [OPTIONS] <COMMAND>

Commands:
  taxonomy  Download taxonomy files from NCBI (alias: tax)
  genomes   Download genomes data from NCBI (alias: gen)
  help      Print this message or the help of the given subcommand(s)

Options:
  -d, --download-dir <DOWNLOAD_DIR>  Directory to store downloaded files [default: lib]
  -n, --num-threads <NUM_THREADS>    Number of threads to use for downloading [default: 20]
  -h, --help                         Print help (see more with '--help')
  -V, --version                      Print version
```
