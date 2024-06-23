This workspace contains two projects: `kr2r` and `ncbi`. The `kr2r` project includes an example that demonstrates how to use the `kun_peng` binary, a tool for processing gene classification, to build a database and process a sample file.

## Get Started

Follow these steps to build the projects and run the example.

### Build the Projects

First, ensure that both projects are built. You can do this by running the following command from the root of the workspace:

```sh
cargo build --release
```

This will build the kr2r and ncbi project in release mode.


### Run the `kun_peng` Example

Next, run the example script that demonstrates how to use the kun_peng binary. Execute the following command from the root of the workspace:

```sh
cargo run --release --example build_and_classify --package kr2r
```

This will run the build_and_classify.rs example located in the kr2r project's examples directory.


Example Output
You should see output similar to the following:

```txt
Executing command: /path/to/workspace/target/release/kun_peng build --download-dir data/ --db test_database
kun_peng build output: [build output here]
kun_peng build error: [any build errors here]

Executing command: /path/to/workspace/target/release/kun_peng direct --db test_database data/COVID_19.fa
kun_peng direct output: [direct output here]
kun_peng direct error: [any direct errors here]
```

This output confirms that the kun_peng commands were executed successfully and the files were processed as expected.


Run the `ncbi` Example
Run the example script in the ncbi project to download the necessary files. Execute the following command from the root of the workspace:

```sh
cargo run --release --example run_download --package ncbi
```

This will run the run_download.rs example located in the ncbi project's examples directory. The script will:

1. Ensure the necessary directories exist.
2. Download the required files using the ncbi binary with the following commands:
  * ./target/release/ncbi -d downloads gen -g archaea
  * ./target/release/ncbi -d downloads tax


Example Output
You should see output similar to the following:

```txt
Executing command: /path/to/workspace/target/release/ncbi -d /path/to/workspace/downloads gen -g archaea
NCBI binary output: [download output here]

Executing command: /path/to/workspace/target/release/ncbi -d /path/to/workspace/downloads tax
NCBI binary output: [download output here]
```


## ncbi tool

The ncbi binary is used to download resources from the NCBI website. Here is the help manual for the ncbi binary:
```sh
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


## kun_peng tool

```sh
Usage: kun_peng <COMMAND>

Commands:
  estimate   estimate capacity
  build      build `k2d` files
  hashshard  split hash file
  splitr     Split fast(q/a) file into ranges
  annotate   annotate a set of sequences
  resolve    resolve taxonomy tree
  classify   Integrates 'splitr', 'annotate', and 'resolve' into a unified workflow for sequence classification. classify a set of sequences
  direct     Directly load all hash tables for classification annotation
  merge-fna  A tool for processing genomic files
  help       Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```


### build database

Build the kun_peng database like Kraken2, specifying the directory for the data files downloaded from NCBI, as well as the database directory.

```sh
./target/release/kun_peng build -h
build database

Usage: kun_peng build [OPTIONS] --download-dir <DOWNLOAD_DIR> --db <DATABASE>

Options:
  -d, --download-dir <DOWNLOAD_DIR>
          Directory to store downloaded files
      --db <DATABASE>
          ncbi library fna database directory
  -k, --k-mer <K_MER>
          Set length of k-mers, k must be positive integer, k=35, k cannot be less than l [default: 35]
  -l, --l-mer <L_MER>
          Set length of minimizers, 1 <= l <= 31 [default: 31]
      --minimizer-spaces <MINIMIZER_SPACES>
          Number of characters in minimizer that are ignored in comparisons [default: 7]
  -T, --toggle-mask <TOGGLE_MASK>
          Minimizer ordering toggle mask [default: 16392584516609989165]
      --min-clear-hash-value <MIN_CLEAR_HASH_VALUE>

  -r, --requested-bits-for-taxid <REQUESTED_BITS_FOR_TAXID>
          Bit storage requested for taxid 0 <= r < 31 [default: 0]
  -p, --threads <THREADS>
          Number of threads [default: 10]
      --cache
          estimate capacity from cache if exists
      --max-n <MAX_N>
          Set maximum qualifying hash code [default: 4]
      --load-factor <LOAD_FACTOR>
          Proportion of the hash table to be populated (build task only; def: 0.7, must be between 0 and 1) [default: 0.7]
  -h, --help
          Print help
  -V, --version
          Print version
```


### classify

The classification process is divided into three modes:

1. Direct Processing Mode:

* Description: In this mode, all database files are loaded simultaneously, which requires a significant amount of memory. Before running this mode, you need to check the total size of hash_*.k2d files in the database directory using the provided script (bash cal_memory.sh out_dir). Ensure that your available memory meets or exceeds this size.
* Characteristics:
    * High memory requirements
    * Fast performance

```sh
./target/release/kun_peng direct -h
Directly load all hash tables for classification annotation

Usage: kun_peng direct [OPTIONS] --db <DATABASE> [INPUT_FILES]...

Arguments:
  [INPUT_FILES]...  A list of input file paths (FASTA/FASTQ) to be processed by the classify program

Options:
      --db <DATABASE>
          database hash chunk directory and other files
  -P, --paired-end-processing
          Enable paired-end processing
  -S, --single-file-pairs
          Process pairs with mates in the same file
  -Q, --minimum-quality-score <MINIMUM_QUALITY_SCORE>
          Minimum quality score for FASTQ data [default: 0]
  -T, --confidence-threshold <CONFIDENCE_THRESHOLD>
          Confidence score threshold [default: 0]
  -K, --report-kmer-data
          In comb. w/ -R, provide minimizer information in report
  -z, --report-zero-counts
          In comb. w/ -R, report taxa w/ 0 count
  -g, --minimum-hit-groups <MINIMUM_HIT_GROUPS>
          The minimum number of hit groups needed for a call [default: 2]
  -p, --num-threads <NUM_THREADS>
          The number of threads to use [default: 10]
      --output-dir <KRAKEN_OUTPUT_DIR>
          File path for outputting normal Kraken output
  -h, --help
          Print help (see more with '--help')
  -V, --version
          Print version
```

2. Chunk Processing Mode:

* Description: This mode processes the sample data in chunks, loading only a small portion of the database files at a time. This reduces the memory requirements, needing a minimum of 4GB of memory plus the size of one pair of sample files.
* Characteristics:
    * Low memory consumption
    * Slower performance compared to Direct Processing Mode


```sh
./target/release/kun_peng classify -h
Integrates 'splitr', 'annotate', and 'resolve' into a unified workflow for sequence classification. classify a set of sequences

Usage: kun_peng classify [OPTIONS] --db <DATABASE> --chunk-dir <CHUNK_DIR> [INPUT_FILES]...

Arguments:
  [INPUT_FILES]...  A list of input file paths (FASTA/FASTQ) to be processed by the classify program

Options:
      --db <DATABASE>

      --chunk-dir <CHUNK_DIR>
          chunk directory
      --output-dir <KRAKEN_OUTPUT_DIR>
          File path for outputting normal Kraken output
  -P, --paired-end-processing
          Enable paired-end processing
  -S, --single-file-pairs
          Process pairs with mates in the same file
  -Q, --minimum-quality-score <MINIMUM_QUALITY_SCORE>
          Minimum quality score for FASTQ data [default: 0]
  -p, --num-threads <NUM_THREADS>
          The number of threads to use [default: 10]
      --batch-size <BATCH_SIZE>
          [default: 16777216]
  -T, --confidence-threshold <CONFIDENCE_THRESHOLD>
          Confidence score threshold [default: 0]
  -g, --minimum-hit-groups <MINIMUM_HIT_GROUPS>
          The minimum number of hit groups needed for a call [default: 2]
      --kraken-db-type
          Enables use of a Kraken 2 compatible shared database
  -K, --report-kmer-data
          In comb. w/ -R, provide minimizer information in report
  -z, --report-zero-counts
          In comb. w/ -R, report taxa w/ 0 count
      --full-output
          output file contains all unclassified sequence
  -h, --help
          Print help (see more with '--help')
  -V, --version
          Print version
```

3. Step-by-Step Processing Mode:

* Description: This mode breaks down the chunk processing mode into individual steps, providing greater flexibility in managing the entire classification process.
* Characteristics:
    * Flexible processing steps
    * Similar memory consumption to Chunk Processing Mode
    * Performance varies based on execution steps
