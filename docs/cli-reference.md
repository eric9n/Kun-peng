# Kun-peng CLI Reference

This page summarizes the `kun_peng` command-line interface by workflow.

Use it as a practical reference for:

- what each subcommand does
- when to use it
- the most important flags
- minimal command examples

For the complete machine-generated help of any command, run:

```bash
kun_peng --help
kun_peng <subcommand> --help
```

## Command Groups

Build and database preparation:

- `build`: run the full build flow from downloaded data
- `merge-fna`: merge downloaded genome files into a library
- `add-library`: add local FASTA files into an existing library
- `estimate`: estimate required hash table capacity
- `build-db`: build final database artifacts from an existing library
- `hashshard`: convert a Kraken 2 database into Kun-peng's sharded format

Classification:

- `classify`: integrated chunk-based workflow
- `direct`: load all hash tables at once for maximum speed
- `splitr`: split inputs into chunk files
- `annotate`: annotate chunk files against the database
- `resolve`: resolve taxonomy assignments and write reports

## Global Entry Point

Basic form:

```bash
kun_peng <subcommand> [options]
```

Examples:

```bash
kun_peng --help
kun_peng classify --help
kun_peng build-db --help
```

## Build and Database Commands

### `build`

Purpose:

- One-command database construction from downloaded input data
- Internally covers library preparation and final database build steps

Use when:

- you have a download directory with taxonomy and genome data
- you want the shortest path from raw downloads to a usable database

Basic usage:

```bash
kun_peng build --download-dir data/ --db test_database --hash-capacity 1G
```

Important options:

- `--download-dir <DIR>`: input download directory
- `--db <DIR>`: target database directory
- `--hash-capacity <SIZE>`: shard capacity, for example `1G`
- `--max-file-size <SIZE>`: maximum temporary library shard size
- `-p, --threads <N>`: build threads
- `--load-factor <FLOAT>`: hash table occupancy target
- `-k`, `-l`, `--minimizer-spaces`: minimizer parameters

Notes:

- `1G` hash capacity corresponds to about a 4 GiB shard file.
- If you already have `library/*.fna`, use `build-db` instead.

### `merge-fna`

Purpose:

- Merge downloaded genomes into the database library layout

Use when:

- you want to prepare `library/*.fna` before running `build-db`
- you want library preparation as a separate step

Basic usage:

```bash
kun_peng merge-fna --download-dir data/ --db test_database
```

Important options:

- `--download-dir <DIR>`: source downloads
- `--db <DIR>`: database directory
- `--max-file-size <SIZE>`: maximum library shard size

### `add-library`

Purpose:

- Add your own FASTA files into an existing database library

Use when:

- you want to extend a database with local sequence files
- you are not starting from an NCBI download layout

Basic usage:

```bash
kun_peng add-library --db test_database -i /path/to/fastas
```

Important options:

- `--db <DIR>`: existing database directory
- `-i, --input-library <PATH>...`: one or more FASTA files or directories
- `--max-file-size <SIZE>`: maximum library shard size

Notes:

- Accepted inputs include `.fa`, `.fna`, `.fasta`, `.fsa`, and `.gz` variants.
- After `add-library`, you must run `build-db` to rebuild `hash_*.k2d`.

### `estimate`

Purpose:

- Estimate the hash table capacity required for a library

Use when:

- you want to inspect sizing before a full build
- you plan to pass a manual slot count to `build-db -c`

Basic usage:

```bash
kun_peng estimate --database test_database
```

Important options:

- `--database <PATH>`: database directory or library path
- `--cache`: reuse cached estimation data when available
- `--load-factor <FLOAT>`: occupancy target
- `-n, --n <N>`: maximum qualifying hash code
- `-p, --threads <N>`: worker threads

### `build-db`

Purpose:

- Build the final database artifacts from an existing library
- Runs estimate, chunk, and build steps

Use when:

- `library/*.fna` already exists
- you used `merge-fna` or `add-library` first

Basic usage:

```bash
kun_peng build-db --db test_database --hash-capacity 1G
```

Important options:

- `--db <DIR>`: database directory
- `--hash-capacity <SIZE>`: shard capacity
- `-c, --required-capacity <SLOTS>`: skip estimation and force exact slot count
- `--cache`: reuse cached estimation data
- `--load-factor <FLOAT>`: occupancy target
- `-p, --threads <N>`: worker threads

Notes:

- `-c` is an advanced option. Too small can fail or slow classification later; too large wastes disk and memory.
- This is the rebuild step you need after `add-library`.

### `hashshard`

Purpose:

- Convert a Kraken 2 database into Kun-peng's sharded hash layout

Use when:

- you already have `hash.k2d`, `opts.k2d`, and `taxo.k2d`
- you want to classify with Kun-peng without rebuilding from source FASTA

Basic usage:

```bash
kun_peng hashshard --db /path/to/kraken_db --hash-capacity 1G
```

Important options:

- `--db <DIR>`: Kraken 2 database directory
- `--hash-capacity <SIZE>`: target shard capacity

Notes:

- If `hash_config.k2d` already exists in the target directory, the command stops to avoid overwriting.
- After conversion, you can use both `classify` and `direct`.

## Classification Commands

### `classify`

Purpose:

- Run the full chunk-based classification workflow
- Internally runs `splitr`, `annotate`, and `resolve`

Use when:

- you want the standard low-memory workflow
- the database is too large to load fully into RAM

Basic usage:

```bash
mkdir -p temp_chunk test_out
kun_peng classify \
  --db test_database \
  --chunk-dir temp_chunk \
  --output-dir test_out \
  data/COVID_19.fa
```

Important options:

- `--db <DIR>`: database directory
- `--chunk-dir <DIR>`: temp working directory
- `--output-dir <DIR>`: Kraken-style output directory
- `-p, --num-threads <N>`: threads
- `--buffer-size <BYTES>`: read/annotation buffering
- `--batch-size <N>`: taxid aggregation batch size
- `-T, --confidence-threshold <FLOAT>`: confidence threshold
- `-g, --minimum-hit-groups <N>`: minimum hit groups
- `-P, --paired-end-processing`: paired-end mode
- `-Q, --minimum-quality-score <N>`: FASTQ quality threshold

Input support:

- FASTA
- FASTQ
- gzipped FASTA/FASTQ
- multiple files
- a single `.txt` file containing one input path per line

Notes:

- `--chunk-dir` must be clean. Leftover `sample_*.k2`, `sample_id*.map`, or `sample_*.bin` files will cause an error.

### `direct`

Purpose:

- Load all hash tables into memory and classify directly

Use when:

- you want maximum throughput
- your machine has enough RAM for the entire database

Basic usage:

```bash
bash cal_memory.sh test_database
kun_peng direct --db test_database data/COVID_19.fa
```

Important options:

- `--db <DIR>`: database directory
- `--output-dir <DIR>`: output directory
- `-p, --num-threads <N>`: threads
- `-T, --confidence-threshold <FLOAT>`: confidence threshold
- `-g, --minimum-hit-groups <N>`: minimum hit groups
- `-P, --paired-end-processing`: paired-end mode
- `-Q, --minimum-quality-score <N>`: FASTQ quality threshold

Notes:

- Required RAM is roughly the sum of all `hash_*.k2d` files.
- If this is too large, switch to `classify`.

### `splitr`

Purpose:

- Split input reads into chunk files for later processing

Use when:

- you want to run the classification pipeline step by step
- you need to inspect or benchmark the chunking stage separately

Basic usage:

```bash
kun_peng splitr --db test_database --chunk-dir temp_chunk data/COVID_19.fa
```

Important options:

- `--db <DIR>`: database directory
- `--chunk-dir <DIR>`: temp working directory
- `-p, --num-threads <N>`: threads
- `-P, --paired-end-processing`: paired-end mode
- `-Q, --minimum-quality-score <N>`: FASTQ quality threshold

### `annotate`

Purpose:

- Annotate previously created chunk files against the database

Use when:

- `splitr` output already exists
- you want to tune the annotation stage independently

Basic usage:

```bash
kun_peng annotate --db test_database --chunk-dir temp_chunk
```

Important options:

- `--db <DIR>`: database directory
- `--chunk-dir <DIR>`: temp working directory
- `--buffer-size <BYTES>`: internal buffer size
- `--batch-size <N>`: taxid aggregation batch size
- `-p, --num-threads <N>`: threads

### `resolve`

Purpose:

- Resolve taxonomy assignments from annotated chunk data
- Write per-read output and report files

Use when:

- `annotate` has already completed
- you want to rerun taxonomy resolution with different thresholds

Basic usage:

```bash
kun_peng resolve --db test_database --chunk-dir temp_chunk --output-dir test_out
```

Important options:

- `--db <DIR>`: database directory
- `--chunk-dir <DIR>`: temp working directory
- `--output-dir <DIR>`: output directory
- `-p, --num-threads <N>`: threads
- `-T, --confidence-threshold <FLOAT>`: confidence threshold
- `-g, --minimum-hit-groups <N>`: minimum hit groups
- `-K, --report-kmer-data`: include minimizer details in report output
- `-z, --report-zero-counts`: include zero-count taxa in reports

## Common Patterns

Full build from downloads:

```bash
kun_peng build --download-dir data/ --db test_database --hash-capacity 1G
```

Build from an existing library:

```bash
kun_peng merge-fna --download-dir data/ --db test_database
kun_peng build-db --db test_database --hash-capacity 1G
```

Add local FASTA files and rebuild:

```bash
kun_peng add-library --db test_database -i /path/to/fastas
kun_peng build-db --db test_database --hash-capacity 1G
```

Integrated low-memory classification:

```bash
mkdir -p temp_chunk test_out
kun_peng classify --db test_database --chunk-dir temp_chunk --output-dir test_out data/COVID_19.fa
```

Stepwise classification:

```bash
kun_peng splitr --db test_database --chunk-dir temp_chunk data/COVID_19.fa
kun_peng annotate --db test_database --chunk-dir temp_chunk
kun_peng resolve --db test_database --chunk-dir temp_chunk --output-dir test_out
```

Convert Kraken 2 database and classify:

```bash
kun_peng hashshard --db /path/to/kraken_db --hash-capacity 1G
kun_peng classify --db /path/to/kraken_db --chunk-dir temp_chunk --output-dir test_out data/COVID_19.fa
```

## Related Docs

- [build-db-demo.md](build-db-demo.md): step-by-step database build guide
- [classify-demo.md](classify-demo.md): detailed classification guide
- [hashshard-demo.md](hashshard-demo.md): Kraken 2 conversion walkthrough
- [../README_en.md](../README_en.md): streamlined project overview
