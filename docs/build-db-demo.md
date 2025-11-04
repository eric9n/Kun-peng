# Kun-peng Database Build: Complete Demo

This guide shows multiple ways to build a Kun-peng database from NCBI downloads or from your own FASTA files, plus a quick classification sanity check. It is written for new users and reviewers who prefer a step‑by‑step walkthrough.

See also: Detailed Classification Demo (docs/classify-demo.md)

## Assumptions

- You have the `kun_peng` binary available (either installed or built at `./target/release/kun_peng`).
- Example directories in this demo:
  - Download directory: `data/`
  - Database directory: `test_database`
  - Example reads: `data/COVID_19.fa`

If your binary is not on PATH, replace `kun_peng` with `./target/release/kun_peng` in the commands below.

## A. One‑Command Build (Recommended)

Runs the full pipeline: merge-fna → estimate → chunk → build.

```bash
kun_peng build --download-dir data/ --db test_database --hash-capacity 1G
```

Notes:
- This command expects that `data/` contains NCBI downloads (taxonomy and genomes), typically prepared with your own download workflow. It will:
  - Merge genomes into `test_database/library/*.fna`
  - Create/update `seqid2taxid.map`, `taxonomy/`, `taxo.k2d`
  - Estimate hash capacity, generate chunk files, and build `hash_*.k2d`
- Useful options:
  - `--hash-capacity 1G` sets the number of hash slots (defaults to `1G`, producing ~4 GiB hash shards); raise or lower to fit your memory/disk budget.
  - `--max-file-size 2G` controls library shard size
  - `-k 35 -l 31 --minimizer-spaces 7` control KLMT parameters
  - `--load-factor 0.7`, `--max-n 4` control capacity estimation details

Expected log highlights: “merge fna start…”, “estimate start…”, “chunk db took: …”, “build k2 db took: …”.

## B. Build Only From an Existing Library

Use this when you already have `test_database/library/*.fna` and `seqid2taxid.map`, or you want to add more sequences and then rebuild the index.

### Step B1: Prepare/Update the Library

Choose one of the following:

- Option 1 (merge downloaded genomes):
  ```bash
  kun_peng merge-fna --download-dir data/ --db test_database --max-file-size 2G
  ```

- Option 2 (add your own FASTA/FASTA.GZ):
  ```bash
  kun_peng add-library --db test_database -i /path/to/fasta_or_folder
  ```
  - Accepts files with extensions: `.fa`, `.fna`, `.fasta`, `.fsa`, and their `.gz` variants
  - Updates `library/*.fna` and appends new entries to `seqid2taxid.map`
  - After adding, you must run `build-db` to rebuild hash tables

### Step B2: Build the Database Artifacts (estimate + chunk + build)

- Automatic estimation:
  ```bash
  kun_peng build-db --db test_database --hash-capacity 1G
  ```

- Skip estimation with a manually chosen capacity (advanced users only):
  ```bash
  kun_peng build-db --db test_database --hash-capacity 1G -c <EXACT_SLOT_COUNT>
  ```
  Caution: an undersized value may fail to build or result in a high load factor (hurting speed); an oversized value wastes disk and memory.

`--hash-capacity` works the same way here as in the all-in-one `build` command: leave it at `1G` for ~4 GiB hash shards, or pick the size that best matches your storage and load-time constraints.

KLMT, threads, and load factor options are the same as in section A.

## What Files Should Exist After a Successful Build?

In `test_database/` you should see:
- `seqid2taxid.map`
- `taxonomy/` and `taxo.k2d`
- `opts.k2d`
- `hash_config.k2d`
- One or more `hash_*.k2d` files

Quick check:
```bash
ls -lh test_database
```

Optional: estimate memory for direct loading (sum of `hash_*.k2d`):
```bash
bash cal_memory.sh test_database
```

## Quick Classification Sanity Check

Prepare output and temp directories, then run the integrated classification flow (splitr → annotate → resolve):

```bash
mkdir -p temp_chunk test_out
kun_peng classify --db test_database --chunk-dir temp_chunk --output-dir test_out data/COVID_19.fa
```

Expected output: timing logs for splitr/annotate/resolve and two files in `test_out/` (`output_*.txt`, `*.kreport2`).

## Direct Mode (Advanced)

If you have sufficient RAM, you can load all hash tables at once for maximal speed:

```bash
# First, see how much memory is needed (≈ total size of hash_*.k2d)
bash cal_memory.sh test_database

# Then run direct classification
kun_peng direct --db test_database data/COVID_19.fa
```

## Practical Tips

- If you run `add-library`, always rebuild with:
  ```bash
  kun_peng build-db --db test_database --hash-capacity 1G
  ```
- `--cache` reuses capacity estimation caches (enabled by default).
- Tune performance and resource usage with:
  - `-p <threads>` number of threads
  - `--max-file-size` shard size for library files (I/O parallelism)
  - `--hash-capacity` to control hash shard sizing (~4× the numeric capacity in bytes)
- About `-c <EXACT_SLOT_COUNT>`: use only if you understand capacity sizing and load factor tradeoffs.

## Minimal Path From Zero to Classification

```bash
# 1) One‑command build
kun_peng build --download-dir data/ --db test_database --hash-capacity 1G

# 2) Classification
mkdir -p temp_chunk test_out
kun_peng classify --db test_database --chunk-dir temp_chunk --output-dir test_out data/COVID_19.fa
```

Adjust `--hash-capacity` in both commands if you need smaller or larger hash shards.

Next: Detailed Classification Demo (docs/classify-demo.md)
