# Kun-peng Classification: Complete Demo

This guide explains how to run Kun-peng classification in two modes: the integrated chunk workflow and the direct loading mode. It also covers inputs, outputs, common options, and practical tips.

See also: Detailed Database Build Demo (docs/build-db-demo.md)

## Prerequisites

- A built database directory, e.g. `test_database`, containing `hash_*.k2d`, `hash_config.k2d`, `taxo.k2d`, `opts.k2d`, and `seqid2taxid.map`.
- One or more input read files in FASTA/FASTQ format (optionally gzipped), e.g. `data/COVID_19.fa`.
- The `kun_peng` binary on PATH (or use `./target/release/kun_peng`).

## A. Integrated Classification Workflow

Runs an end-to-end pipeline: splitr → annotate → resolve. Suitable for typical usage and modest memory machines.

```bash
mkdir -p temp_chunk test_out
kun_peng classify \
  --db test_database \
  --chunk-dir temp_chunk \
  --output-dir test_out \
  data/COVID_19.fa
```

Notes
- `--chunk-dir` must be an empty or clean directory. It must not contain files like `sample_*.k2`, `sample_id*.map`, or `sample_*.bin`.
- `--output-dir` stores the Kraken-style outputs.
- Input supports FASTA/FASTQ and their `.gz` variants. You can pass multiple files.
- You can also pass a single `.txt` file that lists inputs (one path per line).

Useful options
- `-p, --num-threads <N>`: number of threads (default: system CPUs)
- `--buffer-size <BYTES>`: internal buffer sizing (default: 16 MiB)
- `--batch-size <N>`: controls memory when aggregating taxid matches (default: 4 or project default)
- `-T, --confidence-threshold <FLOAT>`: confidence threshold for reporting
- `-g, --minimum-hit-groups <N>`: minimum hit groups for a call (default: 2)
- `-P, --paired-end-processing`: enable paired-end mode
- `-S, --single-file-pairs`: indicates mates are in the same file

Paired-end examples
```bash
# Separate files per mate
kun_peng classify --db test_database --chunk-dir temp_chunk --output-dir test_out \
  -P read1.fq.gz read2.fq.gz

# Single file containing interleaved pairs
kun_peng classify --db test_database --chunk-dir temp_chunk --output-dir test_out \
  -P -S interleaved_pairs.fq.gz
```

## B. Direct Loading Mode (High-Memory, Max Speed)

Loads all hash tables into memory. Fastest, but requires RAM ≥ total size of `hash_*.k2d` files.

```bash
# Estimate memory requirement
bash cal_memory.sh test_database

# Run direct classification
kun_peng direct --db test_database data/COVID_19.fa
```

When to use
- You want maximum throughput.
- Your machine has enough RAM for the full database.

## Outputs

With the integrated workflow (`classify`), you will find in `--output-dir`:
- `output_*.txt`: Standard Kraken output per input file. Each line contains:
  1) C/U (classified or unclassified)
  2) sequence ID
  3) taxid (0 if unclassified)
  4) sequence length(s)
  5) space-delimited LCA mapping for k-mers
- `*.kreport2`: A hierarchical report with percentage, clade counts, direct counts, rank code, taxid, and name.

Timing logs show durations for splitr/annotate/resolve, helpful for performance checks.

## Tips and Common Issues

- Always use a clean `--chunk-dir`; the tool validates this and will error if leftover files are present.
- For large inputs, tune `-p`, `--buffer-size`, and `--batch-size` to balance speed and memory.
- Gzipped input is supported; no manual decompression required.
- You may run multiple inputs in one command; outputs are matched to inputs by index.
- If you recently added new FASTA to the database with `add-library`, run `build-db` before classifying.

## Minimal Examples

Single-end
```bash
mkdir -p temp_chunk test_out
kun_peng classify --db test_database --chunk-dir temp_chunk --output-dir test_out data/COVID_19.fa
```

Paired-end
```bash
mkdir -p temp_chunk test_out
kun_peng classify --db test_database --chunk-dir temp_chunk --output-dir test_out -P read1.fq.gz read2.fq.gz
```

Direct mode
```bash
bash cal_memory.sh test_database
kun_peng direct --db test_database data/COVID_19.fa
```

Back to: Detailed Database Build Demo (docs/build-db-demo.md)
