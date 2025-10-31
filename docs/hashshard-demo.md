# Kun-peng Hashshard: Convert a Kraken2 Database

This demo shows how to convert an existing Kraken 2 database directory into Kun-peng’s sharded hash format using the `hashshard` subcommand.

## When to Use

- You already have a Kraken 2 database folder containing `hash.k2d`, `opts.k2d`, and `taxo.k2d` (standard Kraken 2 files), and you want to use that database with Kun-peng.
- You want to shard the monolithic `hash.k2d` into multiple `hash_*.k2d` chunks to support Kun-peng’s memory- and I/O-efficient workflows.

## Prerequisites

- A Kraken 2 database directory, e.g. `/path/to/kraken_db`, containing:
  - `hash.k2d`
  - `opts.k2d`
  - `taxo.k2d`
- The `kun_peng` binary available (PATH or `./target/release/kun_peng`).

## Basic Conversion

```bash
# Choose a target directory (can be the same as the Kraken 2 DB directory)
DB=/path/to/kraken_db

# Shard size (capacity per shard); index size ≈ 4×capacity
# Example: 1G capacity → each shard file ≈ 4G
kun_peng hashshard --db "$DB" --hash-capacity 1G
```

What happens
- Reads `hash.k2d` and writes:
  - `hash_config.k2d` with partition metadata
  - `hash_1.k2d, hash_2.k2d, …` shards
- Copies `opts.k2d` and `taxo.k2d` into the directory if missing.

Important
- If `hash_config.k2d` already exists in the target directory, the command will abort to avoid accidental overwrites. Use a fresh directory or remove/backup existing `hash_config.k2d` before running.

## Choosing `--hash-capacity`

- Capacity is the number of slots per shard (not bytes). File size per shard is roughly `capacity × 4 bytes`.
- Example sizing:
  - `--hash-capacity 1G` → ~4 GiB per shard file
  - `--hash-capacity 512M` → ~2 GiB per shard file
- More, smaller shards can improve I/O parallelism and reduce per-file memory footprints, with a modest overhead in file count.

## Verifying Outputs

```bash
ls -lh "$DB"
# Expect: hash_config.k2d, hash_*.k2d, opts.k2d, taxo.k2d
```

## Classifying with the Converted Database

Integrated workflow:
```bash
mkdir -p temp_chunk test_out
kun_peng classify --db "$DB" --chunk-dir temp_chunk --output-dir test_out data/COVID_19.fa
```

Direct loading (high memory, fastest):
```bash
bash cal_memory.sh "$DB"   # estimates RAM needed (≈ sum of hash_*.k2d sizes)
kun_peng direct --db "$DB" data/COVID_19.fa
```

## Tips and Pitfalls

- If you see an error about `hash_config.k2d` existing, move or remove it, or choose another output directory.
- `--hash-capacity` affects how many shard files are created (total slots ÷ capacity). Tune it to balance shard size and file count.
- After conversion, you can use all Kun-peng classify modes without rebuilding from source sequences.

## See Also

- Database Build Demo: docs/build-db-demo.md
- Classification Demo: docs/classify-demo.md

