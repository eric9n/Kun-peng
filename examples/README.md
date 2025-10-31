# Example Programs

This directory contains a few runnable demos that mirror common Kun-peng workflows
and show how to interact with the taxonomy data programmatically.

All programs expect the CLI binary to be available. Build it once before running an
example:

```bash
cargo build --release
```

| Example | What it does |
| ------- | ------------ |
| `cargo run --example build_and_classify` | Rebuilds the bundled toy database and runs two direct mode classifications (FASTA + interleaved FASTQ). |
| `cargo run --example classify_pipeline` | Executes the full `kun_peng classify` pipeline and writes reports into `target/examples/`. |
| `cargo run --example taxonomy_inspect` | Loads the prebuilt taxonomy and prints a human-readable lineage for a few reference sequences. |

The shared helpers live in `examples/common/` and are compiled into each example.
They look for the binary under `target/{debug,release}/kun_peng` or via
`CARGO_BIN_EXE_kun_peng`, so custom build locations can be supported by exporting
that environment variable.
