#[path = "common/mod.rs"]
mod common;

use std::io;

fn main() -> io::Result<()> {
    let workspace_root = common::workspace_root();
    let data_dir = workspace_root.join("data");
    let database_dir = workspace_root.join("test_database");

    common::ensure_dir(&data_dir)?;
    common::ensure_dir(&database_dir)?;

    println!("Running full build + direct classification smoke tests\n");

    // Step 1: build database (idempotent if it already exists)
    let build_args = vec![
        "build".to_string(),
        "--download-dir".to_string(),
        common::path_to_string(&data_dir),
        "--db".to_string(),
        common::path_to_string(&database_dir),
    ];
    let build_output = common::run_kun_peng(&build_args)?;
    common::report_command("kun_peng build", &build_output);
    common::require_success("kun_peng build", &build_output)?;

    // Step 2: direct classification on a FASTA genome
    let covid_fasta = data_dir.join("COVID_19.fa");
    if !covid_fasta.exists() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "Example FASTA `{}` is missing — run the download step first.",
                covid_fasta.display()
            ),
        ));
    }

    let covid_args = vec![
        "direct".to_string(),
        "--db".to_string(),
        common::path_to_string(&database_dir),
        common::path_to_string(&covid_fasta),
    ];
    let covid_output = common::run_kun_peng(&covid_args)?;
    common::report_command("kun_peng direct (COVID-19)", &covid_output);
    common::require_success("kun_peng direct (COVID-19)", &covid_output)?;

    // Step 3: direct classification on an interleaved FASTQ file
    let interleaved_fastq = data_dir.join("test_interleaved.fastq");
    if !interleaved_fastq.exists() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "Example FASTQ `{}` is missing — run the download step first.",
                interleaved_fastq.display()
            ),
        ));
    }

    let interleaved_args = vec![
        "direct".to_string(),
        "--db".to_string(),
        common::path_to_string(&database_dir),
        common::path_to_string(&interleaved_fastq),
    ];
    let interleaved_output = common::run_kun_peng(&interleaved_args)?;
    common::report_command("kun_peng direct (interleaved FASTQ)", &interleaved_output);
    common::require_success("kun_peng direct (interleaved FASTQ)", &interleaved_output)?;

    Ok(())
}
