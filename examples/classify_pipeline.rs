#[path = "common/mod.rs"]
mod common;

use std::io;

fn main() -> io::Result<()> {
    let workspace_root = common::workspace_root();
    let database_dir = workspace_root.join("test_database");
    if !database_dir.exists() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "Database directory `{}` not found. Run the `build_and_classify` example first.",
                database_dir.display()
            ),
        ));
    }

    let reads = workspace_root.join("data").join("COVID_19.fa");
    if !reads.exists() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "Example reads `{}` missing. Re-run the build step to download them.",
                reads.display()
            ),
        ));
    }

    let base_examples_dir = workspace_root.join("target").join("examples");
    let chunk_dir = base_examples_dir.join("classify_chunks");
    let output_dir = base_examples_dir.join("classify_output");
    common::recreate_dir(&chunk_dir)?;
    common::recreate_dir(&output_dir)?;

    println!("Running staged classify pipeline against `{}`", reads.display());

    let classify_args = vec![
        "classify".to_string(),
        "--db".to_string(),
        common::path_to_string(&database_dir),
        "--chunk-dir".to_string(),
        common::path_to_string(&chunk_dir),
        "--output-dir".to_string(),
        common::path_to_string(&output_dir),
        common::path_to_string(&reads),
    ];

    let classify_output = common::run_kun_peng(&classify_args)?;
    common::report_command("kun_peng classify", &classify_output);
    common::require_success("kun_peng classify", &classify_output)?;

    println!(
        "Classification results are available under `{}`",
        output_dir.display()
    );

    Ok(())
}
