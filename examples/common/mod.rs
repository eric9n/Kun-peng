#![allow(dead_code)]

use std::env;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// Returns the root of the Cargo workspace.
pub fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Ensures that the given directory exists, creating it if needed.
pub fn ensure_dir(path: &Path) -> io::Result<()> {
    if path.exists() {
        return Ok(());
    }
    fs::create_dir_all(path)
}

/// Removes the directory if it exists and recreates an empty one.
pub fn recreate_dir(path: &Path) -> io::Result<()> {
    if path.exists() {
        fs::remove_dir_all(path)?;
    }
    fs::create_dir_all(path)
}

fn candidate_binary_paths() -> Vec<PathBuf> {
    let root = workspace_root();
    let exe_name = format!("kun_peng{}", env::consts::EXE_SUFFIX);
    vec![
        option_env!("CARGO_BIN_EXE_kun_peng")
            .map(PathBuf::from)
            .unwrap_or_default(),
        root.join("target").join("debug").join(&exe_name),
        root.join("target").join("release").join(&exe_name),
    ]
}

/// Locates the compiled `kun_peng` binary, searching common build targets.
pub fn kun_peng_binary() -> io::Result<PathBuf> {
    for path in candidate_binary_paths() {
        if !path.as_os_str().is_empty() && path.exists() {
            return Ok(path);
        }
    }

    Err(io::Error::new(
        io::ErrorKind::NotFound,
        "Unable to locate `kun_peng` binary. Run `cargo build --release` first.",
    ))
}

/// Converts any iterator of items implementing `ToString` into owned `String`s.
pub fn to_strings<I, T>(iter: I) -> Vec<String>
where
    I: IntoIterator<Item = T>,
    T: ToString,
{
    iter.into_iter().map(|item| item.to_string()).collect()
}

/// Converts a path into an owned `String`, using lossless conversion when possible.
pub fn path_to_string(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

/// Runs the `kun_peng` binary with the provided arguments.
pub fn run_kun_peng(args: &[String]) -> io::Result<Output> {
    let binary = kun_peng_binary()?;
    Command::new(binary).args(args.iter()).output()
}

/// Pretty prints the command output in a consistent format.
pub fn report_command(label: &str, output: &Output) {
    println!("=== {} ===", label);
    println!("status: {}", output.status);
    if !output.stdout.is_empty() {
        println!("stdout:\n{}", String::from_utf8_lossy(&output.stdout));
    }
    if !output.stderr.is_empty() {
        eprintln!("stderr:\n{}", String::from_utf8_lossy(&output.stderr));
    }
    println!();
}

/// Fails early if the command did not exit successfully, returning enriched error context.
pub fn require_success(label: &str, output: &Output) -> io::Result<()> {
    if output.status.success() {
        return Ok(());
    }

    let mut msg = format!("Command `{}` failed with status {}", label, output.status);
    if !output.stderr.is_empty() {
        msg.push_str("\n---- stderr ----\n");
        msg.push_str(&String::from_utf8_lossy(&output.stderr));
    }
    Err(io::Error::new(io::ErrorKind::Other, msg))
}
