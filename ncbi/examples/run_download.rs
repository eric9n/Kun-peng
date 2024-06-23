use std::fs;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let workspace_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf();

    // Run the NCBI binary to download files
    let ncbi_binary = workspace_root.join("target/release/ncbi");
    let download_dir = workspace_root.join("downloads");
    // Ensure the download directory exists
    fs::create_dir_all(&download_dir).expect("Failed to create download directory");

    let args = vec![
        "-d".to_string(),
        download_dir.to_string_lossy().to_string(),
        "gen".to_string(),
        "-g".to_string(),
        "archaea".to_string(),
    ];

    let command_str = format!("{} {}", ncbi_binary.to_string_lossy(), args.join(" "));
    println!("Executing command: {}", command_str);

    // Run the NCBI binary to download files
    let output = Command::new(&ncbi_binary)
        .args(&args)
        .output()
        .expect("Failed to run NCBI binary");
    println!(
        "NCBI binary output: {}",
        String::from_utf8_lossy(&output.stdout)
    );

    let args = vec![
        "-d".to_string(),
        download_dir.to_string_lossy().to_string(),
        "tax".to_string(),
    ];

    let command_str = format!("{} {}", ncbi_binary.to_string_lossy(), args.join(" "));
    println!("Executing command: {}", command_str);

    // Run the NCBI binary to download files
    let output = Command::new(&ncbi_binary)
        .args(&args)
        .output()
        .expect("Failed to run NCBI binary");
    println!(
        "NCBI binary output: {}",
        String::from_utf8_lossy(&output.stdout)
    );
}
