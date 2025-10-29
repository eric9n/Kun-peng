use std::fs;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    // Define the paths and directories
    let workspace_root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).to_path_buf();
    let kr2r_binary = workspace_root.join("target/release/kun_peng");
    let data_dir = workspace_root.join("data");
    let test_dir = workspace_root.join("test_database");

    // Ensure the necessary directories exist
    fs::create_dir_all(&data_dir).expect("Failed to create download directory");
    fs::create_dir_all(&test_dir).expect("Failed to create database directory");

    // Command 1: ./target/release/kun_peng build --download-dir data/ --db test_database
    let build_args = vec![
        "build".to_string(),
        "--download-dir".to_string(),
        data_dir.to_string_lossy().to_string(),
        "--db".to_string(),
        test_dir.to_string_lossy().to_string(),
    ];

    let build_command_str = format!("{} {}", kr2r_binary.to_string_lossy(), build_args.join(" "));
    println!("Executing command: {}", build_command_str);

    let build_output = Command::new(&kr2r_binary)
        .args(&build_args)
        .output()
        .expect("Failed to run kun_peng build command");
    println!(
        "kun_peng build output: {}",
        String::from_utf8_lossy(&build_output.stdout)
    );
    if !build_output.stderr.is_empty() {
        println!(
            "kun_peng build error: {}",
            String::from_utf8_lossy(&build_output.stderr)
        );
    }

    // Command 2: ./target/release/kun_peng direct --db test_database data/COVID_19.fa
    let covid_fa = data_dir.join("COVID_19.fa");
    if !covid_fa.exists() {
        println!(
            "kun_peng error: fasta file {} does not exists",
            covid_fa.to_string_lossy().to_string()
        );
    }
    let direct_args = vec![
        "direct".to_string(),
        "--db".to_string(),
        test_dir.to_string_lossy().to_string(),
        covid_fa.to_string_lossy().to_string(),
    ];

    let direct_command_str = format!(
        "{} {}",
        kr2r_binary.to_string_lossy(),
        direct_args.join(" ")
    );
    println!("Executing command: {}", direct_command_str);

    let direct_output = Command::new(&kr2r_binary)
        .args(&direct_args)
        .output()
        .expect("Failed to run kun_peng direct command");
    println!(
        "kun_peng direct output: {}",
        String::from_utf8_lossy(&direct_output.stdout)
    );
    if !direct_output.stderr.is_empty() {
        println!(
            "kun_peng direct error: {}",
            String::from_utf8_lossy(&direct_output.stderr)
        );
    }

    // --- Test 3: ./target/release/kun_peng direct --db test_database data/test_interleaved.fastq ---
    println!("\n--- Testing Interleaved FASTQ Classification ---");

    let interleaved_fq = data_dir.join("test_interleaved.fastq");
    if !interleaved_fq.exists() {
        println!(
            "kun_peng error: fastq file {} does not exists",
            interleaved_fq.to_string_lossy().to_string()
        );
        // 考虑在这里 panic，如果文件不存在，测试就没有意义
        // panic!("Test FASTQ file not found!"); 
    }
    let direct_fastq_args = vec![
        "direct".to_string(),
        "--db".to_string(),
        test_dir.to_string_lossy().to_string(),
        interleaved_fq.to_string_lossy().to_string(),
    ];

    let direct_fastq_command_str = format!(
        "{} {}",
        kr2r_binary.to_string_lossy(),
        direct_fastq_args.join(" ")
    );
    println!("Executing command: {}", direct_fastq_command_str);

    let direct_fastq_output = Command::new(&kr2r_binary)
        .args(&direct_fastq_args)
        .output()
        .expect("Failed to run kun_peng direct command with interleaved FASTQ");
    
    println!(
        "kun_peng direct (fastq) output: {}",
        String::from_utf8_lossy(&direct_fastq_output.stdout)
    );
    if !direct_fastq_output.stderr.is_empty() {
        println!(
            "kun_peng direct (fastq) error: {}",
            String::from_utf8_lossy(&direct_fastq_output.stderr)
        );
    }
}
