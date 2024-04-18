# kraken2-rust

## 0.Installation Instructions
To install kraken2-rust, follow these steps:

1. Download the appropriate version for your system:

Navigate to the Releases page of the kraken2-rust GitHub repository.
Select the release suitable for your operating system. For example, if you are using CentOS 7, download kraken-rust-${VERSION}-centos7.tar.gz, where ${VERSION} is the version number of the release you wish to install.

2. Extract the downloaded archive:

Open a terminal.
Use the tar command to extract the files from the archive

```bash
tar -xvf kraken-rust-${VERSION}-centos7.tar.gz
```

## 1. NCBI download tool

Downloading Genome Data with NCBI Tool
The ncbi command-line tool offers functionality to download genome data from the NCBI database. This can be done for various groups including archaea, bacteria, viral, fungi, plant, human, protozoa, vertebrate_mammalian, vertebrate_other and invertebrate.


## Key Features
* Resumable Downloads: The tool supports breakpoint resumption, allowing downloads to pause and resume without starting over. This is particularly useful for large files or in conditions of unstable network connections.

* Incremental Download: Users can perform incremental downloads, where only new or updated files in the directory are downloaded. This saves time and bandwidth by avoiding redundant downloads of previously obtained data.

* Automatic MD5 Verification: To ensure data integrity, the tool automatically verifies the MD5 checksum of downloaded files. This step confirms that the files are not corrupted or tampered with during the download process.


### Genomes Command
To download genome data, use the genomes command with the following syntax:

```bash
ncbi genomes [OPTIONS] --group <GROUP> [COMMAND]
```

### Subcommands

* md5: Checks the md5 of the file only.
* fna: Parses genomic files and generates library fna files.
* assembly: Downloads and parses assembly files only.
* url: Downloads genomic files from a specified URL address.
* help: Print this message or the help of the given subcommand(s).

### Options
* --site <SITE>: Choose the NCBI site directory to download from (RefSeq or GenBank). Defaults to refseq. Possible values are:
*genbank*: Download genbank resources.
*refseq*: Download refseq resources.
*all*: Download genbank and refseq resources.

* --asm-level <ASM_LEVEL>: Set the assembly level for the download. Default is `basic`. ["Complete Genome", "Chromosome"]. `all` is ["Complete Genome", "Chromosome", "Scaffold", "Contig"].
* -g, --group <GROUP>: Specifies the category of data to download from the NCBI site. The group can be one or a comma-separated list of the following: archaea, bacteria, viral, fungi, plant, human, protozoa, vertebrate_mammalian, vertebrate_other, invertebrate.
* -h, --help: Print help information (for a summary, use '-h').

### Examples

To download genome data for bacteria from RefSeq:

```bash
ncbi genomes --group bacteria --site refseq
```

To check the md5 of genomic files for fungi:
```bash
ncbi genomes --group fungi md5
```

For more detailed help on a specific command, you can use the help subcommand:
```bash
ncbi help genomes
```

This tool simplifies the process of downloading and processing genome data from NCBI, making it accessible for various research and analysis purposes.


### Generate fna file

```bash
ncbi gen --site all -g archaea fna
```


## 2 Squid Tool

Squid is a versatile command-line tool designed for the efficient processing and classification of biological sequences. With its suite of functionalities, Squid facilitates various tasks related to sequence analysis, taxonomy resolution, and database management, making it an essential utility for bioinformatics workflows.

### Features
Squid offers a wide range of commands, each tailored for specific aspects of sequence data processing:

* estimate: Estimate the capacity requirements for database construction or analysis, aiding in resource planning.
* seqid2taxid: Generate a mapping file from sequence identifiers to taxonomic IDs, facilitating the association of sequences with their respective taxonomic lineage.
* build: Construct a Squid database from a collection of sequences, optimizing it for subsequent analysis tasks.
* hashshard: Divide a hash file into smaller, more manageable shards, improving the efficiency of data processing.
* splitr: Split FASTQ or FASTA files into ranges based on sequence identifiers or other criteria, aiding in the parallel processing of large datasets.
* annotate: Annotate a set of sequences with taxonomic or other relevant information, enriching the dataset for further analysis.
* resolve: Resolve the taxonomic tree for a set of sequences, identifying their positions within the taxonomic hierarchy.
* classify: A comprehensive workflow that integrates splitr, annotate, and resolve into a unified process for the classification of sequence data. This command streamlines the task of * assigning taxonomic classifications to sequences.

### Getting Started
To get started with Squid, you can invoke the tool with the -h or --help option to display detailed help messages for each command:


```bash
./squid -h
Usage: squid <COMMAND>

Commands:
  estimate     estimate capacity
  seqid2taxid  seqid to taxid map file
  build        build database
  hashshard    split hash file
  splitr       Split fast(q/a) file into ranges
  annotate     annotate a set of sequences
  resolve      resolve taxonomy tree
  classify     Integrates 'splitr', 'annotate', and 'resolve' into a unified workflow for sequence classification. classify a set of sequences
  help         Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```

This will provide you with an overview of all available commands and options. For specific information about a subcommand, use:


```bash
./squid <COMMAND> -h
```
Replace <COMMAND> with the name of the subcommand for which you need detailed help, such as estimate, build, or classify.


## 2.1 Seqid2taxid Tool

The seqid2taxid tool is a utility within the kr2r package designed to facilitate the mapping of sequence identifiers (seqid) to taxonomic identifiers (taxid). This tool is essential for researchers and bioinformaticians working with genomic data, enabling them to easily relate sequence data to taxonomic information.

### Usage

```bash
squid seqid2taxid -h

Usage: squid seqid2taxid [OPTIONS] --source <SOURCE>

Options:
      --source <SOURCE>      the database directory
  -f, --map-file <MAP_FILE>  seqid2taxid.map file path, default = $source/seqid2taxid.map
  -h, --help                 Print help
  -V, --version              Print version

```

To use the seqid2taxid tool, execute it with the required and optional arguments as follows:

```bash
squid seqid2taxid [OPTIONS] --source <SOURCE>
```

### Required Options
* --source <SOURCE>: Specifies the database directory containing the sequence and taxonomic data.

### Optional Options
* -f, --map-file <MAP_FILE>: Path to the seqid2taxid.map file. If not specified, the tool defaults to using $source/seqid2taxid.map, where $source is the path provided by the * --source option.
* -h, --help: Displays help information about the tool and its options.
* -V, --version: Prints the version of the seqid2taxid tool.

### Example Command
To run the seqid2taxid tool with a specific source directory:

```bash
squid seqid2taxid --source /path/to/database
```

To specify a custom map file path:

```bash
squid seqid2taxid --source /path/to/database -f /path/to/custom/seqid2taxid.map
```

## 2.2 Estimate Capacity Tool

The estimate_capacity tool is designed for estimating the capacity required for building a database from genomic data. It takes into consideration various parameters related to the genomic data processing, such as k-mer length, minimizer length, and hash table load factor, to provide an efficient estimate of the necessary resources.

### Usage

To use the estimate_capacity tool, execute it from the command line with the desired options:

```bash
squid estimate_capacity [OPTIONS]
```

Options
* --source <SOURCE>: Specifies the build database directory or file. Default is lib.
* --cache: Estimates capacity from cache if exists.
* -k, --k-mer <K_MER>: Sets the length of k-mers. K must be a positive integer (default is 35). K cannot be less than L.
* -l, --l-mer <L_MER>: Sets the length of minimizers. L must be between 1 and 31 (default is 31).
* -n, --n <N>: Sets the maximum qualifying hash code (default is 4).
* --minimizer-spaces <MINIMIZER_SPACES>: Specifies the number of characters in the minimizer that are ignored in comparisons (default is 7).
* -T, --toggle-mask <TOGGLE_MASK>: Defines the minimizer ordering toggle mask.
* --load-factor <LOAD_FACTOR>: Sets the proportion of the hash table to be populated (only for build task; default is 0.7, must be between 0 and 1).
* -p, --threads <THREADS>: Specifies the number of threads to use (default is 4).
* -h, --help: Prints help information (for more details, use '--help').
* -V, --version: Prints the version of the tool.

### Example

```bash
squid estimate_capacity -k 35 -l 31 --source /data/ncbi/path -p 10 --load-factor 0.7
```

### Output

```bash
estimate count: 1213069985, required capacity: 1732968825.0, Estimated hash table requirement: 6.46GB
```


## 2.3 build

```bash
./squid build -h
build database

Usage: squid build [OPTIONS] --source <SOURCE> -H <HASHTABLE_FILENAME> -o <OPTIONS_FILENAME> -t <TAXONOMY_FILENAME> -m <ID_TO_TAXON_MAP_FILENAME> --ncbi-taxonomy-directory <NCBI_TAXONOMY_DIRECTORY> --required-capacity <REQUIRED_CAPACITY> --chunk-dir <CHUNK_DIR>

Options:
      --source <SOURCE>
          build database directory or file
  -H <HASHTABLE_FILENAME>
          Kraken 2 hash table filename
  -o <OPTIONS_FILENAME>
          Kraken 2 options filename
  -k, --k-mer <K_MER>
          Set length of k-mers, k must be positive integer, k=35, k cannot be less than l [default: 35]
  -l, --l-mer <L_MER>
          Set length of minimizers, 1 <= l <= 31 [default: 31]
      --minimizer-spaces <MINIMIZER_SPACES>
          Number of characters in minimizer that are ignored in comparisons [default: 7]
  -T, --toggle-mask <TOGGLE_MASK>
          Minimizer ordering toggle mask [default: 16392584516609989165]
      --min-clear-hash-value <MIN_CLEAR_HASH_VALUE>

  -r, --requested-bits-for-taxid <REQUESTED_BITS_FOR_TAXID>
          Bit storage requested for taxid 0 <= r < 31 [default: 0]
  -p, --threads <THREADS>
          Number of threads [default: 4]
  -t <TAXONOMY_FILENAME>
          Kraken 2 taxonomy filename
  -m <ID_TO_TAXON_MAP_FILENAME>
          Sequence ID to taxon map filename
  -n, --ncbi-taxonomy-directory <NCBI_TAXONOMY_DIRECTORY>
          NCBI taxonomy directory name
  -c, --required-capacity <REQUIRED_CAPACITY>

      --chunk-dir <CHUNK_DIR>
          chunk directory
      --chunk-size <CHUNK_SIZE>
          chunk size 1-4(GB) [default: 1073741824]
      --cache
          estimate capacity from cache if exists
      --max-n <MAX_N>
          Set maximum qualifying hash code [default: 4]
      --load-factor <LOAD_FACTOR>
          Proportion of the hash table to be populated (build task only; def: 0.7, must be between 0 and 1) [default: 0.7]
  -h, --help
          Print help
  -V, --version
          Print version
```

## 2.4 hashshard

```bash
./squid hashshard -h
split hash file

Usage: squid hashshard [OPTIONS] --db <DB>

Options:
      --db <DB>                The database directory for the Kraken 2 index. contains index file(hash.k2d opts.k2d taxo.k2d)
      --hash-dir <HASH_DIR>    database hash chunk directory and other files
      --hash-capacity <HASH_CAPACITY> default: 1073741824(capacity 1G = file size 4G)
  -h, --help                   Print help (see more with '--help')
  -V, --version                Print version
```

## 2.5 splitr

```bash
./squid splitr -h
Split fast(q/a) file into ranges

Usage: squid splitr [OPTIONS] --hash-dir <HASH_DIR> --chunk-dir <CHUNK_DIR> [INPUT_FILES]...

Arguments:
  [INPUT_FILES]...  A list of input file paths (FASTA/FASTQ) to be processed by the classify program

Options:
      --hash-dir <HASH_DIR>
          database hash chunk directory and other files
  -P, --paired-end-processing
          Enable paired-end processing
  -S, --single-file-pairs
          Process pairs with mates in the same file
  -Q, --minimum-quality-score <MINIMUM_QUALITY_SCORE>
          Minimum quality score for FASTQ data, default is 0 [default: 0]
  -p, --num-threads <NUM_THREADS>
          The number of threads to use, default is 1 [default: 10]
      --chunk-dir <CHUNK_DIR>
          chunk directory
  -h, --help
          Print help (see more with '--help')
  -V, --version
          Print version
```

## 2.6 annotate

```bash
annotate a set of sequences

Usage: squid annotate [OPTIONS] --hash-dir <HASH_DIR> --chunk-dir <CHUNK_DIR>

Options:
      --hash-dir <HASH_DIR>      database hash chunk directory and other files
      --chunk-dir <CHUNK_DIR>    The file path for the Kraken 2 options. chunk directory
      --batch-size <BATCH_SIZE>  批量处理大小 default: 8MB [default: 8388608]
  -h, --help                     Print help (see more with '--help')
  -V, --version                  Print version
```


## 2.7 resolve

```bash
resolve taxonomy tree

Usage: squid resolve [OPTIONS] --hash-dir <HASH_DIR> --chunk-dir <CHUNK_DIR>

Options:
      --hash-dir <HASH_DIR>
          database hash chunk directory and other files
      --chunk-dir <CHUNK_DIR>
          chunk directory
      --full-output
          output file contains all unclassified seq
  -T, --confidence-threshold <CONFIDENCE_THRESHOLD>
          Confidence score threshold, default is 0.0 [default: 0]
  -K, --report-kmer-data
          In comb. w/ -R, provide minimizer information in report
  -z, --report-zero-counts
          In comb. w/ -R, report taxa w/ 0 count
  -g, --minimum-hit-groups <MINIMUM_HIT_GROUPS>
          The minimum number of hit groups needed for a call [default: 2]
      --batch-size <BATCH_SIZE>
          批量处理大小 default: 8MB [default: 8388608]
      --output-dir <KRAKEN_OUTPUT_DIR>
          File path for outputting normal Kraken output
  -h, --help
          Print help (see more with '--help')
  -V, --version
          Print version
```


## 2.8 classify

```bash
./squid classify -h
Integrates 'splitr', 'annotate', and 'resolve' into a unified workflow for sequence classification. classify a set of sequences

Usage: squid classify [OPTIONS] --hash-dir <HASH_DIR> --chunk-dir <CHUNK_DIR> [INPUT_FILES]...

Arguments:
  [INPUT_FILES]...  A list of input file paths (FASTA/FASTQ) to be processed by the classify program

Options:
      --hash-dir <HASH_DIR>
          database hash chunk directory and other files
  -P, --paired-end-processing
          Enable paired-end processing
  -S, --single-file-pairs
          Process pairs with mates in the same file
  -Q, --minimum-quality-score <MINIMUM_QUALITY_SCORE>
          Minimum quality score for FASTQ data, default is 0 [default: 0]
  -p, --num-threads <NUM_THREADS>
          The number of threads to use, default is 1 [default: 10]
      --chunk-dir <CHUNK_DIR>
          chunk directory
      --batch-size <BATCH_SIZE>
          批量处理大小 default: 8MB [default: 8388608]
  -T, --confidence-threshold <CONFIDENCE_THRESHOLD>
          Confidence score threshold, default is 0.0 [default: 0]
  -g, --minimum-hit-groups <MINIMUM_HIT_GROUPS>
          The minimum number of hit groups needed for a call [default: 2]
      --output-dir <KRAKEN_OUTPUT_DIR>
          File path for outputting normal Kraken output
  -K, --report-kmer-data
          In comb. w/ -R, provide minimizer information in report
  -z, --report-zero-counts
          In comb. w/ -R, report taxa w/ 0 count
      --full-output
          output file contains all unclassified sequence
  -h, --help
          Print help (see more with '--help')
  -V, --version
          Print version
```

## 3. build_k2_db

The build_k2_db command-line tool facilitates the construction of a Kraken 2 database. It requires specific filenames for the hash table, taxonomy, and the sequence ID to taxon map, among other parameters.


### Introduction
The build_k2_db tool introduces a novel approach to constructing Kraken 2-compatible databases, specifically addressing the challenges associated with the large memory requirements of previous methods. This documentation outlines the process flow, working principles, and the inherent advantages of using the build_k2_db tool for genomic database construction.

The build_k2_db tool revolutionizes the process of building genomic databases for Kraken 2 by introducing a novel, two-step approach to database construction. This method significantly mitigates the challenges associated with the large memory requirements of traditional database building processes, particularly vital for constructing databases like the NCBI RefSeq, which are substantial in size.

### Working Principle
#### Step 1: Preprocessing and Generation of k2 Formatted Files
Initially, the tool preprocesses .fna files to generate intermediary files in a k2 format. This step involves scanning the .fna files to extract relevant k-mer and minimizer information, mapping these to taxonomic IDs, and then hashing these elements to produce indexed intermediary data. These intermediary files are crucial for the next step of the process, as they contain indexed positions and taxonomic IDs necessary for constructing the hash table efficiently.

#### Step 2: Iterative Construction of the Hash Table
In the second phase, the tool iteratively processes the k2 formatted intermediary files to build segments of the hash table. This method involves reading the intermediary files in batches, resolving any taxonomic ID conflicts using a Lowest Common Ancestor (LCA) algorithm, and updating the hash table with the resolved IDs. This step-by-step processing significantly reduces the memory footprint compared to loading the entire hash table into memory at once.

#### Efficiency and Advantages
The build_k2_db tool introduces several advantages over traditional database building methods:

* Memory Efficiency: By generating intermediary files and processing these in chunks, the tool drastically reduces the required memory, enabling the construction of large databases on systems with limited memory capacity.
* Scalability: The approach is highly scalable, thanks to parallel processing and efficient handling of large .fna files, making it suitable for building extensive databases.
Time Efficiency: Despite the intermediary files being substantially larger than the final hash table, the overall time taken to build the database is comparable to methods that process all data at once.
Performance Insights
In a performance test involving the NCBI RefSeq database, approximately 500GB of .fna files were processed to generate 850GB of k2 intermediary files. The final hash table size amounted to 188GB. Utilizing a machine equipped with a 16-core CPU and 32GB of memory, the entire database construction process was completed in just 9 hours and 42 minutes. This showcases the tool's ability to handle large datasets efficiently, both in terms of time and hardware resource requirements.

#### Comparative Analysis with Kraken 2 C++ Version in Fast Mode

In addition to the innovative build_k2_db Rust-based tool, it's informative to compare its performance and resource utilization with that of the traditional Kraken 2 C++ version, particularly in its fast mode operation. Such a comparison underscores the advancements and efficiencies introduced by the Rust implementation.

#### Kraken 2 C++ Version in Fast Mode:
For processing the same dataset from the NCBI RefSeq database (~500GB of .fna files), the Kraken 2 C++ version in fast mode presents the following resource requirements and performance metrics:

CPU and Memory Usage: Requires a machine with a 16-core CPU and 200GB of memory, indicating a significantly higher demand for memory resources compared to the Rust-based build_k2_db tool.
Time Efficiency: Completes the database construction process in approximately 9 hours and 32 minutes. This duration is slightly shorter than that of the build_k2_db tool but at the cost of substantially higher memory requirements.


#### Key Insights and Implications:
Memory Optimization: The build_k2_db tool demonstrates exceptional memory efficiency by requiring only 32GB of memory to process and construct a database from a large genomic dataset. In contrast, the C++ version's fast mode requires 200GB of memory, highlighting the Rust-based tool's optimization in memory usage.
Comparable Time Efficiency: Despite the vast difference in memory consumption, the time taken to build the database is remarkably similar between the two tools, with the Rust version completing the task in 9 hours and 42 minutes versus 9 hours and 32 minutes for the C++ version.
Accessibility and Cost-effectiveness: By drastically reducing the memory requirement, the build_k2_db tool makes the process of building large genomic databases more accessible to researchers and institutions with limited hardware resources. This can significantly lower the computational costs associated with database construction in bioinformatics research.


####  Conclusion
The build_k2_db tool stands out for its innovative approach to genomic database construction, offering a memory-efficient, scalable, and time-effective solution. Its ability to preprocess data into intermediary files before iteratively constructing the hash table addresses the significant challenges of working with large-scale genomic databases, making it an invaluable asset in the field of bioinformatics.

The build_k2_db tool not only matches the Kraken 2 C++ version in terms of processing time but does so with far less memory, making it a highly efficient and accessible option for constructing large genomic databases. Its innovative approach, leveraging Rust's performance and memory management capabilities, offers a more practical solution for the bioinformatics community, particularly when handling extensive datasets like the NCBI RefSeq database.



### Usage
To build the Kraken 2 database, you must specify source, hash table, taxonomy, ID to taxon map filenames, Kraken 2 options filename, NCBI taxonomy directory, required capacity, and chunk directory.

```bash
build_k2_db [OPTIONS] --source <SOURCE> -H <HASHTABLE_FILENAME> -t <TAXONOMY_FILENAME> -m <ID_TO_TAXON_MAP_FILENAME> -o <OPTIONS_FILENAME> --ncbi-taxonomy-directory <NCBI_TAXONOMY_DIRECTORY> --required-capacity <REQUIRED_CAPACITY> --chunk-dir <CHUNK_DIR>
```

### Options

* --source <SOURCE>: Directory or file for database build.
* -H <HASHTABLE_FILENAME>: Filename for the Kraken 2 hash table.
* -t <TAXONOMY_FILENAME>: Filename for the Kraken 2 taxonomy.
* -m <ID_TO_TAXON_MAP_FILENAME>: Filename for the sequence ID to taxon map.
* -o <OPTIONS_FILENAME>: Filename for Kraken 2 options.
* -n, --ncbi-taxonomy-directory <NCBI_TAXONOMY_DIRECTORY>: Directory name for NCBI taxonomy.
* -k, --k-mer <K_MER>: Length of k-mers (default: 35).
* -l, --l-mer <L_MER>: Length of minimizers (default: 31).
* -r, --requested-bits-for-taxid <REQUESTED_BITS_FOR_TAXID>: Bit storage for taxid (default: 0).
* -T, --toggle-mask <TOGGLE_MASK>: Minimizer ordering toggle mask (default: 16392584516609989165).
* --minimizer-spaces <MINIMIZER_SPACES>: Characters in minimizer ignored in comparisons (default: 7).
* -c, --required-capacity <REQUIRED_CAPACITY>: Required capacity for the database.
* -p, --threads <THREADS>: Number of threads (default: 4).
* --chunk-dir <CHUNK_DIR>: Directory for chunks.
* --chunk-size <CHUNK_SIZE>: Size of chunks in GB (default: 1GB).
* --chunk-prefix <CHUNK_PREFIX>: Prefix for chunk files (default: chunk).
* --only-k2: Process only k2 file.
* -h, --help: Prints help information.
* -V, --version: Prints the version of the tool.

### Example

Building a database with custom parameters:

```bash
build_k2_db --source /path/to/source -H hash_table.k2 -t taxonomy.k2 -m id_to_taxon.map -o options.k2 --ncbi-taxonomy-directory /path/to/ncbi/taxonomy --required-capacity 1000000 --chunk-dir /path/to/chunks
```


## 4. classify

The classify tool is a powerful sequence classification program designed for rapid and accurate classification of nucleotide sequences. It leverages the Kraken 2 indexing and taxonomy systems to efficiently assign taxonomic labels to sequences from FASTA/FASTQ files. This document provides a comprehensive guide on how to use the classify tool, including its options and arguments.

### Usage

To classify sequences using the classify tool, execute the command with the required options and input files:

```bash
classify [OPTIONS] --index-filename <INDEX_FILENAME> --taxonomy-filename <TAXONOMY_FILENAME> --options-filename <OPTIONS_FILENAME> [INPUT_FILES]...
```

#### Arguments

* [INPUT_FILES]...: Specifies a list of input file paths. These files should be in FASTA or FASTQ format and contain the sequences to be classified.

#### Options

* -H, --index-filename <INDEX_FILENAME>: Path to the Kraken 2 index file. This file is essential for the classification process.
* -t, --taxonomy-filename <TAXONOMY_FILENAME>: Path to the Kraken 2 taxonomy file. This file contains taxonomic information used for classification.
* -o, --options-filename <OPTIONS_FILENAME>: Path to the Kraken 2 options file. This file includes additional configuration options for Kraken 2.
* -T, --confidence-threshold <CONFIDENCE_THRESHOLD>: Sets the confidence score threshold for classification. Sequences with a confidence score below this threshold will not be * classified. The default value is 0.0.
* -p, --num-threads <NUM_THREADS>: Specifies the number of threads to use for processing. Increasing the number of threads can speed up the classification process. The default is 1.
* -g, --minimum-hit-groups <MINIMUM_HIT_GROUPS>: The minimum number of hit groups required for a classification call. The default is 2.
* -P, --paired-end-processing: Enables processing of paired-end reads. This option should be used if your input files contain paired-end sequence data.
* -S, --single-file-pairs: Indicates that pairs with mates are located in the same file. This option is relevant for paired-end processing.
* -O, --kraken-output-filename <KRAKEN_OUTPUT_FILENAME>: Specifies the file path for outputting the standard Kraken output. This output includes the classification results for * each sequence.
* -Q, --minimum-quality-score <MINIMUM_QUALITY_SCORE>: Sets the minimum quality score for FASTQ data. Sequences with a quality score below this threshold will not be classified. * The default is 0.
* -h, --help: Prints help information, providing a summary of options and usage.
* -V, --version: Displays the version of the classify tool.

#### Example

To classify sequences from a FASTQ file using 4 threads and a confidence threshold of 0.5:

```bash
classify --index-filename path/to/index --taxonomy-filename path/to/taxonomy --options-filename path/to/options -T 0.5 -p 4 input_file.fastq
```


## 5. inspect

The inspect tool is designed for analyzing the content of hash table files used by Kraken 2. It provides insights into the index file, allowing users to verify and understand the structure and statistics of their Kraken 2 databases.

### Usage
To utilize the inspect tool, execute the command with the necessary options:

```bash
inspect [OPTIONS] --index-filename <INDEX_FILENAME>
```

### Options
* -H, --index-filename <INDEX_FILENAME>: Specifies the file path to the Kraken 2 index file. This option is required as it directs the tool to the hash table file that needs to be inspected.
* -t, --taxonomy-filename <TAXONOMY_FILENAME>: Provides the file path to the Kraken 2 taxonomy file. This file contains the taxonomy information that corresponds to the data in the index file. Including this option allows for a more comprehensive inspection that may involve taxonomy data.
* -o, --options-filename <OPTIONS_FILENAME>: Indicates the file path to the Kraken 2 options file. This file can contain various configurations and options used by Kraken 2. * Specifying this option can help understand the configurations under which the index was created or used.
* -v, --value-count: This flag, when set, instructs the tool to iterate through the index file and count the values. It is useful for users who wish to understand the * distribution of data within their Kraken 2 index file.
* -h, --help: Prints out help information, providing a brief summary of all the available options and their usage.
* -V, --version: Displays the version of the inspect tool, helping users to identify the tool's version they are currently using.
