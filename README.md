# kraken2-rust



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
* genbank: Download genbank resources.
* refseq: Download refseq resources.
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
