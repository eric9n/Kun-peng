#!/bin/bash

# === Configuration ===
readonly OUTPUT_FILENAME="assembly_summary_temp.txt"
readonly PLACEHOLDER="na"
readonly NUM_COLS=22 # Expected number of columns for NCBI assembly summary format

# === Function for error logging and exit ===
error_exit() {
    echo "Error: $1" >&2
    exit "${2:-1}" # Default exit code 1
}

# === Usage instructions ===
print_usage() {
    echo "Usage: $0 <added_directory_containing_fna_and_or_fna_gz_files> <target_directory>"
    echo ""
    echo "Important Notes:"
    echo "1. The taxid extracted from FASTA headers must exist in the taxonomy database:"
    echo "   - nodes.dmp: Contains the taxonomic tree structure"
    echo "   - names.dmp: Contains the scientific names of taxa"
    echo "2. If the taxid is not found in these files, you will need to:"
    echo "   - Edit nodes.dmp to add the taxid to the taxonomic tree"
    echo "   - Edit names.dmp to add the scientific name for the taxid"
    echo "3. The script will only process files where the first string after '>' is a valid taxid number"
    echo ""
    echo "Example FASTA file format (example.fna):"
    echo ">2697049 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome"
    echo "ATGTCTGATAATGGACCCCAAAATCAGCGAAATGCACCCCGCATTACGTTTGGTGGACCCTCAGATTCAACTGGCAGTAACCAGAATGGAGAACGC"
    echo "AGTGGGGCGCGATCAAACAACTGCCCGGGGGCACTGGCCTCACCTTCATGGACAGGGTCTTGTCACTCAGGGTCTGGAGCCAGTGCTGGTGTCAAA"
    echo "..."
    echo ""
    echo "IMPORTANT: The first number (2697049 in this example) MUST be a valid taxid."
    echo "If your FASTA header starts with a different identifier (e.g., NC_045512.2),"
    echo "you need to modify it to use the correct taxid before processing."
    echo ""
    echo "Example:"
    echo "  $0 /path/source /path/to/target"
    exit 1
}

# --- Check required tools ---
# Function to check command existence
check_command() {
    if ! command -v "$1" >/dev/null; then
        error_exit "'$1' command not found, required by this script."
    fi
}

check_command head
check_command gzip
check_command grep
check_command basename
check_command find
check_command sed

# === Parameter handling ===
if [ -z "$1" ] || [ -z "$2" ]; then
    print_usage
fi

readonly ADDED_DIR="$1"
readonly TARGET_DIR="$2"

# === Directory validation ===
if [ ! -d "$ADDED_DIR" ]; then
    error_exit "Directory '$ADDED_DIR' does not exist."
fi

if [ ! -d "$TARGET_DIR" ]; then
    error_exit "Directory '$TARGET_DIR' does not exist."
fi

# Create directories if they don't exist
readonly TARGET_ADDED_DIR="$TARGET_DIR/added"
readonly TARGET_TEMP_DIR="$TARGET_ADDED_DIR/temp"
mkdir -p "$TARGET_ADDED_DIR" || error_exit "Cannot create directory '$TARGET_ADDED_DIR'."
mkdir -p "$TARGET_TEMP_DIR" || error_exit "Cannot create directory '$TARGET_TEMP_DIR'."

readonly OUTPUT_FILE="$ADDED_DIR/$OUTPUT_FILENAME"

# === Prepare output file (write header) ===
echo "Generating summary file '$OUTPUT_FILE'..."

# Define header columns
readonly HEADER=(
    "# assembly_accession" "bioproject" "biosample" "wgs_master"
    "refseq_category" "taxid" "species_taxid" "organism_name"
    "infraspecific_name" "isolate" "version_status" "assembly_level"
    "release_type" "genome_rep" "seq_rel_date" "asm_name"
    "submitter" "wgs_project" "assembly_type" "ftp_path"
    "gbrs_paired_asm" "paired_asm_comp"
)

# Check header column count
if [ ${#HEADER[@]} -ne $NUM_COLS ]; then
    echo "Warning: Script's internal header has ${#HEADER[@]} columns, but expected $NUM_COLS columns." >&2
fi

# Write header (overwrite file)
{
    printf "%s\t" "${HEADER[@]}" | sed 's/\t$//'
    echo ""
} > "$OUTPUT_FILE" || error_exit "Cannot write header to '$OUTPUT_FILE'."

# === Process files ===
file_count=0
processed_count=0
declare -a processed_files  # Array to store successfully processed files

# Use find to locate *.fna or *.fna.gz files
# Note the use of \( ... \) to group -o (OR) conditions
while IFS= read -r fna_file; do
    ((file_count++))
    filename_only="${fna_file##*/}"
    first_line=""
    base_name=""
    file_type="" # For logging purposes

    # --- Choose processing method based on file extension ---
    if [[ "$fna_file" == *.fna.gz ]]; then
        file_type="gzipped"
        echo "Processing [$file_type] [${file_count}]: '$filename_only'..."
        # Decompress and read first line
        first_line=$(gzip -dc "$fna_file" 2>/dev/null | head -n 1)
        # Get base name, remove .fna.gz suffix
        base_name=$(basename "$fna_file" .fna.gz)

    elif [[ "$fna_file" == *.fna ]]; then
        file_type="plain"
        echo "Processing [$file_type] [${file_count}]: '$filename_only'..."
        # Read first line directly
        first_line=$(head -n 1 "$fna_file")
        # Get base name, remove .fna suffix
        base_name=$(basename "$fna_file" .fna)
    else
        # This branch should theoretically never execute as find has already filtered file types
        echo "  Warning: Skipping unrecognized file type '$filename_only'." >&2
        continue
    fi

    # --- Common processing logic (after getting first_line and base_name) ---

    # Check if read/decompress was successful or if file is empty
    if [ -z "$first_line" ]; then
        echo "  Warning: File '$filename_only' ($file_type) is empty or first line could not be read. Skipped." >&2
        continue
    fi

    # Check if first line is a FASTA header
    if [[ ! "$first_line" == ">"* ]]; then
        echo "  Warning: File '$filename_only' ($file_type) does not start with '>'. Skipped." >&2
        continue
    fi

    # Extract TaxID
    taxid=$(echo "$first_line" | awk '{sub(/^>/, ""); print $1}')

    if [[ ! "$taxid" =~ ^[1-9][0-9]*$ ]]; then
        echo "  Warning: Invalid or missing taxid in header of '$filename_only' ($file_type). Header: '$first_line'. Skipped." >&2
        continue
    fi

    # Check if base_name was successfully obtained
    if [ -z "$base_name" ]; then
        echo "  Warning: Could not get base name for file '$filename_only' ($file_type). Skipped." >&2
        continue
    fi

    # Build output fields array
    declare -a fields
    for ((i=0; i<NUM_COLS; i++)); do
        fields[i]="$PLACEHOLDER"
    done

    # Fill required fields
    fields[5]="$taxid"       # taxid
    fields[15]="$base_name"  # asm_name
    fields[19]="$base_name"  # ftp_path

    # Append line to output file
    {
        printf "%s\t" "${fields[@]}" | sed 's/\t$//'
        echo ""
    } >> "$OUTPUT_FILE" || echo "  Warning: Could not write entry for '$filename_only' ($file_type) to '$OUTPUT_FILE'." >&2

    echo "  -> Successfully extracted taxid: $taxid. Entry added."
    ((processed_count++))
    processed_files+=("$fna_file")  # Add successfully processed file to array

done < <(find "$ADDED_DIR" -maxdepth 1 \( -name '*.fna' -o -name '*.fna.gz' \) -type f)

# === Copy files to target directory ===
if [ "$processed_count" -gt 0 ]; then
    echo "Copying files to target directories..."
    # Copy assembly summary file to added directory
    cp "$OUTPUT_FILE" "$TARGET_ADDED_DIR/" || error_exit "Failed to copy assembly summary file to '$TARGET_ADDED_DIR'."

    # Copy and process fna files to added/temp directory
    for file in "${processed_files[@]}"; do
        filename=$(basename "$file")
        # Remove .fna or .fna.gz extension while keeping the rest of the name
        base_name="${filename%.fna*}"
        target_file="$TARGET_TEMP_DIR/${base_name}_genomic.fna.gz"

        # If file is already gzipped, just copy and rename
        if [[ "$file" == *.gz ]]; then
            cp "$file" "$target_file" || error_exit "Failed to copy file '$file' to '$target_file'."
        else
            # If file is not gzipped, compress it
            gzip -c "$file" > "$target_file" || error_exit "Failed to compress file '$file' to '$target_file'."
        fi
    done

    echo "Files have been copied:"
    echo "  - Assembly summary file -> '$TARGET_ADDED_DIR'"
    echo "  - FASTA files (renamed to _genomic.fna.gz) -> '$TARGET_TEMP_DIR'"
else
    echo "No files were successfully processed, nothing to copy."
fi

# === Final messages ===
echo "Finished generating '$OUTPUT_FILE'."
echo "Found ${file_count} '.fna' or '.fna.gz' files, successfully processed ${processed_count} entries."

# Output final warnings based on processing results
if [ "$file_count" -gt 0 ] && [ "$processed_count" -eq 0 ]; then
     echo "Warning: Could not process any of the found FASTA files. Please check file format, permissions, and if first line contains valid FASTA header with taxid." >&2
elif [ "$file_count" -ne "$processed_count" ]; then
     echo "Warning: Some FASTA files were skipped. Please check the warnings above for details." >&2
fi

exit 0