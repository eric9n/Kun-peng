#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Run in a directory containing output_*.txt and output_*.kreport2
# Output: output_<N>_name.txt for each sample

for txt in output_*.txt; do
  base="${txt%.txt}"               # e.g., output_1
  kreport2="${base}.kreport2"      # e.g., output_1.kreport2
  out="${base}_name.txt"           # e.g., output_1_name.txt

  [[ -f "$kreport2" ]] || { echo "[WARN] skip $txt (missing $kreport2)"; continue; }

  awk '
  BEGIN{OFS="\t"}
  NR==FNR{
    taxid=$5
    name=$6
    for(i=7;i<=NF;i++) name=name " " $i
    sub(/^[[:space:]]+/, "", name)
    if (taxid ~ /^[0-9]+$/ && name != "") names[taxid]=name
    next
  }
  {
    taxid=$3
    print $0, (taxid in names ? names[taxid] : "NA")
  }' "$kreport2" "$txt" > "$out"

  echo "[OK] $out"
done