#!/bin/bash
# Usage: ./append_taxonomy_names.sh.sh <kreport2> <txt_input> <output>

KREPORT=$1
INPUT=$2
OUTPUT=$3

awk 'BEGIN{OFS="\t"}
NR==FNR{
  # taxid is column 5; name starts from column 6 (may include leading spaces / indentation)
  taxid = $5
  name  = $6
  for(i=7;i<=NF;i++) name = name " " $i
  sub(/^[[:space:]]+/, "", name)

  if (taxid ~ /^[0-9]+$/ && name != "") names[taxid] = name
  next
}
{
  taxid = $3
  print $0, (taxid in names ? names[taxid] : "NA")
}' "$KREPORT" "$INPUT" > "$OUTPUT"