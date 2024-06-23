#!/bin/bash

directory=$1

# Find all hash_*.k2d files and calculate their total size
total_size=$(find "$directory" -name "hash_*.k2d" -exec du -ch {} + | grep total$ | awk '{print $1}')
echo "Total size of hash_*.k2d files: $total_size"
