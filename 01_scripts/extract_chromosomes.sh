#!/bin/bash

# Script to extract chromosome names from a .dict file
# Removes the "SN:" prefix

# Path to the .dict file (adjust if necessary)
DICT_FILE="../03_genome/GCF_949316345.1_Punpun_genome.dict"

# Output file
OUTPUT_FILE="../02_info_files/chromosome_list.txt"

# Check if dict file exists
if [ ! -f "$DICT_FILE" ]; then
    echo "Error: $DICT_FILE not found. Please ensure the .dict file is available."
    exit 1
fi

# Extract chromosome names
grep "^@SQ" "$DICT_FILE" | sed 's/.*SN:\([^[:space:]]*\).*/\1/' > "$OUTPUT_FILE"

echo "Chromosome list extracted to $OUTPUT_FILE"