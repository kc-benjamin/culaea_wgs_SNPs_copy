#!/bin/bash

# Script to rename samples in a VCF file using a mapping file
# Usage: ./rename_vcf_samples.sh <input.vcf> <output.vcf> <name-conversion.txt>

if [ $# -ne 3 ]; then
    echo "Usage: $0 <input_vcf> <output_vcf> <mapping_file>"
    echo "Example: $0 input.vcf output.vcf name-conversion.txt"
    exit 1
fi

INPUT_VCF="$1"
OUTPUT_VCF="$2"
MAPPING_FILE="$3"

# Check if input files exist
if [ ! -f "$INPUT_VCF" ]; then
    echo "Error: Input VCF file '$INPUT_VCF' not found"
    exit 1
fi

if [ ! -f "$MAPPING_FILE" ]; then
    echo "Error: Mapping file '$MAPPING_FILE' not found"
    exit 1
fi

echo "Processing VCF file: $INPUT_VCF"
echo "Using mapping file: $MAPPING_FILE"
echo "Output file: $OUTPUT_VCF"

# Process the VCF file
awk -v mapfile="$MAPPING_FILE" '
BEGIN {
    # Load the mapping into an associative array
    while ((getline line < mapfile) > 0) {
        split(line, parts, "\t")
        mapping[parts[1]] = parts[2]
    }
    close(mapfile)
}
{
    if ($0 ~ /^#CHROM/) {
        # This is the header line with sample names
        # Fields 1-9 are: #CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT
        # Samples start from field 10
        for (i = 1; i <= 9; i++) {
            printf "%s", $i
            if (i < 9) printf "\t"
        }
        
        # Process sample columns
        for (i = 10; i <= NF; i++) {
            printf "\t"
            if ($i in mapping) {
                printf "%s", mapping[$i]
            } else {
                printf "%s", $i
            }
        }
        printf "\n"
    } else {
        # For all other lines (including comments and data), print as-is
        print $0
    }
}
' "$INPUT_VCF" > "$OUTPUT_VCF"

# Check if conversion was successful
if [ -f "$OUTPUT_VCF" ]; then
    echo "SUCCESS: VCF file renamed and saved to '$OUTPUT_VCF'"
    echo "Sample count: $(grep -c '^#CHROM' "$OUTPUT_VCF") (should be 1)"
else
    echo "Error: Failed to create output file"
    exit 1
fi
