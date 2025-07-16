#!/bin/bash

# Usage:
# ./create_toy_dataset.sh <input_gaf.gz> <output_gaf.gz> [fraction]

set -e  # Exit on any error

# Check arguments
if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 <input_gaf.gz> <output_gaf.gz> [fraction]"
    echo "  input_gaf.gz: Large gzipped GAF file to sample from"
    echo "  output_gaf.gz: Output toy dataset file"
    echo "  fraction: Fraction of data lines to sample (default: 0.1 = 10%)"
    echo ""
    echo "This script creates a toy dataset by:"
    echo "  1. Taking a fraction of the input file (default 10%)"
    echo "  2. Removing all IEA annotations (column 7)"
    echo "  3. Preserving all header lines"
    exit 1
fi

INPUT_GAF=$1
OUTPUT_GAF=$2
FRACTION=${3:-0.1}  # Default to 10% if not specified
TEMP_DIR=$(dirname "$OUTPUT_GAF")/temp_dir

echo "Creating toy dataset from $INPUT_GAF"
echo "Sampling fraction: $FRACTION"
echo "Output file: $OUTPUT_GAF"
mkdir -p "$TEMP_DIR"
echo "Temporary directory: $TEMP_DIR"

# Validate fraction is between 0 and 1
if ! awk "BEGIN {exit !($FRACTION > 0 && $FRACTION <= 1)}"; then
    echo "Error: Fraction must be between 0 and 1 (got: $FRACTION)"
    exit 1
fi

# Create AWK script for efficient sampling and filtering
cat > "$TEMP_DIR/sample_filter.awk" << EOF
BEGIN {
    FS = "\t"
    srand()  # Initialize random seed
    fraction = $FRACTION
    header_lines = 0
    data_lines = 0
    sampled_lines = 0
    filtered_lines = 0
}

# Print all header lines (starting with !)
/^!/ {
    print
    header_lines++
    next
}

# For data lines, sample and filter
{
    data_lines++
    
    # Sample based on fraction
    if (rand() <= fraction) {
        sampled_lines++
        
        # Filter out IEA annotations (column 7)
        if (\$7 != "IEA") {
            print
            filtered_lines++
        }
    }
}

END {
    printf "Header lines: %d\n", header_lines > "/dev/stderr"
    printf "Total data lines processed: %d\n", data_lines > "/dev/stderr"
    printf "Lines sampled (%.1f%%): %d\n", fraction * 100, sampled_lines > "/dev/stderr"
    printf "Lines after removing IEA: %d\n", filtered_lines > "/dev/stderr"
    printf "Final reduction: %.2f%% of original\n", (filtered_lines / data_lines) * 100 > "/dev/stderr"
}
EOF

# Process the file with AWK
echo "Processing GAF file..."
echo "Step 1: Sampling $FRACTION of data lines"
echo "Step 2: Removing IEA annotations from column 7"

zcat "$INPUT_GAF" | awk -f "$TEMP_DIR/sample_filter.awk" | gzip > "$OUTPUT_GAF"

# Get file sizes for comparison
INPUT_SIZE=$(du -h "$INPUT_GAF" | cut -f1)
OUTPUT_SIZE=$(du -h "$OUTPUT_GAF" | cut -f1)

echo "Cleaning up temporary files..."
rm -rf "$TEMP_DIR"

echo ""
echo "Done! Toy dataset created successfully."
echo "Input file size:  $INPUT_SIZE"
echo "Output file size: $OUTPUT_SIZE"
echo "Toy dataset saved to: $OUTPUT_GAF"
