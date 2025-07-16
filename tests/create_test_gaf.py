#!/usr/bin/env python3

import os
import gzip
import random
import sys
import argparse

def create_sample_gaf(source_file, output_file, sample_size=10000, seed=42):
    """Create a representative sample of a GAF file for testing"""
    random.seed(seed)
    
    # Determine file types
    is_gzipped = source_file.endswith('.gz')
    open_source = gzip.open if is_gzipped else open
    source_mode = 'rt' if is_gzipped else 'r'
    
    is_output_gzipped = output_file.endswith('.gz')
    open_output = gzip.open if is_output_gzipped else open
    output_mode = 'wt' if is_output_gzipped else 'w'
    
    # Process the file
    header_lines = []
    data_lines = []
    
    with open_source(source_file, source_mode) as f:
        for line in f:
            if line.startswith('!'):
                header_lines.append(line)
            else:
                data_lines.append(line)
                # If we've collected enough data lines, randomly replace some
                # to maintain a manageable in-memory footprint
                if len(data_lines) > sample_size * 1000:
                    break
                #     if random.random() < 0.5:  # 50% chance to replace an existing entry
                #         idx = random.randrange(len(data_lines))
                #         data_lines[idx] = line
    # TODO: sample file all IEA or protein binding lines
    
    # Sample if we have more lines than needed
    if len(data_lines) > sample_size:
        data_lines = random.sample(data_lines, sample_size)
    
    # Write the sample file
    with open_output(output_file, output_mode) as out:
        for header in header_lines:
            out.write(header)
        
        for line in data_lines:
            out.write(line)
    
    print(f"Created sample file with {len(header_lines)} header lines and {len(data_lines)} data records")
    return len(header_lines) + len(data_lines)

def main():
    parser = argparse.ArgumentParser(description='Create test data samples')
    parser.add_argument('--source', '-s', required=True, help='Source GAF file')
    parser.add_argument('--output', '-o', required=True, help='Output sample file')
    parser.add_argument('--size', '-n', type=int, default=10000, help='Number of records to sample')
    
    args = parser.parse_args()
    create_sample_gaf(args.source, args.output, args.size)

if __name__ == "__main__":
    main()