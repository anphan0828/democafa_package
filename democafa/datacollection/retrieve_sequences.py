#!/usr/bin/env python3

"""
Script to process UniProt FASTA file and extract selected sequences from a .txt file

This script:
1. Gunzips the UniProt FASTA file
2. Reads entries from a .txt file and write these sequences to a new FASTA file
"""

import gzip
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import re
import sys
import argparse

# def process_uniprot_fasta(input_fasta, input_terms, output_taxonomy, output_fasta, seq_limit=None):
#     """
#     Process UniProt FASTA file and extract sequences and taxonomy information. 
#     This file only contains sequences that are in train_terms.tsv (proteins labeled with GO terms).
#     """
#     print(f"Processing {input_fasta}...")
#     df = pd.read_csv(input_terms, sep='\t', header=0, usecols=[0], names=['EntryID'])
#     proteins_with_terms = set(df['EntryID'])
    
#     # Compile regex pattern for taxonomy ID extraction
#     tax_pattern = re.compile(r"OX=(\d+)")
    
#     # Counter for progress tracking
#     seq_count = 0
    
#     # Process gzipped file and write both outputs simultaneously
#     with gzip.open(input_fasta, "rt") as gz_file, \
#             open(output_fasta, "w") as fasta_out, \
#             open(output_taxonomy, "w") as mapping_out:
        
#         # Process each sequence record
#         for record in SeqIO.parse(gz_file, "fasta"):
#             seq_count += 1
            
#             # Extract accession (EntryID)
#             entry_id = record.id.split("|")[1] if "|" in record.id else record.id
#             if entry_id not in proteins_with_terms:
#                 continue
#             # Extract taxonomy ID using regex
#             tax_match = tax_pattern.search(record.description)
#             tax_id = tax_match.group(1) if tax_match else "N/A"
            
#             # Write to mapping file
#             mapping_out.write(f"{entry_id}\t{tax_id}\n")
            
#             # Write to FASTA file
#             SeqIO.write(record, fasta_out, "fasta")
            
#             # Print progress every 10 sequences
#             if seq_count % 10000 == 0:
#                 print(f"Processed {seq_count} sequences...")
            
#             # Stop after reaching sequence limit
#             if seq_limit is not None and seq_limit > 0:
#                 if seq_count >= seq_limit:
#                     break
        
#     print(f"Total sequences processed: {seq_count}")
#     print(f"Output files created:")
#     print(f"- FASTA file: {output_fasta}")
#     print(f"- Mapping file: {output_taxonomy}")


def retrieve_sequences(fasta, input, out_fasta):
    """
    Retrieve sequences from a FASTA file based on selected entries in a .txt file.
    
    Args:
        fasta: Path to the input FASTA file (can be gzipped).
        input: Path to a text file containing protein IDs (one per line).
        out_fasta: Path to the output FASTA file with selected proteins.
    """
    print(f"Retrieving sequences from {fasta} based on entries in {input}...")
    
    # Read protein IDs from input file
    with open(input, 'r') as f:
        selected_ids = set(line.strip() for line in f if line.strip())
    
    # Open output FASTA file
    with open(out_fasta, 'w') as out_handle:
        # Process the FASTA file
        with gzip.open(fasta, 'rt') if fasta.endswith('.gz') else open(fasta, 'r') as fasta_handle:
            for record in SeqIO.parse(fasta_handle, 'fasta'):
                entry_id = record.id.split("|")[1] if "|" in record.id else record.id
                if entry_id in selected_ids:
                    SeqIO.write(record, out_handle, 'fasta')
    
    print(f"Sequences written to {out_fasta}")
    
    
def parse_inputs(argv):
    parser = argparse.ArgumentParser(
        description='Retrieve sequences from UniProt FASTA file based on selected entries in a .txt file')

    parser.add_argument('--fasta', '-f', required=True,
                        help='Path to fasta file (can be gzipped)')
    parser.add_argument('--input', '-i', required=True, 
                        help='Path to a text file containing protein IDs (one per line)')
    # parser.add_argument('--out_taxonomy', '-ot', required=True,
    #                     help='File name for output taxonomy mapping')
    parser.add_argument('--out_fasta', '-of', required=True,
                        help='File name for output FASTA file (with selected proteins)')
    # parser.add_argument('--seq_limit', '-s', required=False,
    #                     help='Number of sequences to process (for testing purposes)')
    return parser.parse_args(argv)

def main():
    args = parse_inputs(sys.argv[1:])
    retrieve_sequences(
        fasta = args.fasta,
        input = args.input,
        out_fasta = args.out_fasta
    )
    
if __name__ == "__main__":
    main()