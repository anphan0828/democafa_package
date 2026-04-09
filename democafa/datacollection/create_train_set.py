#!/usr/bin/env python3

"""
Script to create training set (sequences, taxonomy) from SwissProt FASTA file and train_terms.tsv.
"""

import time
import os
import sys
import requests
import argparse
import re
import gzip
import logging
from datetime import datetime
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Create a specific logger for this module (not the root logger)
logger = logging.getLogger('create_train_set')
logger.setLevel(logging.INFO)

# Prevent messages from propagating to the root logger (so multiple loggers can coexist)
logger.propagate = False

# Create file handler
log_dir = 'logs'
os.makedirs(log_dir, exist_ok=True)  # Create logs directory if it doesn't exist
log_filename = os.path.join(log_dir, datetime.now().strftime('create_train_set_%Y%m%d_%H%M%S.log'))
file_handler = logging.FileHandler(log_filename)
file_handler.setLevel(logging.INFO)
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)


def get_proteins_with_all_aspects(terms_file):
    """
    Get proteins that have all three GO aspects (molecular_function, 
    biological_process, cellular_component).
    
    Args:
        df: DataFrame containing GO terms data
        
    Returns:
        Set of protein IDs with all three aspects
    """
    df = pd.read_csv(terms_file, sep='\t')
    # Group by protein ID and get unique aspects for each protein
    protein_aspects = df.groupby('EntryID')['aspect'].unique()
    all_proteins = set(df['EntryID'].unique())
    
    # Find proteins with all three aspects
    complete_proteins = {
        protein for protein, aspects in protein_aspects.items() if len(aspects) == 3
    }
    
    return all_proteins, complete_proteins


def create_train_sequences(proteins_with_terms, sequences_gzfile, train_out_fasta, train_out_taxonomy):
    """
    Process UniProt FASTA file and extract sequences and taxonomy information. 
    This file only contains sequences that are in train_terms.tsv (proteins labeled with GO terms).
    """
    
    # Compile regex pattern for taxonomy ID extraction
    tax_pattern = re.compile(r"OX=(\d+)")
    
    # Counter for progress tracking
    seq_count = 0
    all_taxid = set()
    
    # Process gzipped file and write both outputs simultaneously
    with gzip.open(sequences_gzfile, "rt") as gz_file, \
            open(train_out_fasta, "w") as fasta_out, \
            open(train_out_taxonomy, "w") as mapping_out:
        
        # Process sequences in SwissProt
        for record in SeqIO.parse(gz_file, "fasta"):
            # Extract accession (EntryID)
            entry_id = record.id.split("|")[1] if "|" in record.id else record.id
            if entry_id not in proteins_with_terms:
                continue
            seq_count += 1
            
            # Extract taxonomy ID using regex
            tax_match = tax_pattern.search(record.description)
            tax_id = tax_match.group(1) if tax_match else "N/A"
            if tax_id != "N/A":
                all_taxid.add(tax_id)
            mapping_out.write(f"{entry_id}\t{tax_id}\n")
            
            SeqIO.write(record, fasta_out, "fasta")
            
            # Print progress every 10 sequences
            if seq_count % 10000 == 0:
                logger.debug(f"Processed {seq_count} sequences in SwissProt...")

    logger.info(f"Total SwissProt sequences in training data: {seq_count} proteins in {len(all_taxid)} taxa")

def parse_inputs(args):
    parser = argparse.ArgumentParser(
        description='Create training set from train_terms and all sequences.'
    )
    parser.add_argument('--terms', '-t', required=True, 
                        help='Tab-separated file with UniProtKB accessions and GO terms and GO aspects with header')
    parser.add_argument('--fasta_gz', '-f', required=True,
                        help='Path to gzipped SwissProt FASTA file')
    parser.add_argument('--train_out_fasta', '-tf', required=True,
                        help='Path to training FASTA file with annotated proteins')
    parser.add_argument('--train_out_taxonomy', '-tt', required=True,
                        help='Path to training taxonomy mapping file')
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], 
                        default='INFO', help='Set the logging level (default: INFO)')
 
    args = parser.parse_args(args)
    
    # Configure logging level based on argument
    logger.setLevel(getattr(logging, args.log_level))
    for handler in logger.handlers:
        handler.setLevel(getattr(logging, args.log_level))

    return args

def main():
    args = parse_inputs(sys.argv[1:])
    
    logger.info(f"Arguments: {vars(args)}")
    
    # Read GO terms data
    logger.info("Reading GO terms data...")
    all_proteins, complete_proteins = get_proteins_with_all_aspects(args.terms)
    
    # Create training sequences and taxonomy files
    create_train_sequences(
        proteins_with_terms=all_proteins,
        sequences_gzfile=args.fasta_gz,
        train_out_fasta=args.train_out_fasta,
        train_out_taxonomy=args.train_out_taxonomy
    )

if __name__ == "__main__":
    main()
    