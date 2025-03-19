#!/usr/bin/env python3

"""
Script to create test superset from train_terms and train_sequences.

This script:
1. Extracts proteins missing GO aspects from train_terms and gets their sequences from train_sequences
2. Collects TrEMBL proteins annotated with GO terms but not in train_sequences and adds them to the test set by using the UniProt API
3. Creates a test superset FASTA file
4. Creates a test superset taxonomy file
"""

import sys
import requests
import argparse
import re
import gzip
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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


def read_fasta_proteins(fasta_file: str):
    """
    Read proteins from FASTA file using Biopython.
    Filter for taxonomy (for demo purposes, only human proteins are considered).
    
    Args:
        fasta_file: Path to the FASTA file
        
    Returns:
        Dictionary mapping protein IDs to their sequences
    """
    tax_pattern = re.compile(r"OX=(\d+)")
    species_pattern = re.compile(r"OS=([^O]+)")
    fasta_proteins = {}
    species = {}
    
    seq_count = 0
    
    # Process gzipped fasta file to get all SwissProt proteins
    with gzip.open(fasta_file, "rt") as gz_file:
        for record in SeqIO.parse(gz_file, 'fasta'):
            # Extract accession (EntryID)
            entry_id = record.id.split("|")[1] if "|" in record.id else record.id
            
            # Extract taxonomy ID using regex
            tax_match = tax_pattern.search(record.description)
            tax_id = tax_match.group(1) if tax_match else "N/A"
                
            species_match = species_pattern.search(record.description)
            species_name = species_match.group(1) if species_match else "N/A"
            if tax_id == '9606': # TODO: human only for demo
                seq_count += 1
                fasta_proteins[entry_id] = [str(record.id), str(record.description), str(record.seq), tax_id]
            if species_name != "N/A":
                species[species_name] = tax_id
    print(f"Total sequences processed: {seq_count}")
            
    return fasta_proteins, species


def write_missing_aspects_fasta(fasta_proteins, complete_proteins, trembl_proteins, uniprot_api_version, output_file):
    """
    Write FASTA file containing proteins missing GO aspects.
    
    Args:
        fasta_proteins: Dictionary of protein sequences
        complete_proteins: Set of proteins with all aspects
        trembl_proteins: Set of TrEMBL proteins with missing aspects
        output_file: Path to output FASTA file
    """
    # Get proteins missing aspects
    missing_aspects = set(fasta_proteins.keys()) - complete_proteins
         
    with open(output_file, 'w') as fasta_out:
        for protein_id in missing_aspects:
            seq_record = SeqRecord(Seq(fasta_proteins[protein_id][2]), id=fasta_proteins[protein_id][0], description=fasta_proteins[protein_id][1])
            # values = fasta_proteins[protein_id]
            # f.write(f'>{protein_id}\t{values[1]}\n{values[0]}\n')
            SeqIO.write(seq_record, fasta_out, "fasta")
            
    for protein_id in trembl_proteins:
        header, sequence = get_fasta_from_API(protein_id, uniprot_api_version)
        
        if header and sequence:
            seq_record = SeqRecord(Seq(sequence), id=header, description="")
            with open(output_file, "a") as handle:
                SeqIO.write(seq_record, handle, "fasta")
        else:
            print(f"Release version {uniprot_api_version} not found.")


def write_fasta(fasta_proteins, trembl_proteins, uniprot_api_version, output_file):
    """
    Write FASTA file all proteins with or without missing GO aspects (for partial knowledge evaluation).
    
    Args:
        fasta_proteins: Dictionary of protein sequences
        trembl_proteins: Set of TrEMBL proteins with missing aspects
        uniprot_api_version: UniProt release version to retrieve TrEMBL sequences
        output_file: Path to output FASTA file
    """
    with open(output_file, 'w') as fasta_out:
        for protein_id in fasta_proteins:
            seq_record = SeqRecord(Seq(fasta_proteins[protein_id][2]), id=fasta_proteins[protein_id][0], description=fasta_proteins[protein_id][1])
            # values = fasta_proteins[protein_id]
            # f.write(f'>{protein_id}\t{values[1]}\n{values[0]}\n')
            SeqIO.write(seq_record, fasta_out, "fasta")
    
    print(f"Appending {len(trembl_proteins)} TrEMBL proteins to {output_file}...")
    for protein_id in trembl_proteins:
        header, sequence = get_fasta_from_API(protein_id, uniprot_api_version)
        
        if header and sequence:
            seq_record = SeqRecord(Seq(sequence), id=header, description="")
            with open(output_file, "a") as handle:
                SeqIO.write(seq_record, handle, "fasta")
    
           
def get_fasta_from_API(accession, uniprot_api_version):
    url = f"https://rest.uniprot.org/unisave/{accession}?format=fasta"
    
    response = requests.get(url)
    response.raise_for_status()
    fasta_content = response.text
    
    header = fasta_content.splitlines()[0][1:]
    sequence = ""
    found_release = False
    
    for line in fasta_content.splitlines():
        if line.startswith(">"):
            if re.search(uniprot_api_version, line):
                found_release = True
            elif found_release:
                break
        elif found_release:
            sequence += line
    
    if sequence:
        return header, sequence
    else:
        return None, None
            

def write_species(species, output_file):
    """
    Write text file containing species names and taxonomy IDs.
    
    Args:
        species: Dictionary of species names and taxonomy IDs
        output_file: Path to output text file
    """
    with open(output_file, 'w') as f:
        f.write('ID\tSpecies\n')
        for name, tax_id in species.items():
            f.write(f'{tax_id}\t{name}\n')
            
            
def create_test_set(terms_file, sequences_gzfile, out_fasta, uniprot_api_version, out_taxonomy, include_all=False):
    # Read GO terms data
    print("Reading GO terms data...")
    all_proteins, complete_proteins = get_proteins_with_all_aspects(terms_file)
 
    # Read FASTA file
    print("Reading FASTA sequences...")
    fasta_proteins, species = read_fasta_proteins(sequences_gzfile)
    
    # Get TrEMBL proteins with missing aspects
    trembl_proteins = all_proteins - set(fasta_proteins.keys())
    print(f"Found {len(trembl_proteins)} TrEMBL proteins with missing aspects.")
    
    # Write output files
    if include_all: # including partial knowledge proteins, only append TrEMBL proteins to fasta file
        write_fasta(fasta_proteins, trembl_proteins, uniprot_api_version, out_fasta)
    else: # only include proteins missing aspects (no knowledge and limited knowledge)
        write_missing_aspects_fasta(fasta_proteins, complete_proteins, trembl_proteins, uniprot_api_version, out_fasta)
    write_species(species, out_taxonomy)
    

def parse_inputs(args):
    parser = argparse.ArgumentParser(
        description='Create test set from train_terms and train_sequences.'
    )
    parser.add_argument('--terms', '-t', required=True, 
                        help='Tab-separated file with UniProtKB accessions and GO terms and GO aspects with header')
    parser.add_argument('--fasta_gz', '-f', required=True,
                        help='Path to gzipped SwissProt FASTA file')
    parser.add_argument('--out_fasta', '-o', required=True,
                        help='Path to test superset FASTA file')
    parser.add_argument('--uniprot_api', '-u', required=True,
                        help='UniProt API version to retrieve TrEMBL sequences')
    parser.add_argument('--out_taxonomy', required=True)
    parser.add_argument('--include_all', action='store_true', default=False)
    return parser.parse_args(args)


def main():
    args = parse_inputs(sys.argv[1:])
    create_test_set(
        terms_file=args.terms,
        sequences_gzfile=args.fasta_gz,
        out_fasta=args.out_fasta,
        uniprot_api_version=args.uniprot_api,
        out_taxonomy=args.out_taxonomy,
        include_all=True
    )

if __name__ == "__main__":
    main()
    