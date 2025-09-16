#! /usr/bin/env python3
"""
Script to download FASTA sequences from UniProt for a given taxon ID using UniProt REST API
"""

import sys
import requests
import re
import pandas as pd
import gzip
import argparse
from requests.adapters import HTTPAdapter, Retry

re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)                


def download_fasta_by_taxon(taxon, output_file):
    """
    Download FASTA sequences from UniProt for a given taxon ID and save to output_file
    """
    taxon_id = get_taxon_id_from_file(taxon)
    taxon_query = f'%28taxonomy_id%3A{"+OR+taxonomy_id%3A".join(taxon_id)}%29'
    # Get plain text response, then write to gzip file
    url = f"https://rest.uniprot.org/uniprotkb/search?format=fasta&query={taxon_query}+AND+%28reviewed%3Atrue%29&size=500"
    print(f"Downloading FASTA sequences for taxon ID {taxon_id}...")

    # Bug: 2 sequences merged into 1 header, total number does not match expected total
    # progress = 0
    # with gzip.open(output_file, 'wt') as f:
    #     for batch, total in get_batch(url):
    #         lines = batch.text.splitlines()
    #         if not progress:
    #             print(lines[0], file=f)
    #         for line in lines[1:]:
    #             if line.startswith('>'):
    #                 progress += 1
    #             print(line, file=f)
    #         print(f'{progress} / {total}')
    # Fix: Make progress cumulative
    total_progress = 0
    expected_total = None
    
    with gzip.open(output_file, 'wt') as f:
        for batch_num, (batch, total) in enumerate(get_batch(url), 1):
            if expected_total is None:
                expected_total = int(total)
            
            lines = batch.text.splitlines()
            batch_sequences = 0
                      
            for line in lines:
                if line.startswith('>'):
                    batch_sequences += 1
                    total_progress += 1
                print(line, file=f)
            
            print(f'Batch {batch_num}: +{batch_sequences} sequences | Total: {total_progress:,} / {expected_total:,}')

def download_tsv_by_taxon(taxon_id, output_file):
    """
    Download FASTA sequences from UniProt for a given taxon ID and save to output_file
    """
    taxon_query = f'%28taxonomy_id%3A{"+OR+taxonomy_id%3A".join(taxon_id)}%29'
    # Get plain text response, then write to gzip file
    url = f"https://rest.uniprot.org/uniprotkb/search?&fields=accession%2Ccc_function&format=tsv&query={taxon_query}+AND+%28reviewed%3Atrue%29&size=500"
    print(f"Downloading TSV file for taxon ID {taxon_id}...")
    
    progress = 0
    with gzip.open(output_file, 'wt') as f:
        for batch, total in get_batch(url):
            lines = batch.text.splitlines()
            if not progress:
                print(lines[0], file=f)
            for line in lines[1:]:
                # if line.startswith('>'):
                #     progress += 1
                print(line, file=f)
            progress += len(lines)
            print(f'{progress} / {total}')
            
    # get number of proteins with a Function summary
    data = pd.read_csv(output_file, sep='\t', compression='gzip')
    print(data['Function [CC]'].notnull().sum())


def get_taxon_id_from_file(taxon_path):
    """
    Read taxon ID from a text file
    """
    taxon_file =  pd.read_csv(taxon_path, header=0, sep='\t', encoding='ISO-8859-1')
    selected_taxon = [str(taxon_id) for taxon_id in set(taxon_file.iloc[:,0].tolist())]
    
    return selected_taxon


def parse_args(args):
    parser = argparse.ArgumentParser(description='Download FASTA sequences from UniProt for a given taxon ID')
    parser.add_argument('--taxon', '-t', required=True, 
                        help='Path to the file with taxon IDs')
    
    parser.add_argument('--output', '-o', required=True, 
                        help='Output file path for the downloaded FASTA sequences')
    
    return parser.parse_args(args)
    
    
def main():
    args = parse_args(sys.argv[1:])
    download_fasta_by_taxon(args.taxon, args.output)
    
if __name__ == "__main__":
    main()