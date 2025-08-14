#! /usr/bin/env python3
"""
Script to download FASTA sequences from UniProt for a given taxon ID using UniProt REST API
"""

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


def download_fasta_by_taxon(taxon_id, output_file):
    """
    Download FASTA sequences from UniProt for a given taxon ID and save to output_file
    """
    taxon_query = f'%28taxonomy_id%3A{"+OR+taxonomy_id%3A".join(taxon_id)}%29'
    # Get plain text response, then write to gzip file
    url = f"https://rest.uniprot.org/uniprotkb/search?format=fasta&query={taxon_query}+AND+%28reviewed%3Atrue%29&size=500"
    print(f"Downloading FASTA sequences for taxon ID {taxon_id}...")

    progress = 0
    with gzip.open(output_file, 'wt') as f:
        for batch, total in get_batch(url):
            lines = batch.text.splitlines()
            if not progress:
                print(lines[0], file=f)
            for line in lines[1:]:
                if line.startswith('>'):
                    progress += 1
                print(line, file=f)
            print(f'{progress} / {total}')
    
    
def get_taxon_id_from_file(taxon_path):
    """
    Read taxon ID from a text file
    """
    taxon_file =  pd.read_csv(taxon_path, header=0, sep='\t', encoding='ISO-8859-1')
    selected_taxon = [str(taxon_id) for taxon_id in set(taxon_file.iloc[:,0].tolist())]
    
    return selected_taxon


def main():
    parser = argparse.ArgumentParser(description="Download FASTA sequences from UniProt for a given taxon ID")
    parser.add_argument('--taxon', type=str, required=True, help='Path to the file with taxon IDs')
    parser.add_argument('--output', type=str, required=True, help='Output file path for the downloaded FASTA sequences')
    
    args = parser.parse_args()
    
    taxon_ids = get_taxon_id_from_file(args.taxon)
    download_fasta_by_taxon(taxon_ids, args.output)
    
if __name__ == "__main__":
    main()