#!/usr/bin/env python3
"""
UniProt non-experimental annotations baseline
Makes predictions based on the non-experimental GO terms that the protein has

"""

# TODO: find a way to import these modules from the package
from democafa.datacollection.retrieve_terms import wrapper_retrieve_terms
from democafa.config import GO_CODES, RAW_FILE_PATHS
from democafa.utils.dask_write import write_dask_dataframe_to_gzipped_tsv
import sys
import os
import argparse
import pandas as pd
import dask.dataframe as dd

def create_predictions(terms_file, query_ids, output_baseline):
    terms_df = pd.read_csv(terms_file, sep='\t', header=0, names=['EntryID', 'term', 'aspect'])
    terms_df = terms_df[terms_df['EntryID'].isin(query_ids)]
    del terms_df['aspect']
    terms_df['value'] = [1] * len(terms_df)
    
    # Use dask to write to gzipped TSV
    # terms_df.to_csv(output_baseline, sep='\t', index=False, header=False)
    dask_df = dd.from_pandas(terms_df, npartitions=16)
    write_dask_dataframe_to_gzipped_tsv(dask_df, output_baseline) # this function is defined in utils/write_dask.py
    print(f"Predictions for {len(set(terms_df['EntryID']))} proteins written to {output_baseline}")
    
def goa_nonexp_predict(annot_file, selected_go, query_ids, output_baseline):
    config_go_codes = GO_CODES

    wrapper_retrieve_terms(
        annot_file=annot_file,
        filetype='goa',
        go_codes=config_go_codes,
        selected_go_codes=selected_go,
        graph=RAW_FILE_PATHS['obo'],
        output_tsv='data/processed/nonexp_terms.tsv' # just a temporary file
    )
    
    create_predictions('data/processed/nonexp_terms.tsv', query_ids, output_baseline)
    os.remove('data/processed/nonexp_terms.tsv')
    
def parse_args(argv):
    parser = argparse.ArgumentParser(description='UniProt non-experimental annotations baseline')
    parser.add_argument('--annot_file',  help='Path to the UniProt GOA file', required=True)
    parser.add_argument('--selected_go',  help='Selected GO codes', required=True)
    parser.add_argument('--query_ids',  help='List of query sequence IDs', required=True)
    parser.add_argument('--output_baseline',  help='Path to the output file', required=True)
    return parser.parse_args(argv)    

def main():
    args = parse_args(sys.argv[1:])
    goa_nonexp_predict(
        annot_file=args.annot_file,
        selected_go=args.selected_go,
        query_ids=args.query_ids,
        output_baseline=args.output_baseline
    )
    
if __name__ == '__main__':
    main()
