#!/usr/bin/env python3

import os
from democafa.config import VERSIONS, GO_CODES, RAW_FILE_PATHS, PROCESSED_PATHS
from datacollection.retrieve_terms import wrapper_retrieve_terms
from datacollection.retrieve_sequences import process_uniprot_fasta
from datacollection.create_test_set import create_test_set
from utils.ontology import propagate_and_ia

from predictors.blast import blast_predict
from predictors.naive import naive_predict
from predictors.goa_nonexp import goa_nonexp_predict

def main():
    # 0. Load configuration
    config_go_codes = GO_CODES
    
    # 1. Collecting data for release
    if not os.path.exists('data/processed/train_terms.tsv'):
        wrapper_retrieve_terms(
            annot_file=RAW_FILE_PATHS['uniprot_goa'],
            filetype='goa',
            go_codes=config_go_codes,
            selected_go_codes='Experimental,IC,TAS',
            graph=RAW_FILE_PATHS['obo'],
            output_tsv=PROCESSED_PATHS['train_terms']
        )
    else:
        print("train_terms.tsv already exists, moving to sequences retrieval")
    
    if not os.path.exists('data/processed/train_sequences.fasta'): 
        process_uniprot_fasta(
            input_fasta = RAW_FILE_PATHS['uniprot_sprot'],
            input_terms = PROCESSED_PATHS['train_terms'], # has to match output_tsv from wrapper_retrieve_terms
            output_taxonomy = 'data/processed/train_taxonomy.tsv',
            output_fasta = PROCESSED_PATHS['train_sequences'], # only annotated proteins 
            seq_limit=None
        )
    else:
        print("train_sequences.fasta already exists, moving to test set creation")
    
    if not os.path.exists('data/processed/test_superset.fasta'):
        create_test_set(
            terms_file=PROCESSED_PATHS['train_terms'],
            sequences_gzfile=RAW_FILE_PATHS['uniprot_sprot'], # full SwissProt fasta file 
            out_fasta='data/processed/test_superset.fasta',
            uniprot_api_version=VERSIONS["UniProt_version"],
            out_taxonomy='data/processed/testsuperset-taxon-list.tsv',
            include_all=True
        )
    else:    
        print("test_superset.fasta already exists, all data collected")
        
    propagate_and_ia(
        terms_file=PROCESSED_PATHS['train_terms'],
        graph=RAW_FILE_PATHS['obo'],
        matrix_propagated='data/processed/train_terms_propagated.npz', 
        output_tsv='data/processed/ia.tsv' # change to None if don't need IA calculation
    )
    
    # Use os module to copy files from processed to release folder
    os.system('cp data/processed/* data/release/')
    os.system(f'cp {RAW_FILE_PATHS["obo"]} data/release/')
    print(f"Data for release ready with {len(os.listdir('data/release'))} files")
    
    # 2. Create baseline predictors
    goa_nonexp_predict(
        annot_file=RAW_FILE_PATHS['uniprot_goa'],
        selected_go='Computational,Phylogenetical,Electronic,ND,NAS',
        query_ids='data/release/test_superset.fasta',
        output_baseline='predictors/goa_nonexp_baseline.tsv'
    )
    
if __name__ == "__main__":
    main()
