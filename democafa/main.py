#!/usr/bin/env python3

import os
from config import VERSIONS, GO_CODES, RAW_FILE_PATHS, PROCESSED_PATHS, RELEASE_PATHS, EXTERNAL_TOOLS
from datacollection.retrieve_terms import wrapper_retrieve_terms
from datacollection.retrieve_sequences import process_uniprot_fasta
from datacollection.create_test_set import create_test_set
from utils.ontology import propagate_and_ia

from baselines.blast import blast_predict
from baselines.naive import naive_predict
from baselines.goa_nonexp import goa_nonexp_predict


def main():
    # 0. Load configuration
    config_go_codes = GO_CODES
    DATA_DIR = os.environ.get('DATA_DIR', 'data')
    # 1. Collecting data for release
    if not os.path.exists(PROCESSED_PATHS['train_terms']):
        wrapper_retrieve_terms(
            annot_file=RAW_FILE_PATHS['uniprot_goa'],
            filetype='goa',
            go_codes=config_go_codes,
            selected_go_codes='Experimental,IC,TAS',
            graph=RAW_FILE_PATHS['obo'],
            output_tsv=PROCESSED_PATHS['train_terms']
        )
    else:
        print(f"{PROCESSED_PATHS['train_terms']} already exists, moving to sequences retrieval")
    
    if not os.path.exists(PROCESSED_PATHS['train_sequences']): 
        process_uniprot_fasta(
            input_fasta = RAW_FILE_PATHS['uniprot_sprot'],
            input_terms = PROCESSED_PATHS['train_terms'], # has to match output_tsv from wrapper_retrieve_terms
            output_taxonomy = PROCESSED_PATHS['train_taxonomy'],
            output_fasta = PROCESSED_PATHS['train_sequences'], # only annotated proteins 
            seq_limit=None
        )
    else:
        print(f"{PROCESSED_PATHS['train_sequences']} already exists, moving to test set creation")
    
    if not os.path.exists(PROCESSED_PATHS['test_sequences_all']):
        create_test_set(
            terms_file=PROCESSED_PATHS['train_terms'],
            sequences_gzfile=RAW_FILE_PATHS['uniprot_sprot'], # full SwissProt fasta file 
            out_fasta=PROCESSED_PATHS['test_sequences'],
            uniprot_api_version=VERSIONS["UniProt_version"],
            out_taxonomy=PROCESSED_PATHS['test_taxonomy'],
            include_all=True # include all proteins in test superset (not just proteins with missing aspects)
        )
    else:    
        print(f"{PROCESSED_PATHS['test_sequences_all']} already exists, all data collected")
    if not os.path.exists(PROCESSED_PATHS['train_matrix']):    
        propagate_and_ia(
            terms_file=PROCESSED_PATHS['train_terms'],
            graph=RAW_FILE_PATHS['obo'],
            matrix_propagated=PROCESSED_PATHS['train_matrix'],
            matrix_indices=PROCESSED_PATHS['matrix_indices'], 
            output_tsv=None # change to None if don't need IA calculation
        )
    
    # Use os module to copy files from processed to release folder
    # os.system('cp data/processed/* data/release/')
    # os.system(f'cp {RAW_FILE_PATHS["obo"]} data/release/')
    # print(f"Data for release ready with {len(os.listdir('data/release'))} files")
    
    # 2. Create baseline predictors
    if not os.path.exists(PROCESSED_PATHS['baseline_goa_nonexp']):
        goa_nonexp_predict(
            annot_file=RAW_FILE_PATHS['uniprot_goa'],
            selected_go='Computational,Phylogenetical,Electronic,ND,NAS',
            query_file=RELEASE_PATHS['test_sequences'],
            output_baseline=PROCESSED_PATHS['baseline_goa_nonexp']
        )
    
    if not os.path.exists(PROCESSED_PATHS['baseline_naive']):
        naive_predict(
            annotations=PROCESSED_PATHS['train_matrix'],
            query_file=RELEASE_PATHS['test_sequences'],
            indices=PROCESSED_PATHS['matrix_indices'],
            output_baseline=PROCESSED_PATHS['baseline_naive'] 
        )
    
    # os.system(f'mkdir -p {DATA_DIR}/processed/blast_db && cp {PROCESSED_PATHS["train_sequences"]} {DATA_DIR}/processed/blast_db/')
    # os.system(f'sbatch {EXTERNAL_TOOLS["blast"]}') # TODO: consider a separate shell script
    if not os.path.exists(PROCESSED_PATHS['baseline_blast']):
        blast_predict(
            annotations=PROCESSED_PATHS['train_matrix'],
            query_file=RELEASE_PATHS['test_sequences'],
            indices=PROCESSED_PATHS['matrix_indices'],
            blast_results=PROCESSED_PATHS['output_blast'],
            output_baseline=PROCESSED_PATHS['baseline_blast'],
            keep_self_hits=False,
            use_rscore=False
        )
    
    os.system(f'mkdir -p {DATA_DIR}/processed/prott5')
    os.system(f'sbatch {EXTERNAL_TOOLS["prott5"]} -q {RELEASE_PATHS["test_sequences"]} -d {PROCESSED_PATHS["train_sequences"]} -m {VERSIONS["prott5_model"]} -o {PROCESSED_PATHS["output_prott5"]}')
    
    # 3. Evaluation
if __name__ == "__main__":
    main()
