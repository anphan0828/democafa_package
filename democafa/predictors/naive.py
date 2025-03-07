#!/usr/bin/env python3
"""
Naive predictor implementation for CAFA2 evaluation.
Makes predictions based on term frequency in training data.

This is a direct port of the MATLAB pfp_naive.m function.
"""

import numpy as np
from scipy import sparse
from typing import Dict
import pandas as pd
import pickle as cp
import dask.dataframe as dd
from utils.dask_write import write_dask_dataframe_to_gzipped_tsv
from Bio import SeqIO   

def naive_predict(query_ids: Dict, annotations: sparse.csr_matrix, terms: Dict) -> pd.DataFrame:
    """
    Make predictions based on term frequencies in training annotations.
    
    Args:
        query_ids: List of query sequence IDs
        annotations: Training annotations sparse matrix (n_sequences x n_terms)
        
    Returns:
        Dictionary containing:
            - object: Query sequence IDs
            - score: Sparse matrix of prediction scores
            - date: Timestamp
            
    Note:
        Assigns the same term frequency-based scores to all query sequences.
    """
    # Calculate term frequencies from annotations
    term_frequencies = np.array(annotations.sum(axis=0)).flatten() # sum along rows (within columns), not affected by number of columns (terms)
    
    # Normalize scores by maximum frequency
    max_freq = np.max(term_frequencies)
    if max_freq > 0:
        normalized_scores = term_frequencies / max_freq
    else:
        normalized_scores = term_frequencies
        
    # Create output matrix with same scores for all queries
    n_queries = len(query_ids)
    n_terms = annotations.shape[1]
    scores = sparse.lil_matrix((n_queries, n_terms), dtype=np.float32)
    for i in range(n_queries):
        scores[i] = normalized_scores # TODO: this can be optimized further
    scores = scores.tocsr()
        
    # Turn sparse matrix into dataframe format
    scores_df = pd.DataFrame(scores.toarray(), index=query_ids)
    scores_df.columns = terms
    melted_df = scores_df.reset_index().melt(id_vars = 'index',var_name = 'term', value_name = 'value')
    melted_df = melted_df.rename(columns ={'index':'EntryID'})
    print(melted_df.head())
    # melted_df.to_csv('naive_predictions.tsv.gz', sep="\t", index=False, header=False, compression='gzip')
    # melted_df.to_parquet('naive_predictions.parquet.gz', compression='gzip', engine='pyarrow')
    # dask_df.to_csv('naive_predictions.tsv.gz', sep="\t", index=False, header=False, compression='gzip')
    
    # Create Dask DataFrame
    dask_df = dd.from_pandas(melted_df, npartitions=16)
    # Define the output file path
    output_file = 'data/processed/naive_predictions.tsv.gz'
    # Write the DataFrame to a gzipped TSV file using the Dask function
    write_dask_dataframe_to_gzipped_tsv(dask_df, output_file) # this function is defined in utils/write_dask.py
    
    
    return 

if __name__ == '__main__':
    # Load annotation matrix
    annotations = sparse.load_npz('S_matrix.npz')
    
    # Load protein and term indices
    with open('row_col_idx.pkl', 'rb') as f:
        proteins, terms = cp.load(f)
    
    # For naive baseline predictor, only use terms, don't need training proteins    
    query_ids = []
    for record in SeqIO.parse('/work/idoerg/ahphan/democafa/A_data/data/test_superset.fasta','fasta'):
        query_ids.append(record.id)
        
    # Make predictions
    naive_predict(query_ids, annotations, terms)