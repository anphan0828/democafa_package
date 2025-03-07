#!/usr/bin/env python3
"""
BLAST-based predictor implementation for CAFA2 evaluation.
Makes predictions based on BLAST sequence similarity.

This is a direct port of the MATLAB pfp_blast.m function.
"""

import numpy as np
import pandas as pd
from scipy import sparse
from datetime import datetime
import pickle as cp
from Bio import SeqIO
from typing import Dict, List, Optional
import dask.dataframe as dd
from utils.dask_write import write_dask_dataframe_to_gzipped_tsv

def calculate_rscore(evalue: float, max_rscore: float = 500.0) -> float:
    """
    Calculate R-score from BLAST E-value.
    R-score = -log(E-value) + 2, capped at max_rscore
    
    Args:
        evalue: BLAST E-value
        max_rscore: Maximum allowed R-score (default: 500)
        
    Returns:
        R-score value
    """
    if evalue == 0:
        return max_rscore
    rscore = -np.log10(evalue) + 2
    return min(rscore, max_rscore)

def blast_predict(query_ids: Dict, 
                 blast_results: pd.DataFrame,
                 annotations: sparse.csr_matrix,
                 proteins: Dict,
                 terms: Dict,
                 keep_self_hits: bool = False,
                 use_rscore: bool = False) -> pd.DataFrame:
    """
    Make predictions for query sequences based on BLAST hits.
    
    Args:
        query_ids: List of query sequence IDs
        blast_results: Dictionary mapping query IDs to hit information
        annotations: Training annotations sparse matrix (n_sequences x n_terms)
        use_rscore: Whether to use R-score (True) or sequence identity (False)
                   for weighting hits
                   
    Returns:
        Dictionary containing:
            - object: Query sequence IDs
            - score: Sparse matrix of prediction scores
            - date: Timestamp
            
    Note:
        For each query sequence, finds BLAST hits in training data,
        weights their annotations by sequence similarity, and takes
        the maximum score per term across all hits.
    """
    n_queries = len(query_ids)
    n_terms = len(terms)
    scores = sparse.lil_matrix((n_queries, n_terms), dtype=np.float32)
    
    # Remove self-hits
    blast_results['sseqid_acc'] = blast_results['sseqid'].apply(lambda x: x.split('|')[1])
    if not keep_self_hits:
        blast_results = blast_results[blast_results['sseqid_acc'] != blast_results['qseqid']]
    
    proteins_with_hits = set(blast_results['qseqid'])
    # Process each query sequence
    for i, qid in enumerate(query_ids):
        if qid not in proteins_with_hits:
            continue
            
        hits = blast_results[blast_results['qseqid'] == qid]
        hits = hits.sort_values('evalue', ascending=True)
        hits_unique = hits.drop_duplicates('sseqid_acc', keep='first') # keep only the first hit for each sseqid_acc result (lowest evalue)
        # hit_sseqids = set(hits['sseqid_acc'].unique())
        
        # Calculate weights for hits, make sure the order is retained
        weights = {}
        for _, hit in hits_unique.iterrows():
            if use_rscore:
                weights[hit['qseqid']] = calculate_rscore(hit['evalue'])
            else:
                weights[hit['sseqid_acc']] = hit['pident']/100.0
            
        # Get annotations for hit sequences and weight them
        query_scores = sparse.lil_matrix((len(weights), n_terms), dtype=np.float32)
        for i, (hit_seq, weight) in enumerate(weights.items()):
            if hit_seq in proteins:
                hit_annots = annotations[proteins[hit_seq]].multiply(weight) # assign same weight (pident) to all terms from this hit
                query_scores[i] = hit_annots
        query_scores = query_scores.tocsr()
        max_scores = query_scores.max(axis=0)
                
        scores[i] = max_scores
        
    # Add row and column indices
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
    output_file = 'data/processed/blast_predictions.tsv.gz'
    # Write the DataFrame to a gzipped TSV file using the Dask function
    write_dask_dataframe_to_gzipped_tsv(dask_df, output_file) # this function is defined in utils/write_dask.py
    return 
        
        
    
if __name__ == '__main__':
    # Load annotation matrix
    annotations = sparse.load_npz('S_matrix.npz')
    
    # Load protein and term indices
    with open('row_col_idx.pkl', 'rb') as f:
        proteins, terms = cp.load(f)
    
    # For blast baseline predictor, use terms and training proteins
    query_ids = []
    for record in SeqIO.parse('/work/idoerg/ahphan/democafa/A_data/data/test_superset.fasta','fasta'):
        query_ids.append(record.id)
    
    # Load BLAST results
    blast_results = pd.read_csv('/work/idoerg/ahphan/democafa/B_predictors/1.1_blast_database/blastp3.out', sep='\t', header=None, 
                                names=['qseqid', 'sseqid', 'evalue', 'length', 'pident', 'nident'])
    
    # Make predictions
    blast_predict(query_ids, blast_results, annotations, proteins, terms)