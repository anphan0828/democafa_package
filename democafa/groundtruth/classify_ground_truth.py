#!/usr/bin/env python3

"""
Script to classify proteins in evaluation set into 3 classes No Knowledge, Limited Knowledge, and Partial Knowledge

This script:
1. Reads training data and ground truth data on proteins and their GO terms
2. Compares annotation status of proteins in the evaluation set with the training data
3. Creates 3 tsv files with 2 columns: EntryID, term, each file for one type of protein class

This script can be run from command line with arguments, or imported as a module.
"""

import sys
import obonet
import argparse
import pandas as pd
import numpy as np
import networkx as nx
from Bio import SeqIO
import gzip
import scipy.sparse 
from democafa.utils.ontology import clean_ontology_edges, fetch_aspect, propagate_terms, approach3_optimized


def parse_inputs(argv):
    parser = argparse.ArgumentParser(
        description='Compare two annotation files and classify proteins into No Knowledge, Limited Knowledge, and Partial Knowledge')
    
    parser.add_argument('--annot', '-a', required=True,
                        help='Path to first annotation file (can be gzipped). This file is BEFORE.')
    parser.add_argument('--annot2', '-a2', required=True,
                        help='Path to second annotation file (can be gzipped). This file is AFTER.')
    parser.add_argument('--query_file', '-1', required=True,
                        help='Path to target set of proteins, either in .fasta format or .tsv with protein IDs')
    parser.add_argument('--filetype', '-t', required=True, choices=['goa', 'dat'], 
                        help='Input file type')
    parser.add_argument('--graph', '-g', default=None, 
                        help='Path to OBO ontology graph.')
    parser.add_argument('--out_prefix', default='groundtruth.tsv',
                        help='Prefix for 3 output files')
    return parser.parse_args(argv)


def wrapper_ground_truth(annot, annot2, query_file, filetype, current_graph, pivot_graph, out_prefix):
    # Collect all target proteins
    query_ids = []
    if query_file.endswith('.fasta'):
        print("Reading query IDs from FASTA file")
        with open(query_file, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                entry_id = record.id.split("|")[1] if "|" in record.id else record.id
                query_ids.append(entry_id)
    elif query_file.endswith('.txt'):
        print("Reading query IDs from text file")
        with open(query_file, 'r') as handle:
            query_ids = [line.strip() for line in handle]
    else:
        print("Please provide a fasta file or a text file with query IDs")
        sys.exit(1)
    
    df1 = pd.read_csv(annot, sep='\t', header=0)
    df2 = pd.read_csv(annot2, sep='\t', header=0)
    
    # NK: proteins in query_ids that are in df2 but not in df1
    NK = set.intersection(set.difference(set(df2['EntryID']), set(df1['EntryID'])), set(query_ids))
    NK_df = df2[df2['EntryID'].isin(NK)].reset_index(drop=True)
    NK_t1 = filter_terms_given_obo(NK_df, current_graph, pivot_graph)
    
    
    # LK: proteins in query_ids that gained new aspect in df2
    remaining_ids = set(query_ids) - set(NK)
    aspects_df1 = df1[df1['EntryID'].isin(remaining_ids)].groupby('EntryID')['aspect'].apply(set).reset_index()
    aspects_df2 = df2[df2['EntryID'].isin(remaining_ids)].groupby('EntryID')['aspect'].apply(set).reset_index()
    LK_dict = {}
    compare_df = pd.merge(aspects_df1, aspects_df2, on='EntryID', suffixes=('_before', '_after'))
    for _, row in compare_df.iterrows():
        if len(row['aspect_before']) > len(row['aspect_after']): # from 2 to 1 or 0, from 1 to 0
            print(f"{row['EntryID']} lost aspects {row['aspect_before'] - row['aspect_after']}")
        elif len(row['aspect_before']) <= len(row['aspect_after']): # from 1 to 2 or 3, from 2 to 3; or {'F'} to {'P'}, {'F','P'} to {'F','C'}
            diff = row['aspect_after'] - row['aspect_before']
            if diff:
                LK_dict[row['EntryID']] = diff
                # print(f"{row['EntryID']} gained aspects {diff} from {row['aspect_before']}")
    LK_df = pd.DataFrame()
    for protein, aspects in LK_dict.items():
        df = df2[(df2['EntryID'] == protein) & (df2['aspect'].isin(aspects))]
        LK_df = pd.concat([LK_df, df], ignore_index=True)
    LK_t1 = filter_terms_given_obo(LK_df, current_graph, pivot_graph)
    
    
    # PK: proteins in query_ids that have the same aspects in df1 and df2
    compare_df['aspect_common'] = compare_df.apply(lambda row: row['aspect_before'].intersection(row['aspect_after']), axis=1)
    # temp_PK_dict = {p: aspects for p,aspects in zip(compare_df['EntryID'], compare_df['aspect_common']) if aspects != set()}
    temp_PK_aspects = pd.DataFrame([(p,aspect) for p, aspects in zip(compare_df['EntryID'], compare_df['aspect_common']) if aspects != set() for aspect in aspects], columns=['EntryID', 'aspect'])
    # no protein in NK can be in PK or LK
    assert not set(NK).intersection(set(temp_PK_aspects['EntryID'])) and not set(NK).intersection(set(LK_dict.keys())), "NK proteins should not be in PK or LK"
    
    propagate_and_compare(temp_PK_aspects, df1, df2, current_graph, pivot_graph, query_ids)
    

def filter_terms_given_obo(terms_df, current_graph, pivot_graph):
    """Remove terms on a future obo that is not in the chosen pivot graph"""
    
    # Propagate using current graph
    ontology_graph = clean_ontology_edges(obonet.read_obo(current_graph))
    roots = {'P': 'GO:0008150', 'C': 'GO:0005575', 'F': 'GO:0003674'}
    subontologies = {aspect: fetch_aspect(ontology_graph, roots[aspect]) for aspect in roots} 
    
    prop_terms_df = propagate_terms(terms_df, subontologies)
    before_length = len(prop_terms_df)
    
    # Compare with pivot graph, remove terms not in pivot graph
    ontology_graph = clean_ontology_edges(obonet.read_obo(pivot_graph))
    
    prop_terms_df = prop_terms_df[prop_terms_df['term'].isin(ontology_graph.nodes())].reset_index(drop=True)
    after_length = len(prop_terms_df)
    print(f"Filtered terms from {before_length} to {after_length} using ontology graph {pivot_graph.split('/')[-1]}")
    
    # Propagate terms using the chosen pivot graph
    subontologies = {aspect: fetch_aspect(ontology_graph, roots[aspect]) for aspect in roots} 
    annotation_df = propagate_terms(prop_terms_df, subontologies)
    
    return annotation_df

    
def propagate_and_compare(temp_PK_aspects, df1, df2, current_graph, pivot_graph, query_ids):
    
    # Get only terms of proteins in PK
    filter1 = pd.merge(df1,temp_PK_aspects, on=['EntryID', 'aspect'])
    filter2 = pd.merge(df2,temp_PK_aspects, on=['EntryID', 'aspect'])
    
    # Propagate filter1 using pivot graph
    roots = {'P': 'GO:0008150', 'C': 'GO:0005575', 'F': 'GO:0003674'}
    ontology_graph = clean_ontology_edges(obonet.read_obo(pivot_graph))
    subontologies = {aspect: fetch_aspect(ontology_graph, roots[aspect]) for aspect in roots}
    annotation_df1 = propagate_terms(filter1, subontologies)
    # Filter terms then propagate of filter2 with two graphs (do not filter terms in filter1 because it was generated from pivot graph)
    annotation_df2 = filter_terms_given_obo(filter2, current_graph, pivot_graph)
    
    ## Approach 3: create csr directly from coordinates
    matrix1, pidx1, tidx1, _ = approach3_optimized(annotation_df1)
    matrix2, pidx2, tidx2, _ = approach3_optimized(annotation_df2)
    
    # TODO: compare these matrices with different shapes
    diff_matrix, union_pidx, union_tidx = efficient_matrix_comparison(matrix1, matrix2, pidx1, tidx1, pidx2, tidx2)
    
    
def efficient_matrix_comparison(matrix1, matrix2, pidx1, tidx1, pidx2, tidx2):
    """More efficient implementation using COO format for construction"""
    # Create union of indices
    all_proteins = sorted(set(pidx1.keys()) | set(pidx2.keys()))
    all_terms = sorted(set(tidx1.keys()) | set(tidx2.keys()))
    
    union_pidx = {pid: i for i, pid in enumerate(all_proteins)}
    union_tidx = {tid: i for i, tid in enumerate(all_terms)}
    new_shape = (len(all_proteins), len(all_terms))
    
    # Build COO data for matrix1
    rows1, cols1, data1 = [], [], []
    for pid in pidx1:
        for tid in tidx1:
            val = matrix1[pidx1[pid], tidx1[tid]]
            if val != 0:
                rows1.append(union_pidx[pid])
                cols1.append(union_tidx[tid])
                data1.append(val)
    
    # Build COO data for matrix2
    rows2, cols2, data2 = [], [], []
    for pid in pidx2:
        for tid in tidx2:
            val = matrix2[pidx2[pid], tidx2[tid]]
            if val != 0:
                rows2.append(union_pidx[pid])
                cols2.append(union_tidx[tid])
                data2.append(val)
    
    # Create sparse matrices
    new_m1 = scipy.sparse.csr_matrix((data1, (rows1, cols1)), shape=new_shape)
    new_m2 = scipy.sparse.csr_matrix((data2, (rows2, cols2)), shape=new_shape)
    
    diff_matrix = new_m2 - new_m1
    
    return diff_matrix, union_pidx, union_tidx
    
    
def main():
    args = parse_inputs(sys.argv[1:])
    wrapper_ground_truth(
        annot=args.annot,
        annot2=args.annot2,
        query_file=args.query_file,
        filetype=args.filetype,
        graph=args.graph,
        out_prefix=args.out_prefix        
    )
    
if __name__ == "__main__":
    main()