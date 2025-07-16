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
from Bio import SeqIO
import gzip
import scipy.sparse 
from democafa.utils.ontology import clean_ontology_edges, fetch_aspect, propagate_terms, filter_terms_given_obo


def wrapper_ground_truth(annot_known, annot2, query_file, graph, graph2, out_prefix):
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
    
    print(f"Number of proteins in test superset: {len(query_ids)}")
    # df1 = pd.read_csv(annot, sep='\t', header=0)
    dfk = pd.read_csv(annot_known, sep='\t', header=0) # we don't need df1 since it's replaced by dfk
    df2 = pd.read_csv(annot2, sep='\t', header=0)
    
    # Remove known terms from dfk out of df2 (it removes all the intermediate terms though)
    # common_pairs = dfk.merge(df2, on=['EntryID', 'term', 'aspect'], how='inner')
    # remove common pairs from df2
    dfk_df2 = dfk.merge(df2, on=df2.columns.to_list(), how='outer', indicator=True)
    print(dfk_df2['_merge'].value_counts())
    # print(f"Number of proteins in df2 after removing known terms: {len(df2_gain['EntryID'].unique())}")
    
    # During re-propagation with frozen graph, some obsolete terms (at t1) reappeared, so remove them
    # Creating "terms-of-interest" set (use this for GOslim too)
    toi = set.intersection(set(obonet.read_obo(graph).nodes()), set(obonet.read_obo(graph2).nodes()))
    dfk_df2_toi = dfk_df2[dfk_df2['term'].isin(toi)]
    dfk = dfk[dfk['term'].isin(toi)]
    df2 = df2[df2['term'].isin(toi)]
    
    df2_gain = dfk_df2_toi.loc[dfk_df2._merge=='right_only',dfk_df2.columns!='_merge']
    
    # NK: proteins in query_ids that are in df2 but not in dfk  
    print("Processing No Knowledge...")
    NK = set.intersection(set.difference(set(df2['EntryID']), set(dfk['EntryID'])), set(query_ids))
    NK_df = df2[df2['EntryID'].isin(NK)].reset_index(drop=True) # propagated already 
    # NK_t1 = filter_terms_given_obo(NK_df, current_graph=graph2, pivot_graph=graph)
    NK_df.to_csv(f'{out_prefix.replace(".tsv","_NK.tsv")}', sep="\t", index=False, header=True)
    
    # LK: proteins in query_ids that gained new aspect in df2
    print("Processing Limited Knowledge...")
    remaining_ids = set(query_ids) - set(NK)
    aspects_dfk = dfk[dfk['EntryID'].isin(remaining_ids)].groupby('EntryID')['aspect'].apply(set).reset_index()
    aspects_df2 = df2[df2['EntryID'].isin(remaining_ids)].groupby('EntryID')['aspect'].apply(set).reset_index()
    LK_dict = {}
    compare_df = pd.merge(aspects_dfk, aspects_df2, on='EntryID', suffixes=('_before', '_after'))
    for _, row in compare_df.iterrows():
        # if len(row['aspect_before']) > len(row['aspect_after']): # from 2 to 1 or 0, from 1 to 0
        #     # print(f"{row['EntryID']} lost aspects {row['aspect_before'] - row['aspect_after']}")
        #     pass
        # Note: edge case, protein lost 2 aspects but gained the third aspect, still LK
        # elif len(row['aspect_before']) <= len(row['aspect_after']): # from 1 to 2 or 3, from 2 to 3; or {'F'} to {'P'}, {'F','P'} to {'F','C'}
        diff = row['aspect_after'] - row['aspect_before']
        if diff:
            LK_dict[row['EntryID']] = diff
                # print(f"{row['EntryID']} gained aspects {diff} from {row['aspect_before']}")
    LK_df = pd.DataFrame()
    for protein, aspects in LK_dict.items():
        df = df2[(df2['EntryID'] == protein) & (df2['aspect'].isin(aspects))]
        LK_df = pd.concat([LK_df, df], ignore_index=True) # leaf only
    # LK_t1 = filter_terms_given_obo(LK_df, current_graph=graph2, pivot_graph=graph)
    LK_df.to_csv(f'{out_prefix.replace(".tsv","_LK.tsv")}', sep="\t", index=False, header=True) # new terms propagated

    # TODO: this PK part very slow    
    temp = df2_gain.groupby('EntryID')['aspect'].apply(set).reset_index()
    PK_dict = {}
    for _, row in temp.iterrows():
        if row['EntryID'] in LK_dict:
            gained_aspects = row['aspect'] - LK_dict[row['EntryID']]
            if gained_aspects:
                PK_dict[row['EntryID']] = gained_aspects
        else:
            if row['EntryID'] in remaining_ids: 
                PK_dict[row['EntryID']] = row['aspect']
    PK_df = pd.DataFrame()
    for protein, aspects in PK_dict.items():
        df = df2_gain[(df2_gain['EntryID'] == protein) & (df2_gain['aspect'].isin(aspects))]
        PK_df = pd.concat([PK_df, df], ignore_index=True)
    PK_df.to_csv(f'{out_prefix.replace(".tsv","_PK.tsv")}', sep="\t", index=False, header=True) # new terms propagated
    
    # TODO: protein-specific toi file for evaluation
    # # PK: proteins in query_ids that have the same aspects in df1 and df2, appear in df2_gain
    # print("Processing Partial Knowledge...")
    # compare_df['aspect_common'] = compare_df.apply(lambda row: row['aspect_before'].intersection(row['aspect_after']), axis=1)
    # temp_PK_aspects = pd.DataFrame([(p,aspect) for p, aspects in zip(compare_df['EntryID'], compare_df['aspect_common']) if aspects != set() for aspect in aspects], columns=['EntryID', 'aspect'])
    
    # # no protein in NK can be in PK or LK
    # assert not set(NK).intersection(set(temp_PK_aspects['EntryID'])) and not set(NK).intersection(set(LK_dict.keys())), "NK proteins should not be in PK or LK"
    
    # # Get only terms of proteins in PK
    # filter1 = pd.merge(dfk,temp_PK_aspects, on=['EntryID', 'aspect'])
    # filter2 = pd.merge(df2,temp_PK_aspects, on=['EntryID', 'aspect'])
    
    # # # Propagate filter1 using pivot graph
    # # roots = {'P': 'GO:0008150', 'C': 'GO:0005575', 'F': 'GO:0003674'}
    # # ontology_graph = clean_ontology_edges(obonet.read_obo(graph))
    # # subontologies = {aspect: fetch_aspect(ontology_graph, roots[aspect]) for aspect in roots}
    # # annotation_df1 = propagate_terms(filter1, subontologies)
    
    # # # Filter terms then propagate of filter2 with two graphs (do not filter terms in filter1 because it was generated from pivot graph)
    # # annotation_df2 = filter_terms_given_obo(filter2, current_graph=graph2, pivot_graph=graph)
    
    # # Dataframe comparison
    # temp_PK_df = pd.merge(temp_PK_aspects, df2_gain, on=['EntryID', 'aspect'], how='left')
    # # temp_PK_df['term_before'] = temp_PK_df['term']
    # # temp_PK_df = pd.merge(temp_PK_df, df2, on=['EntryID', 'aspect'], how='left')
    # # if term is in df2 but not in df1, term_before is NaN
    # sum(temp_PK_df['term'].isna()) # NaN means these proteins did not gain any new terms in the aspect 
    # PK_df = temp_PK_df[~temp_PK_df['term'].isna()][['EntryID', 'aspect', 'term']] # new terms only
    
    # # group by EntryID and aspect, if any term_before is NaN, then it is true, create new column for the boolean value
    # temp_PK_df['gained_PK'] = temp_PK_df.groupby(['EntryID', 'aspect'])['term_before'].transform(lambda x: x.isna().any())
    # PK_t1 = temp_PK_df[temp_PK_df['gained_PK'] == True][['EntryID','aspect', 'term', 'term_before']] # new terms and their old parents
    # PK_aspects = {p: aspect for p, aspect in zip(PK_t1['EntryID'], PK_t1['aspect'])}
    # assert len(PK_df['EntryID'].unique()) == len(PK_aspects), "PK_df should have unique EntryID"
    # PK_df.to_csv(f'{out_prefix.replace(".tsv","_PK.tsv")}', sep="\t", index=False, header=True) # write out leaf terms only
    
    print(f"Number of proteins in NK: {len(NK)}, with {len(NK_df)} terms")
    print(f"Number of proteins in LK: {len(LK_dict)}, with {len(LK_df)} terms")
    print(f"Number of proteins in PK: {len(PK_dict)}, with {len(PK_df)} new child terms gained")


# def efficient_matrix_comparison(matrix1, matrix2, pidx1, tidx1, pidx2, tidx2):
#     """Highly optimized implementation for sparse matrix comparison"""
#     # Create union of indices
#     all_proteins = sorted(set(pidx1.keys()) | set(pidx2.keys()))
#     all_terms = sorted(set(tidx1.keys()) | set(tidx2.keys()))
    
#     # Create reverse mappings for faster lookups
#     rev_pidx1 = {v: k for k, v in pidx1.items()}
#     rev_tidx1 = {v: k for k, v in tidx1.items()}
#     rev_pidx2 = {v: k for k, v in pidx2.items()}
#     rev_tidx2 = {v: k for k, v in tidx2.items()}
    
#     union_pidx = {pid: i for i, pid in enumerate(all_proteins)}
#     union_tidx = {tid: i for i, tid in enumerate(all_terms)}
#     new_shape = (len(all_proteins), len(all_terms))
    
#     # Create sparse matrices directly from nonzero entries
#     rows1, cols1, data1 = [], [], []
#     for i, j in zip(*matrix1.nonzero()):
#         pid = rev_pidx1[i]
#         tid = rev_tidx1[j]
#         rows1.append(union_pidx[pid])
#         cols1.append(union_tidx[tid])
#         data1.append(matrix1[i, j])
    
#     rows2, cols2, data2 = [], [], []
#     for i, j in zip(*matrix2.nonzero()):
#         pid = rev_pidx2[i]
#         tid = rev_tidx2[j]
#         rows2.append(union_pidx[pid])
#         cols2.append(union_tidx[tid])
#         data2.append(matrix2[i, j])
    
#     # Create sparse matrices
#     new_m1 = scipy.sparse.csr_matrix((data1, (rows1, cols1)), shape=new_shape)
#     new_m2 = scipy.sparse.csr_matrix((data2, (rows2, cols2)), shape=new_shape)
    
#     # Compute difference
#     diff_matrix = new_m2 - new_m1
    
#     # Get information about differences
#     diff_coords = diff_matrix.nonzero()
#     rev_union_pidx = {v: k for k, v in union_pidx.items()}
#     rev_union_tidx = {v: k for k, v in union_tidx.items()}
    
#     # Optional: create a more useful result format (protein_id, term_id, diff_value)
#     diff_info = []
#     for i, j in zip(*diff_coords):
#         diff_value = diff_matrix[i, j]
#         protein_id = rev_union_pidx[i]
#         term_id = rev_union_tidx[j]
#         diff_info.append((protein_id, term_id, diff_value))
    
#     return diff_matrix, union_pidx, union_tidx, diff_info
    

def parse_inputs(argv):
    parser = argparse.ArgumentParser(
        description='Compare two annotation files and classify proteins into No Knowledge, Limited Knowledge, and Partial Knowledge')
    
    # parser.add_argument('--annot', '-a', required=True,
    #                     help='Path to first annotation file (can be gzipped). This file is BEFORE.')
    parser.add_argument('--annot_known', '-ak', required=True,
                        help='Path to first annotation file (can be gzipped). This file is t0.')                        
    parser.add_argument('--annot2', '-a2', required=True,
                        help='Path to second annotation file (can be gzipped). This file is t1.')
    parser.add_argument('--query_file', '-q', required=True,
                        help='Path to target set of proteins, either in .fasta format or .tsv with protein IDs')
    parser.add_argument('--graph', '-g', default=None,
                        help='Path to OBO ontology graph at pivot timepoint (frozen graph).')
    parser.add_argument('--graph2', '-g2', default=None, 
                        help='Path to OBO ontology graph at a ground truth timepoint.')
    parser.add_argument('--out_prefix', default='groundtruth.tsv',
                        help='Prefix for 3 output files')
    return parser.parse_args(argv)

    # python3 -m democafa.groundtruth.classify_ground_truth -ak data/processed/cafa5/t0_terms.tsv -a2 data/processed/cafa5/t1_terms.tsv 
    # -q data/processed/cafa5/test_superset_all.fasta -g data/raw/go-basic.obo -g2 data/raw/t1_go-basic.obo --out_prefix data/processed/cafa5/t1_groundtruth.tsv  
    
    # python3 -m democafa.groundtruth.classify_ground_truth -ak data/processed/cafa5/t0_terms.tsv -a2 data/processed/cafa5/t1_terms.tsv 
    # -q data/processed/cafa5/test_superset_all.fasta -g data/raw/go-basic.obo -g2 ../robot_obo/results/t1_iba-slim-roots.obo --out_prefix data/processed/cafa5/goslim_t1_groundtruth.tsv
    
def main():
    args = parse_inputs(sys.argv[1:])
    wrapper_ground_truth(
        annot_known=args.annot_known,
        annot2=args.annot2,
        query_file=args.query_file,
        graph=args.graph,
        graph2=args.graph2,
        out_prefix=args.out_prefix        
    )
    
if __name__ == "__main__":
    main()