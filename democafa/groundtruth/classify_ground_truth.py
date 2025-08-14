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
    is_gzipped = query_file.endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    if '.fasta' in query_file:
        print("Reading query IDs from FASTA file")
        with open_func(query_file, mode) as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                entry_id = record.id.split("|")[1] if "|" in record.id else record.id
                query_ids.append(entry_id)
    elif '.txt' in query_file:
        print("Reading query IDs from text file")
        with open(query_file, 'r') as handle:
            query_ids = [line.strip() for line in handle]
    else:
        print("Please provide a fasta file or a text file with query IDs")
        sys.exit(1)
    
    print(f"Number of proteins in query file: {len(query_ids)}")
    dfk = pd.read_csv(annot_known, sep='\t', header=0) 
    df2 = pd.read_csv(annot2, sep='\t', header=0)
    
    # During re-propagation with frozen graph, some obsolete terms (at t1) reappeared, so remove them
    # Creating "terms-of-interest" set (use this for GOslim too)
    toi = set.intersection(set(obonet.read_obo(graph).nodes()), set(obonet.read_obo(graph2).nodes()))
    # dfk_df2_toi = dfk_df2[dfk_df2['term'].isin(toi)]
    dfk = dfk[dfk['term'].isin(toi)]
    df2 = df2[df2['term'].isin(toi)]
    
    
    # Compare dfk and df2 to find proteins that gained new terms
    dfk_df2_toi = dfk.merge(df2, on=df2.columns.to_list(), how='outer', indicator=True)
    print(dfk_df2_toi['_merge'].value_counts())
    df2_gain = dfk_df2_toi.loc[dfk_df2_toi._merge=='right_only',dfk_df2_toi.columns!='_merge'] # new terms
    print(f"Number of terms gained in each aspect: {df2_gain['aspect'].value_counts()}")
    
    p_asp_pairs = df2_gain[['EntryID', 'aspect']].drop_duplicates().reset_index(drop=True)
    print(f"Number of proteins gaining terms: {len(df2_gain['EntryID'].unique())}, with {len(p_asp_pairs)} protein-aspect pairs")
    
    # Remove already existing annotations in df2 (biocuration redundancy)
    redundant_proteins = set.difference(set(df2['EntryID']), set(df2_gain['EntryID']))
    if len(redundant_proteins) > 0:
        print(f"There were {len(redundant_proteins)} proteins whose new annotations already existed: {redundant_proteins}")
        
    # NK: proteins in query_ids that are in df2 but not in dfk  
    print("Processing No Knowledge...")
    NK = set.intersection(set.difference(set(df2['EntryID']), set(dfk['EntryID'])), set(query_ids))
    NK_df = df2[df2['EntryID'].isin(NK)].reset_index(drop=True) # propagated already 
    NK_asp_pairs = NK_df[['EntryID', 'aspect']].drop_duplicates().reset_index(drop=True)
    # NK_t1 = filter_terms_given_obo(NK_df, current_graph=graph2, pivot_graph=graph)
    NK_df.to_csv(f'{out_prefix.replace(".tsv","_NK.tsv")}', sep="\t", index=False, header=True)
    
    # LK: proteins in query_ids that gained new aspect in df2
    print("Processing Limited Knowledge and Partial Knowledge...")
    remaining_ids = set.intersection(set(query_ids), set(df2_gain['EntryID'])) - set(NK) # proteins in query_ids that gained new terms
    dfk_asp_pairs = dfk[dfk['EntryID'].isin(remaining_ids)].groupby('EntryID')['aspect'].apply(set).reset_index() # protein-aspect pairs in t0
    
    # Check each protein-aspect pair in gained dataframe
    # If aspect is in dfk then it is Limited Knowledge, otherwise Partial Knowledge
    LK_dict = {}
    PK_dict = {}
    for _, row in p_asp_pairs.iterrows():
        if row['EntryID'] not in remaining_ids:
            continue # already classified as NK
        dfk_asp = dfk_asp_pairs[dfk_asp_pairs['EntryID'] == row['EntryID']]['aspect'].values[0]
        if row['aspect'] in dfk_asp:
            if row['EntryID'] not in PK_dict:
                PK_dict[row['EntryID']] = set(row['aspect'])
            else:
                PK_dict[row['EntryID']].add(row['aspect'])
        else:
            if row['EntryID'] not in LK_dict:
                LK_dict[row['EntryID']] = set(row['aspect'])
            else:
                LK_dict[row['EntryID']].add(row['aspect'])
    
    LK_df = pd.DataFrame()
    for protein, aspects in LK_dict.items():
        df = df2_gain[(df2_gain['EntryID'] == protein) & (df2_gain['aspect'].isin(aspects))]
        LK_df = pd.concat([LK_df, df], ignore_index=True) # leaf only
    LK_asp_pairs = LK_df[['EntryID', 'aspect']].drop_duplicates()
    LK_df.to_csv(f'{out_prefix.replace(".tsv","_LK.tsv")}', sep="\t", index=False, header=True) # new terms propagated

    PK_df = pd.DataFrame()
    for protein, aspects in PK_dict.items():
        df = df2_gain[(df2_gain['EntryID'] == protein) & (df2_gain['aspect'].isin(aspects))]
        PK_df = pd.concat([PK_df, df], ignore_index=True)
    PK_asp_pairs = PK_df[['EntryID', 'aspect']].drop_duplicates().reset_index(drop=True)
    PK_df.to_csv(f'{out_prefix.replace(".tsv","_PK.tsv")}', sep="\t", index=False, header=True) # new terms propagated
    
    assert len(NK_asp_pairs) + len(LK_asp_pairs) + len(PK_asp_pairs) == len(p_asp_pairs), "Sum of aspect pairs in NK, LK, PK should equal total gained aspect pairs"
    # TODO: protein-specific toi file for evaluation
        
    print(f"Number of proteins in NK: {len(NK)}, with {len(NK_df)} terms")
    print(f"Number of proteins in LK: {len(LK_dict)}, with {len(LK_df)} terms")
    print(f"Number of proteins in PK: {len(PK_dict)}, with {len(PK_df)} new child terms gained")
    
    all_targets = set.union(set(NK), set(LK_dict.keys()), set(PK_dict.keys()))
    print(f"Total number of target proteins: {len(all_targets)}")
    with open(f'{out_prefix.replace(".tsv","_targets.tsv")}', 'w') as f:
        for target in all_targets:
            f.write(f"{target}\n")
            


def parse_inputs(argv):
    parser = argparse.ArgumentParser(
        description='Compare two annotation files and classify proteins into No Knowledge, Limited Knowledge, and Partial Knowledge')
    
    
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