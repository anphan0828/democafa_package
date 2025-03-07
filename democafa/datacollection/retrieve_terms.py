#!/usr/bin/env python3

"""
Script to extract protein and GO terms labels from a GOA file (or a .dat file).

This script:
1. Gunzips the GOA gzipped file 
2. Process NOT annotations and selected evidence codes
3. Creates a tsv file with 3 columns: EntryID, term, aspect

This script can be run from command line with arguments, or imported as a module.
"""


import os
import sys
import argparse
import gzip
import obonet
import pandas as pd
import networkx as nx
import numpy as np
from Bio.UniProt import GOA
from Bio import SwissProt as sp

def parse_inputs(argv):
    parser = argparse.ArgumentParser(
        description='Retrieve terms annotated to UniProtKB proteins from GOA file and save to TSV file')
    
    parser.add_argument('--annot', '-a', required=True,
                        help='Path to first annotation file (can be gzipped)')
    parser.add_argument('--evidence', '-e', required=False, default='Experimental,IC,TAS',
                        help='Comma-separated list of evidence codes to include in the analysis')
    parser.add_argument('--filetype', '-t', required=True, choices=['goa', 'dat'], 
                        help='Input file type')
    parser.add_argument('--graph', '-g', default=None, 
                        help='Path to OBO ontology graph file if local. If empty (default) current OBO structure at run-time will be downloaded from http://purl.obolibrary.org/obo/go/go-basic.obo')
    parser.add_argument('--tsv', default='data/release/train_terms.tsv',
                        help='Path to save annotations in TSV format')
    return parser.parse_args(argv)


def process_gaf_file(gaf_file):
    '''
    function : given a file handle, find the !gaf-version line and return all content from there
    input    : file path
    output   : content starting from !gaf-version line
    '''
    # Check if file is gzipped
    is_gzipped = gaf_file.endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    with open_func(gaf_file, mode) as f:
        content = f.read()
        
    # Find the position of !gaf-version line
    gaf_version_pos = content.find('!gaf-version')
    if gaf_version_pos == -1:
        return content  # If not found, return all content
    return content[gaf_version_pos:]  # Return content from !gaf-version onwards


def clean_ontology_edges(ontology):
    """
    Remove all ontology edges except types "is_a" and "part_of" and ensure there are no inter-ontology edges
    :param ontology: Ontology stucture (networkx DiGraph or MultiDiGraph)
    """
    
    # keep only "is_a" and "part_of" edges (All the "regulates" edges are in BPO)
    remove_edges = [(i, j, k) for i, j, k in ontology.edges if not(k=="is_a" or k=="part_of")]
    
    ontology.remove_edges_from(remove_edges)
    
    # There should not be any cross-ontology edges, but we verify here
    crossont_edges = [(i, j, k) for i, j, k in ontology.edges if
                      ontology.nodes[i]['namespace']!= ontology.nodes[j]['namespace']]
    if len(crossont_edges)>0:
        ontology.remove_edges_from(crossont_edges)
    
    return ontology
    

def fetch_aspect(ontology, root:str):
    """
    Return a subgraph of an ontology starting at node <root>
    
    :param ontology: Ontology stucture (networkx DiGraph or MultiDiGraph)
    :param root: node name (GO term) to start subgraph
    """
    namespace = ontology.nodes[root]['namespace']
    aspect_nodes = [n for n,v in ontology.nodes(data=True) 
                    if v['namespace']==namespace]
    subont_ = ontology.subgraph(aspect_nodes)
    return subont_

def filter_evidence_codes(go_codes, selected='Experimental,IC,TAS'):
    ALL_LISTS = go_codes
    ALL_CODES = set(code for codes in ALL_LISTS.values() for code in codes)
    input_codes = [code.strip().upper() for code in selected.split(',')]
    accepted_codes = set()
    for code in input_codes:
        if code in ALL_LISTS:
            accepted_codes.update(ALL_LISTS[code])
        elif code in ALL_CODES:
            accepted_codes.add(code)
        else:
            print(f"Warning: '{code}' is not a recognized evidence code.")
    return {'Evidence': set(accepted_codes)}

def read_gaf(handle, go_codes, selected):
    """
    Read and process a GAF file (gzipped or plain text)
    """
    name = os.path.splitext(os.path.basename(handle))[0]
    if name.endswith('.gaf'):
        name = os.path.splitext(name)[0]
    
    all_protein_name = set()
    selected_codes = filter_evidence_codes(go_codes, selected)
    data = []
    
    # Get content starting from !gaf-version
    content = process_gaf_file(handle)
    
    # Process annotations using StringIO to create file-like object
    from io import StringIO
    gaf_file = StringIO(content)
    
    # Process annotations
    for rec in GOA.gafiterator(gaf_file):
        # Remove NOT annotations
        if 'NOT' in rec['Qualifier']:
            continue
        all_protein_name.add(rec['DB_Object_ID'])
        # Add the ancestral terms to the dictionary
        if GOA.record_has(rec, selected_codes) and rec['DB'] == 'UniProtKB':
            data.append({'EntryID': rec['DB_Object_ID'], 'term': rec['GO_ID'], 'aspect': rec['Aspect']})
    
    df = pd.DataFrame(data)
    df = df.drop_duplicates()
    return name, df, all_protein_name


def process_go_from_dat(file_path, go_codes, selected):
    # TODO: Add support for other evidence codes
    entries = []
    with open(file_path, 'r') as file:
        for record in sp.parse(file):
            if not record.taxonomy_id:
                continue
            if '9606' not in record.taxonomy_id:  # Only human proteins
                continue
            current_id = record.accessions[0]
            for dr in record.cross_references:    #dr -> db cross refernce
                if dr[0] == 'GO' and len(dr) >= 4:
                    go_id = dr[1]
                    aspect = dr[2][0]  # Getting only the first letter (its either P/ C/ F)
                    # aspect_description = dr[2][2:]  # Getting the rest of the description
                    evidence = dr[3] if len(dr) >= 4 else ''
                    evidence_code = evidence[:3]  # First 3 letters of evidence
                    if evidence_code not in selected:
                        continue
                    # evidence_source = evidence[4:] if len(evidence) > 4 else ''
                    entries.append({
                        "EntryID": current_id,
                        "term": go_id,
                        "aspect": aspect
                        #"Aspect_Description": aspect_description.strip(),
                        # "evidence": evidence_code.strip(),
                        # #"Evidence Source": evidence_source.strip()
                    })
    df = pd.DataFrame(entries)
    df = df.drop_duplicates()
    return df

# TODO: combine train_terms.tsv and train_taxonomy.tsv (maybe not for blast reasons)
def replace_alternate_GO_terms(df, ontology_graph):
    """
    Replace alternate GO terms with their main GO terms using vectorized operations.
    
    Args:
        df: pandas DataFrame containing GO terms in 'term' column
        graph: path to the ontology graph file (obo)
        
    Returns:
        DataFrame with updated GO terms where alternate IDs are replaced with main IDs
    """
    # Create a dictionary to map alternate GO terms to main GO terms
    alt_id_to_id = {}
    for id_, data in ontology_graph.nodes(data=True):
        alt_ids = data.get('alt_id')
        if alt_ids:
            for alt_id in alt_ids:
                alt_id_to_id[alt_id] = id_
    
    # Use map with fillna to replace terms - more efficient than apply
    df['term'] = df['term'].map(lambda x: alt_id_to_id.get(x, x))
    return df

def wrapper_retrieve_terms(annot_file, filetype, go_codes, selected_go_codes, graph, output_tsv):
    # load ontology graph and GO terms. obonet doesn't store OBSOLETE terms
    ontology_graph = clean_ontology_edges(obonet.read_obo(graph))
    # roots = {'P': 'GO:0008150', 'C': 'GO:0005575', 'F': 'GO:0003674'}
    # subontologies = {aspect: fetch_aspect(ontology_graph, roots[aspect]) for aspect in roots}
    
    if filetype == 'goa':
        name,annotation_df,all_protein = read_gaf(annot_file, go_codes, selected_go_codes)
    elif filetype == 'dat':
        annotation_df = process_go_from_dat(annot_file, go_codes, selected_go_codes)
    obsolete_terms = set(annotation_df['term']) - set(ontology_graph.nodes())
    if obsolete_terms:
        print(f"Warning: {len(obsolete_terms)} obsolete terms ({obsolete_terms}) found in the annotation file.")
        annotation_df = annotation_df[~annotation_df['term'].isin(obsolete_terms)]
    annotation_df = replace_alternate_GO_terms(annotation_df, ontology_graph)
    # Write annotations to TSV file
    annotation_df.to_csv(output_tsv, sep='\t', index=False)
    
def main():
    args = parse_inputs(sys.argv[1:])
    wrapper_retrieve_terms(
        annot_file=args.annot,
        filetype=args.filetype,
        go_codes=args.go_codes,
        selected_go_codes=args.evidence,
        graph=args.graph,
        output_tsv=args.tsv
    )
    
if __name__ == "__main__":
    main()
    
