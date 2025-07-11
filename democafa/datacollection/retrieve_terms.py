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
from democafa.utils.ontology import clean_ontology_edges, fetch_aspect, propagate_terms, filter_terms_given_obo
from democafa.config import GO_CODES


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


def read_gaf(file_path, selected_codes):
    """
    Read and process a GAF file (gzipped or plain text)
    Takes 1h45m for 22Gb uniprot goa
    """
    #selected_codes = filter_evidence_codes(go_codes, selected) # handle this before passing to this function
    data = []
    
    # Get content starting from !gaf-version
    # content = process_gaf_file(handle)
    
    # from io import StringIO
    # gaf_file = StringIO(content)
    
    # Process annotations
    is_gzipped = file_path.endswith('.gz')
    if is_gzipped:
        handle = gzip.open(file_path, 'rt')
    else:
        handle = open(file_path, 'r')
    for rec in GOA.gafiterator(handle):
        if 'NOT' in rec['Qualifier']:
            continue
        if GOA.record_has(rec, selected_codes) and rec['DB'] == 'UniProtKB':
            data.append({'EntryID': rec['DB_Object_ID'], 'term': rec['GO_ID'], 'aspect': rec['Aspect']})
    
    df = pd.DataFrame(data)
    df = df.drop_duplicates()
    return df
    

def process_chunk_file(chunk_file, selected_codes):
    chunk_data = []
    with open(chunk_file, 'r') as f:
        for rec in GOA.gafiterator(f):
            if 'NOT' in rec['Qualifier']:
                continue
            if GOA.record_has(rec, selected_codes) and rec['DB'] == 'UniProtKB':
                chunk_data.append({
                    'EntryID': rec['DB_Object_ID'], 
                    'term': rec['GO_ID'], 
                    'aspect': rec['Aspect']
                })
    return chunk_data


def read_gaf_mp(file_path, taxon, selected_codes, use_mp=False, num_processes=None, chunk_size=10000000):
    """
    Read and process a GAF file (gzipped or plain text) with multiprocessing.
    Takes 30m for 22Gb uniprot goa for 15 processes.
    """
    if not use_mp:
        return read_gaf(file_path, selected_codes)
    else:
        import multiprocessing
        from functools import partial
        import tempfile
        if num_processes is None:
            num_processes = min(int(os.environ.get('SLURM_CPUS_PER_TASK', multiprocessing.cpu_count())) - 1, 16)
            num_processes = max(1, num_processes)
            print(f"Using {num_processes} processes for parallel computation")
        
        
        
        temp_dir = tempfile.mkdtemp(dir=os.path.dirname(file_path))
        chunk_files = []
        try:
            # Read and chunk the file
            print(f"Splitting GAF file into chunks of {chunk_size} records...")
            chunk_count = 0
            record_count = 0
            is_gzipped = file_path.endswith('.gz')
            open_func = gzip.open if is_gzipped else open
            mode = 'rt' if is_gzipped else 'r'
            
            with open_func(file_path, mode) as handle:
                current_chunk_lines = []
                headers = []
                line = handle.readline()
                while line.startswith('!'):
                    headers.append(line)
                    line = handle.readline()
                
                # Process first non-header line
                if line:
                    current_chunk_lines.append(line)
                
                # Process remaining lines
                for line in handle:
                    current_chunk_lines.append(line)
                    record_count += 1
                    
                    if len(current_chunk_lines) >= chunk_size:
                        # Write chunk to temporary file
                        chunk_file = os.path.join(temp_dir, f"chunk_{chunk_count}.gaf")
                        with open(chunk_file, 'w') as f:
                            # Write headers to each chunk
                            f.writelines(headers)
                            f.writelines(current_chunk_lines)
                        
                        chunk_files.append(chunk_file)
                        chunk_count += 1
                        current_chunk_lines = []
                        if chunk_count % 10 == 0:
                            print(f"Created chunk {chunk_count} with {chunk_size} records")
                # Write any remaining lines to the final chunk
                if current_chunk_lines:
                    chunk_file = os.path.join(temp_dir, f"chunk_{chunk_count}.gaf")
                    with open(chunk_file, 'w') as f:
                        # Write headers to each chunk
                        f.writelines(headers)
                        f.writelines(current_chunk_lines)
                    
                    chunk_files.append(chunk_file)
                    chunk_count += 1
                    print(f"Created final chunk {chunk_count} with {len(current_chunk_lines)} records")
                    
            # Process chunks in parallel
            print(f"Processing {len(chunk_files)} chunks with {num_processes} processes")
            process_func = partial(process_chunk_file, selected_codes=selected_codes)
            with multiprocessing.Pool(processes=num_processes) as pool:
                results = pool.map(process_func, chunk_files)
            
            all_data = [item for sublist in results for item in sublist]
            df = pd.DataFrame(all_data)
            df = df.drop_duplicates()
            return df
        finally:
            for chunk_file in chunk_files:
                if os.path.exists(chunk_file):
                    os.remove(chunk_file)
            if os.path.exists(temp_dir):
                os.rmdir(temp_dir)


def process_go_from_dat(file_path, selected_taxon, selected_codes):
    entries = []
    # selected_codes = filter_evidence_codes(GO_CODES, selected).get('Evidence') # handle this before passing to this function
    is_gzipped = file_path.endswith('.gz')
    if is_gzipped:
        handle = gzip.open(file_path, 'rt')
    else:
        handle = open(file_path, 'r')
        
    for record in sp.parse(handle):
        if not record.taxonomy_id:
            continue
        if selected_taxon is not None: # only filter if selected_taxon is not None
            if f'taxon:{record.taxonomy_id[0]}' not in selected_taxon: 
                continue
        current_id = record.accessions[0]
        for dr in record.cross_references:    #dr -> db cross refernce
            if dr[0] == 'GO' and len(dr) >= 4:
                go_id = dr[1]
                aspect = dr[2][0]  # Getting only the first letter (its either P/ C/ F)
                # aspect_description = dr[2][2:]  # Getting the rest of the description
                evidence = dr[3] if len(dr) >= 4 else ''
                evidence_code = evidence[:3]  # First 3 letters of evidence
                if evidence_code not in selected_codes:
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
    
    replaced = df['term'].isin(alt_id_to_id.keys())
    if replaced.any():
        df.loc[replaced, 'term'] = df.loc[replaced, 'term'].map(alt_id_to_id)
    return df


def wrapper_retrieve_terms(annot_file, selected_go_codes, graph, add_graph, output_tsv):
    # Load annotations from GOA or DAT file
    if 'gaf' in annot_file:
        filetype = 'gaf'
    elif 'dat' in annot_file:
        filetype = 'dat'
    else:
        raise ValueError("Unsupported annotation file type. Please provide a GAF or DAT file.")
    # if taxon is None:
    #     selected_taxon = None
    # elif taxon.isdigit(): 
    #     # If taxon is a single ID, convert it to a list
    #     selected_taxon = [f'taxon:{taxon}']
    # elif isinstance(taxon, str) and os.path.exists(taxon):
    #     taxon_file = pd.read_csv(taxon, header=0, sep='\t', encoding='ISO-8859-1')
    #     selected_taxon = [f'taxon:{taxon_id}' for taxon_id in set(taxon_file.iloc[:,0].tolist())]
    # print(f"Selected taxon: {selected_taxon}")
    
    if filetype == 'gaf':
        selected_codes = filter_evidence_codes(GO_CODES, selected_go_codes)
        annotation_df = read_gaf_mp(annot_file, selected_codes, use_mp=False, chunk_size=10000000)
    elif filetype == 'dat':
        selected_codes = filter_evidence_codes(GO_CODES, selected_go_codes).get('Evidence')
        annotation_df = process_go_from_dat(annot_file, selected_taxon, selected_codes)
    if annotation_df.empty:
        print(f"No annotations found for the given evidence codes {selected_go_codes}.")
        return
    
    # load ontology graph and GO terms. obonet doesn't store OBSOLETE terms
    if add_graph is None:
        ontology_graph = clean_ontology_edges(obonet.read_obo(graph))
    else:
        ontology_graph = clean_ontology_edges(obonet.read_obo(add_graph))
    annotation_df = replace_alternate_GO_terms(annotation_df, ontology_graph)
    
    obsolete_terms = set(annotation_df['term']) - set(ontology_graph.nodes())
    if obsolete_terms:
        print(f"Warning: {len(obsolete_terms)} obsolete terms ({obsolete_terms}) found in the annotation file.")
        print(f"These terms will not appear in terms file.")
        annotation_df = annotation_df[~annotation_df['term'].isin(obsolete_terms)]
    
    # Remove terms that are not in the frozen graph in 3 steps 
    # (propagate using graph2, intersect with graph terms, propagate again with graph)
    if add_graph is not None:
        annotation_df_filtered = filter_terms_given_obo(annotation_df, current_graph=add_graph, pivot_graph=graph)
        annotation_df_filtered = annotation_df_filtered.drop_duplicates()
        annotation_df_filtered.to_csv(output_tsv, sep='\t', index=False)
    else:
        annotation_df_filtered = annotation_df.drop_duplicates()
        annotation_df_filtered.to_csv(output_tsv, sep='\t', index=False)
    print(f"Annotations saved to {output_tsv} with {len(annotation_df_filtered)} annotations for {len(set(annotation_df_filtered['EntryID']))} proteins.")
       

def parse_inputs(argv):
    parser = argparse.ArgumentParser(
        description='Retrieve terms annotated to UniProtKB proteins from GOA file and save to TSV file')
    
    parser.add_argument('--annot', '-a', required=True,
                        help='Path to first annotation file (can be gzipped)')
    # parser.add_argument('--taxon', '-t', default='9606', required=False,
    #                     help='Taxon ID file to filter proteins (default: 9606 for human)')
    parser.add_argument('--selected_go_codes', '-sgc', required=False, default='Experimental,IC,TAS',
                        help='Comma-separated list of evidence codes to include in the analysis')
    parser.add_argument('--graph', '-g', default=None, 
                        help='Path to OBO ontology graph file if local. If empty (default) current OBO structure at run-time will be downloaded from http://purl.obolibrary.org/obo/go/go-basic.obo')
    parser.add_argument('--add_graph', '-ag', default=None,
                        help='Path to OBO ontology graph of a later timepoint. Provide this graph to remove terms that are not in frozen graph.')
    parser.add_argument('--tsv', default='data/release/train_terms.tsv',
                        help='Path to save annotations in TSV format')
    return parser.parse_args(argv)

    # python3 -m democafa.datacollection.retrieve_terms --annot data/processed/cafa5/goa_uniprot_filtered_mp.gaf.213.gz
    # --taxon data/raw/testsuperset-taxon-list.tsv -sgc 'Experimental,IC,TAS' -g data/raw/go-basic.obo --tsv data/processed/cafa5/train_terms.tsv
    
    # python3 -m democafa.datacollection.retrieve_terms --annot data/processed/cafa5/goa_uniprot_filtered_mp.gaf.216.gz
    # --taxon data/raw/testsuperset-taxon-list.tsv -sgc 'Experimental,IC,TAS' -g data/raw/go-basic.obo -ag data/raw/t0_go-basic.obo --tsv data/processed/cafa5/t0_terms.tsv
    
    # python3 -m democafa.datacollection.retrieve_terms --annot data/processed/cafa5/goa_uniprot_filtered_mp.gaf.224.gz
    # --taxon data/raw/testsuperset-taxon-list.tsv -sgc 'Experimental,IC,TAS' -g data/raw/go-basic.obo -ag data/raw/t1_go-basic.obo --tsv data/processed/cafa5/t1_terms.tsv
    
    # TODO: test this without taxon filtering and mp
def main():
    args = parse_inputs(sys.argv[1:])
    wrapper_retrieve_terms(
        annot_file=args.annot,
        # taxon=args.taxon,
        selected_go_codes=args.selected_go_codes,
        graph=args.graph,
        add_graph=args.add_graph,
        output_tsv=args.tsv
    )
    
if __name__ == "__main__":
    main()
    
