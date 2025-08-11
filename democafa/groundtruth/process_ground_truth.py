#!/usr/bin/env python3

import sys
import argparse
import requests
import obonet
import numpy as np
import pandas as pd
import networkx as nx
import pickle as cp
import time
from collections import Counter
from democafa.utils.ontology import add_aspect_column, clean_ontology_edges, fetch_aspect, propagate_terms
from democafa.datacollection.retrieve_terms import filter_evidence_codes
from democafa.config import GO_CODES

input_file = '/work/idoerg/ahphan/democafa_package/data/raw/GO_annotations_CAFA_20230610.tsv'
ontology_graph = '/work/idoerg/ahphan/democafa_package/data/raw/go-basic-20250601.obo'
ontology_slim = '/work/idoerg/ahphan/robot_obo/results/go-slim-20250601.obo'
taxon_file = '/work/idoerg/ahphan/democafa_package/data/raw/testsuperset-taxon-list.tsv'
output_file = '/work/idoerg/ahphan/democafa_package/data/processed/cafa6/ground_truth_20250610.tsv'

def process_tsv(input_file, ontology_graph, output_file):
    """
    Process a TSV file to add an 'aspect' column based on GO terms.
    
    Args:
        input_file: Path to the input TSV file.
        output_file: Path to the output TSV file with added 'aspect' column.
    """
    ontology_graph = clean_ontology_edges(obonet.read_obo(ontology_graph))
    roots = {'P': 'GO:0008150', 'C': 'GO:0005575', 'F': 'GO:0003674'}
    subontologies = {aspect: fetch_aspect(ontology_graph, roots[aspect]) for aspect in roots}
    
    annotation_df = pd.read_csv(input_file, sep='\t', header=0)
    print(annotation_df.head())
    annotation_df = annotation_df[['ENTITY_ID', 'GO_ID', 'GO_EVIDENCE']]
    annotation_df.columns = ['EntryID', 'term', 'evidence']
    annotation_df = add_aspect_column(annotation_df, subontologies)
    
    # Filter for EXPERIMENTAL,IC,TAS
    selected_codes = filter_evidence_codes(GO_CODES, selected='Experimental,IC,TAS').get('Evidence')
    annotation_df = annotation_df[annotation_df['evidence'].isin(selected_codes)]
    
    obsolete_terms = set(annotation_df['term']) - set(ontology_graph.nodes())
    if obsolete_terms:
        print(f"Warning: {len(obsolete_terms)} obsolete terms ({obsolete_terms}) found in the annotation file.")
        print(f"These terms will not appear in terms file.")
        annotation_df = annotation_df[~annotation_df['term'].isin(obsolete_terms)]
    
    annotation_df = propagate_terms(annotation_df, subontologies) 
    print(annotation_df.describe())
    
    # TODO: Remove 3 binding terms of proteins even though they are annotated in other aspects
    binding_terms = set(nx.descendants(subontologies['F'],'GO:0005515'))
    binding_terms.add('GO:0005515')
    binding_only = annotation_df[annotation_df['term'].isin(binding_terms)]
    binding_only = binding_only.groupby(['EntryID','aspect']).size().reset_index(name='count')
    
    # print(f"Removed {len(annotation_df['EntryID'].unique()) - len(nonbinding_df['EntryID'].unique())} proteins with only protein-binding terms.")
    # annotation_df = annotation_df[annotation_df['EntryID'].isin(nonbinding_df['EntryID'])]
    # TODO: remove trembl proteins
    
    
    # Save to 3-column TSV
    annotation_df[['EntryID', 'term', 'aspect']].to_csv(output_file, sep='\t', index=False, header=True)
    
    # (Solved) filter for new species to be added to training and test sets: Using all SwissProt species now
    old_taxon = pd.read_csv(taxon_file, header=0, sep='\t', encoding='ISO-8859-1')
    old_taxid = [str(taxon_id) for taxon_id in set(old_taxon.iloc[:,0].tolist())]
    
    batch = list(annotation_df['EntryID'].unique())
    url = "https://rest.uniprot.org/idmapping/run"
    taxon = dict()
    data = {
            'from': 'UniProtKB_AC-ID',
            'to': 'UniProtKB',
            'ids': f"{','.join(batch)}"
        }
    review_status = dict()
    # Submit job
    response = requests.post(url, data=data)
    job_id = response.json()['jobId']
    results_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    while True:
        status_response = requests.get(results_url)
        status = status_response.json()
        if 'results' in status:
            # Download results
            stream_url = f"https://rest.uniprot.org/idmapping/uniprotkb/results/stream/{job_id}?fields=accession%2Creviewed%2Corganism_name%2Corganism_id&format=tsv"
            tsv_response = requests.get(stream_url)
            tsv_content = tsv_response.text
            
            for line in tsv_content.strip().split('\n')[1:]:
                accession, _, reviewed, organism_name, taxid = line.split('\t')
                review_status[accession] = reviewed
                taxon[accession] = (organism_name, taxid)
            break
        time.sleep(1)

    # new_taxon = {name: id for (name, id) in taxon.values() if id not in old_taxid}

    
    # Intersect with GOslim
    ontology_slim = clean_ontology_edges(obonet.read_obo(ontology_slim))
    slim_df = annotation_df[(annotation_df['term'].isin(ontology_slim.nodes())) & 
                            (~annotation_df['term'].isin(roots.values()))] # proteins with only roots (after slim) will be removed
    print(slim_df.describe())                            
    
    
    