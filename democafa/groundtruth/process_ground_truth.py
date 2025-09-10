#!/usr/bin/env python3

import sys
import argparse
import os
import gzip
import pandas as pd
import re
from Bio import SeqIO
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
from Bio import Entrez
Entrez.email = "ahphan@iastate.edu"

input_file = '/work/idoerg/ahphan/democafa_package/data/cafa6/raw/GO_annotations_CAFA_20230610.tsv'
ontology_graph = '/work/idoerg/ahphan/democafa_package/data/cafa6/raw/go-basic-20250601.obo'
# ontology_slim = '/work/idoerg/ahphan/robot_obo/results/go-slim-20250601.obo'
taxon_path = '/work/idoerg/ahphan/democafa_package/data/raw/testsuperset-taxon-list.tsv'
output_file = '/work/idoerg/ahphan/democafa_package/data/cafa6/processed/terms_20250610.tsv'
swissprot_fasta = '/work/idoerg/ahphan/democafa_package/data/cafa6/raw/uniprot_sprot.fasta.2025.03.gz'
holdout_groundtruth = '/work/idoerg/ahphan/democafa_package/data/cafa6/processed/leaf_20250610_targets.tsv'

def process_tsv(input_file, ontology_graph, output_file, testsuperset):
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
    
    # Union-ize isoform annotations
    annotation_df['EntryID'] = annotation_df['EntryID'].str.split('-').str[0]
    annotation_df = annotation_df.drop_duplicates(subset=['EntryID', 'term', 'aspect'])
    
    print(annotation_df.describe())
    
    # Remove 3 binding terms of proteins even though they are annotated in other aspects
    # annotation_df = propagate_terms(annotation_df, subontologies) 
    binding_terms = set(nx.descendants(subontologies['F'],'GO:0005515'))
    binding_terms.add('GO:0005515')
    no_binding = annotation_df[~annotation_df['term'].isin(binding_terms)]
    binding_only = set.difference(set(annotation_df['EntryID']), set(no_binding['EntryID']))
    print(f"Removed {len(binding_only)} proteins with only protein-binding terms.")
    annotation_df = annotation_df[~annotation_df['EntryID'].isin(binding_only)]
    
    # Remove TrEMBL proteins
    fasta_proteins = read_fasta_proteins(swissprot_fasta)
    trembl_gain_proteins = set.difference(set(annotation_df['EntryID']), set(fasta_proteins.keys()))
    print(f"Removed {len(trembl_gain_proteins)} TrEMBL proteins: {trembl_gain_proteins}.")
    annotation_df = annotation_df[~annotation_df['EntryID'].isin(trembl_gain_proteins)]
    
    # Save to 3-column TSV
    annotation_df[['EntryID', 'term', 'aspect']].to_csv(output_file, sep='\t', index=False, header=True)
    print(annotation_df['aspect'].value_counts())
    # # Intersect with GOslim
    # ontology_slim = clean_ontology_edges(obonet.read_obo(ontology_slim))
    # slim_df = annotation_df[(annotation_df['term'].isin(ontology_slim.nodes())) & 
    #                         (~annotation_df['term'].isin(roots.values()))] # proteins with only roots (after slim) will be removed
    # print(slim_df.describe())   
    
    
    batch = list(annotation_df['EntryID'].unique())
    return batch
    
    
def fetch_taxonomy(batch):
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
        time.sleep(0.5)
    
    return taxon


def read_fasta_proteins(fasta_file: str):
    """
    Read proteins from FASTA file using Biopython.
    Filter for taxonomy (for demo purposes, only human proteins are considered).
    
    Args:
        fasta_file: Path to the FASTA file
        
    Returns:
        Dictionary mapping protein IDs to their sequences
    """
    # Read in taxon file
    # if taxon is None:
    #     selected_taxon = None
    # elif taxon.isdigit(): 
    #     # If taxon is a single ID, convert it to a list
    #     selected_taxon = [str(taxon)]
    # elif isinstance(taxon, str) and os.path.exists(taxon):
    #     taxon_file = pd.read_csv(taxon, header=0, sep='\t', encoding='ISO-8859-1')
    #     selected_taxon = [str(taxon_id) for taxon_id in set(taxon_file.iloc[:,0].tolist())]
    
    tax_pattern = re.compile(r"OX=(\d+)")
    species_pattern = re.compile(r"OS=([^O]+)")
    fasta_proteins = {}
    missing_taxon = {}
    # species = {}
    
    seq_count = 0
    
    # Process gzipped fasta file to get all SwissProt proteins
    with gzip.open(fasta_file, "rt") as gz_file:
        for record in SeqIO.parse(gz_file, 'fasta'):
            # Extract accession (EntryID)
            entry_id = record.id.split("|")[1] if "|" in record.id else record.id
            
            # Extract taxonomy ID using regex
            tax_match = tax_pattern.search(record.description)
            tax_id = tax_match.group(1) if tax_match else "N/A"
                
            # species_match = species_pattern.search(record.description)
            # species_name = species_match.group(1) if species_match else "N/A"
            # if selected_taxon is None or tax_id in selected_taxon:
                # seq_count += 1
            fasta_proteins[entry_id] = tax_id
            # elif tax_id not in selected_taxon:
                # missing_taxon[entry_id] = tax_id
            # if species_name != "N/A":
                # species[species_name] = tax_id
    
    return fasta_proteins

# gained_entries = process_tsv(input_file, ontology_graph, output_file)
# groundtruth_df = pd.read_csv('/work/idoerg/ahphan/democafa_package/data/temp/groundtruth_targets.tsv', sep='\t', header=None)


# Comparing with 90-species set
taxon_file =  pd.read_csv(taxon_path, header=0, sep='\t', encoding='ISO-8859-1')
selected_taxon = [str(taxon_id) for taxon_id in set(taxon_file.iloc[:,0].tolist())]
fasta_proteins = read_fasta_proteins('/work/idoerg/ahphan/democafa_package/data/cafa6/raw/uniprot_sprot_taxid_2759.fasta.gz')
# len(set.difference(set(selected_taxon), set(fasta_proteins.values())))

# Getting taxonomy names from taxon_id, compare with CAFA5 and write to file
taxon_new_file = '/work/idoerg/ahphan/democafa_package/data/cafa6/raw/selected-taxon-list.tsv'
taxon_new_set = set()
print("Fetching scientific names for proteins that are not in the Eukaryotes but in 90-species")
with open(taxon_new_file,'w') as f:
    f.write("ID\tSpecies\n")
    f.write("2759\tEukaryota\n") # add all eukaryotes
    taxon_new_set.add('2759')
    for taxid in set.difference(set(selected_taxon), set(fasta_proteins.values())):
        s = Entrez.esummary(db="taxonomy", id=taxid, retmode="xml")
        result = Entrez.read(s)
        print(taxid, result[0]['ScientificName'])
        f.write(f"{taxid}\t{result[0]['ScientificName']}\n")
        taxon_new_set.add(taxid)

# Checking hold out ground truth for any additional taxon
groundtruth_df = pd.read_csv(holdout_groundtruth, sep='\t', header=None)
gained_entries = groundtruth_df.iloc[:,0].tolist()
gained_taxon = fetch_taxonomy(gained_entries)
gained_taxon_count = Counter(gained_taxon.values())
gained_taxon_not_eukaryotes = {(name,id):count for (name,id),count in gained_taxon_count.items() if id not in set(fasta_proteins.values())}
with open(taxon_new_file,'a') as f:
    for (name,taxid), count in gained_taxon_not_eukaryotes.items():
        s = Entrez.esummary(db="taxonomy", id=taxid, retmode="xml")
        result = Entrez.read(s)
        print(taxid, result[0]['ScientificName'])
        if taxid not in taxon_new_set:
            f.write(f"{taxid}\t{result[0]['ScientificName']}\n")
            taxon_new_set.add(taxid)

# For getting most gained taxa, use groundtruth of 2-year rather than hold out ground truth
groundtruth_df = pd.read_csv('/work/idoerg/ahphan/democafa_package/data/temp/groundtruth_targets.tsv', sep='\t', header=None)
gained_entries = groundtruth_df.iloc[:,0].tolist()
gained_taxon = fetch_taxonomy(gained_entries)
gained_taxon_count = Counter(gained_taxon.values())
gained_taxon_not_eukaryotes = {(name,id):count for (name,id),count in gained_taxon_count.items() if id not in set(fasta_proteins.values())}
new_taxon = {name: id for (name, id) in gained_taxon.values() if id not in selected_taxon}
gained_taxon_not_eukaryotes = {(name,id):count for (name,id),count in gained_taxon_count.items() if id not in set(fasta_proteins.values())}
most_common_gained_taxon = Counter(gained_taxon_not_eukaryotes).most_common(20)
df = []
for (name, id), count in gained_taxon_not_eukaryotes.items():
    df.append({'count': count, 'taxon_id': id, 'species_name': name})
df = pd.DataFrame(df)
df.sort_values(by='count', ascending=False, inplace=True)
df.to_csv('/work/idoerg/ahphan/democafa_package/data/temp/gained_taxon_not_eukaryotes.tsv', sep='\t', index=False)