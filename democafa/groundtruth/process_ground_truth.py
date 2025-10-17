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
import pandas as pd
import networkx as nx
import time
from collections import Counter
from democafa.utils.ontology import add_aspect_column, clean_ontology_edges, fetch_aspect
from democafa.datacollection.retrieve_terms import filter_evidence_codes
from democafa.config import GO_CODES
from Bio import Entrez
Entrez.email = "ahphan@iastate.edu"


def process_tsv(input_file, ontology_graph):
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
        print(f"Warning: {len(obsolete_terms)} obsolete (or new) terms ({obsolete_terms}) found in annotation file but not in graph.")
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
    # no_binding = annotation_df[(~annotation_df['term'].isin(binding_terms)) & (annotation_df['aspect']=='F')]
    # binding_only = set.difference(set(annotation_df['EntryID']), set(no_binding['EntryID']))
    # print(f"Removed {len(binding_only)} proteins with only protein-binding terms.")
    # annotation_df = annotation_df[~annotation_df['EntryID'].isin(binding_only)]
    
    annotated_F_proteins = set(annotation_df[annotation_df['aspect']=='F']['EntryID'])
    not_binding_only_proteins = set(annotation_df[(annotation_df['EntryID'].isin(annotated_F_proteins)) & 
                                                  (~annotation_df['term'].isin(binding_terms)) & 
                                                  (annotation_df['aspect']=='F')]['EntryID']) # these proteins have other MFO annotations
    binding_df = annotation_df[(annotation_df['EntryID'].isin(annotated_F_proteins - not_binding_only_proteins)) & 
                               (annotation_df['aspect']=='F')] # these annotations are protein-binding only
    # Remove rows in binding_df from annotation_df
    annotation_df = annotation_df[~annotation_df.index.isin(binding_df.index)]
                               
    # # Remove TrEMBL proteins
    # fasta_proteins = read_fasta_proteins(swissprot_fasta)
    # trembl_gain_proteins = set.difference(set(annotation_df['EntryID']), set(fasta_proteins.keys()))
    # print(f"Removed {len(trembl_gain_proteins)} TrEMBL proteins: {trembl_gain_proteins}.")
    # annotation_df = annotation_df[~annotation_df['EntryID'].isin(trembl_gain_proteins)]
    
    # Save to 3-column TSV
    
    print(annotation_df['aspect'].value_counts())
    print(f"Number of unique proteins: {len(annotation_df['EntryID'].unique())}")
    # # Intersect with GOslim
    # ontology_slim = clean_ontology_edges(obonet.read_obo(ontology_slim))
    # slim_df = annotation_df[(annotation_df['term'].isin(ontology_slim.nodes())) & 
    #                         (~annotation_df['term'].isin(roots.values()))] # proteins with only roots (after slim) will be removed
    # print(slim_df.describe())   
        
    batch = list(annotation_df['EntryID'].unique())
    return annotation_df
    
    
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


def read_fasta_proteins(fasta_file: str, entries):
    """
    Read proteins from FASTA file using Biopython.
    Filter for taxonomy (for demo purposes, only human proteins are considered).
    
    Args:
        fasta_file: Path to the FASTA file
        
    Returns:
        Dictionary mapping protein IDs to their sequences
    """
    tax_pattern = re.compile(r"OX=(\d+)")
    species_pattern = re.compile(r"OS=([^O]+)")
    fasta_proteins = {}
    # species = {}
    
    
    # Process gzipped fasta file to get all SwissProt proteins
    with gzip.open(fasta_file, "rt") as gz_file:
        for record in SeqIO.parse(gz_file, 'fasta'):
            # Extract accession (EntryID)
            entry_id = record.id.split("|")[1] if "|" in record.id else record.id
            if entry_id not in entries:
                continue
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


def fetch_eukaryote_taxa():
    """Fetch eukaryote taxonomy data from UniProt"""
    url = "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Corganism_name%2Corganism_id&format=tsv&query=%28%28taxonomy_id%3A2759%29+AND+%28reviewed%3Atrue%29%29"
    
    response = requests.get(url)
    response.raise_for_status()
    
    # Decompress the content
    decompressed_data = gzip.decompress(response.content).decode('utf-8')
    
    # Parse into lines
    lines = decompressed_data.strip().split('\n')
    
    # Extract taxonomy IDs
    eukaryote_taxa_ids = set()
    for line in lines[1:]:  # Skip header
        if line.strip():
            fields = line.split('\t')
            if len(fields) >= 3:
                tax_id = fields[2]  # Organism ID column
                eukaryote_taxa_ids.add(tax_id)
    
    return eukaryote_taxa_ids


def add_additional_taxon(old_taxon_path, new_taxon_path, holdout_groundtruth, fasta_file):
    """
    Compare old taxon list with new taxon list and add any additional taxon found in holdout ground truth.
    
    Args:
        old_taxon_path: Path to the old taxon list file (CAFA5).
        new_taxon_path: Path to the new taxon list file (CAFA6).
        holdout_groundtruth: Path to the holdout ground truth file.
    """
    groundtruth_df = pd.read_csv(holdout_groundtruth, sep='\t', header=0)
    gained_entries = set(groundtruth_df.iloc[:,0].tolist())
    print(f"Holdout ground truth has {len(gained_entries)} entries.")
    # gained_taxon = fetch_taxonomy(gained_entries) # TODO: fetch_taxonomy in batches is not optimized
    gained_protein_taxon = read_fasta_proteins(fasta_file, gained_entries)
    gained_taxon_counts = {id:count for id,count in Counter(gained_protein_taxon.values()).items()}
    print(gained_taxon_counts)
    gained_taxon_ids = set(gained_taxon_counts.keys())
    
    # species of interest in CAFA5: 90 species
    taxon_file =  pd.read_csv(old_taxon_path, header=0, sep='\t', encoding='ISO-8859-1')
    old_taxon = [str(taxon_id) for taxon_id in set(taxon_file.iloc[:,0].tolist())]
    print(f"Old taxon list has {len(old_taxon)} entries.")
    
    # CAFA6: all Eukaryotes (id=2759)
    taxon_new_set = set()
    eukaryote_taxa = fetch_eukaryote_taxa()
    
    old_prokaryote_taxa = set(old_taxon) - eukaryote_taxa
    if old_prokaryote_taxa:
        print(f"Adding back {len(old_prokaryote_taxa)} prokaryote taxa from old taxon list.")
        taxon_new_set.update(old_prokaryote_taxa)
    
    new_prokaryote_taxa = (set(gained_taxon_ids) - eukaryote_taxa) - old_prokaryote_taxa
    if new_prokaryote_taxa:
        print(f"Adding {len(new_prokaryote_taxa)} prokaryote taxa from holdout ground truth.")
        taxon_new_set.update(new_prokaryote_taxa)
    
    # write to new taxon file
    with open(new_taxon_path,'w') as f:
        f.write("ID\tSpecies\n")
        f.write("2759\tEukaryota\n") # write Eukaryota first
        for taxid in taxon_new_set:
            s = Entrez.esummary(db="taxonomy", id=taxid, retmode="xml")
            result = Entrez.read(s)
            print(taxid, result[0]['ScientificName'])
            f.write(f"{taxid}\t{result[0]['ScientificName']}\n")
            taxon_new_set.add(taxid)
    print(f"New taxon list has {len(taxon_new_set)} entries.")


def parse_args(argv):
    parser = argparse.ArgumentParser(description="Process ground truth data for CAFA challenge.")
    parser.add_argument('--input_file', type=str, required=True, help='Path to the input TSV file with GO annotations.')
    parser.add_argument('--graph', type=str, required=True, help='Path to the ontology OBO file.')
    parser.add_argument('--output_terms', type=str, required=True, help='Path to the output TSV file with processed terms.')
    parser.add_argument('--fasta_file', type=str, required=True, help='Path to the FASTA file with protein sequences.')
    
    return parser.parse_args(argv)


def process_ground_truth(input_file, graph, output_terms):
    annotation_df = process_tsv(input_file, graph)
    annotation_df[['EntryID', 'term', 'aspect']].to_csv(output_terms, sep='\t', index=False, header=True)


def main():
    args = parse_args(sys.argv[1:])
    process_ground_truth(args.input_file, args.graph, args.output_terms)
    # add_additional_taxon('data/raw/testsuperset-taxon-list.tsv', 
    #                      'data/cafa6/raw/selected-taxon-list.tsv', 
    #                      args.output_terms,
    #                      args.fasta_file)
    
if __name__ == "__main__":
    main()