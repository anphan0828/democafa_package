#!/usr/bin/env python3
"""Filter a GAF file by taxonomy, optinally using multiprocessing."""

import os
import sys
import argparse
import gzip
import logging
import obonet
import pandas as pd
import networkx as nx
import numpy as np
from Bio.UniProt import GOA
from Bio import SwissProt as sp
from Bio import SeqIO
from democafa.utils.ontology import clean_ontology_edges, fetch_aspect, propagate_terms, filter_terms_given_obo
from democafa.config import GO_CODES

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def process_chunk_file(chunk_file, taxon, entries):
    chunk_data = []
    processed_count = 0
    with open(chunk_file, 'r') as f:
        for rec in GOA.gafiterator(f):
            processed_count += 1
            if entries is not None:
                if rec['DB_Object_ID'][0] not in entries:
                    continue
            if taxon is not None:
                if rec['Taxon_ID'][0] not in taxon:
                    continue
            chunk_data.append(rec)
    
    logger.debug(f"Chunk {os.path.basename(chunk_file)}: processed {processed_count} records, kept {len(chunk_data)} records")
    return chunk_data
    

def get_entries_from_file(query_file):
    """
    Extract entries from a FASTA file or a text file.
    Returns a set of entry IDs.
    """
    query_ids = set()
    if 'fasta' in query_file:
        is_gzipped = query_file.endswith('.gz')
        open_func = gzip.open if is_gzipped else open
        mode = 'rt' if is_gzipped else 'r'
        logger.info(f"Reading query IDs from FASTA file: {query_file}")
        with open_func(query_file, mode) as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                entry_id = record.id.split("|")[1] if "|" in record.id else record.id
                query_ids.add(entry_id)
    elif query_file.endswith('.txt'):
        logger.info(f"Reading query IDs from text file: {query_file}")
        with open(query_file, 'r') as handle:
            query_ids = {line.strip() for line in handle}
    else:
        logger.error("Please provide a fasta file or a text file with query IDs")
        sys.exit(1)
    logger.info(f"Extracted {len(query_ids)} entries from input query file")
    
    return query_ids
    

def get_taxon_from_file(taxon):
    if taxon is None:
        selected_taxon = None
    elif taxon.isdigit(): 
        # If taxon is a single ID, convert it to a list
        selected_taxon = [f'taxon:{taxon}']
    elif isinstance(taxon, str) and os.path.exists(taxon):
        logger.info(f"Reading taxon IDs from file: {taxon}")
        taxon_file = pd.read_csv(taxon, header=0, sep='\t', encoding='ISO-8859-1')
        selected_taxon = [f'taxon:{taxon_id}' for taxon_id in set(taxon_file.iloc[:,0].tolist())]
    logger.info(f"Selected {len(selected_taxon) if selected_taxon else 0} taxa: {selected_taxon[:5] if selected_taxon and len(selected_taxon) > 5 else selected_taxon}")
    return selected_taxon
    
        
def filter_gaf(file_path, taxon, entries, output, use_mp=True, num_processes=None, chunk_size=10000000):
    """
    Read and process a GAF file (gzipped or plain text) with multiprocessing.
    Takes 30m for 22Gb uniprot goa for 15 processes.
    """
    logger.info("Starting GAF filtering process")
    logger.info(f"Input file: {file_path}")
    logger.info(f"Output file: {output}")
    logger.info(f"Using multiprocessing: {use_mp}")
    
    is_gzipped = file_path.endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    handle = open_func(file_path, mode)
    mode_out = 'wt' if output.endswith('.gz') else 'w'
    
    if not use_mp: 
        logger.info("Processing file without multiprocessing")
        with open_func(output, mode_out) as out_handle:
            for rec in GOA.gafiterator(handle):
                if entries is not None:
                    if rec['DB_Object_ID'][0] not in entries:
                        continue
                if taxon is not None:
                    if rec['Taxon_ID'][0] not in taxon:
                        continue
                GOA.writerec(rec,out_handle)
    else:
        import multiprocessing
        from functools import partial
        import tempfile
        if num_processes is None:
            num_processes = min(int(os.environ.get('SLURM_CPUS_PER_TASK', multiprocessing.cpu_count())) - 1, 16)
            num_processes = max(1, num_processes)
        logger.info(f"Using {num_processes} processes for parallel computation")
        
        temp_dir = tempfile.mkdtemp(dir=os.path.dirname(file_path))
        logger.debug(f"Created temporary directory: {temp_dir}")
        chunk_files = []
        try:
            # Read and chunk the file
            logger.info(f"Splitting GAF file into chunks of {chunk_size} records...")
            chunk_count = 0
            record_count = 0
            
            with open_func(file_path, mode) as handle:
                current_chunk_lines = []
                headers = []
                line = handle.readline()
                while line.startswith('!'):
                    headers.append(line)
                    line = handle.readline()
                
                logger.debug(f"Found {len(headers)} header lines")
                
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
                            logger.info(f"Created chunk {chunk_count} with {chunk_size} records")
                            
                # Write any remaining lines to the final chunk
                if current_chunk_lines:
                    chunk_file = os.path.join(temp_dir, f"chunk_{chunk_count}.gaf")
                    with open(chunk_file, 'w') as f:
                        # Write headers to each chunk
                        f.writelines(headers)
                        f.writelines(current_chunk_lines)
                    
                    chunk_files.append(chunk_file)
                    chunk_count += 1
                    logger.info(f"Created final chunk {chunk_count} with {len(current_chunk_lines)} records")
                    
            # Process chunks in parallel
            logger.info(f"Processing {len(chunk_files)} chunks with {num_processes} processes")
            process_func = partial(process_chunk_file, taxon=taxon, entries=entries)
            with multiprocessing.Pool(processes=num_processes) as pool:
                results = pool.map(process_func, chunk_files)
            
            all_data = [item for sublist in results for item in sublist]
            logger.info(f"Filtered data contains {len(all_data)} records")

            logger.info(f"Writing filtered data to output file: {output}")
            with open_func(output, mode_out) as out_handle:
                out_handle.writelines(headers)  # Write headers to the output file
                for rec in all_data:
                    GOA.writerec(rec, out_handle)
                    
            logger.info("Successfully completed GAF filtering")
            
        except Exception as e:
            logger.error(f"Error during GAF filtering: {str(e)}")
            raise
        finally:
            # Cleanup temporary files
            logger.debug("Cleaning up temporary files")
            for chunk_file in chunk_files:
                if os.path.exists(chunk_file):
                    os.remove(chunk_file)
            if os.path.exists(temp_dir):
                os.rmdir(temp_dir)
                

def parse_inputs(argv):
    parser = argparse.ArgumentParser(
        description='Filter a GAF file and write a new GAF gzipped file with selected taxon')
    
    parser.add_argument('--annot', '-a', required=True,
                        help='Path to first annotation file (can be gzipped)')
    parser.add_argument('--taxon', '-t', required=False,
                        help='Taxon ID file to filter proteins (default: None)')
    parser.add_argument('--query', '-q', required=False,
                        help='Query file with protein IDs (or FASTA file, can be gzipped) to filter annotations (default: None)')
    parser.add_argument('--output', '-o', default='data/processed/goa_uniprot_filtered.gaf.gz',
                        help='Path to save annotations in TSV format')
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], 
                        default='INFO', help='Set the logging level (default: INFO)')
    
    args = parser.parse_args(argv)
    
    # Configure logging level based on argument
    logging.getLogger().setLevel(getattr(logging, args.log_level))
    
    return args

    
def main():
    args = parse_inputs(sys.argv[1:])
    
    logger.info("Starting GAF filtering script")
    logger.info(f"Arguments: {vars(args)}")
    
    entries = None
    taxon = None
    
    if args.query is not None:
        entries = get_entries_from_file(args.query)
    if args.taxon is not None:
        taxon = get_taxon_from_file(args.taxon)
        
    filter_gaf(args.annot, taxon, entries, args.output, use_mp=True, num_processes=None)
    
    logger.info("GAF filtering script completed successfully")
    
    # python3 -m democafa.datacollection.filter_taxonomy -a data/raw/goa_uniprot_all.gaf.216.gz -t data/raw/testsuperset-taxon-list.tsv -o data/processed/cafa5/goa_uniprot_filtered_mp.gaf.216.gz
    # python3 -m democafa.datacollection.filter_taxonomy -a data/raw/goa_uniprot_all.gaf.224.gz -t data/raw/testsuperset-taxon-list.tsv -o data/processed/cafa5/goa_uniprot_filtered_mp.gaf.224.gz
    # python3 -m democafa.datacollection.filter_taxonomy -a data/raw/goa_uniprot_all.gaf.226.gz -q data/raw/uniprot_sprot.fasta.gz -o data/processed/cafa5/goa_uniprot_filtered_mp.gaf.226.gz
    # TODO: test this script with the new filtering and logging setup
if __name__ == "__main__":
    main()