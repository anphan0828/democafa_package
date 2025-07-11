#!/usr/bin/env python3

import sys
import argparse
from democafa.utils.ontology import propagate_and_ia

def parse_args(args):
    parser = argparse.ArgumentParser(description='Propagate and Compute Information Accretion of GO annotations')
    parser.add_argument('--terms', '-t', required=True, 
                        help='Path to annotation file')
    
    parser.add_argument('--graph', '-g', default=None, 
                        help='Path to OBO ontology graph file if local. If empty (default) current OBO structure at run-time will be downloaded from http://purl.obolibrary.org/obo/go/go-basic.obo')
    
    parser.add_argument('--matrix_propagated', '-mp', default=None, 
                        help='Path to save propagated term counts')
    
    parser.add_argument('--matrix_indices', '-mi', default=None,
                        help='Path to save protein and term indices')
    
    parser.add_argument('--output_tsv', '-ot', default=None, 
                        help='Path to save computed IA for each term in the GO. If empty, will be saved to ./IA.txt')  
    
    return parser.parse_args(args)

    
def main():    
    args = parse_args(sys.argv[1:])
    propagate_and_ia(
        terms_file = args.terms,
        graph = args.graph,
        matrix_propagated = args.matrix_propagated,
        matrix_indices = args.matrix_indices,
        output_tsv = args.output_tsv
    )

        
if __name__ == "__main__":
    main()  