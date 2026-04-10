#!/usr/bin/env bash

# Get first 1000 lines and filter not IEA from uniprot_goa.gaf
zcat '/work/idoerg/ahphan/democafa_package/data/cafa6/processed/goa_uniprot_selected.gaf.gz' | head -n 1000 | awk '$7 !~ /IEA/' > tests/test_data/goa_uniprot_selected.gaf

# Get testsuperset.fasta
python3 -m democafa.datacollection.create_test_set -t tests/test_data/train_terms_selected.tsv -f data/cafa6/raw/uniprot_sprot.fasta.2025.03.gz -o tests/test_data/test_superset_all.2025.03.fasta --train_out_fasta tests/test_data/train_sequences.2025.03.fasta --train_out_taxonomy tests/test_data/train_taxonomy.2025.03.tsv --include_all -n 100

# Propagate train_terms and calculate IA
python3 -m democafa.datacollection.propagate_and_ia -t '/work/idoerg/ahphan/democafa_package/tests/test_data/train_terms_selected.tsv' -g data/cafa6/raw/go-basic-20250601.obo -tp tests/test_data/train_terms_selected_propagated.tsv --output_tsv tests/test_data/IA.txt

# Generate baseline predictions at t-1
# Create baseline prediction from propagated terms file is better (some ancestor terms get higher score compared to propagating with 'fill' option in cafaeval)
# python3 -m democafa.baselines.naive -a '/work/idoerg/ahphan/democafa_package/tests/test_data/train_terms_selected.tsv' --graph data/cafa6/raw/go-basic-20250601.obo -q '/work/idoerg/ahphan/democafa_package/tests/test_data/test_superset_all.2025.03.fasta' -o tests/test_predictions/naive_predictions_from_leaf.tsv
python3 -m democafa.baselines.naive -a '/work/idoerg/ahphan/democafa_package/tests/test_data/train_terms_selected_propagated.tsv' --graph data/cafa6/raw/go-basic-20250601.obo -q '/work/idoerg/ahphan/democafa_package/tests/test_data/test_superset_all.2025.03.fasta' -o tests/test_predictions/naive_predictions_from_tsv.tsv
python3 -m democafa.baselines.goa_nonexp --annot_file '/work/idoerg/ahphan/democafa_package/tests/test_data/goa_uniprot_selected.gaf' --query_file '/work/idoerg/ahphan/democafa_package/tests/test_data/test_superset_all.2025.03.fasta' --graph data/cafa6/raw/go-basic-20250601.obo --output_baseline tests/test_predictions/goa_nonexp_predictions.tsv --selected_go 'Computational,Phylogenetical,Electronic,ND,NAS'
cd democafa/baselines/blast_container
python3 blast_main.py --annot_file '/work/idoerg/ahphan/democafa_package/tests/test_data/train_terms_selected_propagated.tsv' --query_file '/work/idoerg/ahphan/democafa_package/tests/test_data/test_superset_all.2025.03.fasta' --graph /work/idoerg/ahphan/democafa_package/data/cafa6/raw/go-basic-20250601.obo --train_sequences '/work/idoerg/ahphan/democafa_package/tests/test_data/train_sequences.2025.03.fasta' --train_taxonomy '/work/idoerg/ahphan/democafa_package/tests/test_data/train_taxonomy.2025.03.tsv' --output_baseline '/work/idoerg/ahphan/democafa_package/tests/test_predictions/blast_predictions.tsv'
cd democafa/baselines/prott5_container
micromamba activate prott5-env
python3 prott5_main.py --annot_file '/work/idoerg/ahphan/democafa_package/tests/test_data/train_terms_selected_propagated.tsv' --query_file '/work/idoerg/ahphan/democafa_package/tests/test_data/test_superset_all.2025.03.fasta' --graph /work/idoerg/ahphan/democafa_package/data/cafa6/raw/go-basic-20250601.obo --train_sequences '/work/idoerg/ahphan/democafa_package/tests/test_data/train_sequences.2025.03.fasta' --output_baseline '/work/idoerg/ahphan/democafa_package/tests/test_predictions/prott5_predictions.tsv' --model_dir /work/idoerg/ahphan/.cache/huggingface