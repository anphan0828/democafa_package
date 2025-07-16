#!/bin/bash

#SBATCH --time=06:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 8 processor core(s) per node 
#SBATCH --mem=128G   # maximum memory per node
#SBATCH --job-name="blast"
#SBATCH --mail-user=ahphan@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load micromamba
module load gcc/14.2.0-cuda12-vx6uhdf # for dask module
eval "$(micromamba shell hook --shell=bash)"
micromamba activate democafaenv
cd /work/idoerg/ahphan/democafa_package


# Build local blast-able database with train_sequences.fasta and train_taxonomy.tsv (https://www.ncbi.nlm.nih.gov/books/NBK569841/)
makeblastdb -in data/processed/cafa6/train_sequences.fasta -parse_seqids -blastdb_version 5 -taxid_map data/processed/cafa6/train_taxonomy.tsv -title 'demo' -dbtype prot
# Run blastp (on 16 cores)
blastp -query data/processed/cafa6/test_superset_all.fasta -db data/processed/cafa6/train_sequences.fasta -outfmt "6 qseqid sseqid evalue length pident nident" -mt_mode 1 -num_threads 16 -out data/processed/cafa6/blast_results.tsv