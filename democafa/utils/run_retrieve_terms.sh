#!/bin/bash

#SBATCH --time=03:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 8 processor core(s) per node 
#SBATCH --mem=128G   # maximum memory per node
#SBATCH --job-name="retrieve_terms"
#SBATCH --mail-user=ahphan@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load micromamba
module load gcc/14.2.0-cuda12-vx6uhdf # for dask module
eval "$(micromamba shell hook --shell=bash)"
micromamba activate democafaenv
cd /work/idoerg/ahphan/democafa_package

echo "Running filter taxonomy"
python3 -m democafa.datacollection.filter_gaf -a data/raw/goa_uniprot_all.gaf.226.gz -q data/raw/uniprot_sprot.fasta.2025.03.gz -o data/processed/cafa6/goa_uniprot_filtered_mp.gaf.226.gz
echo "Running retrieve terms"
python3 -m democafa.datacollection.retrieve_terms --annot data/processed/cafa6/goa_uniprot_filtered_mp.gaf.226.gz -sgc 'Experimental,IC,TAS' -g data/raw/go-basic.obo --tsv data/processed/cafa6/train_terms.tsv 

