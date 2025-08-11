#!/bin/bash


#SBATCH --time=05:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=8   # 8 processor core(s) per node 
#SBATCH --mem=128G   # maximum memory per node
#SBATCH --job-name="groundtruth"
#SBATCH --mail-user=ahphan@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
# TODO: Find a way to store scripts and output data paths
module load micromamba
module load gcc/14.2.0-cuda12-vx6uhdf # for dask module
eval "$(micromamba shell hook --shell=bash)"
micromamba activate democafaenv
cd /work/idoerg/ahphan/democafa_package

# Read arguments
while getopts "f:g:h:i:q:o:" opt; do
  case $opt in
    f) gaf1="$OPTARG"
    ;;
    g) gaf2="$OPTARG"
    ;;
    h) obo1="$OPTARG"
    ;;
    i) obo2="$OPTARG"
    ;;
    q) query="$OPTARG"
    ;;
    o) output="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&6
        exit 1
    ;;
  esac
done

gaf_filtered1=${gaf1/.gz/_filtered.gz}
gaf_filtered2=${gaf2/.gz/_filtered.gz}

# Filter full gaf files for swissprot proteins
python3 -m democafa.datacollection.filter_gaf -a "$gaf1" -q "$query" -o "$gaf_filtered1"
python3 -m democafa.datacollection.filter_gaf -a "$gaf2" -q "$query" -o "$gaf_filtered2"

# Retrieve leaf terms from two gaf.gz files
python3 -m democafa.datacollection.retrieve_terms -a "$gaf_filtered1" -sgc "Experimental,IC,TAS" -g "$obo1" --tsv "data/processed/terms1.tsv" 
python3 -m democafa.datacollection.retrieve_terms -a "$gaf_filtered2" -sgc "Experimental,IC,TAS" -g "$obo1" -ag "$obo2" --tsv "data/processed/terms2.tsv" 

# Propagate terms to the root
python3 -m democafa.datacollection.propagate_and_ia -t "data/processed/terms1.tsv" -g "$obo1" -tp "data/processed/terms1_propagated.tsv" -ot "data/processed/terms1_IA.tsv"
python3 -m democafa.datacollection.propagate_and_ia -t "data/processed/terms2.tsv" -g "$obo2" -tp "data/processed/terms2_propagated.tsv"

# Classify ground truth based on two propagated term files
python3 -m democafa.groundtruth.classify_ground_truth -ak "data/processed/terms1_propagated.tsv" -a2 "data/processed/terms2_propagated.tsv" -g "$obo1" -g2 "$obo2" -q "$query" -o "$output"
# TODO: add a function in classify_ground_truth to get gain queries (in fasta format)
# TODO: need train_sequences and train_taxonomy for blast predictor
gain_query="${query/.fasta/_gain.fasta}"

# Run evaluations
cd /work/idoerg/ahphan/CAFA_forever
mkdir -p current/
mkdir -p current/predictions

module load singularity
singularity exec --pwd /app --bind /work/idoerg/ahphan/democafa_package/data:/app/data --bind /work/idoerg/ahphan/CAFA_forever/current/predictions:/app/output \
  test_blast_latest.sif python3 blast_main.py \
  --annot_file /app/data/processed/terms1_propagated.tsv --query_file "$gain_query" --graph "$obo1" \
  --train_sequences /app/data/processed/train_sequences.2025.02.fasta --train_taxonomy /app/data/processed/train_taxonomy.2025.02.tsv \
  --output_baseline /app/output/blast_predictions.tsv.gz
singularity exec --pwd /app   --bind /work/idoerg/ahphan/democafa_package/data:/app/data   --bind /work/idoerg/ahphan/CAFA-forever/AprJun/predictions:/app/output   test_naive_latest.sif   python3 naive.py --annot_file /app/data/processed/cafa6/goa_uniprot_filtered_mp.gaf.225.gz --query_file /work/idoerg/ahphan/CAFA-forever/AprJun/targets2.txt --graph /app/data/raw/go-basic-20250316.obo --output_baseline /app/output/naive_predictions.tsv.gz 
singularity exec --pwd /app   --bind /work/idoerg/ahphan/democafa_package/data:/app/data   --bind /work/idoerg/ahphan/CAFA-forever/AprJun/predictions:/app/output   test_goa_nonexp_latest.sif   python3 goa_nonexp.py --annot_file /app/data/processed/cafa6/goa_uniprot_filtered_mp.gaf.225.gz --query_file /work/idoerg/ahphan/CAFA-forever/AprJun/targets2.txt --graph /app/data/raw/go-basic-20250316.obo --output_baseline /app/output/goa_nonexp_predictions.tsv.gz --selected_go 'Computational,Phylogenetical,Electronic,ND,NAS'