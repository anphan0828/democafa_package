#!/bin/bash


#SBATCH --time=05:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --cpus-per-task=16 
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
# python3 -m democafa.datacollection.propagate_and_ia -t "data/processed/terms2.tsv" -g "$obo1" -tp "data/processed/terms2_propagated.tsv" # no need to propagate

# Classify ground truth based on two propagated term files
python3 -m democafa.groundtruth.classify_ground_truth -ak "data/processed/terms1_propagated.tsv" -a2 "data/processed/terms2.tsv" -g "$obo1" -g2 "$obo2" -q "$query" --out_prefix "$output"
# TODO: need train_sequences and train_taxonomy for blast predictor
gain_query="${output/.tsv/_targets.tsv}"
gain_query_fasta="${gain_query/.tsv/.fasta}"

python3 -m democafa.datacollection.retrieve_sequences -i "$gain_query" -f "$query" -of "$gain_query_fasta"
# # Run evaluations
cd /work/idoerg/ahphan/CAFA_forever
mkdir -p current/
mkdir -p current/predictions

module load singularity
singularity exec --pwd /app --bind /work/idoerg/ahphan/democafa_package/data:/app/data --bind /work/idoerg/ahphan/CAFA_forever/current/predictions:/app/output \
  test_blast_latest.sif python3 blast_main.py \
  --annot_file /app/data/processed/terms1_propagated.tsv --query_file "$gain_query" --graph "$obo1" \
  --train_sequences /app/data/processed/train_sequences.2025.02.fasta --train_taxonomy /app/data/processed/train_taxonomy.2025.02.tsv \
  --output_baseline /app/output/blast_predictions.tsv.gz
singularity exec --pwd /app --bind /work/idoerg/ahphan/democafa_package/data:/app/data --bind /work/idoerg/ahphan/CAFA-forever/current/predictions:/app/output \
  test_naive_latest.sif python3 naive.py \
  --annot_file /app/data/processed/terms1_propagated.tsv --query_file "$gain_query" --graph "$obo1" \
  --output_baseline /app/output/naive_predictions.tsv.gz 
singularity exec --pwd /app --bind /work/idoerg/ahphan/democafa_package/data:/app/data --bind /work/idoerg/ahphan/CAFA-forever/current/predictions:/app/output \
  test_goa_nonexp_latest.sif python3 goa_nonexp.py \
  --annot_file /app/data/processed/terms1_propagated.tsv --query_file "$gain_query" --graph "$obo1" \
  --selected_go 'Computational,Phylogenetical,Electronic,ND,NAS' \
  --output_baseline /app/output/goa_nonexp_predictions.tsv.gz 
