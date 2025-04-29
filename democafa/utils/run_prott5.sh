#!/bin/bash


#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=8   # 8 processor core(s) per node 
#SBATCH --mem=122G   # maximum memory per node
#SBATCH --gres=gpu:a100:1
#SBATCH --job-name="prott5-baseline"
#SBATCH --mail-user=ahphan@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
# TODO: Find a way to store scripts and output data paths
module load micromamba
module load gcc/14.2.0-cuda12-vx6uhdf # for dask module
eval "$(micromamba shell hook --shell=bash)"
micromamba activate prott5-env
cd /work/idoerg/ahphan/democafa_package/democafa

while getopts "q:d:m:o:" opt; do
  case $opt in
    q) query="$OPTARG"
    ;;
    d) database="$OPTARG"
    ;;
    m) model="$OPTARG"
    ;;
    o) output="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&4
        exit 1
    ;;
  esac
done

out_evalset="../data/processed/prott5/evalset_embeddings.h5"
out_dbset="../data/processed/prott5/blast_db_embeddings.h5"
echo "Running prott5 embeddings for evaluation set"
python utils/prott5-baseline/prott5_embedder.py --input $query --output $out_evalset --per_protein 1 --model $model

echo "Running prott5 embeddings for database set"
python utils/prott5-baseline/prott5_embedder.py --input $database --output $out_dbset --per_protein 1 --model $model

echo "Processing embeddings to $output"
python utils/prott5-baseline/process_embeddings_gpu.py $out_evalset $out_dbset $output

python utils/prott5-baseline/normalize_embeddings.py $output

# replace .tsv with _norm.tsv from the output file name
output_norm=${output/.tsv/_norm.tsv}
# Format clearning
sed -ri 's/sp\|//g' $output_norm
sed -ri 's/tr\|//g' $output_norm
sed -ri 's/\|[^\t]*\t/\t/g' $output_norm
