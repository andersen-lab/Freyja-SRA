#!/bin/bash
#SBATCH --job-name=freyja-sra
#Resource request
#SBATCH --time=24:00:00
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

BATCH_SIZE=100
python scripts/get_accession_list_wendy.py $BATCH_SIZE
nextflow run main.nf -entry sra -profile singularity --accession_list data/accession_list.csv --num_samples $BATCH_SIZE
rm -rf work
python scripts/update_sample_status.py $BATCH_SIZE