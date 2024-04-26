#!/bin/bash
#SBATCH --job-name=freyja-sra
#Resource request
#SBATCH --time=24:00:00
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

python scripts/get_accession_list_wendy.py 500
nextflow run main.nf -entry sra -profile singularity --accession_list data/accession_list.csv --num_samples 500
rm -rf work
python scripts/update_sample_status.py 500