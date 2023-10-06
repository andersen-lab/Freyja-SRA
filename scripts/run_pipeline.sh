#!/bin/bash
set -x
nextflow run main.nf \
    --input data/new_samples.csv \
    --num_samples 5 \
    -profile singularity \
    -entry fetch_sra

# Clean up
rm -rf work