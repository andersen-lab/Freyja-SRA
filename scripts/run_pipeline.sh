#!/bin/bash
set -x

BATCH_SIZE=3
for i in {1..10}
do
    nextflow run main.nf \
        --input data/new_samples.csv \
        --num_samples $BATCH_SIZE \
        -profile singularity \
        -entry fetch_sra &
    BACK_PID=$!
    wait $BACK_PID
    rm -rf work

    # Remove the first 3 data lines from the input file
    sed -i '2,4d' data/new_samples.csv
done
