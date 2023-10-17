#!/bin/bash
set -x

for i in {1..10}
do
    nextflow run main.nf \
        --input data/new_samples.csv \
        --num_samples 10 \
        -profile docker \
        -entry fetch_sra &
    BACK_PID=$!
    wait $BACK_PID
    rm -rf work

    # Remove the first 10 lines from the input file
    sed -i '1,10d' data/new_samples.csv
done
