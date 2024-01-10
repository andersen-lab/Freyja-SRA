#!/bin/bash
set -x

BATCH_SIZE=50


# for i in {1..10}
# do
nextflow run main.nf \
    --accession_list data/samples_to_run.csv \
    --num_samples $BATCH_SIZE \
    -profile docker \
    -entry sra &
BACK_PID=$!
wait $BACK_PID
rm -rf work

sed -i '2,51d' data/samples_to_run.csv
# done


