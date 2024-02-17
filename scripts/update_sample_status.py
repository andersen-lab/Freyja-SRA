import os
import sys
import pandas as pd

metadata = pd.read_csv('data/all_metadata.csv')
demix = os.listdir('outputs/demix')
variants = os.listdir('outputs/variants')
batch_size = int(sys.argv[1])

accessions_processed = []
with open('data/samples_to_run.csv', 'r') as f:
    accessions_processed = f.readlines()
    accessions_processed = [x.strip() for x in accessions_processed]
    accessions_processed = accessions_processed[:batch_size]

for accession in accessions_processed:
    if f'{accession}.variants.tsv' not in variants:
        metadata.loc[metadata['accession'] == accession, 'sample_status'] = 'fastq_error'
for file in demix:
    if file.endswith('demix.tsv'):
        sample_id = file.split('.')[0]
        metadata.loc[metadata['accession'] == sample_id, 'sample_status'] = 'completed'

for file in variants:
    if file.endswith('variants.tsv'):
        sample_id = file.split('.')[0]
        if f'{sample_id}.demix.tsv' not in demix:
            metadata.loc[metadata['accession'] == sample_id, 'sample_status'] = 'demix_error'
        
metadata.to_csv('data/all_metadata.csv', index=False)