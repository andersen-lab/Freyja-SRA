import pandas as pd
import os

metadata = pd.read_csv('data/all_metadata.csv')
original_sample_status = metadata['sample_status']

# Add a new column to the metadata dataframe called 'sample_status' and set it to 'to_run'
metadata['sample_status'] = 'to_run'

for file in os.listdir('outputs/variants'):
    if file.endswith('variants.tsv'):
        sample_id = file.split('.')[0]
        metadata.loc[metadata['accession'] == sample_id, 'sample_status'] = 'demix_error'

for file in os.listdir('outputs/demix'):
    if file.endswith('demix.tsv'):
        sample_id = file.split('.')[0]
        metadata.loc[metadata['accession'] == sample_id, 'sample_status'] = 'completed'

# If the original sample status was 'fastq_error', set the new sample status to 'fastq_error'
metadata['sample_status'] = metadata.apply(lambda x: 'fastq_error' if x['sample_status'] == 'fastq_error' else x['sample_status'], axis=1)

print(metadata['sample_status'].value_counts())
metadata.to_csv('data/all_metadata.csv', index=False)