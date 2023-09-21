import pandas as pd
import json

metadata = pd.read_csv('data/wastewater_ncbi.csv')
print('Total samples: ', len(metadata))
# Select only samples where 'amplicon_PCR_primer_scheme' is not empty
metadata = metadata[metadata['amplicon_PCR_primer_scheme'].notna()]
print('Samples with amplicon_PCR_primer_scheme: ', len(metadata))
current_samples = pd.read_csv('outputs/aggregate/aggregate_demix.tsv', sep='\t')['Unnamed: 0'].apply(lambda x: x.split('.')[0]).values

print('Current samples: ', len(current_samples))

new_samples = metadata[~metadata['Unnamed: 0'].isin(current_samples)]
print ('New samples: ', len(new_samples))
print ('Total samples: ', len(current_samples)+len(new_samples))

new_samples.to_csv('data/new_samples.csv', index=False)