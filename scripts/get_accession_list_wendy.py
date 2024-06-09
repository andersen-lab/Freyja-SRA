import sys
import pandas as pd

batch_size = int(sys.argv[1])
metadata = pd.read_csv('data/all_metadata.csv')

metadata = metadata[metadata['sample_status'] == 'to_run'] 

accessions = metadata['accession'][:batch_size] # Get first batch_size accessions

accessions.to_csv('data/accession_list.csv', index=False, header=False)
