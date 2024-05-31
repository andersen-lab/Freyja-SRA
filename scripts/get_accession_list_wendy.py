import sys
import pandas as pd

batch_size = int(sys.argv[1])
metadata = pd.read_csv('data/all_metadata.csv')

metadata = metadata[metadata['sample_status'] == 'to_run'] 
accessions = metadata['accession'][-1 * batch_size:] # Get last batch_size accessions

accessions.to_csv('data/accession_list.csv', index=False, header=False)
