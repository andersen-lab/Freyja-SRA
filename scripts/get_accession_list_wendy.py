import sys
import pandas as pd

batch_size = int(sys.argv[1])
metadata = pd.read_csv('data/all_metadata.csv')

## Temporary: Restrict to Freyja-global period (USA, 04-01-2022 - 10-01-2023)
metadata = metadata[metadata['collection_date'] >= '2022-04-01']
metadata = metadata[metadata['collection_date'] < '2023-10-01']
metadata = metadata[metadata['geo_loc_name'].str.contains('USA', na=False)]
######################################################################

metadata = metadata[metadata['sample_status'] == 'to_run']
accessions = metadata['accession'][-1 * batch_size:] # Get last batch_size accessions

accessions.to_csv('data/accession_list.csv', index=False, header=False)
