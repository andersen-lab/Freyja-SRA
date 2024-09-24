import sys
import pandas as pd

batch_size = int(sys.argv[1])
metadata = pd.read_csv('data/all_metadata.csv')

metadata = metadata[metadata['sample_status'] == 'fastq_error'] # Only run with fasterq_dump if aspera download failed

accessions = metadata['accession'][:batch_size] # Get first batch_size accessions

if len(accessions) == 0:
    print('No accessions to run')
    sys.exit(0)

accessions.to_csv('data/accession_list.csv', index=False, header=False)