import pandas as pd

metadata = pd.read_csv('data/all_metadata.csv', index_col=0)
metadata['SRA_published_date'] = pd.to_datetime(metadata['SRA_published_date'])

# Get accessions collected in the week 6 months ago
six_months_ago = pd.Timestamp.now() - pd.DateOffset(months=6)

print(six_months_ago)

# Select samples from 6 months ago
metadata = metadata[(metadata['SRA_published_date'] > six_months_ago) &
                    (metadata['SRA_published_date'] < six_months_ago + pd.DateOffset(days=1))]

metadata = metadata[metadata['sample_status'] == 'completed']
# save index to txt
metadata.index.to_series().to_csv('data/samples_to_rerun.txt', index=False, header=False)