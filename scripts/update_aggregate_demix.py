import pandas as pd

# If samples_to_rerun.txt is not empty
# Read samples_to_rerun.txt
try:
    samples_to_rerun = pd.read_csv('data/samples_to_rerun.txt', header=None)
except pd.errors.EmptyDataError:
    print('No samples to rerun.')
    exit()

existing = pd.read_json('outputs/aggregate/aggregate_demix.json', lines=True, orient='records')
new = pd.read_json('outputs/aggregate/aggregate_demix_new.json', lines=True, orient='records')

# Delete rows containing 'sra_accesion' in existing that are in new
existing = existing[~existing['sra_accession'].isin(new['sra_accession'])]

# Concatenate existing and new
df = pd.concat([existing, new], ignore_index=True)

# Save to json
df.to_json('outputs/aggregate/aggregate_demix.json', orient='records', lines=True)