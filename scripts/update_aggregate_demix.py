import pandas as pd

existing = pd.read_json('outputs/aggregate/aggregate_demix.json', lines=True, orient='records')
new = pd.read_json('outputs/aggregate/aggregate_demix_new.json', lines=True, orient='records')

# Delete rows containing 'sra_accesion' in existing that are in new
existing = existing[~existing['sra_accession'].isin(new['sra_accession'])]

# Concatenate existing and new
df = pd.concat([existing, new], ignore_index=True)

# Save to json
df.to_json('outputs/aggregate/aggregate_demix.json', orient='records', lines=True)