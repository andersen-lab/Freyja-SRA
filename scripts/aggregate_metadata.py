import os
import pandas as pd

paths_list = [entry.path for entry in os.scandir('outputs/variants') if 'variants' in entry.name]

accessions = [p.split('/')[-1].split('.')[0] for p in paths_list]

# Create metadata json
metadata = pd.read_csv('data/all_metadata.csv', index_col=None)
metadata = metadata[metadata['accession'].isin(accessions)]
metadata = metadata[['accession', 'collection_date', 'geo_loc_country', 'geo_loc_region', 'ww_population', 'ww_surv_target_1_conc', 'collection_site_id']]
metadata = metadata.rename(columns={'accession':'sra_accession', 'ww_surv_target_1_conc':'viral_load'})

os.makedirs('outputs/aggregate', exist_ok=True)
metadata.to_json('outputs/aggregate/aggregate_metadata_new.json', orient='records', lines=True)