import os
import pandas as pd
import numpy as np

def get_intervals(accession):
    try:
        df_depth = pd.read_csv(f'outputs/variants/{accession}.depths.tsv', sep='\t', header=None, index_col=1)[3]
    except:
        return []
    vec = df_depth >= 10
    vec = vec.astype(int)
    if len(vec)==0:
        return []
    elif not isinstance(vec, np.ndarray):
        vec = np.array(vec)

    edges, = np.nonzero(np.diff((vec==0)*1))
    edge_vec = [edges+1]
    if vec[0] != 0:
        edge_vec.insert(0, [0])
    if vec[-1] != 0:
        edge_vec.append([len(vec)])
    edges = np.concatenate(edge_vec)
    return np.dstack((edges[::2], edges[1::2]))

paths_list = [entry.path for entry in os.scandir('outputs/variants') if 'variants' in entry.name]
depths_list = [entry.path for entry in os.scandir('outputs/variants') if 'depths' in entry.name]
demix_success = [entry.path.split('.')[0].split('/')[-1] for entry in os.scandir('outputs/demix') if 'demix' in entry.name]

accessions = [p.split('/')[-1].split('.')[0] for p in paths_list]

# Create metadata json
metadata = pd.read_csv('data/all_metadata.csv', index_col=None, low_memory=False)
metadata = metadata[metadata['accession'].isin(accessions)]
metadata = metadata[['accession', 'collection_date', 'geo_loc_country', 'geo_loc_region', 'ww_population', 'ww_surv_target_1_conc', 'collection_site_id']]
metadata = metadata.rename(columns={'accession':'sra_accession', 'ww_surv_target_1_conc':'viral_load'})
metadata['ww_population'] = metadata['ww_population'].fillna(-1.0)
metadata = metadata.drop_duplicates(subset='sra_accession', keep='first')


# Check if demix output exists and has coverage > 0
metadata['demix_success'] = metadata['sra_accession'].isin(demix_success)
agg_demix = pd.read_json('outputs/aggregate/aggregate_demix_new.json', orient='records', lines=True).drop_duplicates(subset='sra_accession', keep='first')

try:
    metadata['demix_success'] = metadata['sra_accession'].isin(agg_demix['sra_accession']) & (metadata['sra_accession'].isin(demix_success) | metadata['sra_accession'].isin(agg_demix[agg_demix['coverage'] == 0.0]['sra_accession']))
except:
    metadata['demix_success'] = False

agg_variants = pd.read_json('outputs/aggregate/aggregate_variants_new.json', orient='records', lines=True).drop_duplicates(subset='sra_accession', keep='first')

try:
    metadata['variants_success'] = metadata['sra_accession'].isin(agg_variants['sra_accession'])
except:
    metadata['variants_success'] = False

metadata['coverage_intervals'] = metadata['sra_accession'].apply(get_intervals)

os.makedirs('outputs/aggregate', exist_ok=True)
metadata.to_json('outputs/aggregate/aggregate_metadata_new.json', orient='records', lines=True)