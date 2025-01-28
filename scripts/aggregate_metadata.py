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
    return np.dstack((edges[::2], edges[1::2]))[0]

def format_intervals(intervals):
    return [{'start': int(i[0]), 'end': int(i[1])} for i in intervals]

paths_list = [entry.path for entry in os.scandir('outputs/variants') if 'variants' in entry.name]
depths_list = [entry.path for entry in os.scandir('outputs/variants') if 'depths' in entry.name]
demix_success = [entry.path.split('.')[0].split('/')[-1] for entry in os.scandir('outputs/demix') if 'demix' in entry.name]

for acc in demix_success:
    with open(f'outputs/demix/{acc}.demix.tsv', 'r') as f:
        lines = f.readlines()
        if len(lines[2].split(' ')) == 1:
            demix_success.remove(acc)

accessions = [p.split('/')[-1].split('.')[0] for p in paths_list]

# Create metadata json
metadata = pd.read_csv('data/all_metadata.csv', index_col=None, low_memory=False)
metadata = metadata[metadata['accession'].isin(accessions)]
metadata = metadata[['accession', 'collection_date', 'geo_loc_country', 'geo_loc_region', 'ww_population', 'collected_by', 'ww_surv_target_1_conc','ww_surv_target_1_conc_unit', 'collection_site_id']]
metadata = metadata.rename(columns={'accession':'sra_accession', 'ww_surv_target_1_conc':'viral_load', 'ww_surv_target_1_conc_unit':'viral_load_unit'})
metadata['ww_population'] = metadata['ww_population'].fillna(-1.0)
metadata = metadata.drop_duplicates(subset='sra_accession', keep='first')


# Check if demix output exists and has coverage > 0
metadata['demix_success'] = metadata['sra_accession'].isin(demix_success)
agg_demix = pd.read_json('outputs/aggregate/aggregate_demix_new.json', orient='records', lines=True).drop_duplicates(subset='sra_accession', keep='first')

try:
    metadata['demix_success'] = metadata['sra_accession'].isin(agg_demix['sra_accession']) & (metadata['sra_accession'].isin(demix_success))
except:
    metadata['demix_success'] = False

agg_variants = pd.read_json('outputs/aggregate/aggregate_variants_new.json', orient='records', lines=True).drop_duplicates(subset='sra_accession', keep='first')


try:
    metadata['variants_success'] = metadata['sra_accession'].isin(agg_variants['sra_accession'])
except KeyError:
    metadata['variants_success'] = False

# if variants_success is false, set demix_success to false
metadata['demix_success'] = np.where(metadata['variants_success']==False, False, metadata['demix_success'])

metadata['coverage_intervals'] = metadata['sra_accession'].apply(get_intervals)
metadata['coverage_intervals'] = metadata['coverage_intervals'].apply(format_intervals)

os.makedirs('outputs/aggregate', exist_ok=True)
metadata.to_json('outputs/aggregate/aggregate_metadata_new.json', orient='records', lines=True)