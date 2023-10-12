import json
import pandas as pd
import shortuuid

def merge_collapsed(lin_dict):
    new_dict = {}
    for k in lin_dict.keys():
        if k.endswith('-like'):
            true_lin = k.split('-')[0]
            if true_lin in lin_dict:
                new_dict[true_lin] = lin_dict[k] + lin_dict[true_lin]
            else:
                new_dict[true_lin] = lin_dict[k]
        elif 'Misc' not in k and 'like' not in k:
            new_dict[k] = lin_dict[k]
    return new_dict

# Merge collapsed lineages (*-like)

agg_demix = pd.read_csv('outputs/aggregate/aggregate_demix.tsv', sep='\t')


agg_demix['lin_dict'] = [dict(zip(row['lineages'].split(' '), map(float, row['abundances'].split(' ')))) for _, row in agg_demix.iterrows()]
agg_demix['lin_dict'].apply(merge_collapsed)
agg_demix['lineages'] = agg_demix['lin_dict'].apply(lambda x: ' '.join(list(x.keys())))
agg_demix['abundances'] = agg_demix['lin_dict'].apply(lambda x: ' '.join([str(v) for v in list(x.values())]))
agg_demix.drop('lin_dict', axis=1, inplace=True)
metadata = pd.read_csv('data/wastewater_ncbi.csv')

columns = ['accession', 'lineages', 'abundances', 'collection_date', 'geo_loc_name', 'ww_population', 'ww_surv_target_1_conc', 'collection_site_id']

df = pd.DataFrame(columns=columns)

agg_demix['Unnamed: 0'] = agg_demix['Unnamed: 0'].apply(lambda x: x.split('.')[0])

# Drop samples that are not in the metadata
agg_demix = agg_demix[agg_demix['Unnamed: 0'].isin(metadata['Unnamed: 0'])]

df['accession'] = agg_demix['Unnamed: 0']
df['lineages'] = agg_demix['lineages']
df['abundances'] = agg_demix['abundances']

for col in ['collection_date', 'geo_loc_name', 'ww_population','ww_surv_target_1_conc', 'collection_site_id']:
    df[col] = [metadata[metadata['Unnamed: 0'] == x][col] for x in df['accession']]


df['ww_surv_target_1_conc'] = df['ww_surv_target_1_conc'].astype(float)
df = df.rename(columns={'ww_surv_target_1_conc':'viral_load'})

merged = df['geo_loc_name']+df['ww_population'].astype(str)
merged = merged.apply(lambda x: str(x))
merged = merged.apply(lambda x: shortuuid.uuid(x)[0:12])

df['collection_site_id'] = df['collection_site_id'].astype(str).str.split('\n').apply(lambda x: x[0]).str.split(' ').apply(lambda x: x[4])

# iteratively replace NaN values in collection_site_id with values from merged
site_id_list = []
for row in df.iterrows():
    if row[1]['collection_site_id'] == 'NaN':
        site_id_list.append(merged[row[0]])
    else:
        site_id_list.append(row[1]['collection_site_id'])

df['site_id'] = pd.Series(site_id_list, index=df.index)

df.set_index('accession', inplace=True)

df['viral_load'] = df['viral_load'].fillna(-1.0)
df['ww_population'] = df['ww_population'].fillna(-1.0)

with open('outputs/aggregate/aggregate_demix.json', 'w') as f:
    for row in df.iterrows():
        json_row = {
            'sample_id': row[0],
            'lineages': [
                {'name': lineage, 'abundance': float(abundance)} for lineage, abundance in zip(row[1]['lineages'].split(' '), row[1]['abundances'].split(' '))
            ],
            'collection_date': row[1]['collection_date'].values[0],
            'geo_loc_name': row[1]['geo_loc_name'].values[0],
            'ww_population': float(row[1]['ww_population'].values[0]),
            'viral_load': row[1]['viral_load'],
            'site_id': row[1]['site_id']
        }
        
        json_row = json.dumps(json_row)
        f.write(json_row+'\n')
