import json
import pandas as pd
import shortuuid

agg_demix = pd.read_csv('aggregate/aggregate_demix.tsv', sep='\t')
metadata = pd.read_csv('data/wastewater_ncbi.csv')

columns = ['accession', 'lineages', 'abundances', 'collection_date', 'geo_loc_name', 'ww_population', 'ww_surv_target_1_conc', 'collection_site_id']

df = pd.DataFrame(columns=columns)

agg_demix['Unnamed: 0'] = agg_demix['Unnamed: 0'].apply(lambda x: x.split('.')[0])

df['accession'] = agg_demix['Unnamed: 0']
df['lineages'] = agg_demix['lineages']
df['abundances'] = agg_demix['abundances']


for col in ['collection_date', 'geo_loc_name', 'ww_population','ww_surv_target_1_conc', 'collection_site_id']:
    df[col] = [metadata[metadata['Unnamed: 0'] == x][col].values[0] for x in df['accession']]

df['collection_date'] = pd.to_datetime(df['collection_date'].apply(lambda x: x.split('/')[0] if '/' in x and len(x.split('/')[0])>2 else x))
df['ww_population'] = df['ww_population'].astype(float).astype(int)
df['ww_surv_target_1_conc'] = df['ww_surv_target_1_conc'].astype(float)
df = df.rename(columns={'ww_surv_target_1_conc':'viral_load'})

merged = df['geo_loc_name']+df['ww_population'].fillna('').astype(str)
merged = merged.apply(lambda x:shortuuid.uuid(x)[0:12])
df['site_id'] = df['collection_site_id'].combine_first(merged)
df.drop('collection_site_id', axis=1, inplace=True)

df.set_index('accession', inplace=True)

with open('aggregate/aggregate_demix.json', 'a') as f:
    for row in df.iterrows():
        json_row = {
            'sample_id': row[0],
            'lineages': [
                {'name': lineage, 'abundance': float(abundance)} for lineage, abundance in zip(row[1]['lineages'].split(' '), row[1]['abundances'].split(' '))
            ],
            'collection_date': row[1]['collection_date'].strftime('%Y-%m-%d'),
            'geo_loc_name': row[1]['geo_loc_name'],
            'ww_population': row[1]['ww_population'],
            'viral_load': row[1]['viral_load'],
            'site_id': row[1]['site_id']
        }
        
        json_row = json.dumps(json_row)
        f.write(json_row+'\n')
