#!/usr/bin/env python3

import hashlib
import os
import json
import pandas as pd 

def md5hash(s: str): 
    return hashlib.md5(s.encode('utf-8')).hexdigest()[:64]

paths_list = []
for file in os.listdir('outputs/variants'):
    if 'variants' in file:
        paths_list.append(os.path.join('outputs/variants', file))

# Aggregate variants to one dataframe
depth_thresh = 20
freq_thresh = 0.01
j=0
for var_path in paths_list:
    try:
        df = pd.read_csv(var_path,sep='\t')
    except:
        continue

    df = df[df['ALT_DP']>=depth_thresh]
    df = df[df['ALT_FREQ']>=freq_thresh]

    #drop frame shifts
    sname=var_path.split('/')[-1].split('.')[0]
    df = df[df['ALT'].apply(lambda x: True if (('+' not in x and '-' not in x) or ((len(x)-1)%3==0)) else False)]
    if len(df)==0:
        continue
    df['mutName'] = df['REF'] + df['POS'].astype(str) + df['ALT']
    df['sample'] = sname
    df = df[['mutName','sample','ALT_FREQ','ALT_DP']]
    if j==0:
        variants = df.copy()
    else:
        variants = pd.concat((variants,df), axis=0)
    j+=1

# Load in metadata
metadata = pd.read_csv('data/all_metadata.csv')
metadata.set_index('accession', inplace=True)

# Remove samples that are not in metadata
variants = variants[variants['sample'].isin(metadata.index)]

variants = variants.rename(columns={'ALT_FREQ':'frequency', 'ALT_DP':'depth'})

acc_df = variants.groupby('sample').agg({'mutName': lambda x: ' '.join(x), 'frequency': lambda x: ' '.join([str(v) for v in x]), 'depth': lambda x: ' '.join([str(v) for v in x])}).reset_index()

for col in ['collection_date', 'geo_loc_country', 'geo_loc_region', 'ww_population', 'ww_surv_target_1_conc', 'site_id']:
    acc_df[col] = acc_df['sample'].map(metadata[col])
    variants[col] = variants['sample'].map(metadata[col])


variants_by_acc_path = 'outputs/aggregate/aggregate_variants_by_acc.json'

with open(variants_by_acc_path, 'w') as output_file:
    for _, row in acc_df.iterrows():
        mut_names = row['mutName'].split(' ')
        frequencies = [float(i) for i in row['frequency'].split(' ')]
        depths = [float(i) for i in row['depth'].split(' ')]

        mutations = [{'mut_name': mut_name, 'frequency': frequency, 'depth': depth}
                     for mut_name, frequency, depth in zip(mut_names, frequencies, depths)]
        
        row_dict = {
            'sra_accession': row['sample'],
            'mutations': mutations,
            'collection_date': row['collection_date'],
            'geo_loc_country': row['geo_loc_country'],
            'geo_loc_region': row['geo_loc_region'],
            'ww_population': float(row['ww_population']),
            'viral_load': float(row['ww_surv_target_1_conc']),
            'site_id': row['site_id']
        }

        output_file.write(f'{json.dumps(row_dict)}\n')

variants_by_mut_path = 'outputs/aggregate/aggregate_variants_by_mut.json'

mut_df = variants.groupby('mutName').agg({'sample': lambda x: ' '.join(x), 'frequency': lambda x: ' '.join([str(v) for v in x]), 'depth': lambda x: ' '.join([str(v) for v in x]), 'collection_date': lambda x: ' '.join([str(v) for v in x]), 'geo_loc_country': lambda x: ' '.join([str(v) for v in x]), 'geo_loc_region': lambda x: ':'.join([str(v) for v in x]), 'ww_population': lambda x: ' '.join([str(v) for v in x]), 'ww_surv_target_1_conc': lambda x: ' '.join([str(v) for v in x]), 'site_id': lambda x: ' '.join([str(v) for v in x])}).reset_index()
mut_df['mut_hash'] = mut_df['mutName'].apply(md5hash)

with open(variants_by_mut_path, 'w') as output_file:
    for _, row in mut_df.iterrows():
        sample_names = row['sample'].split(' ')
        frequencies = [float(i) for i in row['frequency'].split(' ')]
        depths = [float(i) for i in row['depth'].split(' ')]
        collection_dates = row['collection_date'].split(' ')
        geo_loc_countries = row['geo_loc_country'].split(' ')
        geo_loc_regions = row['geo_loc_region'].split(':')
        ww_populations = [float(i) for i in row['ww_population'].split(' ')]
        viral_loads = [float(i) for i in row['ww_surv_target_1_conc'].split(' ')]
        site_ids = row['site_id'].split(' ')

        samples = [
            {
                'sra_accession': sample_name,
                'frequency': frequency,
                'depth': depth,
                'collection_date': collection_date, 
                'geo_loc_country': geo_loc_country,
                'geo_loc_region': geo_loc_region,
                'ww_population': ww_population,
                'viral_load': viral_load,
                'site_id': site_id
            }
            for sample_name, frequency, depth, collection_date, geo_loc_country, geo_loc_region, ww_population, viral_load, site_id 
            in zip(sample_names, frequencies, depths, collection_dates, geo_loc_countries, geo_loc_regions, ww_populations, viral_loads, site_ids)
        ]

        row_dict = {
            'mut_hash': row['mut_hash'],
            'mut_name': row['mutName'],
            'samples': samples,
        }

        output_file.write(f'{json.dumps(row_dict)}\n')
