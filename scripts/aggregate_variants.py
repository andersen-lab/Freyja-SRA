import hashlib
import os
import json
import pandas as pd 

def md5hash(s: str): 
    return hashlib.md5(s.encode('utf-8')).hexdigest()[:64]

paths_list = [entry.path for entry in os.scandir('outputs/variants') if 'variants' in entry.name]

variants_list = []
depth_thresh = 20
for var_path in paths_list:
    try:
        df = pd.read_csv(var_path,sep='\t')
    except:
        continue

    avg_qual = df['ALT_QUAL'].mean()
    freq_thresh = 10 ** (-avg_qual/10)

    df = df[(df['ALT_FREQ']>freq_thresh) & (df['ALT_DP']>=depth_thresh)]
    df = df[(~df['ALT'].str.contains('[+-]')) | ((df['ALT'].str.len()-1)%3==0)]
    if df.empty:
        continue

    sname=var_path.split('/')[-1].split('.')[0]
    df['mutName'] = df['REF'] + df['POS'].astype(str) + df['ALT']
    df['sample'] = sname
    df = df[['mutName','sample','ALT_FREQ','ALT_DP']]
    variants_list.append(df)

if not variants_list:
    sys.exit()

variants = pd.concat(variants_list, axis=0)

metadata = pd.read_csv('data/all_metadata.csv')
metadata.set_index('accession', inplace=True)

variants = variants[variants['sample'].isin(metadata.index)]
variants = variants.rename(columns={'ALT_FREQ':'frequency', 'ALT_DP':'depth'})

acc_df = variants.groupby('sample').agg({'mutName': lambda x: ' '.join(x), 'frequency': lambda x: ' '.join([str(v) for v in x]), 'depth': lambda x: ' '.join([str(v) for v in x])}).reset_index()

for col in ['collection_date', 'geo_loc_country', 'geo_loc_region', 'ww_population', 'ww_surv_target_1_conc', 'site_id']:
    acc_df[col] = acc_df['sample'].map(metadata[col])
    variants[col] = variants['sample'].map(metadata[col])

os.makedirs('outputs/aggregate', exist_ok=True)


with open('outputs/aggregate/aggregate_variants_by_acc_new.json', 'w') as output_file:
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


mut_df = variants.groupby('mutName').agg({'sample': lambda x: ' '.join(x), 'frequency': lambda x: ' '.join([str(v) for v in x]), 'depth': lambda x: ' '.join([str(v) for v in x]), 'collection_date': lambda x: ' '.join([str(v) for v in x]), 'geo_loc_country': lambda x: ' '.join([str(v) for v in x]), 'geo_loc_region': lambda x: ':'.join([str(v) for v in x]), 'ww_population': lambda x: ' '.join([str(v) for v in x]), 'ww_surv_target_1_conc': lambda x: ' '.join([str(v) for v in x]), 'site_id': lambda x: ' '.join([str(v) for v in x])}).reset_index()
mut_df['mut_key'] = mut_df['mutName'].apply(md5hash)

# Since variants_by_mut is keyed by mutation name, we need to append the new data to the existing keys
with open('outputs/aggregate/aggregate_variants_by_mut.json', 'r') as f:
    data = []
    for line in f:
        data.append(json.loads(line))
    existing_data = pd.DataFrame(data) 

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
        'mut_key': row['mut_key'],
        'mut_name': row['mutName'],
        'samples': samples
    }

    # If the mutation is not in the existing data, add it
    if row['mut_key'] not in existing_data['mut_key'].values:
        existing_data = pd.concat([existing_data, pd.DataFrame([row_dict])])
        continue

    # Otherwise, append the new data to the existing data
    existing_data.loc[existing_data['mut_key'] == row['mut_key'], 'samples'] = existing_data.loc[existing_data['mut_key'] == row['mut_key'], 'samples'].apply(lambda x: x + samples)

    
existing_data.to_json('outputs/aggregate/aggregate_variants_by_mut.json', orient='records', lines=True)