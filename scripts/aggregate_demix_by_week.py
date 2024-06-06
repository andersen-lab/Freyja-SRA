import pandas as pd
import os
from epiweeks import Week
import datetime

# Get existing aggregate outputs
demix = pd.read_json('outputs/aggregate/aggregate_demix.json', lines=True)
metadata = pd.read_json('outputs/aggregate/aggregate_metadata.json', lines=True)


demix = demix[demix['sra_accession'].isin(metadata['sra_accession'])]


# Remove samples where the sum of the lineage prevalence is greater than 1
demix = demix[~demix.groupby('sra_accession')['prevalence'].transform('sum').gt(1)]

# Add metadata to demix
df_agg = pd.merge(demix, metadata, on='sra_accession', how='left')


# Calculate population-weighted prevalence for each lineage
df_agg['prevalence'] = pd.to_numeric(df_agg['prevalence'])
df_agg['ww_population'] = pd.to_numeric(df_agg['ww_population'])
df_agg['pop_weighted_prevalence'] = df_agg['prevalence'] * df_agg['ww_population']

# Find the total population for each week and region
df_agg['collection_date'] = pd.to_datetime(df_agg['collection_date'])
df_agg['epiweek'] = df_agg['collection_date'].apply(lambda x: Week.fromdate(x))

# Create a dictionary mapping epiweek and region to the total prevalence
total_lineage_prevalence = df_agg.groupby(['epiweek', 'geo_loc_region']).agg({
    'pop_weighted_prevalence': 'sum',
}).reset_index()

total_prev_dict = {total_lineage_prevalence['epiweek'].astype(str)[i] + '_' + total_lineage_prevalence['geo_loc_region'][i]: total_lineage_prevalence['pop_weighted_prevalence'][i] for i in range(len(total_lineage_prevalence))}

df_agg_weekly = df_agg.groupby(['epiweek', 'geo_loc_region', 'name']).agg({
    'pop_weighted_prevalence': 'sum',
    'collection_site_id': 'nunique',
    'sra_accession': 'nunique',
}).reset_index().rename(columns={
    'collection_site_id': 'num_sites',
    'sra_accession': 'num_samples',
})

df_agg_weekly['id'] = df_agg_weekly['epiweek'].astype(str) + '_' + df_agg_weekly['geo_loc_region']
df_agg_weekly['total_lineage_prevalence'] = df_agg_weekly['id'].map(total_prev_dict)
df_agg_weekly['mean_lineage_prevalence'] = df_agg_weekly['pop_weighted_prevalence'] / df_agg_weekly['total_lineage_prevalence']


df_agg_weekly['id'] = df_agg_weekly['id'] + '_' + df_agg_weekly['name']


df_agg_weekly['week_start'] = df_agg_weekly['epiweek'].apply(lambda x: x.startdate()).astype(str)
df_agg_weekly['week_end'] = df_agg_weekly['epiweek'].apply(lambda x: x.enddate()).astype(str)

df_agg_weekly = df_agg_weekly[['id', 'epiweek', 'week_start', 'week_end', 'geo_loc_region','num_sites', 'num_samples', 'name', 'mean_lineage_prevalence']]

df_agg_weekly.to_csv('outputs/aggregate/aggregate_demix_by_week.csv', index=False)
df_out = pd.read_csv('outputs/aggregate/aggregate_demix_by_week.csv')
os.remove('outputs/aggregate/aggregate_demix_by_week.csv')
df_out.to_json('outputs/aggregate/aggregate_demix_weekly.json', orient='records', lines=True)
