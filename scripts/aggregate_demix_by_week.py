import pandas as pd
from epiweeks import Week
import datetime

# Get existing aggregate outputs
demix = pd.read_json('outputs/aggregate/aggregate_demix.json', lines=True)
metadata = pd.read_json('outputs/aggregate/aggregate_metadata.json', lines=True)


demix = demix[demix['sra_accession'].isin(metadata['sra_accession'])]


# Remove samples where the sum of the lineage prevalence is greater than 1
print('accessions_before', demix['sra_accession'].nunique())
demix = demix[~demix.groupby('sra_accession')['prevalence'].transform('sum').gt(1)]
print('accessions_after', demix['sra_accession'].nunique())

# Add metadata to demix
df_agg = pd.merge(demix, metadata, on='sra_accession', how='left')


# Calculate population-weighted prevalence for each lineage
df_agg['prevalence'] = pd.to_numeric(df_agg['prevalence'])
df_agg['ww_population'] = pd.to_numeric(df_agg['ww_population'])
df_agg['pop_weighted_prevalence'] = df_agg['prevalence'] * df_agg['ww_population']

# Find the total population for each week and region
df_agg['collection_date'] = pd.to_datetime(df_agg['collection_date'])
df_agg['epiweek'] = df_agg['collection_date'].apply(lambda x: Week.fromdate(x))

df_agg_weekly = df_agg.groupby(['epiweek', 'name', 'geo_loc_region']).agg({
    'pop_weighted_prevalence': 'sum',
    'ww_population': 'sum',
    'sra_accession': 'nunique',
    'collection_date': 'first',
    'collection_site_id': 'nunique'
}).reset_index().rename(columns={
    'sra_accession': 'num_samples',
    'ww_population': 'total_population',
    'collection_site_id': 'num_collection_sites'
    })

# get mean prevalence for each lineage
df_agg_weekly['mean_lineage_prevalence'] = df_agg_weekly['pop_weighted_prevalence'] / df_agg_weekly['total_population']


df_agg_weekly['week_start'] = df_agg_weekly['epiweek'].apply(lambda x: x.startdate()).astype(str)
df_agg_weekly['week_end'] = df_agg_weekly['epiweek'].apply(lambda x: x.enddate()).astype(str)

df_agg_weekly = df_agg_weekly[['epiweek','week_start', 'week_end', 'total_population', 'num_collection_sites', 'num_samples', 'geo_loc_region', 'name', 'mean_lineage_prevalence']]
df_agg_weekly['id'] = df_agg_weekly['epiweek'].astype(str) + '_' + df_agg_weekly['geo_loc_region'].str.replace(' ', '_') + '_' + df_agg_weekly['name']


# Get number of duplicate ids
print('Number of duplicate ids:', df_agg_weekly['id'].duplicated().sum())

# Check that the prevalence for each state and week sums to 1
print('Sum of prevalence:', df_agg_weekly.groupby(['epiweek', 'geo_loc_region'])['mean_lineage_prevalence'].sum())

print(df_agg_weekly)
df_agg_weekly.to_json('outputs/aggregate/aggregate_demix_weekly.json', orient='records', lines=True)