import pandas as pd

# Get existing aggregate outputs
agg_demix = pd.read_json('outputs/aggregate/aggregate_demix.json', lines=True)
agg_metadata = pd.read_json('outputs/aggregate/aggregate_metadata.json', lines=True)
df_agg = pd.merge(agg_demix, agg_metadata, on='sra_accession', how='left')

# Calculate population-weighted prevalence for each lineage
df_agg['prevalence'] = pd.to_numeric(df_agg['prevalence'], errors='coerce')
df_agg['ww_population'] = pd.to_numeric(df_agg['ww_population'], errors='coerce')
df_agg['pop_weighted_prevalence'] = df_agg['prevalence'] * df_agg['ww_population']

# Aggregate by week
df_agg['week'] = pd.to_datetime(df_agg['collection_date']).dt.strftime('%U')
df_agg_weekly = df_agg.groupby(['week', 'name', 'geo_loc_region']).agg({
    'pop_weighted_prevalence': 'sum',
    'ww_population': 'sum',
    'sra_accession': 'nunique',
    'collection_date': 'first',
    'collection_site_id': 'nunique'
}).reset_index()

df_agg_weekly['date_start'] = pd.to_datetime(df_agg_weekly['collection_date']).dt.to_period('W').dt.start_time.dt.strftime('%Y-%m-%d')
df_agg_weekly['date_end'] = pd.to_datetime(df_agg_weekly['collection_date']).dt.to_period('W').dt.end_time.dt.strftime('%Y-%m-%d')

df_agg_weekly = df_agg_weekly.rename(columns={
    'sra_accession': 'num_samples',
    'ww_population': 'total_population',
    'collection_site_id': 'num_collection_sites'
    })

# Divide by total population to get population-weighted mean prevalence
df_agg_weekly['mean_lineage_prevalence'] = df_agg_weekly['pop_weighted_prevalence'] / df_agg_weekly['total_population']
df_agg_weekly = df_agg_weekly.drop(columns=['pop_weighted_prevalence', 'collection_date', 'week'])

df_agg_weekly = df_agg_weekly[['date_start', 'date_end', 'total_population', 'num_collection_sites', 'num_samples', 'geo_loc_region', 'name', 'mean_lineage_prevalence']]
df_agg_weekly.to_json('outputs/aggregate/aggregate_demix_weekly.json', orient='records', lines=True)