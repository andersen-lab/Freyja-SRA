import pandas as pd
import os
from epiweeks import Week

# USA census regions
CENSUS_REGIONS = {
    'Northeast': ['Connecticut', 'Maine', 'Massachusetts', 'New Hampshire', 'Rhode Island', 'Vermont', 'New Jersey', 'New York', 'Pennsylvania'],
    'Midwest': ['Illinois', 'Indiana', 'Michigan', 'Ohio', 'Wisconsin', 'Iowa', 'Kansas', 'Minnesota', 'Missouri', 'Nebraska', 'North Dakota', 'South Dakota'],
    'South': ['Delaware', 'Maryland', 'Florida', 'Georgia', 'North Carolina', 'South Carolina', 'Virginia', 'District of Columbia', 'West Virginia', 'Alabama', 'Kentucky', 'Mississippi', 'Tennessee', 'Arkansas', 'Louisiana', 'Oklahoma', 'Texas'],
    'West': ['Arizona', 'Colorado', 'Idaho', 'Montana', 'Nevada', 'New Mexico', 'Utah', 'Wyoming', 'Alaska', 'California', 'Hawaii', 'Oregon', 'Washington']
}   

# Create state to region mapping for easier reference
STATE_TO_REGION = {}
for region, states in CENSUS_REGIONS.items():
    for state in states:
        STATE_TO_REGION[state] = region

# Get existing aggregate outputs
demix = pd.read_json('outputs/aggregate/aggregate_demix.json', lines=True)
metadata = pd.read_json('outputs/aggregate/aggregate_metadata.json', lines=True)

metadata = metadata[metadata['geo_loc_country'] == 'USA'] # USA only for now

demix = demix[demix['sra_accession'].isin(metadata['sra_accession'])]

crumbs = {k : v for k, v in zip(demix['name'], demix['crumbs'])}

# Remove samples where the sum of the lineage prevalence is greater than 1
demix = demix[~demix.groupby('sra_accession')['prevalence'].transform('sum').gt(1)]

# Add metadata to demix
df_agg = pd.merge(demix, metadata, on='sra_accession', how='left')

# Add census region to data
df_agg['census_region'] = df_agg['geo_loc_region'].map(STATE_TO_REGION)

# Calculate population-weighted lineage prevalence for each lineage (prevalence * sampled population)
df_agg['prevalence'] = pd.to_numeric(df_agg['prevalence'])
df_agg['ww_population'] = pd.to_numeric(df_agg['ww_population'])
df_agg['pop_weighted_prevalence'] = df_agg['prevalence'] * df_agg['ww_population']

# Find the total population for each week and region
df_agg['collection_date'] = pd.to_datetime(df_agg['collection_date'])
df_agg['epiweek'] = df_agg['collection_date'].apply(lambda x: Week.fromdate(x))

# First aggregate by state
total_lineage_prevalence_state = df_agg.groupby(['epiweek', 'geo_loc_region']).agg({
    'pop_weighted_prevalence': 'sum',
}).reset_index()

total_prev_state_dict = {total_lineage_prevalence_state['epiweek'].astype(str)[i] + '_' + total_lineage_prevalence_state['geo_loc_region'][i]: total_lineage_prevalence_state['pop_weighted_prevalence'][i] for i in range(len(total_lineage_prevalence_state))}

# Also calculate totals by census region
total_lineage_prevalence_region = df_agg.groupby(['epiweek', 'census_region']).agg({
    'pop_weighted_prevalence': 'sum',
}).reset_index()

total_prev_region_dict = {total_lineage_prevalence_region['epiweek'].astype(str)[i] + '_' + total_lineage_prevalence_region['census_region'][i]: total_lineage_prevalence_region['pop_weighted_prevalence'][i] for i in range(len(total_lineage_prevalence_region))}

# Calculate totals by nation (USA)
total_lineage_prevalence_nation = df_agg.groupby(['epiweek']).agg({
    'pop_weighted_prevalence': 'sum',
}).reset_index()

total_prev_nation_dict = {str(total_lineage_prevalence_nation['epiweek'][i]) + '_USA': total_lineage_prevalence_nation['pop_weighted_prevalence'][i] for i in range(len(total_lineage_prevalence_nation))}

# Aggregate by state
df_agg_weekly = df_agg.groupby(['epiweek', 'geo_loc_region', 'name']).agg({
    'pop_weighted_prevalence': 'sum',
    'collection_site_id': 'nunique',
    'sra_accession': 'nunique',
    'ww_population': 'mean',
}).reset_index().rename(columns={
    'collection_site_id': 'num_sites',
    'sra_accession': 'num_samples',
    'ww_population': 'total_population'
})

df_agg_weekly['id'] = df_agg_weekly['epiweek'].astype(str) + '_' + df_agg_weekly['geo_loc_region']
df_agg_weekly['total_lineage_prevalence'] = df_agg_weekly['id'].map(total_prev_state_dict)
df_agg_weekly['mean_lineage_prevalence'] = df_agg_weekly['pop_weighted_prevalence'] / df_agg_weekly['total_lineage_prevalence']

# Now aggregate by census region
df_agg_census = []
for region, states in CENSUS_REGIONS.items():
    # Filter data for this region
    df_region_data = df_agg[df_agg['census_region'] == region]
    
    # Aggregate by epiweek and lineage
    df_region = df_region_data.groupby(['epiweek', 'name']).agg({
        'pop_weighted_prevalence': 'sum',
        'collection_site_id': 'nunique',
        'sra_accession': 'nunique',
        'ww_population': 'mean',
    }).reset_index().rename(columns={
        'collection_site_id': 'num_sites',
        'sra_accession': 'num_samples',
        'ww_population': 'total_population'
    })
    
    # Calculate proper weighted mean prevalence for region
    df_region['geo_loc_region'] = region
    df_region['id'] = df_region['epiweek'].astype(str) + '_' + df_region['geo_loc_region']
    df_region['total_lineage_prevalence'] = df_region['id'].map(total_prev_region_dict)
    df_region['mean_lineage_prevalence'] = df_region['pop_weighted_prevalence'] / df_region['total_lineage_prevalence']
    
    df_agg_census.append(df_region)

# Aggregate by nation (USA)
df_nation = df_agg.groupby(['epiweek', 'name']).agg({
    'pop_weighted_prevalence': 'sum',
    'collection_site_id': 'nunique',
    'sra_accession': 'nunique',
    'ww_population': 'mean',
}).reset_index().rename(columns={
    'collection_site_id': 'num_sites',
    'sra_accession': 'num_samples',
    'ww_population': 'total_population'
})

# Calculate proper weighted mean prevalence for USA
df_nation['geo_loc_region'] = 'USA'
df_nation['id'] = df_nation['epiweek'].astype(str) + '_USA'
df_nation['total_lineage_prevalence'] = df_nation['id'].map(total_prev_nation_dict)
df_nation['mean_lineage_prevalence'] = df_nation['pop_weighted_prevalence'] / df_nation['total_lineage_prevalence']

# Combine all census regions with state data and national data
df_region_combined = pd.concat(df_agg_census)
df_agg_weekly = pd.concat([df_agg_weekly, df_region_combined, df_nation])
df_agg_weekly['total_population'] = df_agg_weekly.groupby(['epiweek', 'geo_loc_region'])['total_population'].transform('mean') # Ensure total population is consistent across lineages in the same region

df_agg_weekly['id'] = df_agg_weekly['epiweek'].astype(str) + '_' + df_agg_weekly['geo_loc_region'] + '_' + df_agg_weekly['name']
df_agg_weekly['crumbs'] = df_agg_weekly['name'].map(crumbs)
df_agg_weekly['week_start'] = df_agg_weekly['epiweek'].apply(lambda x: x.startdate()).astype(str)
df_agg_weekly['week_end'] = df_agg_weekly['epiweek'].apply(lambda x: x.enddate()).astype(str)

df_agg_weekly = df_agg_weekly[['id', 'epiweek', 'week_start', 'week_end', 'geo_loc_region', 'total_population', 'num_sites', 'num_samples', 'name', 'mean_lineage_prevalence', 'crumbs']]

# Workaround to save to json
df_agg_weekly.to_csv('outputs/aggregate/aggregate_demix_by_week.csv', index=False)
df_out = pd.read_csv('outputs/aggregate/aggregate_demix_by_week.csv')
os.remove('outputs/aggregate/aggregate_demix_by_week.csv')
df_out.to_json('outputs/aggregate/aggregate_demix_weekly.json', orient='records', lines=True)

