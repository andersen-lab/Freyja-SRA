process AGGREGATE_VARIANTS {
    publishDir "${params.output}/aggregate", mode: 'copy'
    input:
    val variants_outputs
    path baseDir

    output:
    path "aggregate_variants.tsv"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd 
    import os

    paths_string = "${variants_outputs}"
    if not paths_string.startswith('['):
        paths_string = f'[{paths_string}]'
    
    paths_list = paths_string[1:-1].split(", ")

    ct_thresh = 20
    j=0
    for var_path in paths_list:
        df = pd.read_csv(var_path,sep='\t')
        df = df[df['ALT_DP']>=ct_thresh]
        #drop frame shifts
        sname=var_path.split('/')[-1].split('.')[0]
        df = df[df['ALT'].apply(lambda x: True if (('+' not in x and '-' not in x) or ((len(x)-1)%3==0)) else False)]
        if len(df)==0:
            print('No data in ',var_path)
            continue
        df['mutName'] = df['REF'] + df['POS'].astype(str) + df['ALT']
        df['sample'] = sname
        df = df[['mutName','sample','ALT_FREQ','ALT_DP']]
        if j==0:
            dfAll = df.copy()
        else:
            dfAll = pd.concat((dfAll,df), axis=0)
        j+=1

    dfAll = dfAll.set_index(['mutName','sample'])
    
    dfAll.to_csv('aggregate_variants.tsv',sep='\t')
    """
}

process AGGREGATE_DEMIX {
    publishDir "${params.output}/aggregate", mode: 'copy'
    input:
    val demix_outputs
    path baseDir

    output:
    path "aggregate_demix.tsv"
    path "aggregate_demix.json"

    script:
    """
    #!/usr/bin/env python3
    import subprocess
    import json
    import pandas as pd
    import shortuuid

    paths_string = "${demix_outputs}"
    paths_list = paths_string[1:-1].split(", ")

    subprocess.run(["mkdir", "aggregate_dir"])
    for file in paths_list:
       subprocess.run(["cp", file, "aggregate_dir"])

    # Save to tsv
    subprocess.run(["freyja", "aggregate", "aggregate_dir/", "--output", "aggregate_demix.tsv"])
    
    # Save to json
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
    """
}

process AGGREGATE_COVARIANTS {
    publishDir "${params.output}/aggregate", mode: 'copy'
    input:
    val covariants_outputs
    path baseDir

    output:
    path "aggregate_covariants.tsv"

    script:
    """
    #!/usr/bin/env python3
    import subprocess
    import pandas as pd

    paths_string = "${covariants_outputs}"
    paths_list = paths_string[1:-1].split(", ")

    agg_df = pd.DataFrame(columns=['Covariants','Sample', 'Count', 'Max_count', 'Freq', 'Coverage_start', 'Coverage_end'])

    for covar_path in paths_list:
        df = pd.read_csv(covar_path,sep='\t')
        sample_name = covar_path.split('/')[-1].split('.')[0]
        df['Sample'] = sample_name
        agg_df = pd.concat((agg_df,df), axis=0)

    agg_df = agg_df.set_index(['Covariants','Sample'])
    agg_df.to_csv('aggregate_covariants.tsv',sep='\t')
    """
}

process DEMIX_TO_JSON {
    publishDir "${params.output}/aggregate", mode: 'copy'

    input:
    path demix_tsv

    output:
    path "aggregate_demix.json"

    script:
    """
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

    """

}