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

    for file in os.listdir('${baseDir}/outputs/variants'):
        paths_list.append(os.path.join('${baseDir}/outputs/variants', file))

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
    import os
    import pandas as pd
    import shortuuid

    paths_string = "${demix_outputs}"
    paths_list = paths_string[1:-1].split(", ")
    for file in os.listdir('${baseDir}/outputs/demix'):
        paths_list.append(os.path.join('${baseDir}/outputs/demix', file))

    subprocess.run(["mkdir", "aggregate_dir"])
    for file in paths_list:
       subprocess.run(["cp", file, "aggregate_dir"])

    # Save to tsv
    subprocess.run(["freyja", "aggregate", "aggregate_dir/", "--output", "aggregate_demix.tsv"])

    # Save to json
    agg_demix = pd.read_csv('aggregate_demix.tsv', sep='\\t')
    metadata = pd.read_csv('${metadata}')

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


    df['ww_population'] = df['ww_population'].astype(float).astype(int)
    df['ww_surv_target_1_conc'] = df['ww_surv_target_1_conc'].astype(float)
    df = df.rename(columns={'ww_surv_target_1_conc':'viral_load'})

    merged = df['geo_loc_name']+df['ww_population'].fillna('').astype(str)
    merged = merged.apply(lambda x: str(x))

    merged = merged.apply(lambda x: shortuuid.uuid(x)[0:12])
    df['site_id'] = df['collection_site_id'].combine_first(merged)

    df.set_index('accession', inplace=True)

    with open('outputs/aggregate/aggregate_demix.json', 'w') as f:
        for row in df.iterrows():
            json_row = {
                'sample_id': row[0],
                'lineages': [
                    {'name': lineage, 'abundance': float(abundance)} for lineage, abundance in zip(row[1]['lineages'].split(' '), row[1]['abundances'].split(' '))
                ],
                'collection_date': row[1]['collection_date'].values[0],
                'geo_loc_name': row[1]['geo_loc_name'].values[0],
                'ww_population': row[1]['ww_population'],
                'viral_load': row[1]['viral_load'],
                'site_id': row[1]['site_id'].values[0]
            }
            
            json_row = json.dumps(json_row)
            f.write(json_row+'\\n')

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
    import os

    paths_string = "${covariants_outputs}"
    paths_list = paths_string[1:-1].split(", ")

    for file in os.listdir('${baseDir}/outputs/covariants'):
        paths_list.append(os.path.join('${baseDir}/outputs/covariants', file))

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