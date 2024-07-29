import os
import sys
import yaml
import numpy as np
import pandas as pd
from datetime import datetime


def agg(results):
    allResults = [pd.read_csv(fn, skipinitialspace=True, sep='\t',
                              index_col=0) for fn in results]
    df_demix = pd.DataFrame()
    try:
        df_demix = pd.concat(allResults, axis=1).T

    except ValueError as e:
        print(e)
        print('No valid demix results found')
        # Create empty json
        os.makedirs('outputs/aggregate', exist_ok=True)
        empty_df = pd.DataFrame(columns=['sra_accession', 'name', 'prevalence', 'coverage', 'spike_coverage'])
        empty_df.to_json('outputs/aggregate/aggregate_demix_new.json', orient='records', lines=True)
        sys.exit()

    df_demix.index = [x.split('/')[-1] for x in df_demix.index]
    df_demix['sra_accession'] = [x.split('.')[0] for x in df_demix.index]
    df_demix = df_demix[~df_demix['lineages'].isna()]
    df_demix = df_demix.drop(columns=['summarized', 'resid'])    

    # Explode lineages and prevalences
    df_demix['lineages'] = df_demix['lineages'].str.split(' ')
    # Rename abundances col to prevalence
    df_demix = df_demix.rename(columns={'abundances': 'prevalence'})
    df_demix['prevalence'] = df_demix['prevalence'].str.split(' ')
    df_demix = df_demix.explode(['lineages', 'prevalence'], ignore_index=True)
    df_demix = df_demix.rename(columns={'lineages': 'name'})
    df_demix = df_demix[['sra_accession', 'name', 'prevalence', 'coverage']]

    return df_demix

def merge_collapsed(df_agg):
    # if -like is in the name, merge with the true lineage, if present
    df_agg['name'] = df_agg['name'].str.replace('-like', '')

    # Sum prevalences for the same lineage in the same sample
    df_agg = df_agg.groupby(['sra_accession', 'name']).agg({'prevalence': 'sum', 'coverage': 'first'}).reset_index()
    return df_agg


def get_alias_key(lineages_yml='data/lineages.yml'):
    with open(lineages_yml, 'r') as alias_key:
        lineage_key = yaml.load(alias_key, Loader=yaml.Loader)
    alias_key = dict([(lin['name'], lin['parent']) for lin in lineage_key if 'parent' in lin])
    alias_key.update([(lin['name'], lin['alias']) for lin in lineage_key if lin['name'] != lin['alias']])
    alias_key.update([r for lin in lineage_key for r in \
        [(lin['name'], lin['name'].split('.')[0]), (lin['name'].split('.')[0], lin['alias'])] \
        if (lin['name'] != lin['alias']) and len(lin['name'].split('.')) == 2 ])
    for n in range(4):
        alias_key.update([(alias, '.'.join(alias.split('.')[:-1])) for name, alias in alias_key.items() if not alias in alias_key and len(alias.split('.')) > 1])
    alias_key.update({'A.1': 'A', 'B.1': 'B'})
    return alias_key


def _crumbs(lin, alias_key):
    return [lin] + ( _crumbs(alias_key[lin], alias_key) if lin in alias_key else [])


def crumbs(lin, alias_key):
    lin = lin.upper()
    return _crumbs(lin, alias_key) if lin in alias_key else crumbs(lin[:-1], alias_key) if len(lin.split('.')) > 1 else []


def process_depth_file(acc):
    SPIKE_START = 21563
    SPIKE_END = 25384

    depth_file = f'{acc}.depths.tsv'
    if 'depth' not in depth_file:
        return 0

    if not os.path.exists(f'outputs/variants/{depth_file}'):
        return 0

    df_depth = pd.read_csv(f'outputs/variants/{depth_file}', sep='\t', names=['chromosome', 'position', 'base', 'depth'], header=None)
    df_depth['position'] = pd.to_numeric(df_depth['position'], errors='coerce')

    df_depth = df_depth[(df_depth['position'] >= SPIKE_START) & (df_depth['position'] <= SPIKE_END)]
    spike_coverage = len(df_depth[df_depth['depth'] >= 10]) / len(df_depth) if len(df_depth) > 0 else 0

    return spike_coverage * 100

def get_spike_coverage(df_agg):
    df_agg['spike_coverage'] = df_agg['sra_accession'].apply(process_depth_file)
    return df_agg


def main():
    # Load demix results
    os.makedirs('outputs/demix', exist_ok=True)
    demix_output = 'outputs/demix/'
    paths = [demix_output + fn for fn in os.listdir(demix_output)]
    df_agg = agg(paths)

    df_agg = merge_collapsed(df_agg)
    df_agg = df_agg[df_agg['coverage'].astype(float) > 0.0]
    
    # Get lineage breadcrumbs
    alias_key = get_alias_key()
    df_agg['crumbs'] = df_agg['name'].apply(lambda x: ';' + ';'.join(crumbs(x, alias_key)[::-1]) + ';')

    # Get spike coverage
    df_agg = get_spike_coverage(df_agg)

    df_agg = df_agg[df_agg['name'] != ''] 

    df_agg = df_agg.drop_duplicates(subset=['sra_accession', 'name'], keep='first')
    df_agg['prevalence'] = pd.to_numeric(df_agg['prevalence'], errors='coerce')
    df_agg = df_agg[~np.isinf(df_agg['prevalence'])]

    barcode_version = open('data/last_barcode_update.txt').readlines()[0].split('-')[0]
    barcode_version_date = datetime.strptime(barcode_version, '%m_%d_%Y')
    df_agg['barcode_version'] = barcode_version_date.strftime('%Y-%m-%d')

    os.makedirs('outputs/aggregate', exist_ok=True)
    df_agg.to_json('outputs/aggregate/aggregate_demix_new.json', orient='records', lines=True)


if __name__ == '__main__':
    main()