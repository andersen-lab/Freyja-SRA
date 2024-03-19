#!/usr/bin/env python3

import subprocess
import json
import os
import sys
import numpy as np
import yaml
import pandas as pd

def agg(results):
    allResults = [pd.read_csv(fn, skipinitialspace=True, sep='\t',
                              index_col=0) for fn in results]
    df_demix = pd.concat(allResults, axis=1).T
    df_demix.index = [x.split('/')[-1] for x in df_demix.index]
    df_demix['accession'] = [x.split('.')[0] for x in df_demix.index]
    return df_demix

def isnumber(x):
    try:
        float(x)
        return True
    except:
        return False

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

def merge_collapsed(lin_dict):
    new_dict = {}
    for k in lin_dict.keys():
        if '-like' in k:
            true_lin = k.split('-')[0]
            if true_lin in lin_dict:
                new_dict[true_lin] = lin_dict[k] + lin_dict[true_lin]
            else:
                new_dict[true_lin] = lin_dict[k]
        elif 'Misc' not in k:
            new_dict[k] = lin_dict[k]
    return new_dict

def main():

    # Load demix results
    results = 'outputs/demix/'
    results_ = [results + fn for fn in os.listdir(results)]
    df = agg(results_)

    df['lin_dict'] = [dict(zip(row['lineages'].split(' '), map(float, row['abundances'].split(' ')))) for _, row in df.iterrows()]
    df['lin_dict'] = df['lin_dict'].apply(merge_collapsed)
    df['lineages'] = df['lin_dict'].apply(lambda x: ' '.join(list(x.keys())))
    df['abundances'] = df['lin_dict'].apply(lambda x: ' '.join([str(v) for v in list(x.values())]))
    df.drop('lin_dict', axis=1, inplace=True)

    df.index.rename('accession', inplace=True)
    df = df.rename(columns={'Unnamed: 0': 'accession'})

    df['accession'] = df['accession'].apply(lambda x: x.split('.')[0])
   
    # Get lineage breadcrumbs
    alias_key = get_alias_key()
    df['crumbs'] = df['lineages'].apply(lambda x: [';' + ';'.join(crumbs(lin, alias_key)[::-1]) + ';' for lin in x.split(' ')])

    df = df.rename(columns={'ww_surv_target_1_conc':'viral_load'})
    df.set_index('accession', inplace=True)
    df = df[df['lineages'] != ''] 

    df = df[~df.index.duplicated(keep='first')]


    os.makedirs('outputs/aggregate', exist_ok=True)
    with open('outputs/aggregate/aggregate_demix_new.json', 'w') as f:
        for row in df.iterrows():

            json_row = {
                'sra_accession': row[0],
                'lineages': [
                    {'name': lineage, 'abundance': float(abundance), 'crumbs': crumbs} for lineage, abundance, crumbs in zip(row[1]['lineages'].split(' '), row[1]['abundances'].split(' '), row[1]['crumbs'])
                ],
                'coverage': row[1]['coverage']
            }
            if not np.isfinite([lin['abundance'] for lin in json_row['lineages']]).all():
                continue

            json_row = json.dumps(json_row)
            f.write(json_row+'\n')

if __name__ == '__main__':
    main()