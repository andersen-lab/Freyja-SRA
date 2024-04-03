import os
import yaml
import numpy as np
import pandas as pd


def agg(results):
    allResults = [pd.read_csv(fn, skipinitialspace=True, sep='\t',
                              index_col=0) for fn in results]
    df_demix = pd.concat(allResults, axis=1).T
    df_demix.index = [x.split('/')[-1] for x in df_demix.index]
    df_demix['sra_accession'] = [x.split('.')[0] for x in df_demix.index]
    df_demix = df_demix[~df_demix['lineages'].isna()]
    df_demix = df_demix.drop(columns=['summarized', 'resid'])

    # Explode lineages and abundances
    df_demix['lineages'] = df_demix['lineages'].str.split(' ')
    df_demix['abundances'] = df_demix['abundances'].str.split(' ')
    df_demix = df_demix.explode(['lineages', 'abundances'], ignore_index=True)
    df_demix = df_demix.rename(columns={'lineages': 'name', 'abundances': 'abundance'})
    df_demix = df_demix[['sra_accession', 'name', 'abundance']]
    return df_demix


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
    demix_output = 'outputs/demix/'
    paths = [demix_output + fn for fn in os.listdir(demix_output)][:5]
    df = agg(paths)
    
    # Get lineage breadcrumbs
    alias_key = get_alias_key()

    # Can store in dict if needed
    df['crumbs'] = df['name'].apply(lambda x: ';' + ';'.join(crumbs(x, alias_key)[::-1]) + ';')

    df = df[df['name'] != ''] 

    df = df.drop_duplicates(subset=['sra_accession', 'name'], keep='first')
    df['abundance'] = pd.to_numeric(df['abundance'], errors='coerce')
    df = df[~np.isinf(df['abundance'])]

    os.makedirs('outputs/aggregate', exist_ok=True)
    df.to_json('outputs/aggregate/aggregate_demix_new.json', orient='records', lines=True)


if __name__ == '__main__':
    main()