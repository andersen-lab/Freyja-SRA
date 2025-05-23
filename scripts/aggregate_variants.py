import os
import sys
import pandas as pd 

def handle_alt_base(mut):
    if '+' in mut:
        return mut[mut.index('+'):]
    elif '-' in mut:
        return mut[mut.index('-'):]
    else:
        return mut[-1]

# Create variants json (by accession * site)
os.makedirs('outputs/variants', exist_ok=True)
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

    df['sra_accession'] = var_path.split('/')[-1].split('.')[0]
    df['mutName'] = df['REF'] + df['POS'].astype(str) + df['ALT']
    df = df[['mutName','sra_accession','ALT_FREQ','ALT_DP']]
    variants_list.append(df)

if not variants_list:
    print('No valid variants found')
    # Create empty json
    empty_df = pd.DataFrame(columns=['sra_accession', 'site', 'ref_base', 'alt_base', 'depth', 'prevalence'])

    empty_df.to_json('outputs/aggregate/aggregate_variants_new.json', orient='records', lines=True)
    sys.exit()

variants = pd.concat(variants_list, axis=0)
variants = variants.rename(columns={'ALT_FREQ':'prevalence', 'ALT_DP':'depth'})

variants['ref_base'] = variants['mutName'].apply(lambda x: x[0])

variants['site'] = variants['mutName'].apply(lambda x: x[1:-1])
variants['site'] = variants['site'].apply(lambda x: x[:x.index('+')] if '+' in x else x)
variants['site'] = variants['site'].apply(lambda x: x[:x.index('-')] if '-' in x else x)

variants['alt_base'] = variants['mutName'].apply(handle_alt_base)

variants = variants.drop(columns=['mutName'])
variants = variants.rename(columns={0:'variants'})
variants = variants.drop_duplicates(subset=['sra_accession', 'site', 'alt_base'], keep='first')

variants['site'] = variants['site'].apply(lambda x: x if str(x).isdigit() else None)
variants['ref_base'] = variants['ref_base'].apply(lambda x: x if len(x)==1 and not x.isdigit() else None)
variants['alt_base'] = variants['alt_base'].apply(lambda x: x if len(x)==1 and not x.isdigit() else None)

variants.to_csv('outputs/aggregate/aggregate_variants.csv', index=False)
variants = variants.dropna(subset=['site', 'depth', 'ref_base', 'alt_base'])
variants = variants.drop_duplicates(subset=['sra_accession', 'site', 'ref_base', 'alt_base'], keep='first')

variants = variants[['sra_accession', 'site', 'ref_base', 'alt_base', 'depth', 'prevalence']]
os.makedirs('outputs/aggregate', exist_ok=True)

variants.to_json('outputs/aggregate/aggregate_variants_new.json', orient='records', lines=True)
