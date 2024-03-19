import pandas as pd
import os

files = ['demix', 'variants', 'metadata']

for file in files:
    with open(f'outputs/aggregate/aggregate_{file}.json', 'a') as outfile:
        with open(f'outputs/aggregate/aggregate_{file}_new.json', 'r') as infile:
            outfile.write(infile.read())

    if file == 'variants':
        df = pd.read_json(f'outputs/aggregate/aggregate_{file}.json', lines=True).drop_duplicates(subset=['sra_accession', 'site'], keep='first')
    else:
        df = pd.read_json(f'outputs/aggregate/aggregate_{file}.json', lines=True).drop_duplicates(subset='sra_accession', keep='first')
        
    df.to_json(f'outputs/aggregate/aggregate_{file}_dedup.json', orient='records', lines=True)

    os.remove(f'outputs/aggregate/aggregate_{file}.json')
    os.remove(f'outputs/aggregate/aggregate_{file}_new.json')
    os.rename(f'outputs/aggregate/aggregate_{file}_dedup.json', f'outputs/aggregate/aggregate_{file}.json')