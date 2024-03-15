import os

files = ['demix', 'variants', 'metadata']

for file in files:
    with open(f'outputs/aggregate/aggregate_{file}.json', 'a') as outfile:
        with open(f'outputs/aggregate/aggregate_{file}_new.json', 'r') as infile:
            outfile.write(infile.read())

    os.system(f'awk \'!seen[$0]++\' outputs/aggregate/aggregate_{file}.json > outputs/aggregate/aggregate_{file}_dedup.json')        
    os.remove(f'outputs/aggregate/aggregate_{file}.json')
    os.remove(f'outputs/aggregate/aggregate_{file}_new.json')
    os.rename(f'outputs/aggregate/aggregate_{file}_dedup.json', f'outputs/aggregate/aggregate_{file}.json')