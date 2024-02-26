import os

files = ['demix', 'variants', 'metadata']

for file in files:
    with open(f'outputs/aggregate/aggregate_{file}.json', 'a') as outfile:
        with open(f'outputs/aggregate/aggregate_{file}_new.json', 'r') as infile:
            outfile.write(infile.read())

    os.remove(f'outputs/aggregate/aggregate_{file}_new.json')
