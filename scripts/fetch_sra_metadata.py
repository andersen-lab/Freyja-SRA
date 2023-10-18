#!/usr/bin/env python

import argparse
import os
import pandas as pd
from Bio import Entrez
from Bio import SeqIO
import shortuuid

argparser = argparse.ArgumentParser(description='Fetch most recent SRA metadata')

def main():
    Entrez.email = "jolevy@scripps.edu"
    handle = Entrez.esearch(db="sra", idtype='acc', retmax=4000,
                            sort='recently_added',
                            term="((wastewater metagenome[Organism] OR wastewater metagenome[All Fields]) AND SARS-CoV-2))") 
    record = Entrez.read(handle)
    handle.close()

    handle = Entrez.efetch(db="sra", id=record['IdList'], rettype="gb",retmode='text')
    string= handle.read()
    handle.close()

    returned_meta=str(string,'UTF-8')

    with open("data/NCBI_metadata.xml", "w") as f:
        f.write(returned_meta)
    # #parse xml
    import xml.etree.ElementTree as ET
    root = ET.fromstring(returned_meta)
    allDictVals = {}

    for root0 in root:
        ### pull all sample attributes
        vals = [r.text for r in root0.findall('.//SAMPLE_ATTRIBUTE/')]
        sampExp = [r.text for r in root0.findall('.//EXPERIMENT/IDENTIFIERS/PRIMARY_ID')]
        seq_meta = [r.text for r in root0.findall('.//RUN_SET/RUN/RUN_ATTRIBUTES/RUN_ATTRIBUTE/')]
        sampID =  [r.text for r in root0.findall('.//RUN_SET/RUN/IDENTIFIERS/PRIMARY_ID')]
        if len(sampID)>1:
            print('more than one experiment... add funcs')
        else:
            sampID = sampID[0]
        ## write to dictionary form
        dictVals = {vals[i].replace(' ','_'):vals[i+1] for i in range(0,len(vals),2)}
        for i in range(0,len(seq_meta),2):
            dictVals[seq_meta[i].replace(' ','')] = seq_meta[i+1]
        dictVals['experiment_id'] = sampID
        dictVals['SRA_id'] = root0[0].attrib['accession']
        allDictVals[sampID] =dictVals

    metadata = pd.DataFrame(allDictVals).T

    metadata.columns = metadata.columns.str.replace(' ','_')
    # Remove samples where the collection date is not of the form 20XX-XX-XX
    metadata  = metadata[metadata['collection_date'].str.contains('20[0-9]{2}-[0-9]{2}-[0-9]{2}')]
    
    metadata['collection_date'] = pd.to_datetime(metadata['collection_date'].apply(lambda x: x.split('/')[0] if '/' in x else x))
    metadata = metadata.sort_values(by='collection_date',ascending=True)

    merged = metadata['geo_loc_name']+metadata['ww_population'].fillna('').astype(str)
    merged = merged.apply(lambda x:shortuuid.uuid(x)[0:12])
    metadata['site_id'] = metadata['collection_site_id'].combine_first(merged)

    metadata = metadata[metadata['collection_date'] >='2023-02-01']
    current_metadata = pd.read_csv('data/wastewater_ncbi.csv', index_col=0)

    # Get the number of newly added samples

    # Add samples that are not in the current metadata
    metadata = metadata[~metadata.index.isin(current_metadata.index)]
    num_new_samples = len(metadata)

    # Concatenate the two metadata files
    metadata = pd.concat([current_metadata, metadata], axis=0)

    print('All samples: ', len(metadata))
    current_samples = pd.read_csv('outputs/aggregate/aggregate_variants.tsv', sep='\t')['sample'].values
    current_samples = set(current_samples)
    print('Newly added samples: ', num_new_samples)
    print('Processed samples: ', len(current_samples))

    new_samples = metadata[~metadata.index.isin(current_samples)]

    print('Samples to run: ', len(new_samples))

    metadata.to_csv('data/wastewater_ncbi.csv')
    new_samples.to_csv('data/new_samples.csv', index=True, header=True)

if __name__ == "__main__":
    main()