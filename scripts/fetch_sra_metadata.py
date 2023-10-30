#!/usr/bin/env python

import argparse
import os
import pandas as pd
from Bio import Entrez
from Bio import SeqIO
import xml.etree.ElementTree as ET
import shortuuid

argparser = argparse.ArgumentParser(description='Fetch most recent SRA metadata')

def get_metadata():
    Entrez.email = "jolevy@scripps.edu"
    handle = Entrez.esearch(db="sra", idtype='acc', retmax=4500,
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
        elif len(sampID)==0:
            continue
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

    return metadata
def main():
    metadata = get_metadata()

    # Convert collection date to datetime
    metadata  = metadata[metadata['collection_date'].str.contains('20[0-9]{2}-[0-9]{2}-[0-9]{2}')]
    metadata['collection_date'] = pd.to_datetime(metadata['collection_date'].apply(lambda x: x.split('/')[0] if '/' in x else x))
    metadata = metadata.sort_values(by='collection_date',ascending=False)

    # Filter to USA samples
    metadata = metadata[metadata['geo_loc_name'].str.contains('USA')]
    
    # For samples with no site id, hash the location and population to create a unique id
    merged = metadata['geo_loc_name']+metadata['ww_population'].fillna('').astype(str)
    merged = merged.apply(lambda x:shortuuid.uuid(x)[0:12])
    metadata['site_id'] = metadata['collection_site_id'].combine_first(merged)

    # Since Entrez returns the most recent samples, we need to concatenate the new metadata with the old metadata
    current_metadata = pd.read_csv('data/all_metadata.csv', index_col=0)
    new_metadata = metadata[~metadata.index.isin(current_metadata.index)]
    all_metadata = pd.concat([current_metadata, metadata], axis=0)

    demixed_samples = pd.read_csv('outputs/aggregate/aggregate_demix.tsv', sep='\t')['Unnamed: 0'].apply(lambda x: x.split('.')[0]).values

    
    # Failed samples will produce variants output but fail in the demixing step
    failed_samples = [file.split('.')[0] for file in os.listdir('outputs/variants') if f'{file.split(".")[0]}.demix.tsv' not in os.listdir('outputs/demix')]

    
    samples_to_run = all_metadata.copy()
    samples_to_run = samples_to_run[~samples_to_run.index.isin(failed_samples)]
    samples_to_run = samples_to_run[~samples_to_run.index.isin(demixed_samples)]
    samples_to_run = samples_to_run[~samples_to_run['ww_surv_target_1_conc'].isna()]
    samples_to_run = samples_to_run[~samples_to_run['ww_population'].str.contains('<')]
    samples_to_run = samples_to_run[~samples_to_run['ww_population'].str.contains('>')]
    samples_to_run = samples_to_run[~samples_to_run['ww_population'].isna()]

    samples_to_run = samples_to_run[samples_to_run['collection_date'] >='2022-04-01']
    samples_to_run = samples_to_run[samples_to_run['collection_date'] <='2023-10-01']

    # 10 samples per collection date
    samples_to_run = samples_to_run.groupby('collection_date').head(10)

    print('All samples: ', len(metadata))
    print('Newly added samples: ', len(new_metadata))
    print('Processed samples: ', len(demixed_samples))
    print('Samples to run: ', len(samples_to_run))

    all_metadata.to_csv('data/all_metadata.csv')
    samples_to_run.to_csv('data/samples_to_run.csv', index=True, header=True)

if __name__ == "__main__":
    main()