import pandas as pd

import argparse
import os
import pandas as pd
from Bio import Entrez
from Bio import SeqIO
import xml.etree.ElementTree as ET
import shortuuid



Entrez.email = "jolevy@scripps.edu"

finished_samples = pd.read_csv('outputs/aggregate/aggregate_demix.tsv', sep='\t')['Unnamed: 0'].apply(lambda x: x.split('.')[0]).values
# TODO: Handle failed samples (samples that have variants output but not demix output)

# Get the SRA ids of the finished samples
sra_ids = [file.split('.')[0] for file in os.listdir('outputs/variants') if f'{file.split(".")[0]}.demix.tsv' in os.listdir('outputs/demix')]

all_metadata = pd.read_csv('data/all_metadata.csv', index_col=0)

# Get the SRA ids of the finished samples that are missing from the metadata
missing_sra_ids = set(sra_ids) - set(all_metadata.index)

print('Missing SRA ids: ', len(missing_sra_ids))

missing_data_df = pd.DataFrame()

for sra_id in missing_sra_ids:
    print(sra_id)

    Entrez.email = "jolevy@scripps.edu"
    handle = Entrez.esearch(db="sra", idtype='acc', retmax=5,
                            sort='recently_added',
                            term=sra_id) 
    try:
        record = Entrez.read(handle)
    except RuntimeError:
        continue
    
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

    # Remove samples where the collection date is not of the form 20XX-XX-XX
    metadata  = metadata[metadata['collection_date'].str.contains('20[0-9]{2}-[0-9]{2}-[0-9]{2}')]
    metadata['collection_date'] = pd.to_datetime(metadata['collection_date'].apply(lambda x: x.split('/')[0] if '/' in x else x))
    metadata = metadata.sort_values(by='collection_date',ascending=False)

    # Filter to USA samples
    metadata = metadata[metadata['geo_loc_name'].str.contains('USA')]
    
    # For samples with no site id, hash the location and population to create a unique id
    merged = metadata['geo_loc_name']+metadata['ww_population'].fillna('').astype(str)
    merged = merged.apply(lambda x:shortuuid.uuid(x)[0:12])

    if 'collection_site_id' not in metadata.columns:
        metadata['site_id'] = merged
    else:
        metadata['site_id'] = metadata['collection_site_id']
    
    missing_data_df = pd.concat([missing_data_df, metadata], axis=0)
    
    
all_metadata = pd.concat([all_metadata, missing_data_df], axis=0)
all_metadata.to_csv('data/all_metadata.csv')
