/*
 * Handle SRA fetching and preprocessing
 */

process GET_NCBI_METADATA {
    input:
    //path acc_in_db
    path baseDir

    output: 
    path "wastewater_ncbi.csv"

    script:
    """
    #!/usr/bin/env python3

    import os
    import pandas as pd
    from Bio import Entrez
    from Bio import SeqIO

    Entrez.email = "jolevy@scripps.edu"
    # handle = Entrez.esearch(db="sra", idtype='acc', retmax=10000, term="((wastewater metagenome[Organism] OR (wastewater metagenome[Organism] OR wastewater metagenome[All Fields])) AND (Severe acute respiratory syndrome coronavirus 2[Organism] OR (Severe acute respiratory syndrome coronavirus 2[Organism] OR sars-cov-2[All Fields]))) AND strategy wgs[Properties]")#, idtype="acc")
    handle = Entrez.esearch(db="sra", idtype='acc', retmax=2000,
                            sort='recently_added',
                            term="((wastewater metagenome[Organism] OR wastewater metagenome[All Fields]) AND SARS-CoV-2))") 
    record = Entrez.read(handle)
    handle.close()

    handle = Entrez.efetch(db="sra", id=record['IdList'], rettype="gb",retmode='text')
    string= handle.read()
    handle.close()

    returned_meta=str(string,'UTF-8')

    with open("NCBI_metadata.xml", "w") as f:
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
            asdfa
        else:
            sampID = sampID[0]
        ## write to dictionary form
        dictVals = {vals[i].replace(' ','_'):vals[i+1] for i in range(0,len(vals),2)}
        for i in range(0,len(seq_meta),2):
            dictVals[seq_meta[i].replace(' ','')] = seq_meta[i+1]
        dictVals['experiment_id'] = sampID
        dictVals['SRA_id'] = root0[0].attrib['accession']
        allDictVals[sampID] =dictVals

    // acc_in_db = pd.read_csv('${acc_in_db}',header=None).index
    // print('acc in db', acc_in_db)

    df = pd.DataFrame(allDictVals).T

    df.columns = df.columns.str.replace(' ','_')
    df = df[df['collection_date'].str.startswith('20')]
    df['collection_date'] = pd.to_datetime(df['collection_date'].apply(lambda x: x.split('/')[0] if '/' in x else x))
    df = df.sort_values(by='collection_date',ascending=False)
    ## add last 


    df = df[df['collection_date'] >='2023-07-20']

    
    df.to_csv('wastewater_ncbi.csv')
    """
}

process GET_ACCESSIONS {
    input:
    file sra_data

    output:
    path "acc_list.csv"
    
    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    
    df = pd.read_csv('${sra_data}',index_col=0)
    pd.Series(df.index).to_csv('acc_list.csv', index=False, header=False)
    """ 
}

process GET_AMPLICON_SCHEME {

    input:
    val sra_accession
    file sra_data    

    output:
    val sra_accession
    path('primer_scheme.txt')


    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    
    df = pd.read_csv('${sra_data}',index_col=0)
    scheme = df.loc['${sra_accession}','amplicon_PCR_primer_scheme']
    
    scheme = str(scheme)

    if scheme == 'nan':
        primer_scheme = 'unknown'
    elif 'QIAseq' in scheme or 'v3' in scheme:
        primer_scheme = 'ARTICv3'
    elif 'V5.3' in scheme:
        primer_scheme = 'ARTICv5.3.2'
    elif 'V4.1' or 'v4.1' in scheme:
        primer_scheme = 'ARTICv4.1'
    elif 'SNAP' in scheme:
        primer_scheme = 'snap_primers'
    else:
        primer_scheme = 'unknown'

    with open('primer_scheme.txt', 'w') as f:
            f.write(primer_scheme)
    """
}


process FASTERQ_DUMP {
    container { params.profile == "docker" ? "ncbi/sra-tools" : "docker://ncbi/sra-tools" }
    input:
    val accession
    path primer_scheme

    output:
    tuple val(accession), path("${accession}_1.fastq"), path("${accession}_2.fastq"), path(primer_scheme)

    script:
    """
    #!/bin/sh
    fasterq-dump --split-files ${accession} 
    """
}