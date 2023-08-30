/*
 * Handle SRA data
 */

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
    val sample_id
    file sra_data

    output:
    path 'primer_scheme.txt'

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    
    df = pd.read_csv('${sra_data}',index_col=0)
    scheme = df.loc['${sample_id}','amplicon_PCR_primer_scheme']
    if 'QIAseq' in scheme or 'v3' in scheme:
        primer_scheme = 'ARTICv3'
    elif 'V5.3' in scheme:
        primer_scheme = 'ARTICv5.3.2'
    elif 'V4.1' in scheme:
        primer_scheme = 'ARTICv4.1'
    elif 'SNAP' in scheme:
        primer_scheme = 'snap_primers'
    else:
        primer_scheme = 'default'

    with open('primer_scheme.txt', 'w') as f:
        f.write(primer_scheme)
    """
}


process FASTERQ_DUMP {
    container { params.profile == "docker" ? "ncbi/sra-tools" : "docker://ncbi/sra-tools" }
    input:
    val accession

    output:
    tuple val(accession), path("${accession}_1.fastq"), path("${accession}_2.fastq")

    script:
    """
    #!/bin/sh
    fasterq-dump --split-files ${accession} 
    """
}