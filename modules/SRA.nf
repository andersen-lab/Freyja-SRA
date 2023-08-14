/*
 * Handle SRA data
 */

process GET_ACCESSIONS {
    publishDir "data/SRA/", mode: 'move'
    input:
    file sra_data

    output:
    path "${sra_data.baseName}_accessions.csv"
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    df = pd.read_csv('${sra_data}')
    df['acc'].to_csv('${sra_data.baseName}_accessions.csv', index=False, header=False)
    """ 
}

process FASTERQ_DUMP {
    container 'ncbi/sra-tools'
    input:
    val accession

    output:
    tuple val(accession), path("${accession}_1.fastq"), path("${accession}_2.fastq")

    script:
    """
    #!/bin/sh
    fasterq-dump -e 2 -p ${accession}
    """
}