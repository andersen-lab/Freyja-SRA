/*
 * Handle SRA data
 */

process GET_ACCESSIONS {
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
    input:
    val accession

    output:
    tuple val(accession), path("${accession}.fastq.gz")

    script:
    """
    fasterq-dump --split-files ${accession}
    """    

}