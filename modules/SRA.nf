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

process FETCH_NGS {
    input:
    path accession_list

    output:
    path "fastq_dir"

    script:
    """
    nextflow run nf-core/fetchngs --input ${accession_list} --outdir fastq_dir
    """ 
}