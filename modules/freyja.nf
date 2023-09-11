process FREYJA_VARIANTS {
    publishDir "${params.output}/variants", mode: 'copy'

    input:
    val sra_accession
    path input_bam
    path bam_index
    path ref

    output:
    tuple val(sra_accession), path("${sra_accession}.variants.tsv"), path("${sra_accession}.depths.tsv")

    script:
    """
    freyja variants ${input_bam} --variants ${sra_accession}.variants.tsv --depths ${sra_accession}.depths.tsv --ref ${ref}
    """
}

process FREYJA_DEMIX {
    publishDir "${params.output}/demix", mode: 'copy'

    input:
    tuple val(sample_id), path(variants), path(depths)
    val eps
    val depthCutoff

    output:
    path "${sample_id}.demix.tsv"

    script:
    """
    freyja demix ${variants} ${depths} --eps ${eps} --depthcutoff ${depthCutoff} --output ${sample_id}.demix.tsv
    """
}

process FREYJA_COVARIANTS {
    publishDir "${params.output}/covariants", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
    val sra_accession
    path input_bam
    path bam_index
    path ref
    path annot

    output:
    path "${sra_accession}.covariants.tsv"

    script:
    """
    freyja covariants ${input_bam} ${params.min_site} ${params.max_site} --output ${sra_accession}.covariants.tsv --ref-genome ${ref} --gff-file ${annot}
    """
}

process AGGREGATE_VARIANTS {
    publishDir "${params.output}/aggregate", mode: 'copy'
    input:
    val variants_outputs
    path baseDir

    output:
    path "aggregate_variants.tsv"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd 
    import os

    paths_string = "${variants_outputs}"
    if not paths_string.startswith('['):
        paths_string = f'[{paths_string}]'
    
    paths_list = paths_string[1:-1].split(", ")

    ct_thresh = 20
    j=0
    for var_path in paths_list:
        df = pd.read_csv(var_path,sep='\t')
        df = df[df['ALT_DP']>=ct_thresh]
        #drop frame shifts
        sname=var_path.split('/')[-1].split('.')[0]
        df = df[df['ALT'].apply(lambda x: True if (('+' not in x and '-' not in x) or ((len(x)-1)%3==0)) else False)]
        if len(df)==0:
            print('No data in ',var_path)
            continue
        df['mutName'] = df['REF'] + df['POS'].astype(str) + df['ALT']
        df['sample'] = sname
        df = df[['mutName','sample','ALT_FREQ','ALT_DP']]
        if j==0:
            dfAll = df.copy()
        else:
            dfAll = pd.concat((dfAll,df), axis=0)
        j+=1

    dfAll = dfAll.set_index(['mutName','sample'])
    
    dfAll.to_csv('aggregate_variants.tsv',sep='\t')
    """
}

process AGGREGATE_DEMIX {
    publishDir "${params.output}/aggregate", mode: 'copy'
    input:
    val demix_outputs
    path baseDir

    output:
    path "aggregate_demix.tsv"

    script:
    """
    #!/usr/bin/env python3
    import subprocess

    paths_string = "${demix_outputs}"
    paths_list = paths_string[1:-1].split(", ")

    subprocess.run(["mkdir", "aggregate_dir"])
    for file in paths_list:
       subprocess.run(["cp", file, "aggregate_dir"])

    subprocess.run(["freyja", "aggregate", "aggregate_dir/", "--output", "aggregate_demix.tsv"])
    """
}

process AGGREGATE_COVARIANTS {
    publishDir "${params.output}/aggregate", mode: 'copy'
    input:
    val covariants_outputs
    path baseDir

    output:
    path "aggregate_covariants.tsv"

    script:
    """
    #!/usr/bin/env python3
    import subprocess
    import pandas as pd

    paths_string = "${covariants_outputs}"
    paths_list = paths_string[1:-1].split(", ")

    agg_df = pd.DataFrame(columns=['Covariants','Sample', 'Count', 'Max_count', 'Freq', 'Coverage_start', 'Coverage_end'])

    for covar_path in paths_list:
        df = pd.read_csv(covar_path,sep='\t')
        sample_name = covar_path.split('/')[-1].split('.')[0]
        df['Sample'] = sample_name
        agg_df = pd.concat((agg_df,df), axis=0)

    agg_df = agg_df.set_index(['Covariants','Sample'])
    agg_df.to_csv('aggregate_covariants.tsv',sep='\t')
    """
}