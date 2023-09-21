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
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample_id), path(variants), path(depths)
    val eps
    val depthCutoff

    output:
    path "${sample_id}.demix.tsv"

    script:
    """
    freyja demix ${variants} ${depths} --eps ${eps} --depthcutoff ${depthCutoff + (5*task.attempt)} --output ${sample_id}.demix.tsv
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