process FREYJA_VARIANTS {
    input:
    path input_bam

    output:
    tuple val(input_bam.baseName), path("${input_bam.baseName}_variants.tsv"), path("${input_bam.baseName}_depths.tsv")

    script:
    """
    samtools sort -o ${input_bam.baseName}.sorted.bam ${input_bam}
    freyja variants ${input_bam}.sorted.bam --variants ${input_bam.baseName}_variants.tsv --depths ${input_bam.baseName}_depths.tsv --ref ${params.ref}
    """
}

process FREYJA_DEMIX {
    input:
    tuple val(sample_id), path(variants), path(depths)

    output:
    path "${sample_id}_demix.tsv"

    script:
    """
    freyja demix ${variants} ${depths} --variants ${variants} --depths ${depths} --output ${sample_id}_demix.tsv
    """
}