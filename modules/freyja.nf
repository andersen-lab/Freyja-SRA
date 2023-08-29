process FREYJA_VARIANTS {
    input:
    path input_bam
    path ref

    output:
    tuple val(input_bam.baseName), path("${input_bam.baseName}_variants.tsv"), path("${input_bam.baseName}_depths.tsv")

    script:
    """
    samtools sort -o ${input_bam.baseName}.sorted.bam ${input_bam}
    freyja variants ${input_bam.baseName}.sorted.bam --variants ${input_bam.baseName}_variants.tsv --depths ${input_bam.baseName}_depths.tsv --ref ${ref}
    """
}

process FREYJA_DEMIX {
    publishDir params.output, mode: 'copy'

    input:
    tuple val(sample_id), path(variants), path(depths)

    output:
    path "${sample_id}.demix.tsv"

    script:
    """
    freyja demix ${variants} ${depths} --output ${sample_id}.demix.tsv
    """
}

process FREYJA_COVARIANTS {
    input:
    input_bam
    path ref
    path annot

    output:
    path "${input_bam.baseName}.covariants.tsv"

    script:
    """
    samtools sort -o ${input_bam.baseName}.sorted.bam ${input_bam}
    freyja covariants ${input_bam}.sorted.bam --variants ${input_bam.baseName}.covariants.tsv --ref ${ref} --annot ${annot}
    """
}