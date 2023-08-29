process FREYJA_VARIANTS {
    publishDir "${params.output}/variants", mode: 'copy'

    input:
    path input_bam
    path bam_index
    path ref

    output:
    tuple val(input_bam.baseName), path("${input_bam.baseName}.variants.tsv"), path("${input_bam.baseName}.depths.tsv")

    script:
    """
    samtools sort -o ${input_bam.baseName}.sorted.bam ${input_bam}
    freyja variants ${input_bam.baseName}.sorted.bam --variants ${input_bam.baseName}.variants.tsv --depths ${input_bam.baseName}.depths.tsv --ref ${ref}
    """
}

process FREYJA_DEMIX {
    publishDir "${params.output}/demix", mode: 'copy'

    input:
    tuple val(sample_id), path(variants_output)

    output:
    path "${sample_id}.demix.tsv"

    script:
    """
    freyja demix ${variants_output[1]} ${variants_output[0]} --output ${sample_id}.demix.tsv
    """
}

process FREYJA_AGGREGATE {
    publishDir params.output, mode: 'copy'

    input:
    path demix

    output:
    path "demix/aggregate.tsv"

    script:
    """
    freyja aggregate ${demix} --output demix/aggregate.tsv
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