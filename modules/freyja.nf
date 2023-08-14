process FREYJA_VARIANTS {
    input:
    path input_bam

    output:
    tuple path("${input_bam.baseName}_variants.tsv"), path("${input_bam.baseName}_depths.tsv")

    script:
    """
    freyja variants ${input_bam} --variants ${input_bam.baseName}_variants.tsv --depths ${input_bam.baseName}_depths.tsv --ref ${params.ref}
    """
}

process FREYJA_BOOT {
    input:
    tuple path(variants), path(depths)

    output:
    tuple path("${variants.baseName}_boot.tsv"), path("${depths.baseName}_boot.tsv")

    script:
    """
    freyja boot ${variants} ${depths} --variants ${variants.baseName}_boot.tsv --depths ${depths.baseName}_boot.tsv --output_basename ${variants.baseName}
    """
}