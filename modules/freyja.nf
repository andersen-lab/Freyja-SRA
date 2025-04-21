process FREYJA_VARIANTS {
    publishDir "${params.output}/variants", mode: 'copy'

    input:
    tuple val(meta), path(bam)
    path reference

    output:
    tuple val(meta), path("${meta.id}.variants.tsv"), path("${meta.id}.depths.tsv")

    script:
    """
    freyja variants ${bam} --variants ${meta.id}.variants.tsv --depths ${meta.id}.depths.tsv --ref ${reference}
    """
}

process FREYJA_DEMIX {
    publishDir "${params.output}/demix", mode: 'copy'
    errorStrategy 'ignore' 

    input:
    tuple val(meta), path(variants), path(depths)
    path barcodes

    output:
    path "*.demix.tsv"

    script:
    """
    freyja demix \\
        $variants \\
        $depths \\
        --barcodes $barcodes \\
        --output ${meta.id}.demix.tsv
    """
}