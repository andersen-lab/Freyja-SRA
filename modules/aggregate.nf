process AGGREGATE_DEMIX {
    publishDir "${params.output}/aggregate", mode: 'copy'

    input:
    path demix_ch
    path baseDir

    output:
    path "aggregate_demix_new.json"

    script:
    """
    python ${baseDir}/scripts/aggregate_demix.py ${baseDir}
    """
}

process AGGREGATE_VARIANTS {
    publishDir "${params.output}/aggregate", mode: 'copy'

    input:
    path variants_ch
    path baseDir

    output:
    path "aggregate_variants_by_acc_new.json"
    path "aggregate_variants_by_mut_new.json"

    script:
    """
    python ${baseDir}/scripts/aggregate_variants.py ${baseDir}
    """
}