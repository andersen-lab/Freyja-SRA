process AGGREGATE_DEMIX {
    publishDir "${params.output}/aggregate", mode: 'copy'

    input:
    path baseDir
    path demix_ch

    output:
    path "aggregate_demix_new.json"

    script:
    """
    python ${baseDir}/scripts/aggregate_demix.py ${baseDir}
    """
}
