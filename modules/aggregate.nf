process AGGREGATE_DEMIX {
    publishDir "${params.output}/aggregate", mode: 'copy'
    input:
    path baseDir
    path demix_ch

    script:
    """
    python ${baseDir}/scripts/aggregate_demix.py ${baseDir}/outputs/demix
    """
}
