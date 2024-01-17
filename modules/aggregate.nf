process AGGREGATE_DEMIX {
    publishDir "${params.output}/aggregate", mode: 'copy'

    script:
    """
    python !{baseDir}/scripts/aggregate_demix.py 
    """
}
