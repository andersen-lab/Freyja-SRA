process AGGREGATE_DEMIX {
    publishDir "${params.output}/aggregate", mode: 'copy'
    input:
    val baseDir

    script:
    """
    python !{baseDir}/scripts/aggregate_demix.py 
    """
}
