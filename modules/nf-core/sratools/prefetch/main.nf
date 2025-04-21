process SRATOOLS_PREFETCH {
    errorStrategy 'ignore'
    tag "$id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), val(id)

    output:
    tuple val(meta), path(id), emit: sra
    path 'versions.yml'      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    args = '5 1 100'  // <num retries> <base delay in seconds> <max delay in seconds>
    template 'retry_with_backoff.sh'
}