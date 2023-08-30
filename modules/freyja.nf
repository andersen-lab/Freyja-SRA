process FREYJA_VARIANTS {
    publishDir "${params.output}/variants", mode: 'copy'

    input:
    val sra_accession
    path input_bam
    path bam_index
    path ref

    output:
    tuple val(sra_accession), path("${sra_accession}.variants.tsv"), path("${sra_accession}.depths.tsv")

    script:
    """
    freyja variants ${input_bam} --variants ${sra_accession}.variants.tsv --depths ${sra_accession}.depths.tsv --ref ${ref}
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
    publishDir "${params.output}/aggregate", mode: 'copy'

    input:
    val demix_outputs
    path baseDir

    output:
    path "aggregate.tsv"

    script:
    """
    #!/usr/bin/env python3
    import subprocess

    paths_string = "${demix_outputs}"
    paths_list = paths_string[1:-1].split(", ")

    subprocess.run(["mkdir", "aggregate_dir"])
    for file in paths_list:
       subprocess.run(["cp", file, "aggregate_dir"])

    subprocess.run(["freyja", "aggregate", "aggregate_dir/", "--output", "aggregate.tsv"])
    """
}

process FREYJA_PLOT {
    publishDir "${params.output}/plot", mode: 'copy'

    input:
    path aggregate_output

    output:
    path "mix_plot.pdf"

    script:
    """
    freyja plot ${aggregate_output} --output mix_plot.pdf
    """
}
process FREYJA_COVARIANTS {
    publishDir "${params.output}/covariants", mode: 'copy'

    input:
    val sra_accession
    path input_bam
    path bam_index
    path ref
    path annot

    output:
    path "${sra_accession}.covariants.tsv"

    script:
    """
    freyja covariants ${input_bam} ${params.min_site} ${params.max_site} --output ${sra_accession}.covariants.tsv --ref-genome ${ref} --gff-file ${annot}
    """
}