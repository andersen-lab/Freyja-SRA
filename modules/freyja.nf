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
    publishDir "${params.output}/demix", mode: 'copy'

    input:
    val demix_outputs

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