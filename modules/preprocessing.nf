/*
 *  NGS processing
 */

process MINIMAP2 {
    input:
    tuple val(sample_id), path(read1), path(read2)
    path ref

    output:
    path "${sample_id}.bam"

    script:
    """
    minimap2 -ax sr ${ref} ${read1} ${read2} | samtools view -bS - > ${sample_id}.bam
    """
}

process GET_PRIMER_BED {
    input:
    val(sample_id)

    output:
    path "${sample_id}_primer.bed"

    script:
    """
    #!/usr/bin/env python3

    from get_ncbi_metadata import *;pull_sample_bed('${sample}')
    """
}
process IVAR_TRIM {
    input:
    path input_bam
    path primer_bed

    output:
    path "${input_bam.baseName}.trimmed.bam"

    script:
    """
    samtools sort -o ${input_bam.baseName}.sorted.bam ${input_bam}
    samtools index ${input_bam.baseName}.sorted.bam
    ivar trim -x 5 -e -i ${input_bam.baseName}.sorted.bam -b ${primer_bed} -p ${input_bam.baseName}.trimmed.bam
    """
}