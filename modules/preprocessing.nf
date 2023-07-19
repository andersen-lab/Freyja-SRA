/*
 *  Preprocessing steps
 */

process MINIMAP2 {
    input:
    tuple val(sample_id), path(reads)
    path ref

    output:
    path "${sample_id}.bam"

    script:
    """
    minimap2 -ax sr ${ref} ${reads[0]} ${reads[1]} | samtools view -bS - > ${sample_id}.bam
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