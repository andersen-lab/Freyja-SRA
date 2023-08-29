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
    samtools sort -o ${sample_id}.bam ${sample_id}.bam
    """
}


process IVAR_TRIM {
    errorStrategy 'ignore'
    
    input:
    path input_bam
    val primer_scheme
    path bedfiles

    output:
    path "${input_bam.baseName}.trimmed.bam"

    script:
    """
    samtools sort -o ${input_bam.baseName}.sorted.bam ${input_bam}
    ivar trim -x 4 -e -m 80 -i ${input_bam.baseName}.sorted.bam -b ${bedfiles}/${primer_scheme}.bed -p ${input_bam.baseName}.trimmed.bam
    """
}