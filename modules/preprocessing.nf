/*
 *  Alignment and primer trimming
 */

process CUTADAPT_TRIM {
    input:
    tuple val(sample_id), path(primer_scheme), path(reads)

    output:
    tuple val(sample_id), val(primer_scheme), path("*_trimmed.fastq")

    script:
    def read1 = reads.first()
    def read2 = reads.last()
    """
    if ${read1} == ${read2}
    then
        cutadapt -l +30 -l -30 -o ${sample_id}_trimmed.fastq ${read1}
    else
        cutadapt -l +30 -l -30 -o ${sample_id}_1_trimmed.fastq -p ${sample_id}_2_trimmed.fastq ${read1} ${read2}
    fi
    """
}

process MINIMAP2 {
    input:
    tuple val(sample_id), val(primer_scheme), path(reads)
    path ref

    output:
    path "${sample_id}.bam"

    script:
    """
    minimap2 -ax sr ${ref} ${reads} | samtools view -bS - > ${sample_id}.bam   
    """
}

process SAMTOOLS_1 {
    input:
    path bamfile

    output:
    val bamfile.baseName
    path "${bamfile.baseName}.sorted.bam"
    path "${bamfile.baseName}.sorted.bam.bai"

    script:
    """
    samtools sort -o ${bamfile.baseName}.sorted.bam ${bamfile}
    samtools index ${bamfile.baseName}.sorted.bam
    """
}

process IVAR_TRIM {
    input:
    val sra_accession
    path sorted_bam
    path bam_index
    val primer_scheme
    path bedfiles

    output:
    val sra_accession
    path "${sra_accession}.trimmed.bam"

    script:
    """
    ivar trim -x 4 -e -m 80 -i ${sorted_bam} -b ${bedfiles}/${primer_scheme}.bed -p ${sra_accession}.trimmed.bam
    """
}

process SAMTOOLS_2 {
    input:
    val sra_accession
    path bamfile

    output:
    val sra_accession
    path "${sra_accession}.sorted.bam"
    path "${sra_accession}.sorted.bam.bai"

    script:
    """
    samtools sort -o ${sra_accession}.sorted.bam ${bamfile}
    samtools index ${sra_accession}.sorted.bam
    """
}
