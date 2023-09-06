/*
 *  Alignment and primer trimming
 */

process BBDUK_TRIM {
    input:
    tuple val(sample_id), path(read1), path(read2), val(primer_scheme)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_1.fastq"), path("${sample_id}_trimmed_2.fastq"), val(primer_scheme)

    script:
    """
    bbduk.sh in=${read1} in2=${read2} out=${sample_id}_trimmed_1.fastq out2=${sample_id}_trimmed_2.fastq ftl=10 ftr=139
    """
}

process MINIMAP2 {
    input:
    tuple val(sample_id), path(read1), path(read2), val(primer_scheme)
    path ref

    output:
    path "${sample_id}.bam"

    script:
    """
    minimap2 -ax sr ${ref} ${read1} ${read2} | samtools view -bS - > ${sample_id}.bam
    """
}

process MINIMAP2_unknown_primer {
    input:
    tuple val(sample_id), path(read1), path(read2), val(primer_scheme)
    path ref

    output:
    path "${sample_id}.bam"

    script:
    """
    minimap2 -ax sr ${ref} ${read1} ${read2} | samtools view -bS - > ${sample_id}.bam
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

process SAMTOOLS_1_unknown_primer {
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
