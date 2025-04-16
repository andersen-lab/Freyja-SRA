#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Automated pipeline for Freyja analysis of SARS-CoV-2 wastewater sequencing data
 */

sra_metadata = file(params.sra_metadata)
ref = file(params.ref)
bedfiles = file(params.bedfiles)
baseDir = file("$baseDir")
barcodes = file("${baseDir}/data/usher_barcodes.feather")

// Import modules

include { SRATOOLS_FASTERQDUMP    } from './modules/nf-core/sratools/fasterqdump/main'
include { SRATOOLS_PREFETCH       } from './modules/nf-core/sratools/prefetch/main'

include {
    CUTADAPT_TRIM;
    MINIMAP2;
    MINIMAP2_UNKNOWN_PRIMER;
    SAMTOOLS_1;
    SAMTOOLS_2;
    IVAR_TRIM;
} from "./modules/preprocessing.nf"

include {
    FREYJA_VARIANTS;
    FREYJA_DEMIX;
} from "./modules/freyja.nf"

workflow fetch {
    Channel.fromPath(sra_metadata)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            meta = [
                id: row.accession.toString(),
                sample_status: row.sample_status.toString(),
                primer_scheme: row.amplicon_PCR_primer_scheme.toString(),
            ]
            [meta, row.accession.toString()]
        }
        .take(5)
        .filter { it[0].sample_status == 'to_run' }
        .set { samples_ch }

    samples_ch.view { meta, acc -> "Sample: ${meta.id}, status: ${meta.sample_status}, primer scheme: ${meta.primer_scheme}" }

    SRATOOLS_PREFETCH(samples_ch)
    SRATOOLS_FASTERQDUMP(SRATOOLS_PREFETCH.out.sra)

    SRATOOLS_FASTERQDUMP.out.reads
        .branch { meta, fastq ->
            known_primer: meta.primer_scheme != 'unknown'
            unknown_primer: meta.primer_scheme == 'unknown'
        }
        .set { fq_ch }
    
    // View the branched channels with labels
    fq_ch.known_primer.view { meta, fastq -> "Known primer: ${meta.id}, scheme: ${meta.primer_scheme}" }
    fq_ch.unknown_primer.view { meta, fastq -> "Unknown primer: ${meta.id}, scheme: ${meta.primer_scheme}" }

    process_unknown_primer(fq_ch.unknown_primer)
    process_known_primer(fq_ch.known_primer)
}

workflow process_unknown_primer {
    take: unknown_primer_fastq_ch

    main:
    CUTADAPT_TRIM(unknown_primer_fastq_ch)
    MINIMAP2_UNKNOWN_PRIMER(CUTADAPT_TRIM.out, ref)
    SAMTOOLS_1(MINIMAP2_UNKNOWN_PRIMER.out)

    freyja(SAMTOOLS_1.out)
}

workflow process_known_primer {
    take: known_primer_fastq_ch

    main:
    known_primer_fastq_ch
        .map { it[1].text }
        .set { primer_ch }
    MINIMAP2(known_primer_fastq_ch, ref)
    SAMTOOLS_1(MINIMAP2.out)
    IVAR_TRIM(SAMTOOLS_1.out, primer_ch, bedfiles)
    SAMTOOLS_2(IVAR_TRIM.out)

    freyja(SAMTOOLS_2.out)
}

workflow freyja {
    take:
    sra_accession
    input_bam
    bam_index
    
    main:
    FREYJA_VARIANTS(sra_accession, input_bam, bam_index, ref)
        .collect()
        .map { it[1] }
        .set { variants_ch }

    FREYJA_DEMIX(FREYJA_VARIANTS.out, params.eps, barcodes)
        .collect()
        .set { demix_ch }   
}


workflow rerun_demix {
    Channel
        .fromFilePairs("${params.variants_dir}/SRR*{variants,depths}.tsv")
        .map { k, v -> tuple(k, v[1], v[0]) }
        .set { variants_ch }

    FREYJA_DEMIX(variants_ch, params.eps, barcodes)
        .collect()
        .set { demix_ch }
}
