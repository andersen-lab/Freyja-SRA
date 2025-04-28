#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Automated pipeline for Freyja analysis of SARS-CoV-2 wastewater sequencing data
 */

include { SRATOOLS_FASTERQDUMP    } from './modules/nf-core/sratools/fasterqdump/main'
include { SRATOOLS_PREFETCH       } from './modules/nf-core/sratools/prefetch/main'

include {
    CUTADAPT_TRIM;
    MINIMAP2;
    IVAR_TRIM;
} from "./modules/preprocessing.nf"

include {
    FREYJA_VARIANTS;
    FREYJA_DEMIX;
} from "./modules/freyja.nf"

workflow fetch {
    Channel.fromPath(params.sra_metadata)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            meta = [
                id: row.accession.toString(),
                sample_status: row.sample_status.toString(),
                primer_scheme: row.amplicon_PCR_primer_scheme.toString(),
            ]
            [meta, row.accession.toString()]
        }
        .filter { it[0].sample_status == 'to_run' }
        .take(params.num_samples)
        .set { samples_ch }


    samples_ch
        .map { meta, accession -> accession }
        .collectFile(name: "accession_list.txt", storeDir: "./data") { accession -> 
            accession + "\n" 
        }

    SRATOOLS_PREFETCH(samples_ch)
    SRATOOLS_FASTERQDUMP(SRATOOLS_PREFETCH.out.sra)

    SRATOOLS_FASTERQDUMP.out.reads
        .filter { it[0].sample_status == 'to_run' }
        .branch { 
            unknown_primer: it[0].primer_scheme == 'unknown'
            known_primer: it[0].primer_scheme != 'unknown'
        }
        .set { fq_ch }

    process_unknown_primer(fq_ch.unknown_primer)
    process_known_primer(fq_ch.known_primer)

}

workflow process_unknown_primer {
    take: unknown_primer_fastq_ch

    main:
    CUTADAPT_TRIM(unknown_primer_fastq_ch)
    MINIMAP2(CUTADAPT_TRIM.out, params.reference)

    freyja(MINIMAP2.out)
}

workflow process_known_primer {
    take: known_primer_fastq_ch

    main:
    MINIMAP2(known_primer_fastq_ch, params.reference)
    IVAR_TRIM(MINIMAP2.out, params.bedfiles)

    freyja(IVAR_TRIM.out)
}

workflow freyja {
    take:
    bam_ch
    
    main:
    FREYJA_VARIANTS(bam_ch, params.reference)
    FREYJA_DEMIX(FREYJA_VARIANTS.out, params.barcodes)
}

workflow rerun_demix {
    Channel
        .fromFilePairs("${params.variants_dir}/SRR*{variants,depths}.tsv")
        .map { k, v -> tuple(k, v[1], v[0]) }
        .set { variants_ch }

    FREYJA_DEMIX(variants_ch, params.barcodes)
}
