#!/usr/bin/env nextflow

/*
 * Automated pipeline for Freyja analysis of SRA data
 */

ref = file(params.ref)
bedfiles = file(params.bedfiles)
baseDir = file("$baseDir")
annot = file(params.annot)

// Import modules
include {
    GET_NCBI_METADATA;
    GET_ACCESSIONS;
    GET_AMPLICON_SCHEME;
    FASTERQ_DUMP;
} from "./modules/sra.nf"

include {
    MINIMAP2;
    SAMTOOLS_1;
    SAMTOOLS_2;
    IVAR_TRIM;
} from "./modules/preprocessing.nf"

include {
    FREYJA_VARIANTS;
    FREYJA_DEMIX;
    FREYJA_COVARIANTS;
    // AGGREGATE_VARIANTS
    // AGGREGATE_DEMIX;
    // AGGREGATE_COVARIANTS;
} from "./modules/freyja.nf"

workflow fetch {

    GET_NCBI_METADATA(baseDir)
        .set { input_ch }

    GET_ACCESSIONS(input_ch)
        .splitCsv()
        .map { line -> line.join('') }
        .take(1)
        .set { acc_ch }

    GET_AMPLICON_SCHEME(acc_ch, input_ch)
        .map { it.text }
        .set { primer_scheme_ch }

    FASTERQ_DUMP(acc_ch)
    MINIMAP2(FASTERQ_DUMP.out, ref)
    SAMTOOLS_1(MINIMAP2.out)
    IVAR_TRIM(SAMTOOLS_1.out, primer_scheme_ch, bedfiles)
    SAMTOOLS_2(IVAR_TRIM.out)

    FREYJA_VARIANTS(SAMTOOLS_2.out, ref)
        .collect()
        .set { variants_ch }

    FREYJA_DEMIX(FREYJA_VARIANTS.out)
        .collect()
        .set { demix_ch }

    FREYJA_COVARIANTS(SAMTOOLS_2.out, ref, annot)
        .collect()
        .set { covariants_ch }

    // AGGREGATE_VARIANTS()
    // AGGREGATE_DEMIX()
    // AGGREGATE_COVARIANTS()
}

workflow rerun_demix {
    Channel
        .fromFilePairs("${params.output}/variants/SRR*{variants,depths}.tsv")
        .set { variants_ch }

    FREYJA_DEMIX(variants_ch)
        .collect()
        .set { demix_ch }

    FREYJA_AGGREGATE(demix_ch, baseDir)
    FREYJA_PLOT(FREYJA_AGGREGATE.out)
}
