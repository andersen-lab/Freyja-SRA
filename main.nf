#!/usr/bin/env nextflow

/*
 * Automated pipeline for Freyja analysis of SRA data
 */

ref = file(params.ref)
bedfiles = file(params.bedfiles)
baseDir = file("$baseDir")
annot = file(params.annot)

host = params.es_host
user = params.es_user
password = params.es_pass

// Import modules
include {
    GET_NCBI_METADATA;
    GET_ACCESSIONS;
    GET_AMPLICON_SCHEME;
    FASTERQ_DUMP;
} from "./modules/sra.nf"

include {
    BBDUK_TRIM;
    MINIMAP2;
    MINIMAP2_unknown_primer;
    SAMTOOLS_1;
    SAMTOOLS_1_unknown_primer;
    SAMTOOLS_2;
    IVAR_TRIM;
} from "./modules/preprocessing.nf"

include {
    FREYJA_VARIANTS;
    FREYJA_DEMIX;
    FREYJA_COVARIANTS;
    AGGREGATE_VARIANTS
    AGGREGATE_DEMIX;
    // AGGREGATE_COVARIANTS;
} from "./modules/freyja.nf"

include {
    CHECK_SAMPLES_IN_ES;
} from "./modules/elasticsearch.nf"

params.input = "./output/metadata/wastewater_ncbi.csv"
input = file(params.input)
workflow fetch {
    // Channel
    //     .fromPath(input)
    //     .set { input_ch }

    GET_NCBI_METADATA(baseDir)
        .set { input_ch }

    GET_ACCESSIONS(input_ch)
        .splitCsv()
        .map { line -> line.join('') }
        .take(10)
        .set { acc_ch }

    GET_AMPLICON_SCHEME(acc_ch, input_ch)
        .set { primer_scheme_ch }

    FASTERQ_DUMP(primer_scheme_ch)
        .branch {
            unknown_primer: it[3].text == 'unknown'
            known_primer: it[3].text != 'unknown'
        }
        .set { fq_ch }

    BBDUK_TRIM(fq_ch.unknown_primer)
    MINIMAP2_unknown_primer(BBDUK_TRIM.out, ref)
        .set { unknown_primer_bam_ch }
    SAMTOOLS_1_unknown_primer(unknown_primer_bam_ch)
        .set { unknown_primer_sorted_trimmed_ch }
    
    MINIMAP2(fq_ch.known_primer, ref)
    SAMTOOLS_1(MINIMAP2.out)
    IVAR_TRIM(SAMTOOLS_1.out, bedfiles)
    SAMTOOLS_2(IVAR_TRIM.out)
        .set { known_primer_sorted_trimmed_ch }

    Channel
        .from(unknown_primer_sorted_trimmed_ch)
        .concat(known_primer_sorted_trimmed_ch)
        .set { sorted_trimmed_ch }

    FREYJA_VARIANTS(known_primer_sorted_trimmed_ch, ref)
        .collect()
        .map { it[1] }
        .set { variants_ch }

    // FREYJA_DEMIX(FREYJA_VARIANTS.out)
    //     .collect()
    //     .set { demix_ch }

    // FREYJA_COVARIANTS(SAMTOOLS_2.out, ref, annot)
    //     .collect()
    //     .set { covariants_ch }

    // AGGREGATE_VARIANTS(variants_ch, baseDir)
    // AGGREGATE_DEMIX(demix_ch, baseDir)
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
