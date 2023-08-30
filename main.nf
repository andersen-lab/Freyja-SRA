#!/usr/bin/env nextflow

/*
 * Automated pipeline for Freyja analysis of SRA data
 */


// SARS-CoV-2 by default, but can be changed to any other pathogen.
params.ref = "$baseDir/data/NC_045512_Hu-1.fasta"
params.bedfiles = "$baseDir/data/bedfiles"
params.workDir = "$baseDir/work"
params.output = "$baseDir/output"

// Freyja covariants specific parameters
params.annot = "$baseDir/data/NC_045512_Hu-1.gff"
params.min_site = 21563
params.max_site = 25384

ref = file(params.ref)
bedfiles = file(params.bedfiles)
workDir = file(params.workDir)
annot = file(params.annot)

// Import modules
include { 
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
    FREYJA_AGGREGATE;
    FREYJA_PLOT;
} from "./modules/freyja.nf"

workflow preprocessing {
    Channel
        .fromPath(params.input)
        .set { input_ch } 

    GET_ACCESSIONS(input_ch)
        .splitCsv()
        .map { line -> line.join('') }
        .take(100)
        .set { acc_ch }

    GET_AMPLICON_SCHEME(acc_ch, input_ch.first())
        .map { it.text }
        .set { primer_scheme_ch }

    FASTERQ_DUMP(acc_ch)
    MINIMAP2(FASTERQ_DUMP.out, ref)
    SAMTOOLS_1(MINIMAP2.out)
    IVAR_TRIM(SAMTOOLS_1.out, primer_scheme_ch, bedfiles)
    SAMTOOLS_2(IVAR_TRIM.out)
    FREYJA_VARIANTS(SAMTOOLS_2.out, ref)
}

workflow demix {
    Channel
        .fromFilePairs("${params.output}/variants/SRR*{variants,depths}.tsv")
        .set { variants_ch }

    FREYJA_DEMIX(variants_ch)
        .collect()
        .set { demix_ch }

    FREYJA_AGGREGATE(demix_ch, workDir)
    FREYJA_PLOT(FREYJA_AGGREGATE.out)
}
