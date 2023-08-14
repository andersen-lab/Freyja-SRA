#!/usr/bin/env nextflow

/*
 * Automated pipeline for Freyja analysis of SRA data
 */


// SARS-CoV-2 by default, but can be changed to any other pathogen.
params.ref = "$baseDir/data/preprocessing/NC_045512_Hu-1.fasta"
ref = file(params.ref)

// Primer file for iVar trimming
params.primer_bed = "$baseDir/data/preprocessing/nCov-2019_v3.primer.bed"
primer_bed = file(params.primer_bed)

// Import modules
include { 
    GET_ACCESSIONS;
    FASTERQ_DUMP;
} from "./modules/SRA.nf"

include {
    MINIMAP2;
    IVAR_TRIM;
} from "./modules/preprocessing.nf"

include {
    FREYJA_VARIANTS;
    FREYJA_BOOT;
} from "./modules/freyja.nf"

workflow {
    Channel
        .fromPath(params.input)
        .set { input_ch } 

    GET_ACCESSIONS(input_ch)
        .splitCsv()
        .map { line -> line.join('') }
        .set { acc_ch }

    FASTERQ_DUMP(acc_ch)
}