#!/usr/bin/env nextflow

/*
 * Automated pipeline for Freyja analysis of SRA data
 */


// SARS-CoV-2 by default, but can be changed to any other pathogen.
params.ref = "$baseDir/data/NC_045512_Hu-1.fasta"
params.bedfiles = "$baseDir/data/bedfiles"


ref = file(params.ref)
bedfiles = file(params.bedfiles)

// Import modules
include { 
    GET_ACCESSIONS;
    GET_AMPLICON_SCHEME;
    FASTERQ_DUMP;
} from "./modules/sra.nf"

include {
    MINIMAP2;
    IVAR_TRIM;
} from "./modules/preprocessing.nf"

include {
    FREYJA_VARIANTS;
    FREYJA_DEMIX;
} from "./modules/freyja.nf"

workflow preprocessing {
    Channel
        .fromPath(params.input)
        .set { input_ch } 

    GET_ACCESSIONS(input_ch)
        .splitCsv()
        .map { line -> line.join('') }
        .take(5)
        .set { acc_ch }


    FASTERQ_DUMP(acc_ch)
    MINIMAP2(FASTERQ_DUMP.out, ref)
    
    GET_AMPLICON_SCHEME(acc_ch, input_ch.first())
        .map { it.text }
        .set { primer_scheme_ch }

    
    IVAR_TRIM(MINIMAP2.out, primer_scheme_ch, bedfiles)
    FREYJA_VARIANTS(IVAR_TRIM.out)
    FREYJA_DEMIX(FREYJA_VARIANTS.out)
}

// workflow DEMIX {
//     FREYJA_DEMIX()
// }