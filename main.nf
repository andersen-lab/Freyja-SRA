#!/usr/bin/env nextflow

/*
 * Automated pipeline for Freyja analysis on SRA data
 */

// Enable DSL 2 syntax
nextflow.enable.dsl=2

// Define the input parameters

// SARS-CoV-2 by default, but can be changed to any other pathogen.
params.ref = "$baseDir/data/preprocessing/NC_045512_Hu-1.fasta"
ref = file(params.ref)

params.primer_bed = "$baseDir/data/preprocessing/nCov-2019_v3.primer.bed"
primer_bed = file(params.primer_bed)

// Import modules
include { 
    GET_ACCESSIONS;
    FETCH_NGS;
} from "./modules/SRA.nf"

Channel
    .fromPath(params.input)
    .set { sra_data_ch }

workflow {
    GET_ACCESSIONS(sra_data_ch)
    FETCH_NGS(GET_ACCESSIONS.out)
    
    Channel.fromFilePairs("${FETCH_NGS.out}/*_{1,2}.fastq.gz")
        .set { fastq_ch }
}