#!/usr/bin/env nextflow

/*
 * Automated SRA download and processing pipeline for Freyja analysis.
 */

// Enable DSL 2 syntax
nextflow.enable.dsl=2

// Define the input parameters
params.sra_data = "$baseDir/data/input/*.csv"

// SARS-CoV-2 by default, but can be changed to any other pathogen.
params.refseq = "$baseDir/data/NC_045512_Hu-1.fasta"
params.primer_bed = "$baseDir/data/nCov-2019_v3.primer.bed"

Channel 
    .fromPath(params.sra_data)
    .splitCsv(header: true, sep: ',')
    .map { row -> row.acc}
    .set { sra_accessions_ch }
    // .fromSRA(apiKey:$ncbi_api_key)
sra_accessions_ch.view()