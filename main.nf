#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ───────────────────────────────────────────────────────────
// Basic checks
// ───────────────────────────────────────────────────────────
if ( !file(params.input).exists() )
    error "Cannot find input TSV: ${params.input}"

// ───────────────────────────────────────────────────────────
// Include and invoke the main workflow
// ───────────────────────────────────────────────────────────
include { rna_similarity } from './workflows/rna_similarity.nf'

workflow {
    rna_similarity()
}
