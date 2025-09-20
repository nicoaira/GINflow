#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ───────────────────────────────────────────────────────────
// Basic checks
// ───────────────────────────────────────────────────────────
if ( !file(params.input).exists() )
    error "Cannot find input TSV: ${params.input}"

if ( params.node_embeddings_tsv && !file(params.node_embeddings_tsv).exists() )
    error "Cannot find node embeddings TSV: ${params.node_embeddings_tsv}"

if ( !file(params.queries).exists() )
    error "Cannot find queries CSV: ${params.queries}"

// ───────────────────────────────────────────────────────────
// Include and invoke the main workflow
// ───────────────────────────────────────────────────────────
include { rna_similarity } from './workflows/rna_similarity.nf'

workflow {
    rna_similarity()
}
