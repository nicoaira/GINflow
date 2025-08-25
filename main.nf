#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ───────────────────────────────────────────────────────────
// Basic checks
// ───────────────────────────────────────────────────────────
if ( !file(params.input).exists() )
    error "Cannot find input TSV: ${params.input}"

if ( params.embeddings_tsv && !file(params.embeddings_tsv).exists() )
    error "Cannot find embeddings TSV: ${params.embeddings_tsv}"

if ( (params.faiss_index && !params.faiss_mapping) || (!params.faiss_index && params.faiss_mapping) )
    error "Both --faiss_index and --faiss_mapping must be provided together"

if ( params.faiss_index && !file(params.faiss_index).exists() )
    error "Cannot find FAISS index: ${params.faiss_index}"

if ( params.faiss_mapping && !file(params.faiss_mapping).exists() )
    error "Cannot find FAISS mapping: ${params.faiss_mapping}"

// ───────────────────────────────────────────────────────────
// Include and invoke the main workflow
// ───────────────────────────────────────────────────────────
include { rna_similarity } from './workflows/rna_similarity.nf'

workflow {
    rna_similarity()
}
