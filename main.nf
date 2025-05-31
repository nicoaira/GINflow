#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ───────────────────────────────────────────────────────────
// User-tunable parameters (override on CLI)
// ───────────────────────────────────────────────────────────
params.input                 = "$baseDir/input.tsv"           // TSV with sequences
params.outdir                = "results"
params.model_path            = "$baseDir/model.pth"
params.structure_column_name = 'secondary_structure'
params.structure_column_num  = null
params.id_column             = 'transcript_id'

params.device                = 'cpu'
params.num_workers           = 4
params.subgraphs             = false
params.L                     = null
params.keep_paired_neighbors = false
params.retries               = 0
params.batch_size_embed      = 100       // sequences per embedding task

params.faiss_k               = 1000      // neighbors to retrieve in FAISS

// aggregation parameters
params.alpha1    = 0.25
params.alpha2    = 0.24
params.beta1     = 0.0057
params.beta2     = 1.15
params.gamma     = 0.41
params.percentile= 1
params.top_n     = 10

// plotting flags
params.plot_distances_distribution = true
params.hist_seed                   = 42
params.hist_frac                   = 0.001
params.hist_bins                   = 200

params.plot_score_distribution     = true
params.score_bins                  = 30

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
