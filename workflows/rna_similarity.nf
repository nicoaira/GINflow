#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Note: parameter declarations and validation moved to main.nf

// Include all modules
include { GENERATE_WINDOWS } from '../modules/generate_windows/main'
include { EXTRACT_META_MAP } from '../modules/extract_meta_map/main'
include { PREP_BATCH } from '../modules/prep_batch/main'
include { GENERATE_EMBEDDINGS } from '../modules/generate_embeddings/main'
include { MERGE_EMBEDDINGS } from '../modules/merge_embeddings/main'
include { BUILD_FAISS_INDEX } from '../modules/build_faiss_index/main'
include { QUERY_FAISS_INDEX } from '../modules/query_faiss_index/main'
include { SORT_DISTANCES } from '../modules/sort_distances/main'
include { PLOT_DISTANCES } from '../modules/plot_distances/main'
include { AGGREGATE_SCORE } from '../modules/aggregate_score/main'
include { PLOT_SCORE } from '../modules/plot_score/main'
include { FILTER_TOP_CONTIGS } from '../modules/filter_top_contigs/main'
include { DRAW_CONTIG_SVGS } from '../modules/draw_contig_svgs/main'
include { DRAW_UNAGG_SVGS } from '../modules/draw_unagg_svgs/main'
include { GENERATE_AGGREGATED_REPORT } from '../modules/generate_aggregated_report/main'
include { GENERATE_UNAGGREGATED_REPORT } from '../modules/generate_unaggregated_report/main'

// ───────────────────────────────────────────────────────────
// Workflow wiring
// ───────────────────────────────────────────────────────────
workflow rna_similarity {
    // 0) Build the metadata map once
    def meta = EXTRACT_META_MAP(file(params.input))

    // Split transcripts into batches
    // Ensuring header and separator are correctly used for TSV
    def transcript_batches = Channel.fromPath(params.input)
        .splitCsv(header: params.header, sep:'\t', strip:true, by: params.batch_size_embed)
    def batch_files_ch = PREP_BATCH(transcript_batches).batch_file

    def gen
    if(params.subgraphs) {
        // Each GENERATE_WINDOWS task processes one batch and outputs window files
        // GENERATE_EMBEDDINGS should start as soon as each window file is ready
        def window_files = GENERATE_WINDOWS(batch_files_ch).window_files
        gen = GENERATE_EMBEDDINGS(window_files)
    } else {
        // Direct processing of batch files without windowing
        gen = GENERATE_EMBEDDINGS(batch_files_ch)
    }
    
    // Explicitly collect all batch embeddings before merging
    // This ensures MERGE_EMBEDDINGS only starts when all GENERATE_EMBEDDINGS tasks have completed
    def all_embeddings = gen.batch_embeddings.collect()
    def merged = MERGE_EMBEDDINGS(all_embeddings)
    
    def idx     = BUILD_FAISS_INDEX(merged.embeddings)
    def dists   = QUERY_FAISS_INDEX(merged.embeddings, idx.faiss_idx, idx.faiss_map)
    def sorted  = SORT_DISTANCES(dists.distances)
    PLOT_DISTANCES(sorted.sorted_distances)

    // 7
    def (agg_all, agg_un) = AGGREGATE_SCORE(sorted.sorted_distances, meta.meta_map)
    PLOT_SCORE(agg_all)

    // 8-9-10-11-12-13
    def (top_c, top_u)  = FILTER_TOP_CONTIGS(agg_all, agg_un)
    def contigs_draw    = DRAW_CONTIG_SVGS(top_c)
    def windows_draw    = DRAW_UNAGG_SVGS(top_u)
    GENERATE_AGGREGATED_REPORT(top_c, contigs_draw.contig_individual)
    GENERATE_UNAGGREGATED_REPORT(top_u, windows_draw.window_individual)
}
