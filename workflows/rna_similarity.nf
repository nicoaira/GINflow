#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Note: parameter declarations and validation moved to main.nf

// Include all modules (removed PREP_BATCH and MERGE_EMBEDDINGS)
include { GENERATE_WINDOWS } from '../modules/generate_windows/main'
include { EXTRACT_META_MAP } from '../modules/extract_meta_map/main'
include { GENERATE_EMBEDDINGS } from '../modules/generate_embeddings/main'
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

    // Split input TSV into batches using splitCsv operator
    def input_records = Channel.fromPath(params.input)
        .splitCsv(header: params.header, sep: '\t', strip: true, by: params.split_size)
    
    // Transform each batch of records into a temporary TSV file using collectFile
    def batch_files_ch = input_records
        .collectFile() { records ->
            // Get header from first record
            def header = records[0].keySet().join('\t')
            // Convert records to TSV lines
            def lines = records.collect { record -> record.values().join('\t') }.join('\n')
            ["batch_${records.hashCode().abs()}.tsv", header + '\n' + lines + '\n']
        }

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
    
    // Use collectFile to merge all embedding files into a single TSV
    def merged_embeddings = gen.batch_embeddings
        .collectFile(name: 'embeddings.tsv', keepHeader: true, skip: 1)
    
    def idx     = BUILD_FAISS_INDEX(merged_embeddings)
    def dists   = QUERY_FAISS_INDEX(merged_embeddings, idx.faiss_idx, idx.faiss_map)
    def sorted  = SORT_DISTANCES(dists.distances)
    PLOT_DISTANCES(sorted.sorted_distances)

    // 7
    def (agg_all, agg_un) = AGGREGATE_SCORE(sorted.sorted_distances, meta.meta_map)
    PLOT_SCORE(agg_all)

    // 8-9-10-11-12-13
    def (top_c, top_u)  = FILTER_TOP_CONTIGS(agg_all, agg_un)
    
    // Conditionally run aggregated report processes
    if (params.run_aggregated_report) {
        def contigs_draw = DRAW_CONTIG_SVGS(top_c)
        GENERATE_AGGREGATED_REPORT(top_c, contigs_draw.contig_individual)
    }
    
    // Conditionally run unaggregated report processes
    if (params.run_unaggregated_report) {
        def windows_draw = DRAW_UNAGG_SVGS(top_u)
        GENERATE_UNAGGREGATED_REPORT(top_u, windows_draw.window_individual)
    }
}
