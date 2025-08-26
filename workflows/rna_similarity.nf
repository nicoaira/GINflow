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
include { MERGE_QUERY_RESULTS } from '../modules/merge_query_results/main'

workflow PER_QUERY {
    take:
        query_id
        embeddings
        faiss_idx
        faiss_map
        meta_map
    main:
        def dists  = QUERY_FAISS_INDEX(embeddings, faiss_idx, faiss_map, query_id)
        def sorted = SORT_DISTANCES(dists.distances, query_id)
        PLOT_DISTANCES(sorted.sorted_distances, query_id)

        def agg    = AGGREGATE_SCORE(sorted.sorted_distances, meta_map, query_id)
        PLOT_SCORE(agg.enriched_all, query_id)

        def filtered = FILTER_TOP_CONTIGS(agg.enriched_all, agg.enriched_unagg, query_id)

        if (params.run_aggregated_report) {
            def contigs_draw = DRAW_CONTIG_SVGS(filtered.top_contigs, query_id)
            GENERATE_AGGREGATED_REPORT(filtered.top_contigs, contigs_draw.contig_individual, query_id)
        }

        if (params.run_unaggregated_report) {
            def windows_draw = DRAW_UNAGG_SVGS(filtered.top_contigs_unagg, query_id)
            GENERATE_UNAGGREGATED_REPORT(filtered.top_contigs_unagg, windows_draw.window_individual, query_id)
        }
    emit:
        sorted_distances = sorted.sorted_distances
        enriched_all     = agg.enriched_all
        enriched_unagg   = agg.enriched_unagg
}

// ───────────────────────────────────────────────────────────
// Workflow wiring
// ───────────────────────────────────────────────────────────
workflow rna_similarity {
    // 0) Build the metadata map once
    def meta = EXTRACT_META_MAP(file(params.input))

    def merged_embeddings
    if (params.embeddings_tsv) {
        // Use provided embeddings file and ensure it's available in the output directory
        merged_embeddings = Channel.fromPath(params.embeddings_tsv)
            .map { path ->
                def out = file(params.outdir).resolve('embeddings.tsv')
                out.parent.mkdirs()
                path.copyTo(out)
                out
            }
    } else {
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
        merged_embeddings = gen.batch_embeddings
            .collectFile(name: 'embeddings.tsv', keepHeader: true, skip: 1, storeDir: params.outdir)
    }

    def faiss_idx_ch
    def faiss_map_ch
    if (params.faiss_index && params.faiss_mapping) {
        faiss_idx_ch = Channel.fromPath(params.faiss_index)
        faiss_map_ch = Channel.fromPath(params.faiss_mapping)
    } else {
        def idx = BUILD_FAISS_INDEX(merged_embeddings)
        faiss_idx_ch = idx.faiss_idx
        faiss_map_ch = idx.faiss_map
    }

    def queries = Channel.fromPath(params.queries)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row -> row['id'] }

    // Convert single-value channels to value channels so they can be reused for each query
    def embeddings_val = merged_embeddings.first()
    def faiss_idx_val  = faiss_idx_ch.first()
    def faiss_map_val  = faiss_map_ch.first()
    def meta_map_val   = meta.meta_map.first()

    def per_query = PER_QUERY(queries, embeddings_val, faiss_idx_val, faiss_map_val, meta_map_val)

    MERGE_QUERY_RESULTS(
        per_query.sorted_distances.collect(),
        per_query.enriched_all.collect(),
        per_query.enriched_unagg.collect()
    )
}
