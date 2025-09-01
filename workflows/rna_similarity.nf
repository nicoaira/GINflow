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
include { AGGREGATE_SCORE_RAW } from '../modules/aggregate_score_raw/main'
include { PLOT_SCORE } from '../modules/plot_score/main'
include { PLOT_SCORE_WITH_NULL } from '../modules/plot_score_with_null/main'
include { FILTER_TOP_CONTIGS } from '../modules/filter_top_contigs/main'
include { DRAW_CONTIG_SVGS } from '../modules/draw_contig_svgs/main'
include { DRAW_UNAGG_SVGS } from '../modules/draw_unagg_svgs/main'
include { GENERATE_AGGREGATED_REPORT } from '../modules/generate_aggregated_report/main'
include { GENERATE_UNAGGREGATED_REPORT } from '../modules/generate_unaggregated_report/main'
include { MERGE_QUERY_RESULTS } from '../modules/merge_query_results/main'
include { NORMALIZE_SCORES } from '../modules/normalize_scores/main'
include { SHUFFLE_AND_FOLD } from '../modules/shuffle_and_fold/main'
include { MERGE_NULL_META } from '../modules/merge_null_meta/main'
include { GENERATE_NULL_WINDOWS } from '../modules/generate_null_windows/main'
include { GENERATE_EMBEDDINGS_FROM_WINDOWS_NULL } from '../modules/generate_embeddings_from_windows_null/main'
include { QUERY_FAISS_INDEX_BY_TUPLE } from '../modules/query_faiss_index_by_tuple/main'
include { SORT_DISTANCES_NULL } from '../modules/sort_distances_null/main'
include { MERGE_NULL_SCORES } from '../modules/merge_null_scores/main'
include { AGGREGATE_SCORE_WITH_NULL } from '../modules/aggregate_score_with_null/main'

workflow PER_QUERY {
    take:
        query_id
        embeddings
        faiss_idx
        faiss_map
        meta_map
    main:
        def dists  = QUERY_FAISS_INDEX(embeddings, faiss_idx, faiss_map, query_id)
        def sorted = SORT_DISTANCES(dists.distances)
        PLOT_DISTANCES(sorted.sorted_distances)

        // Aggregation is handled in the top-level workflow
    emit:
        sorted_distances = sorted.sorted_distances
}

workflow NULL_PER_QUERY {
    take:
        pairs
        faiss_idx
        faiss_map
        meta_map
    main:
        def shuf = SHUFFLE_AND_FOLD(pairs, meta_map)
        // Batch null meta TSVs using the same split_size as real data.
        // Use toList + collate to ensure groups are emitted even with small totals.
        // Collect all shuffled batches, then sort deterministically by new_id
        // to keep chunk composition and order stable across runs for caching.
        def grouped = shuf.shuffled
            .toList()
            .map { all_items ->
                all_items.sort { it[0].toString() }
                all_items
            }
            .flatMap { all_items ->
                def chunks = all_items.collate(params.split_size as int)
                chunks.collect { chunk ->
                    tuple( chunk.collect { it[0] }, chunk.collect { it[1] } )
                }
            }
        // Merge per-batch TSVs (keeps header once)
        def merged = MERGE_NULL_META(grouped)
        def id_to_emb
        if (params.subgraphs) {
            // 1) Generate windows per null batch
            def wins   = GENERATE_NULL_WINDOWS(merged.batch)
            // 2) Embed windows per null batch
            def emb    = GENERATE_EMBEDDINGS_FROM_WINDOWS_NULL(wins.window_files_null)
            // 3) Expand (ids, embeddings.tsv) into (id, embeddings.tsv) for querying
            id_to_emb = emb.embeddings_for_query.flatMap { ids, embfile -> ids.collect { q -> [ q, embfile ] } }
        } else {
            // Direct embeddings from batched TSVs (no window generation)
            def emb_direct = GENERATE_EMBEDDINGS(merged.batch.map{ it[1] })
            // Pair back each batch's IDs with the produced embeddings file
            def paired = merged.batch.map{ it[0] }.zip(emb_direct.batch_embeddings)
            id_to_emb  = paired.flatMap { ids, embfile -> ids.collect { q -> [ q, embfile ] } }
        }
        // 4) Query FAISS with correct (query_id, embeddings.tsv)
        def d    = QUERY_FAISS_INDEX_BY_TUPLE(id_to_emb, faiss_idx, faiss_map)
        def srt  = SORT_DISTANCES_NULL(d.distances)
        def scr  = AGGREGATE_SCORE_RAW(srt.sorted_distances)
    emit:
        scores = scr.scores
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

    // Build two independent query streams so both branches receive all IDs
    def queries_real = Channel.fromPath(params.queries)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row -> row['id'] }
    def queries_null = Channel.fromPath(params.queries)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row -> row['id'] }

    // Convert single-value channels to value channels so they can be reused for each query
    def embeddings_val = merged_embeddings.first()
    def faiss_idx_val  = faiss_idx_ch.first()
    def faiss_map_val  = faiss_map_ch.first()
    def meta_map_val   = meta.meta_map.first()

    def per_query = PER_QUERY(queries_real, embeddings_val, faiss_idx_val, faiss_map_val, meta_map_val)

    if (params.null_shuffles as int > 0) {
        // Build null distributions per query
        def null_ids = queries_null.flatMap { q -> (1..(params.null_shuffles as int)).collect { r -> [q, "${q}_null${r}"] } }
        def null_per = NULL_PER_QUERY(null_ids, faiss_idx_val, faiss_map_val, meta_map_val)

        // Map each null ID back to its base query ID (strip _nullN suffix)
        def null_map = null_per.scores.map { nid, f ->
            def m = (nid =~ /(.*)_null\d+$/)
            def base = m ? m[0][1] : nid
            [ base, f ]
        }
        // Group all null replicate scores per base query and merge into a single file
        // Group all null replicate scores per base query and sort file lists
        // to ensure stable ordering for cache-resume.
        def null_grouped = null_map.groupTuple().map { base_id, files ->
            def sorted_files = files.sort { it.toString() }
            tuple(base_id, sorted_files)
        }
        def null_merged  = MERGE_NULL_SCORES(null_grouped)

        // Join real sorted distances with their corresponding null distributions
        def joined = per_query.sorted_distances.join(null_merged.null_scores)

        // Aggregate and normalize per query, then continue plotting and filtering
        def agg_n = AGGREGATE_SCORE_WITH_NULL(joined, meta_map_val)
        // Plot real vs. null distribution per query
        PLOT_SCORE_WITH_NULL(agg_n.enriched_agg.join(null_merged.null_scores))

        def joined_scores = agg_n.enriched_agg.join(agg_n.windows)
        def filtered = FILTER_TOP_CONTIGS(joined_scores)

        if (params.run_aggregated_report && params.draw_contig_svgs) {
            def contigs_draw = DRAW_CONTIG_SVGS(filtered.top_contigs)
            def agg_report_in = filtered.top_contigs.join(contigs_draw.contig_individual)
            GENERATE_AGGREGATED_REPORT(agg_report_in)
        }

        if (params.run_unaggregated_report) {
            def windows_draw = DRAW_UNAGG_SVGS(filtered.top_contigs_unagg)
            def unagg_report_in = filtered.top_contigs_unagg.join(windows_draw.window_individual)
            GENERATE_UNAGGREGATED_REPORT(unagg_report_in)
        }

        // Merge all final per-query results
        MERGE_QUERY_RESULTS(
            per_query.sorted_distances.map{ it[1] }.collect(),
            agg_n.contigs.map{ it[1] }.collect(),
            agg_n.windows.map{ it[1] }.collect(),
            agg_n.enriched_agg.map{ it[1] }.collect()
        )
    } else {
        // No nulls: aggregate per query directly, then plot/filter/report and merge
        def agg = AGGREGATE_SCORE(per_query.sorted_distances, meta_map_val)
        PLOT_SCORE(agg.enriched_agg)

        def joined_scores = agg.enriched_agg.join(agg.windows)
        def filtered = FILTER_TOP_CONTIGS(joined_scores)

        if (params.run_aggregated_report && params.draw_contig_svgs) {
            def contigs_draw = DRAW_CONTIG_SVGS(filtered.top_contigs)
            def agg_report_in = filtered.top_contigs.join(contigs_draw.contig_individual)
            GENERATE_AGGREGATED_REPORT(agg_report_in)
        }

        if (params.run_unaggregated_report) {
            def windows_draw = DRAW_UNAGG_SVGS(filtered.top_contigs_unagg)
            def unagg_report_in = filtered.top_contigs_unagg.join(windows_draw.window_individual)
            GENERATE_UNAGGREGATED_REPORT(unagg_report_in)
        }

        MERGE_QUERY_RESULTS(
            per_query.sorted_distances.map{ it[1] }.collect(),
            agg.contigs.map{ it[1] }.collect(),
            agg.windows.map{ it[1] }.collect(),
            agg.enriched_agg.map{ it[1] }.collect()
        )
    }
}
