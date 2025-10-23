#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { EXTRACT_META_MAP }                  from '../modules/extract_meta_map/main'
include { GENERATE_NODE_EMBEDDINGS }          from '../modules/generate_node_embeddings/main'
include { GENERATE_WINDOW_VECTORS }           from '../modules/generate_window_vectors/main'
include { MERGE_EMBEDDING_CHUNKS }            from '../modules/merge_embedding_chunks/main'
include { MERGE_WINDOW_CHUNKS }               from '../modules/merge_window_chunks/main'
include { CALCULATE_EVD_PARAMS }              from '../modules/calculate_evd_params/main'
include { BUILD_FAISS_INDEX }                 from '../modules/build_faiss_index/main'
include { SPLIT_QUERIES }                     from '../modules/split_queries/main'
include { QUERY_FAISS_INDEX }                 from '../modules/query_faiss_index/main'
include { CLUSTER_SEEDS }                     from '../modules/cluster_seeds/main'
include { ALIGN_CANDIDATES }                  from '../modules/align_candidates/main'
include { DRAW_ALIGNMENT_SVGS_RNARTISTCORE }  from '../modules/draw_alignment_svgs_rnartistcore/main'
include { DRAW_ALIGNMENT_SVGS_R4RNA }         from '../modules/draw_alignment_svgs_r4rna/main'
include { GENERATE_REPORT }                   from '../modules/generate_report/main'

workflow rna_similarity {
    if (!params.input) {
        error "Parameter --input is required"
    }
    if (!params.queries) {
        error "Parameter --queries is required"
    }

    def input_file = file(params.input)
    if (!input_file.exists()) {
        error "Input file not found: ${input_file}"
    }
    def queries_file = file(params.queries)
    if (!queries_file.exists()) {
        error "Queries file not found: ${queries_file}"
    }

    // Extract metadata for later reporting
    def meta = EXTRACT_META_MAP(input_file)
    meta.meta_map
        .multiMap { path ->
            windows: path
            align: path
        }
        .set { meta_split }

    def meta_for_windows = meta_split.windows
    def meta_ch = meta_split.align

    // Generate or reuse node embeddings (batching via splitCsv for scalability)
    def chunked_embeddings_ch
    def embeddings_for_align
    if (params.node_embeddings_tsv) {
        def provided_embeddings = file(params.node_embeddings_tsv)
        if (!provided_embeddings.exists()) {
            error "Provided --node_embeddings_tsv not found: ${provided_embeddings}"
        }

        chunked_embeddings_ch = Channel.of(tuple('batch_000000', provided_embeddings))
        embeddings_for_align = Channel.of(provided_embeddings)
    } else {
        if (!params.header) {
            error "Batched node embedding generation requires --header true to reconstruct TSV chunks"
        }
        def batches = Channel.fromPath(input_file)
            .splitCsv(header: params.header, sep: '\t', strip: true, by: params.split_size)

        def batch_index = new java.util.concurrent.atomic.AtomicInteger(0)
        def batch_files = batches.collectFile { records ->
            def columns = records[0].keySet().toList()
            def headerLine = columns.join('\t')
            def lines = records.collect { row -> columns.collect { col -> row[col] ?: '' }.join('\t') }
            def content = ([headerLine] + lines).join('\n') + '\n'
            def idx = batch_index.getAndIncrement()
            [ String.format('batch_%06d.tsv', idx), content ]
        }

        def chunk_inputs = batch_files.map { path ->
            def chunk_id = path.simpleName
            tuple(chunk_id, path)
        }

        def node_batches = GENERATE_NODE_EMBEDDINGS(chunk_inputs)
        node_batches.node_embeddings
            .multiMap { tuple ->
                windows: tuple
                merge: tuple
            }
            .set { embedding_split }

        chunked_embeddings_ch = embedding_split.windows
        def embedding_chunks_for_merge_ch = embedding_split.merge
        def embedding_chunk_list = embedding_chunks_for_merge_ch.collect(flat: false)
        def merged_embeddings = MERGE_EMBEDDING_CHUNKS(embedding_chunk_list)

        embeddings_for_align = merged_embeddings.node_embeddings
    }

    def window_inputs = chunked_embeddings_ch.map { chunk_id, embeddings_path ->
        tuple(chunk_id, embeddings_path, queries_file)
    }

    // Build sliding window vectors for database and query sets per chunk
    def window_chunks = GENERATE_WINDOW_VECTORS(window_inputs, meta_for_windows)

    def window_chunk_list = window_chunks.window_chunk.collect(flat: false)

    def merged_windows = MERGE_WINDOW_CHUNKS(window_chunk_list)

    merged_windows.database_vectors
        .multiMap { path -> for_index: path; for_query: path }
        .set { db_vec_split }

    def db_vec_a = db_vec_split.for_index
    def db_vec_b = db_vec_split.for_query

    merged_windows.database_metadata
        .multiMap { path -> for_index: path; for_query: path }
        .set { db_meta_split }

    def db_meta_a = db_meta_split.for_index
    def db_meta_b = db_meta_split.for_query

    def q_vec = merged_windows.query_vectors
    def q_meta = merged_windows.query_metadata

    // Calculate EVD parameters once for the database (if E-value calculation is enabled)
    def evd_params_ch = CALCULATE_EVD_PARAMS(embeddings_for_align).evd_params

    // Build FAISS index from database windows
    def index = BUILD_FAISS_INDEX(db_vec_a.combine(db_meta_a))
    def faiss_idx = index.faiss_idx
    def faiss_map = index.faiss_map

    // Split queries into separate files for parallel processing
    def split_queries_result = SPLIT_QUERIES(q_vec, q_meta)

    // Create tuples of (query_id, query_vectors, query_metadata) for each query
    def query_tuples = split_queries_result.query_vectors
        .flatten()
        .map { vec_file ->
            def safe_id = vec_file.name.replaceAll(/^query_(.+)_vectors\.npy$/, '$1')
            tuple(safe_id, vec_file)
        }
        .join(
            split_queries_result.query_metadata
                .flatten()
                .map { meta_file ->
                    def safe_id = meta_file.name.replaceAll(/^query_(.+)_metadata\.tsv$/, '$1')
                    tuple(safe_id, meta_file)
                }
        )
        .map { safe_id, vec_file, meta_file ->
            tuple(safe_id, vec_file, meta_file)
        }

    // For each query, combine with shared database resources and run QUERY_FAISS_INDEX
    def query_inputs = query_tuples
        .combine(faiss_idx)
        .combine(db_vec_b)
        .combine(db_meta_b)
        .map { query_id, q_vec_file, q_meta_file, faiss, db_vec, db_meta ->
            tuple(query_id, faiss, db_vec, db_meta, q_vec_file, q_meta_file)
        }

    // Query FAISS index for each query independently (parallelized)
    def seeds_per_query = QUERY_FAISS_INDEX(
        query_inputs.map { it[1] },  // faiss_idx
        query_inputs.map { it[2] },  // database_vectors
        query_inputs.map { it[3] },  // database_metadata
        query_inputs.map { it[4] },  // query_vectors
        query_inputs.map { it[5] }   // query_metadata
    )

    // Cluster seeds for each query independently (parallelized)
    def clusters_per_query = CLUSTER_SEEDS(seeds_per_query.seeds)

    // Align candidates for each query independently (parallelized)
    // Extract query_id from filenames and join cluster outputs
    def members_with_id = clusters_per_query.cluster_members
        .map { file ->
            def query_id = file.baseName.replaceAll(/^cluster_members_/, '')
            tuple(query_id, file)
        }

    def clusters_with_id = clusters_per_query.clusters
        .map { file ->
            def query_id = file.baseName.replaceAll(/^clusters_/, '')
            tuple(query_id, file)
        }

    // Join by query_id to pair members with clusters, then split for ALIGN_CANDIDATES
    def paired_clusters = members_with_id
        .join(clusters_with_id)
        .multiMap { query_id, members, clusters ->
            cluster_members: members
            cluster_summaries: clusters
        }

    def alignments = ALIGN_CANDIDATES(
        paired_clusters.cluster_members,
        paired_clusters.cluster_summaries,
        embeddings_for_align,
        meta_ch,
        evd_params_ch
    )

    // Conditionally draw SVG visualizations and generate HTML report
    if (params.enable_drawings || params.enable_report) {
        // Choose drawing backend based on parameter
        def alignment_svgs
        if (params.drawing_backend == 'r4rna') {
            alignment_svgs = DRAW_ALIGNMENT_SVGS_R4RNA(alignments.alignments)
        } else {
            // Default to rnartistcore
            alignment_svgs = DRAW_ALIGNMENT_SVGS_RNARTISTCORE(alignments.alignments)
        }

        // Generate HTML report with embedded SVG visualizations
        if (params.enable_report) {
            // Both channels originate from the same source and are in the same order
            // Use merge to pair them by index
            def paired_for_report = alignments.alignments
                .merge(alignment_svgs.alignment_individual)

            GENERATE_REPORT(paired_for_report)
        }
    }
}
