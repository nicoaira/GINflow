#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { EXTRACT_META_MAP }          from '../modules/extract_meta_map/main'
include { GENERATE_NODE_EMBEDDINGS }  from '../modules/generate_node_embeddings/main'
include { GENERATE_WINDOW_VECTORS }   from '../modules/generate_window_vectors/main'
include { MERGE_EMBEDDING_CHUNKS }    from '../modules/merge_embedding_chunks/main'
include { MERGE_WINDOW_CHUNKS }       from '../modules/merge_window_chunks/main'
include { BUILD_FAISS_INDEX }         from '../modules/build_faiss_index/main'
include { QUERY_FAISS_INDEX }         from '../modules/query_faiss_index/main'
include { CLUSTER_SEEDS }             from '../modules/cluster_seeds/main'
include { ALIGN_CANDIDATES }          from '../modules/align_candidates/main'

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
    def meta_ch = meta.meta_map

    // Generate or reuse node embeddings (batching via splitCsv for scalability)
    def chunked_embeddings_ch
    def embedding_chunks_for_merge_ch
    if (params.node_embeddings_tsv) {
        def provided_embeddings = Channel.fromPath(params.node_embeddings_tsv)
            .map { path -> tuple('batch_000000', path) }

        provided_embeddings
            .multiMap { tuple ->
                windows: tuple
                merge: tuple
            }
            .set { embedding_split }

        chunked_embeddings_ch = embedding_split.windows
        embedding_chunks_for_merge_ch = embedding_split.merge
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
        embedding_chunks_for_merge_ch = embedding_split.merge
    }

    def embedding_chunk_list = embedding_chunks_for_merge_ch.collect(flat: false)
    def merged_embeddings = MERGE_EMBEDDING_CHUNKS(embedding_chunk_list)

    def embeddings_for_align = merged_embeddings.node_embeddings

    def window_inputs = chunked_embeddings_ch.map { chunk_id, embeddings_path ->
        tuple(chunk_id, embeddings_path, queries_file)
    }

    // Build sliding window vectors for database and query sets per chunk
    def window_chunks = GENERATE_WINDOW_VECTORS(window_inputs)

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

    // Build FAISS index from database windows
    def index = BUILD_FAISS_INDEX(db_vec_a.combine(db_meta_a))
    def faiss_idx = index.faiss_idx
    def faiss_map = index.faiss_map

    // Query seeds from FAISS index
    def seeds = QUERY_FAISS_INDEX(faiss_idx, db_vec_b, db_meta_b, q_vec, q_meta)

    // Group seeds into diagonal clusters
    def clusters = CLUSTER_SEEDS(seeds.seeds)

    // Align clustered candidates and produce final report TSV
    ALIGN_CANDIDATES(
        clusters.cluster_members,
        clusters.clusters,
        embeddings_for_align,
        meta_ch
    )
}
