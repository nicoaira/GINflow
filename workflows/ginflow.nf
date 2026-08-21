include { PREPARE_WINDOWS as PREPARE_DB }    from './prepare_windows'
include { PREPARE_WINDOWS as PREPARE_QUERY } from './prepare_windows'
include { BUILD_FAISS_INDEX }                from '../modules/build_faiss_index/main'
include { SEARCH_FAISS }                     from '../modules/search_faiss/main'
include { BUILD_CUVS_INDEX }                 from '../modules/build_cuvs_index/main'
include { SEARCH_CUVS }                      from '../modules/search_cuvs/main'
include { BUILD_HNSWLIB_INDEX }             from '../modules/build_hnswlib_index/main'
include { SEARCH_HNSWLIB }                  from '../modules/search_hnswlib/main'
include { BUILD_PQ_CAGRA_INDEX }            from '../modules/build_pq_cagra_index/main'
include { SEARCH_PQ_CAGRA }                 from '../modules/search_pq_cagra/main'
include { RERANK_CANDIDATES }               from '../modules/rerank_candidates/main'
include { CLUSTER_SEEDS }                    from '../modules/cluster_seeds/main'
include { ALIGN_CLUSTERS }                   from '../modules/align_clusters/main'
include { ESTIMATE_EVD as ESTIMATE_EVD_BUILD } from '../modules/estimate_evd/main'
include { ESTIMATE_EVD as ESTIMATE_EVD_QUERY } from '../modules/estimate_evd/main'
include { SPLIT_ALIGNMENTS as SPLIT_CLUSTERS } from '../modules/split_alignments/main'
include { MERGE_ALIGNMENTS }                 from '../modules/merge_alignments/main'
include { DRAW_RNARTISTCORE }                from '../modules/draw_rnartistcore/main'
include { DRAW_R4RNA }                       from '../modules/draw_r4rna/main'
include { DRAW_SW }                          from '../modules/draw_sw/main'
include { WRITE_REPORT }                     from '../modules/write_report/main'

workflow GINFLOW {

    take:
    structures
    queries
    database
    reuse_windows
    evd_existing
    index_library
    quantize_mode

    main:
    ch_versions        = channel.empty()
    ch_graphs          = channel.empty()
    ch_embeddings      = channel.empty()
    ch_windows         = channel.empty()
    ch_built_database  = channel.empty()
    ch_search_database = channel.empty()
    ch_seed_shards     = channel.empty()
    ch_rerank_metrics  = channel.empty()
    ch_seeds           = channel.empty()
    ch_clusters        = channel.empty()
    ch_cluster_members = channel.empty()
    ch_alignments      = channel.empty()
    ch_alignment_text  = channel.empty()
    ch_query_embeddings = channel.empty()
    ch_query_metadata   = channel.empty()
    ch_quantization     = channel.empty()
    ch_quantization_source = channel.empty()
    ch_quantized_windows = channel.empty()
    ch_evd              = channel.empty()
    ch_published_evd    = channel.empty()
    ch_plots_rnartist   = channel.empty()
    ch_plots_r4rna      = channel.empty()
    ch_plots_sw         = channel.empty()
    ch_report           = channel.empty()
    alignment_params    = file("${projectDir}/assets/alignment.json", checkIfExists: true)
    hnsw_bundle         = file("${projectDir}/vendor/hnswlib-0.8.0", checkIfExists: true)
    def pq_codes        = quantize_mode in ['pq', 'opq']
    def sq_vectors      = quantize_mode == 'sq'
    def need_quantize   = quantize_mode != 'none'
    def exact_index     = index_library == 'faiss' && params.faiss_index in ['flatip', 'flatl2']
    def run_rerank      = (params.exact_rerank || params.hnswlib_rerank) && !exact_index

    if (structures) {
        PREPARE_DB(
            structures,
            (queries && !reuse_windows) ? 'db' : '',
            params.shard_size,
            need_quantize,
            null
        )
        ch_versions   = ch_versions.mix(PREPARE_DB.out.versions)
        ch_graphs     = ch_graphs.mix(PREPARE_DB.out.graphs)
        ch_embeddings = ch_embeddings.mix(PREPARE_DB.out.embeddings)
        ch_windows    = ch_windows.mix(PREPARE_DB.out.windows)
        ch_quantization = ch_quantization.mix(PREPARE_DB.out.quantization)
        ch_quantized_windows = ch_quantized_windows.mix(PREPARE_DB.out.quantized_windows)

        PREPARE_DB.out.windows
            .multiMap { meta, npz, manifest ->
                npz: npz
                manifest: manifest
            }
            .set { ch_db_win }

        ch_db_windows   = ch_db_win.npz.collect()
        ch_db_manifests = ch_db_win.manifest.collect()
        ch_db_embeddings = PREPARE_DB.out.embeddings.map { meta, npz, manifest -> npz }.collect()
        ch_db_metadata   = PREPARE_DB.out.graphs.map { meta, tensors, sidecar -> sidecar }.collect()

        if (index_library == 'faiss') {
            def faiss_windows = sq_vectors ? PREPARE_DB.out.quantized_npz.collect() : ch_db_windows
            def faiss_manifests = sq_vectors ? PREPARE_DB.out.quantized_manifests.collect() : ch_db_manifests
            BUILD_FAISS_INDEX(faiss_windows, faiss_manifests, ch_db_embeddings, ch_db_metadata)
            ch_versions       = ch_versions.mix(BUILD_FAISS_INDEX.out.versions)
            ch_built_database = BUILD_FAISS_INDEX.out.database
        }
        else if (index_library == 'cagra' && pq_codes) {
            BUILD_PQ_CAGRA_INDEX(
                PREPARE_DB.out.quantized_windows.map { meta, directory -> directory },
                PREPARE_DB.out.quantization,
                ch_db_embeddings,
                ch_db_metadata
            )
            ch_versions       = ch_versions.mix(BUILD_PQ_CAGRA_INDEX.out.versions)
            ch_built_database = BUILD_PQ_CAGRA_INDEX.out.database
        }
        else if (index_library in ['cagra', 'ivf']) {
            def cuvs_windows = sq_vectors ? PREPARE_DB.out.quantized_npz.collect() : ch_db_windows
            def cuvs_manifests = sq_vectors ? PREPARE_DB.out.quantized_manifests.collect() : ch_db_manifests
            BUILD_CUVS_INDEX(cuvs_windows, cuvs_manifests, ch_db_embeddings, ch_db_metadata)
            ch_versions       = ch_versions.mix(BUILD_CUVS_INDEX.out.versions)
            ch_built_database = BUILD_CUVS_INDEX.out.database
        }
        else if (index_library == 'hnswlib') {
            BUILD_HNSWLIB_INDEX(
                PREPARE_DB.out.quantized_windows.map { meta, directory -> directory },
                PREPARE_DB.out.quantization,
                ch_db_embeddings,
                ch_db_metadata,
                hnsw_bundle
            )
            ch_versions       = ch_versions.mix(BUILD_HNSWLIB_INDEX.out.versions)
            ch_built_database = BUILD_HNSWLIB_INDEX.out.database
        }
        else {
            error "Unsupported index library: ${index_library}"
        }
        ch_search_database = ch_built_database

        ESTIMATE_EVD_BUILD(ch_built_database, alignment_params)
        ch_versions = ch_versions.mix(ESTIMATE_EVD_BUILD.out.versions)
        ch_evd      = ESTIMATE_EVD_BUILD.out.evd
        ch_published_evd = ESTIMATE_EVD_BUILD.out.evd
    }
    else if (database) {
        ch_search_database = channel.fromPath(database, checkIfExists: true, type: 'dir')
        if (evd_existing) {
            ch_evd = channel.fromPath(evd_existing, checkIfExists: true)
        }
        else {
            ESTIMATE_EVD_QUERY(ch_search_database, alignment_params)
            ch_versions = ch_versions.mix(ESTIMATE_EVD_QUERY.out.versions)
            ch_evd      = ESTIMATE_EVD_QUERY.out.evd
            ch_published_evd = ESTIMATE_EVD_QUERY.out.evd
        }
        if (need_quantize && sq_vectors) {
            ch_quantization_source = channel
                .fromPath("${database}/quantization", checkIfExists: true, type: 'dir')
        }
    }

    if (queries) {
        if (reuse_windows) {
            ch_query_windows    = PREPARE_DB.out.windows
            ch_query_embeddings = PREPARE_DB.out.embeddings
            ch_query_metadata   = PREPARE_DB.out.graphs
            ch_quantized_windows = PREPARE_DB.out.quantized_windows
        }
        else {
            def query_quantization = sq_vectors ? (structures ? ch_quantization : ch_quantization_source) : null
            PREPARE_QUERY(
                queries,
                structures ? 'query' : '',
                params.search_shard_size ?: params.shard_size,
                sq_vectors,
                query_quantization
            )
            ch_versions         = ch_versions.mix(PREPARE_QUERY.out.versions)
            ch_graphs           = ch_graphs.mix(PREPARE_QUERY.out.graphs)
            ch_embeddings       = ch_embeddings.mix(PREPARE_QUERY.out.embeddings)
            ch_windows          = ch_windows.mix(PREPARE_QUERY.out.windows)
            ch_query_windows    = PREPARE_QUERY.out.windows
            ch_query_embeddings = PREPARE_QUERY.out.embeddings
            ch_query_metadata   = PREPARE_QUERY.out.graphs
            ch_quantized_windows = PREPARE_QUERY.out.quantized_windows
            if (!structures && sq_vectors) {
                ch_quantization = PREPARE_QUERY.out.quantization
            }
        }

        def search_windows = ch_query_windows
        if (sq_vectors && !reuse_windows) {
            search_windows = PREPARE_QUERY.out.quantized_npz
                .merge(PREPARE_QUERY.out.quantized_manifests)
                .map { npz, manifest -> tuple([id: npz.baseName], npz, manifest) }
        }
        else if (sq_vectors && reuse_windows) {
            search_windows = PREPARE_DB.out.quantized_npz
                .merge(PREPARE_DB.out.quantized_manifests)
                .map { npz, manifest -> tuple([id: npz.baseName], npz, manifest) }
        }

        if (index_library == 'faiss') {
            SEARCH_FAISS(search_windows, ch_search_database.collect())
            ch_versions    = ch_versions.mix(SEARCH_FAISS.out.versions)
            ch_seed_shards = SEARCH_FAISS.out.seeds
        }
        else if (index_library == 'cagra' && pq_codes) {
            SEARCH_PQ_CAGRA(ch_query_windows, ch_search_database.collect())
            ch_versions    = ch_versions.mix(SEARCH_PQ_CAGRA.out.versions)
            ch_seed_shards = SEARCH_PQ_CAGRA.out.seeds
        }
        else if (index_library in ['cagra', 'ivf']) {
            SEARCH_CUVS(search_windows, ch_search_database.collect())
            ch_versions    = ch_versions.mix(SEARCH_CUVS.out.versions)
            ch_seed_shards = SEARCH_CUVS.out.seeds
        }
        else if (index_library == 'hnswlib') {
            SEARCH_HNSWLIB(ch_query_windows, ch_search_database.collect(), hnsw_bundle)
            ch_versions    = ch_versions.mix(SEARCH_HNSWLIB.out.versions)
            ch_seed_shards = SEARCH_HNSWLIB.out.seeds
        }
        else {
            error "Unsupported index library: ${index_library}"
        }
        if (run_rerank) {
            RERANK_CANDIDATES(
                ch_seed_shards,
                ch_search_database.collect(),
                ch_query_windows.map { meta, npz, manifest -> npz }.collect(),
                ch_query_windows.map { meta, npz, manifest -> manifest }.collect()
            )
            ch_versions       = ch_versions.mix(RERANK_CANDIDATES.out.versions)
            ch_rerank_metrics = RERANK_CANDIDATES.out.metrics
            ch_seed_shards    = RERANK_CANDIDATES.out.seeds
        }
        ch_seeds = ch_seed_shards
            .map { meta, tsv -> tsv }
            .collectFile(name: 'seeds.tsv', keepHeader: true, skip: 1, sort: true)

        CLUSTER_SEEDS(ch_seeds)
        ch_versions        = ch_versions.mix(CLUSTER_SEEDS.out.versions)
        ch_clusters        = CLUSTER_SEEDS.out.clusters
        ch_cluster_members = CLUSTER_SEEDS.out.members

        SPLIT_CLUSTERS(CLUSTER_SEEDS.out.clusters)
        ch_versions = ch_versions.mix(SPLIT_CLUSTERS.out.versions)

        ALIGN_CLUSTERS(
            SPLIT_CLUSTERS.out.alignments.flatten(),
            ch_query_embeddings.map { meta, npz, manifest -> npz }.collect(),
            ch_query_metadata.map { meta, tensors, sidecar -> sidecar }.collect(),
            ch_search_database.collect(),
            alignment_params,
            ch_evd.collect()
        )
        ch_versions = ch_versions.mix(ALIGN_CLUSTERS.out.versions)

        MERGE_ALIGNMENTS(
            ALIGN_CLUSTERS.out.alignments.collect().ifEmpty([]),
            ALIGN_CLUSTERS.out.text.collect().ifEmpty([])
        )
        ch_versions       = ch_versions.mix(MERGE_ALIGNMENTS.out.versions)
        ch_alignments     = MERGE_ALIGNMENTS.out.alignments
        ch_alignment_text = MERGE_ALIGNMENTS.out.text

        def plot_rn = params.plot_backend in ['rnartistcore', 'both']
        def plot_r4 = params.plot_backend in ['r4rna', 'both']
        if (plot_rn) {
            DRAW_RNARTISTCORE(ALIGN_CLUSTERS.out.alignments)
            ch_versions       = ch_versions.mix(DRAW_RNARTISTCORE.out.versions)
            ch_plots_rnartist = DRAW_RNARTISTCORE.out.plots.collect()
        }
        if (plot_r4) {
            DRAW_R4RNA(ALIGN_CLUSTERS.out.alignments)
            ch_versions    = ch_versions.mix(DRAW_R4RNA.out.versions)
            ch_plots_r4rna = DRAW_R4RNA.out.plots.collect()
        }
        if (params.plot_sw) {
            DRAW_SW(
                ALIGN_CLUSTERS.out.alignments,
                CLUSTER_SEEDS.out.clusters,
                ch_query_embeddings.map { meta, npz, manifest -> npz }.collect(),
                ch_search_database.collect(),
                alignment_params
            )
            ch_versions = ch_versions.mix(DRAW_SW.out.versions)
            ch_plots_sw = DRAW_SW.out.plots.collect()
        }

        WRITE_REPORT(
            MERGE_ALIGNMENTS.out.alignments,
            MERGE_ALIGNMENTS.out.text,
            ch_evd.collect(),
            CLUSTER_SEEDS.out.clusters,
            ch_seeds,
            ch_plots_rnartist.ifEmpty([]),
            ch_plots_r4rna.ifEmpty([]),
            ch_plots_sw.ifEmpty([])
        )
        ch_versions = ch_versions.mix(WRITE_REPORT.out.versions)
        ch_report   = WRITE_REPORT.out.report
    }

    emit:
    graphs           = ch_graphs
    embeddings       = ch_embeddings
    windows          = ch_windows
    quantization     = ch_quantization
    quantized_windows = ch_quantized_windows.map { meta, directory -> directory }
    database         = ch_built_database
    rerank_metrics   = ch_rerank_metrics
    seeds            = ch_seeds
    clusters         = ch_clusters
    cluster_members  = ch_cluster_members
    alignments       = ch_alignments
    alignment_text   = ch_alignment_text
    evd              = ch_published_evd
    plots_rnartist   = ch_plots_rnartist
    plots_r4rna      = ch_plots_r4rna
    plots_sw         = ch_plots_sw
    report           = ch_report
    versions         = ch_versions
}
