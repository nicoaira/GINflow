include { PREPARE_WINDOWS as PREPARE_DB }    from './prepare_windows'
include { PREPARE_WINDOWS as PREPARE_QUERY } from './prepare_windows'
include { BUILD_FAISS_INDEX }                from '../modules/build_faiss_index/main'
include { SEARCH_FAISS }                     from '../modules/search_faiss/main'
include { CLUSTER_SEEDS }                    from '../modules/cluster_seeds/main'
include { ALIGN_CLUSTERS }                   from '../modules/align_clusters/main'
include { ESTIMATE_EVD as ESTIMATE_EVD_BUILD } from '../modules/estimate_evd/main'
include { ESTIMATE_EVD as ESTIMATE_EVD_QUERY } from '../modules/estimate_evd/main'
include { DRAW_RNARTISTCORE }                from '../modules/draw_rnartistcore/main'
include { DRAW_R4RNA }                       from '../modules/draw_r4rna/main'
include { WRITE_REPORT }                     from '../modules/write_report/main'

workflow GINFLOW {

    take:
    structures
    queries
    database
    reuse_windows
    evd_existing

    main:
    ch_versions        = channel.empty()
    ch_graphs          = channel.empty()
    ch_embeddings      = channel.empty()
    ch_windows         = channel.empty()
    ch_built_database  = channel.empty()
    ch_search_database = channel.empty()
    ch_seeds           = channel.empty()
    ch_clusters        = channel.empty()
    ch_cluster_members = channel.empty()
    ch_alignments      = channel.empty()
    ch_alignment_text  = channel.empty()
    ch_query_embeddings = channel.empty()
    ch_query_metadata   = channel.empty()
    ch_evd              = channel.empty()
    ch_plots_rnartist   = channel.empty()
    ch_plots_r4rna      = channel.empty()
    ch_report           = channel.empty()
    alignment_params    = file("${projectDir}/assets/alignment.json", checkIfExists: true)

    if (structures) {
        PREPARE_DB(structures, (queries && !reuse_windows) ? 'db' : '')
        ch_versions   = ch_versions.mix(PREPARE_DB.out.versions)
        ch_graphs     = ch_graphs.mix(PREPARE_DB.out.graphs)
        ch_embeddings = ch_embeddings.mix(PREPARE_DB.out.embeddings)
        ch_windows    = ch_windows.mix(PREPARE_DB.out.windows)

        PREPARE_DB.out.windows
            .multiMap { meta, npz, manifest ->
                npz: npz
                manifest: manifest
            }
            .set { ch_db_win }

        BUILD_FAISS_INDEX(
            ch_db_win.npz.collect(),
            ch_db_win.manifest.collect(),
            PREPARE_DB.out.embeddings.map { meta, npz, manifest -> npz }.collect(),
            PREPARE_DB.out.graphs.map { meta, tensors, sidecar -> sidecar }.collect()
        )
        ch_versions        = ch_versions.mix(BUILD_FAISS_INDEX.out.versions)
        ch_built_database  = BUILD_FAISS_INDEX.out.database
        ch_search_database = BUILD_FAISS_INDEX.out.database

        ESTIMATE_EVD_BUILD(BUILD_FAISS_INDEX.out.database, alignment_params)
        ch_versions = ch_versions.mix(ESTIMATE_EVD_BUILD.out.versions)
        ch_evd      = ESTIMATE_EVD_BUILD.out.evd
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
        }
    }

    if (queries) {
        if (reuse_windows) {
            ch_query_windows    = PREPARE_DB.out.windows
            ch_query_embeddings = PREPARE_DB.out.embeddings
            ch_query_metadata   = PREPARE_DB.out.graphs
        }
        else {
            PREPARE_QUERY(queries, structures ? 'query' : '')
            ch_versions         = ch_versions.mix(PREPARE_QUERY.out.versions)
            ch_graphs           = ch_graphs.mix(PREPARE_QUERY.out.graphs)
            ch_embeddings       = ch_embeddings.mix(PREPARE_QUERY.out.embeddings)
            ch_windows          = ch_windows.mix(PREPARE_QUERY.out.windows)
            ch_query_windows    = PREPARE_QUERY.out.windows
            ch_query_embeddings = PREPARE_QUERY.out.embeddings
            ch_query_metadata   = PREPARE_QUERY.out.graphs
        }

        SEARCH_FAISS(ch_query_windows, ch_search_database.collect())
        ch_versions = ch_versions.mix(SEARCH_FAISS.out.versions)
        ch_seeds = SEARCH_FAISS.out.seeds
            .map { meta, tsv -> tsv }
            .collectFile(name: 'seeds.tsv', keepHeader: true, skip: 1, sort: true)

        CLUSTER_SEEDS(ch_seeds)
        ch_versions        = ch_versions.mix(CLUSTER_SEEDS.out.versions)
        ch_clusters        = CLUSTER_SEEDS.out.clusters
        ch_cluster_members = CLUSTER_SEEDS.out.members

        ALIGN_CLUSTERS(
            CLUSTER_SEEDS.out.clusters,
            ch_query_embeddings.map { meta, npz, manifest -> npz }.collect(),
            ch_query_metadata.map { meta, tensors, sidecar -> sidecar }.collect(),
            ch_search_database.collect(),
            alignment_params,
            ch_evd.collect()
        )
        ch_versions       = ch_versions.mix(ALIGN_CLUSTERS.out.versions)
        ch_alignments     = ALIGN_CLUSTERS.out.alignments
        ch_alignment_text = ALIGN_CLUSTERS.out.text

        ch_report_rn = channel.fromPath("${projectDir}/assets/no_plots_rnartist", checkIfExists: true, type: 'dir')
        ch_report_r4 = channel.fromPath("${projectDir}/assets/no_plots_r4rna", checkIfExists: true, type: 'dir')
        if (params.plot_backend in ['rnartistcore', 'both']) {
            DRAW_RNARTISTCORE(ALIGN_CLUSTERS.out.alignments)
            ch_versions       = ch_versions.mix(DRAW_RNARTISTCORE.out.versions)
            ch_plots_rnartist = DRAW_RNARTISTCORE.out.plots
            ch_report_rn      = DRAW_RNARTISTCORE.out.plots
        }
        if (params.plot_backend in ['r4rna', 'both']) {
            DRAW_R4RNA(ALIGN_CLUSTERS.out.alignments)
            ch_versions    = ch_versions.mix(DRAW_R4RNA.out.versions)
            ch_plots_r4rna = DRAW_R4RNA.out.plots
            ch_report_r4   = DRAW_R4RNA.out.plots
        }

        WRITE_REPORT(
            ALIGN_CLUSTERS.out.alignments,
            ALIGN_CLUSTERS.out.text,
            ch_evd.collect(),
            CLUSTER_SEEDS.out.clusters,
            ch_seeds,
            ch_report_rn,
            ch_report_r4
        )
        ch_versions = ch_versions.mix(WRITE_REPORT.out.versions)
        ch_report   = WRITE_REPORT.out.report
    }

    emit:
    graphs           = ch_graphs
    embeddings       = ch_embeddings
    windows          = ch_windows
    database         = ch_built_database
    seeds            = ch_seeds
    clusters         = ch_clusters
    cluster_members  = ch_cluster_members
    alignments       = ch_alignments
    alignment_text   = ch_alignment_text
    evd              = ch_evd
    plots_rnartist   = ch_plots_rnartist
    plots_r4rna      = ch_plots_r4rna
    report           = ch_report
    versions         = ch_versions
}
