include { PREPARE_WINDOWS as PREPARE_DB }    from './prepare_windows'
include { PREPARE_WINDOWS as PREPARE_QUERY } from './prepare_windows'
include { BUILD_FAISS_INDEX }                from '../modules/build_faiss_index/main'
include { SEARCH_FAISS }                     from '../modules/search_faiss/main'

workflow GINFLOW {

    take:
    structures
    queries
    database
    reuse_windows

    main:
    ch_versions        = channel.empty()
    ch_graphs          = channel.empty()
    ch_embeddings      = channel.empty()
    ch_windows         = channel.empty()
    ch_built_database  = channel.empty()
    ch_search_database = channel.empty()
    ch_seeds           = channel.empty()

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
            ch_db_win.manifest.collect()
        )
        ch_versions        = ch_versions.mix(BUILD_FAISS_INDEX.out.versions)
        ch_built_database  = BUILD_FAISS_INDEX.out.database
        ch_search_database = BUILD_FAISS_INDEX.out.database
    }
    else if (database) {
        ch_search_database = channel.fromPath(database, checkIfExists: true, type: 'dir')
    }

    if (queries) {
        if (reuse_windows) {
            ch_query_windows = PREPARE_DB.out.windows
        }
        else {
            PREPARE_QUERY(queries, structures ? 'query' : '')
            ch_versions      = ch_versions.mix(PREPARE_QUERY.out.versions)
            ch_graphs        = ch_graphs.mix(PREPARE_QUERY.out.graphs)
            ch_embeddings    = ch_embeddings.mix(PREPARE_QUERY.out.embeddings)
            ch_windows       = ch_windows.mix(PREPARE_QUERY.out.windows)
            ch_query_windows = PREPARE_QUERY.out.windows
        }

        SEARCH_FAISS(ch_query_windows, ch_search_database.collect())
        ch_versions = ch_versions.mix(SEARCH_FAISS.out.versions)
        ch_seeds = SEARCH_FAISS.out.seeds
            .map { meta, tsv -> tsv }
            .collectFile(name: 'seeds.tsv', keepHeader: true, skip: 1, sort: true)
    }

    emit:
    graphs     = ch_graphs
    embeddings = ch_embeddings
    windows    = ch_windows
    database   = ch_built_database
    seeds      = ch_seeds
    versions   = ch_versions
}
