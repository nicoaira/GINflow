#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nicoaira/ginflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nicoaira/ginflow
----------------------------------------------------------------------------------------
*/

include { GINFLOW } from './workflows/ginflow'

def resolve_optional_path(value) {
    value ? file(value).toAbsolutePath().normalize().toString() : null
}

def validate_run_mode() {
    def has_input    = params.input as boolean
    def has_query    = params.query as boolean
    def has_database = params.database as boolean

    if (has_input && has_query && has_database) {
        error "Specify either --input (optionally with --query) or --query --database, not all three."
    }
    if (!has_input && !has_query) {
        error "Provide --input to build a database and/or --query to search."
    }
    if (has_query && !has_input && !has_database) {
        error "--query requires --database (existing FAISS directory) or --input (build then search)."
    }
    if (has_database && !has_query) {
        error "--database requires --query."
    }
}

workflow {
    main:
    validate_run_mode()

    def input_path    = resolve_optional_path(params.input)
    def query_path    = resolve_optional_path(params.query)
    def database_path = resolve_optional_path(params.database)
    def reuse_windows = input_path && query_path && input_path == query_path
    def evd_existing  = null
    if (params.database) {
        def bundled = file("${params.database}/evd.json")
        if (bundled.exists()) {
            evd_existing = bundled.toString()
        }
    }

    if (params.input) {
        file(params.input, checkIfExists: true)
    }
    if (params.query) {
        file(params.query, checkIfExists: true)
    }
    if (params.database) {
        def db = file(params.database, checkIfExists: true)
        if (!db.isDirectory()) {
            error "--database must be a directory containing index.faiss, windows.tsv, and meta.json."
        }
    }

    result = GINFLOW(
        params.input    ?: [],
        params.query    ?: [],
        params.database ?: [],
        reuse_windows,
        evd_existing
    )

    samples_ch = result.graphs
        .join(result.embeddings)
        .join(result.windows)
        .map { meta, tensors, sidecar, npz, manifest, windows, windows_manifest ->
            [
                id               : meta.id,
                graphs           : tensors,
                metadata         : sidecar,
                embeddings       : npz,
                manifest         : manifest,
                windows          : windows,
                windows_manifest : windows_manifest,
            ]
        }

    publish:
    samples          = samples_ch
    database         = result.database
    seeds            = result.seeds
    clusters         = result.clusters
    cluster_members  = result.cluster_members
    alignments       = result.alignments
    alignment_text   = result.alignment_text
    evd              = result.evd
}

output {
    samples {
        path { sample ->
            sample.graphs >> "graphs/${sample.id}/"
            sample.metadata >> "graphs/${sample.id}/"
            sample.embeddings >> "embeddings/${sample.id}/"
            sample.manifest >> "embeddings/${sample.id}/"
            sample.windows >> "windows/${sample.id}/"
            sample.windows_manifest >> "windows/${sample.id}/"
        }
        index {
            path 'samples.csv'
            header true
        }
    }
    database {
        path '.'
    }
    seeds {
        path '.'
    }
    clusters {
        path '.'
    }
    cluster_members {
        path '.'
    }
    alignments {
        path '.'
    }
    alignment_text {
        path '.'
    }
    evd {
        path { evd_file ->
            evd_file >> 'faiss/evd.json'
        }
    }
}
