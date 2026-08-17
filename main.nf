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
        error "--query requires --database (existing index directory) or --input (build then search)."
    }
    if (has_database && !has_query) {
        error "--database requires --query."
    }
}

def cli_has_flag(cli, token) {
    return cli.contains('--' + token + '=') || cli.contains('--' + token + ' ') || cli.endsWith('--' + token)
}

def parameter_on_cli(name) {
    def cli = (workflow.commandLine ?: '') as String
    def dashed = name.replace('_', '-')
    return cli_has_flag(cli, name) || (dashed != name && cli_has_flag(cli, dashed))
}

def parameter_explicit(name, default_value) {
    if (parameter_on_cli(name)) {
        return true
    }
    def value = params[name]
    if (default_value == null) {
        return value != null
    }
    return value != default_value
}

def normalize_index_library(raw) {
    def key = (raw as String)?.trim()?.toLowerCase()
    if (!key || key == 'faiss') {
        return 'faiss'
    }
    if (key in ['scann', 'scan']) {
        return 'scann'
    }
    if (key in ['ngt', 'anng']) {
        return 'ngt'
    }
    if (key == 'cuvs') {
        return 'cuvs'
    }
    error "--index must be faiss, scann, ngt, or cuvs, got '${raw}'. See docs/indexes.md."
}

def normalize_lowercase_choice(name, default_value, choices, aliases = [:]) {
    def raw = params[name]
    def value = raw == null ? default_value : (raw as String).trim().toLowerCase()
    if (!value) {
        value = default_value
    }
    if (aliases.containsKey(value)) {
        value = aliases[value]
    }
    if (!(value in choices)) {
        error "--${name} must be one of: ${choices.join(', ')} (lowercase), got '${raw}'. See docs/indexes.md."
    }
    params[name] = value
    return value
}

def normalize_parameter_values() {
    normalize_lowercase_choice('embed_device', 'cpu', ['cpu', 'cuda'])
    normalize_lowercase_choice('ngt_index', 'ngt', ['ngt', 'qg', 'qbg'], [
        anng: 'ngt',
        quantizedgraph: 'qg',
        'quantized-graph': 'qg',
        quantizedblobgraph: 'qbg',
        'quantized-blob-graph': 'qbg',
    ])
    normalize_lowercase_choice('cuvs_index', 'cagra', ['cagra', 'ivf', 'ivf-pq'], [
        'ivf-flat': 'ivf',
        ivfpq: 'ivf-pq',
    ])
    normalize_lowercase_choice('cuvs_build_algo', 'nn_descent', ['ivf_pq', 'nn_descent', 'iterative_cagra_search', 'ace'])
    if ((params.faiss_index as String)?.trim()?.equalsIgnoreCase('scann')) {
        error "--faiss_index scann is not valid. FAISS and ScaNN are different index libraries. Use --index scann. See docs/indexes.md."
    }
    normalize_lowercase_choice('faiss_index', 'flatip', ['flatip', 'flatl2', 'hnsw', 'ivfflat', 'lsh', 'sq', 'pq', 'ivfsq', 'ivfpq', 'ivfpqr'], [
        flat: 'flatip',
        indexflatip: 'flatip',
        indexflatl2: 'flatl2',
        hnswflat: 'hnsw',
        indexhnswflat: 'hnsw',
        indexivfflat: 'ivfflat',
        indexlsh: 'lsh',
        scalarquantizer: 'sq',
        indexscalarquantizer: 'sq',
        indexpq: 'pq',
        ivfscalarquantizer: 'ivfsq',
        indexivfscalarquantizer: 'ivfsq',
        indexivfpq: 'ivfpq',
        indexivfpqr: 'ivfpqr',
    ])
    normalize_lowercase_choice('faiss_sq_type', '8bit', ['8bit', '6bit', '4bit', 'fp16'], [
        '8': '8bit',
        qt_8bit: '8bit',
        '6': '6bit',
        qt_6bit: '6bit',
        '4': '4bit',
        qt_4bit: '4bit',
        qt_fp16: 'fp16',
    ])
    normalize_lowercase_choice('plot_backend', 'none', ['none', 'rnartistcore', 'r4rna', 'both'])
    normalize_lowercase_choice('report_theme', 'light', ['light', 'dark'])
}

def detect_database_library() {
    if (!params.database) {
        return null
    }
    def meta_file = file("${params.database}/meta.json")
    if (meta_file.exists()) {
        try {
            def meta = new groovy.json.JsonSlurperClassic().parseText(meta_file.text)
            def kind = meta.index_type as String
            def backend = meta.backend as String
            if (kind?.equalsIgnoreCase('ScaNN') || backend?.equalsIgnoreCase('scann')) {
                return 'scann'
            }
            if (backend?.equalsIgnoreCase('ngt') || ['NGT', 'QG', 'QBG'].contains(kind?.toUpperCase())) {
                return 'ngt'
            }
            if (backend?.equalsIgnoreCase('cuvs') || ['CAGRA', 'IVF', 'IVF-PQ'].contains(kind?.toUpperCase())) {
                return 'cuvs'
            }
            if (backend?.equalsIgnoreCase('faiss') || kind) {
                return 'faiss'
            }
        }
        catch (Exception ignored) {
            // Fall through to the on-disk artifact check.
        }
    }
    if (file("${params.database}/scann").exists()) {
        return 'scann'
    }
    if (file("${params.database}/ngt").exists()) {
        return 'ngt'
    }
    if (file("${params.database}/cuvs").exists()) {
        return 'cuvs'
    }
    if (file("${params.database}/index.faiss").exists()) {
        return 'faiss'
    }
    return null
}

def resolve_index_library() {
    def requested = normalize_index_library(params.index)
    def from_db = detect_database_library()
    if (from_db && parameter_on_cli('index') && requested != from_db) {
        error "--index ${requested} does not match the --database (${from_db}). Use --index ${from_db} or omit --index."
    }
    return from_db ?: requested
}

def faiss_index_kind() {
    def kind = (params.faiss_index as String)?.trim()?.toLowerCase()
    if (kind == 'scann') {
        error "--faiss_index scann is not valid. FAISS and ScaNN are different index libraries. Use --index scann. See docs/indexes.md."
    }
    def aliases = [
        flatip: 'FlatIP',
        flatl2: 'FlatL2',
        hnsw: 'HNSW',
        ivfflat: 'IVFFlat',
        lsh: 'LSH',
        sq: 'SQ',
        pq: 'PQ',
        ivfsq: 'IVFSQ',
        ivfpq: 'IVFPQ',
        ivfpqr: 'IVFPQR',
    ]
    return aliases[kind] ?: 'FlatIP'
}

def validate_faiss_gpu(library, kind) {
    if (!params.faiss_gpu) {
        return
    }
    if (library != 'faiss') {
        error "--faiss_gpu applies only to --index faiss, not --index ${library}."
    }
    def gpu_types = ['FlatIP', 'FlatL2', 'IVFFlat', 'IVFPQ', 'IVFSQ'] as Set
    // Search-only runs read the type from meta.json; Python rejects GPU-incompatible indexes.
    if (params.input && !gpu_types.contains(kind)) {
        error "--faiss_gpu is not supported for --faiss_index ${kind.toLowerCase()}. GPU FAISS indexes: ${gpu_types.collect { it.toLowerCase() }.sort().join(', ')}."
    }
    def profiles = workflow.profile.tokenize(',').collect { it.trim() }
    if (!profiles.contains('gpu')) {
        error "--faiss_gpu requires -profile gpu so BUILD_FAISS_INDEX and SEARCH_FAISS get the faiss-gpu image and NVIDIA runtime."
    }
}

def validate_cuvs_gpu(library) {
    if (library != 'cuvs') {
        return
    }
    def profiles = workflow.profile.tokenize(',').collect { it.trim() }
    if (!profiles.contains('gpu')) {
        error "--index cuvs requires -profile gpu so BUILD_CUVS_INDEX and SEARCH_CUVS get the cuVS image and NVIDIA runtime."
    }
}

def collect_if_unused(unused, name, default_value, applies) {
    if (!applies && parameter_explicit(name, default_value)) {
        unused.add('--' + name)
    }
}

def warn_unused_index_params(library, kind) {
    def ivf     = ['IVFFlat', 'IVFSQ', 'IVFPQ', 'IVFPQR'] as Set
    def pq      = ['PQ', 'IVFPQ', 'IVFPQR'] as Set
    def hnsw    = ['HNSW'] as Set
    def lsh     = ['LSH'] as Set
    def sq      = ['SQ', 'IVFSQ'] as Set
    def unused  = []

    if (library == 'scann') {
        collect_if_unused(unused, 'ngt_index', 'ngt', false)
        collect_if_unused(unused, 'cuvs_index', 'cagra', false)
        collect_if_unused(unused, 'cuvs_n_lists', null, false)
        collect_if_unused(unused, 'cuvs_n_probes', null, false)
        collect_if_unused(unused, 'cuvs_pq_bits', 8, false)
        collect_if_unused(unused, 'cuvs_pq_dim', 0, false)
        collect_if_unused(unused, 'cuvs_intermediate_graph_degree', 128, false)
        collect_if_unused(unused, 'cuvs_graph_degree', 64, false)
        collect_if_unused(unused, 'cuvs_build_algo', 'nn_descent', false)
        collect_if_unused(unused, 'cuvs_itopk_size', 64, false)
        collect_if_unused(unused, 'faiss_index', 'flatip', false)
        collect_if_unused(unused, 'faiss_gpu', false, false)
        collect_if_unused(unused, 'faiss_nlist', null, false)
        collect_if_unused(unused, 'faiss_nprobe', null, false)
        collect_if_unused(unused, 'faiss_pq_m', 16, false)
        collect_if_unused(unused, 'faiss_pq_nbits', 8, false)
        collect_if_unused(unused, 'faiss_pq_m_refine', 4, false)
        collect_if_unused(unused, 'faiss_hnsw_m', 32, false)
        collect_if_unused(unused, 'faiss_hnsw_ef_construction', 40, false)
        collect_if_unused(unused, 'faiss_hnsw_ef_search', 16, false)
        collect_if_unused(unused, 'faiss_lsh_nbits', null, false)
        collect_if_unused(unused, 'faiss_sq_type', '8bit', false)
    }
    else if (library == 'ngt') {
        collect_if_unused(unused, 'cuvs_index', 'cagra', false)
        collect_if_unused(unused, 'cuvs_n_lists', null, false)
        collect_if_unused(unused, 'cuvs_n_probes', null, false)
        collect_if_unused(unused, 'cuvs_pq_bits', 8, false)
        collect_if_unused(unused, 'cuvs_pq_dim', 0, false)
        collect_if_unused(unused, 'cuvs_intermediate_graph_degree', 128, false)
        collect_if_unused(unused, 'cuvs_graph_degree', 64, false)
        collect_if_unused(unused, 'cuvs_build_algo', 'nn_descent', false)
        collect_if_unused(unused, 'cuvs_itopk_size', 64, false)
        collect_if_unused(unused, 'faiss_index', 'flatip', false)
        collect_if_unused(unused, 'faiss_gpu', false, false)
        collect_if_unused(unused, 'faiss_nlist', null, false)
        collect_if_unused(unused, 'faiss_nprobe', null, false)
        collect_if_unused(unused, 'faiss_pq_m', 16, false)
        collect_if_unused(unused, 'faiss_pq_nbits', 8, false)
        collect_if_unused(unused, 'faiss_pq_m_refine', 4, false)
        collect_if_unused(unused, 'faiss_hnsw_m', 32, false)
        collect_if_unused(unused, 'faiss_hnsw_ef_construction', 40, false)
        collect_if_unused(unused, 'faiss_hnsw_ef_search', 16, false)
        collect_if_unused(unused, 'faiss_lsh_nbits', null, false)
        collect_if_unused(unused, 'faiss_sq_type', '8bit', false)
        collect_if_unused(unused, 'scann_reorder', 100, false)
        collect_if_unused(unused, 'scann_ah_dim', 2, false)
        collect_if_unused(unused, 'scann_anisotropic', 0.2, false)
        collect_if_unused(unused, 'scann_soar', false, false)
        collect_if_unused(unused, 'scann_leaves', null, false)
        collect_if_unused(unused, 'scann_leaves_to_search', null, false)
    }
    else if (library == 'cuvs') {
        collect_if_unused(unused, 'ngt_index', 'ngt', false)
        collect_if_unused(unused, 'faiss_index', 'flatip', false)
        collect_if_unused(unused, 'faiss_gpu', false, false)
        collect_if_unused(unused, 'faiss_nlist', null, false)
        collect_if_unused(unused, 'faiss_nprobe', null, false)
        collect_if_unused(unused, 'faiss_pq_m', 16, false)
        collect_if_unused(unused, 'faiss_pq_nbits', 8, false)
        collect_if_unused(unused, 'faiss_pq_m_refine', 4, false)
        collect_if_unused(unused, 'faiss_hnsw_m', 32, false)
        collect_if_unused(unused, 'faiss_hnsw_ef_construction', 40, false)
        collect_if_unused(unused, 'faiss_hnsw_ef_search', 16, false)
        collect_if_unused(unused, 'faiss_lsh_nbits', null, false)
        collect_if_unused(unused, 'faiss_sq_type', '8bit', false)
        collect_if_unused(unused, 'scann_reorder', 100, false)
        collect_if_unused(unused, 'scann_ah_dim', 2, false)
        collect_if_unused(unused, 'scann_anisotropic', 0.2, false)
        collect_if_unused(unused, 'scann_soar', false, false)
        collect_if_unused(unused, 'scann_leaves', null, false)
        collect_if_unused(unused, 'scann_leaves_to_search', null, false)
    }
    else {
        collect_if_unused(unused, 'ngt_index', 'ngt', false)
        collect_if_unused(unused, 'cuvs_index', 'cagra', false)
        collect_if_unused(unused, 'cuvs_n_lists', null, false)
        collect_if_unused(unused, 'cuvs_n_probes', null, false)
        collect_if_unused(unused, 'cuvs_pq_bits', 8, false)
        collect_if_unused(unused, 'cuvs_pq_dim', 0, false)
        collect_if_unused(unused, 'cuvs_intermediate_graph_degree', 128, false)
        collect_if_unused(unused, 'cuvs_graph_degree', 64, false)
        collect_if_unused(unused, 'cuvs_build_algo', 'nn_descent', false)
        collect_if_unused(unused, 'cuvs_itopk_size', 64, false)
        collect_if_unused(unused, 'scann_reorder', 100, false)
        collect_if_unused(unused, 'scann_ah_dim', 2, false)
        collect_if_unused(unused, 'scann_anisotropic', 0.2, false)
        collect_if_unused(unused, 'scann_soar', false, false)
        collect_if_unused(unused, 'scann_leaves', null, false)
        collect_if_unused(unused, 'scann_leaves_to_search', null, false)
        collect_if_unused(unused, 'faiss_nlist', null, ivf.contains(kind))
        collect_if_unused(unused, 'faiss_nprobe', null, ivf.contains(kind))
        collect_if_unused(unused, 'faiss_pq_m', 16, pq.contains(kind))
        collect_if_unused(unused, 'faiss_pq_nbits', 8, pq.contains(kind))
        collect_if_unused(unused, 'faiss_pq_m_refine', 4, kind == 'IVFPQR')
        collect_if_unused(unused, 'faiss_hnsw_m', 32, hnsw.contains(kind))
        collect_if_unused(unused, 'faiss_hnsw_ef_construction', 40, hnsw.contains(kind))
        collect_if_unused(unused, 'faiss_hnsw_ef_search', 16, hnsw.contains(kind))
        collect_if_unused(unused, 'faiss_lsh_nbits', null, lsh.contains(kind))
        collect_if_unused(unused, 'faiss_sq_type', '8bit', sq.contains(kind))
    }

    if (!unused.isEmpty()) {
        def suffix = library == 'faiss' ? (' / --faiss_index ' + kind.toLowerCase()) : ''
        def plural = unused.size() == 1 ? '' : 's'
        log.warn "Unused index parameter${plural} for --index ${library}${suffix} (ignored): ${unused.join(', ')}. See docs/indexes.md."
    }
}

workflow {
    main:
    validate_run_mode()
    normalize_parameter_values()
    def library = resolve_index_library()
    def kind    = faiss_index_kind()
    params.index = library
    validate_faiss_gpu(library, kind)
    validate_cuvs_gpu(library)
    warn_unused_index_params(library, kind)

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
        error "--database must be a directory containing windows.tsv, meta.json, and a vector index."
        }
    }

    result = GINFLOW(
        params.input    ?: [],
        params.query    ?: [],
        params.database ?: [],
        reuse_windows,
        evd_existing,
        library
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
    plots_rnartist   = result.plots_rnartist
    plots_r4rna      = result.plots_r4rna
    plots_sw         = result.plots_sw
    report           = result.report
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
    plots_rnartist {
        path 'plots/rnartistcore'
    }
    plots_r4rna {
        path 'plots/r4rna'
    }
    plots_sw {
        path 'plots/sw'
    }
    report {
        path '.'
    }
}
