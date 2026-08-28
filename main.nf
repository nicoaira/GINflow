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

def cli_has_bare_flag(token) {
    return ((workflow.commandLine ?: '') as String).tokenize().contains('--' + token)
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
    if (key in ['cagra', 'cuvs']) {
        return 'cagra'
    }
    if (key == 'ivf') {
        return 'ivf'
    }
    if (key in ['hnswlib', 'hnsw']) {
        return 'hnswlib'
    }
    error "--index must be faiss, cagra, ivf, or hnswlib, got '${raw}'. See docs/indexes.md."
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

def normalize_device_choice(name, default_value, aliases = [:]) {
    def raw = params[name]
    def value = raw == null ? default_value : (raw as String).trim().toLowerCase()
    if (!value) {
        value = default_value
    }
    if (aliases.containsKey(value)) {
        def canonical = aliases[value]
        log.warn "--${name} value '${value}' is deprecated; use '${canonical}'."
        value = canonical
    }
    if (!(value in ['cpu', 'gpu'])) {
        error "--${name} must be one of: cpu, gpu (lowercase), got '${raw}'."
    }
    return value
}

def normalize_boolean_param(name, default_value = false) {
    // Do not write back to params: Nextflow 26 ignores those assignments.
    try {
        return BooleanParam.parse(params[name], default_value)
    }
    catch (IllegalArgumentException ignored) {
        error "--${name} must be true or false, got '${params[name]}'."
    }
}

def normalize_faiss_device() {
    def selected = normalize_device_choice('faiss_device', 'cpu', [cuda: 'gpu'])
    def legacy_gpu = normalize_boolean_param('faiss_gpu', false)
    def legacy_explicit = parameter_on_cli('faiss_gpu') || params.faiss_gpu != false
    if (!legacy_explicit) {
        return selected
    }

    def legacy_selected = legacy_gpu ? 'gpu' : 'cpu'
    def new_explicit = parameter_on_cli('faiss_device') || selected != 'cpu'
    if (new_explicit && selected != legacy_selected) {
        error "--faiss_device ${selected} conflicts with deprecated --faiss_gpu ${legacy_gpu}. Use only --faiss_device cpu|gpu."
    }
    if (!new_explicit) {
        selected = legacy_selected
    }
    log.warn "--faiss_gpu is deprecated; use --faiss_device ${legacy_selected}."
    return selected
}

def faiss_device_is_gpu() {
    def raw = params.faiss_device
    def selected = raw == null ? 'cpu' : (raw as String).trim().toLowerCase()
    selected = selected == 'cuda' ? 'gpu' : selected
    def legacy_gpu = normalize_boolean_param('faiss_gpu', false)
    def legacy_explicit = parameter_on_cli('faiss_gpu') || params.faiss_gpu != false
    (selected == 'gpu') || (legacy_explicit && legacy_gpu)
}

def search_device_is_gpu() {
    def raw = params.search_device
    def selected = raw == null ? 'gpu' : (raw as String).trim().toLowerCase()
    !(selected == 'cpu')
}

def normalize_parameter_values() {
    normalize_device_choice('embed_device', 'cpu', [cuda: 'gpu'])
    normalize_lowercase_choice('quantize', 'none', ['none', 'sq', 'pq', 'opq'])
    normalize_device_choice('search_device', 'gpu', [auto: 'gpu', cuda: 'gpu'])
    normalize_lowercase_choice('cagra_build_algo', 'nn_descent', ['ivf_pq', 'nn_descent', 'iterative_cagra_search', 'ace'])
    normalize_lowercase_choice('faiss_index', 'flatip', ['flatip', 'flatl2', 'hnsw', 'ivfflat'], [
        flat: 'flatip',
        indexflatip: 'flatip',
        indexflatl2: 'flatl2',
        hnswflat: 'hnsw',
        indexhnswflat: 'hnsw',
        indexivfflat: 'ivfflat',
    ])
    normalize_lowercase_choice('plot_backend', 'none', ['none', 'rnartistcore', 'r4rna', 'both'])
    normalize_lowercase_choice('report_theme', 'light', ['light', 'dark'])
    normalize_boolean_param('hnswlib_rerank', false)
    normalize_boolean_param('exact_rerank', true)
    normalize_boolean_param('cagra_to_hnsw', false)
    normalize_boolean_param('save_graphs', false)
    normalize_boolean_param('save_embeddings', false)
    normalize_boolean_param('save_windows', false)
    normalize_boolean_param('save_quantized_windows', false)
    normalize_device_choice('exact_rerank_device', 'cpu', [cuda: 'gpu'])
    normalize_faiss_device()
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
            if (backend?.equalsIgnoreCase('cagra') || kind?.equalsIgnoreCase('PQ_CAGRA') || kind?.equalsIgnoreCase('CAGRA')) {
                return 'cagra'
            }
            if (backend?.equalsIgnoreCase('cuvs') && kind?.equalsIgnoreCase('IVF')) {
                return 'ivf'
            }
            if (backend?.equalsIgnoreCase('cuvs') || kind?.equalsIgnoreCase('CAGRA')) {
                return 'cagra'
            }
            if (backend?.equalsIgnoreCase('hnswlib') || ['HNSWLIB', 'HNSWLIB_PQ'].contains(kind?.toUpperCase())) {
                return 'hnswlib'
            }
            if (backend?.equalsIgnoreCase('faiss') || kind) {
                return 'faiss'
            }
        }
        catch (Exception ignored) {
            // Fall through to the on-disk artifact check.
        }
    }
    if (file("${params.database}/cagra.index").exists()) {
        return 'cagra'
    }
    if (file("${params.database}/cuvs").exists()) {
        def cuvs_meta_file = file("${params.database}/meta.json")
        if (cuvs_meta_file.exists()) {
            try {
                def cuvs_meta = new groovy.json.JsonSlurperClassic().parseText(cuvs_meta_file.text)
                if ((cuvs_meta.index_type as String)?.equalsIgnoreCase('IVF')) {
                    return 'ivf'
                }
            }
            catch (Exception ignored) {
            }
        }
        return 'cagra'
    }
    if (file("${params.database}/index.bin").exists()) {
        return 'hnswlib'
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
    def aliases = [
        flat: 'FlatIP',
        flatip: 'FlatIP',
        flatl2: 'FlatL2',
        indexflatip: 'FlatIP',
        indexflatl2: 'FlatL2',
        hnsw: 'HNSW',
        hnswflat: 'HNSW',
        indexhnswflat: 'HNSW',
        ivfflat: 'IVFFlat',
        indexivfflat: 'IVFFlat',
    ]
    return aliases[kind] ?: 'FlatIP'
}

def validate_faiss_device(library, kind) {
    if (!faiss_device_is_gpu()) {
        return
    }
    if (library != 'faiss') {
        error "--faiss_device gpu applies only to --index faiss, not --index ${library}."
    }
    def gpu_types = ['FlatIP', 'FlatL2'] as Set
    if (params.input && !gpu_types.contains(kind)) {
        def hint = kind == 'IVFFlat' \
            ? ' FAISS IVF is CPU-only; use --index ivf for GPU IVF-Flat, or --index cagra for a GPU graph.' \
            : ' FAISS HNSW is CPU-only; use --index cagra for a GPU graph.'
        error "--faiss_device gpu is not supported for --faiss_index ${kind.toLowerCase()}. GPU FAISS indexes: flatip, flatl2.${hint}"
    }
}

def validate_cagra_device(library) {
    def converted = BooleanParam.parse(params.cagra_to_hnsw, false)
    def cpu_search = !search_device_is_gpu() || converted
    if (library == 'ivf' && params.query && cpu_search) {
        error "--index ivf search requires --search_device gpu; cuVS IVF has no CPU search path."
    }
}

def validate_quantize_index(library, quantize) {
    def pq = quantize in ['pq', 'opq']
    def standard = quantize in ['none', 'sq']
    if (pq && !(library in ['cagra', 'hnswlib'])) {
        error "--quantize ${quantize} cannot be used with --index ${library}. Use --index cagra (GPU graph, recommended) or --index hnswlib (CPU-only build). See docs/indexes.md."
    }
    if (standard && library == 'hnswlib') {
        error "--index hnswlib is the custom-distance PQ/OPQ CPU path. For uncompressed or SQ windows use --index faiss or --index cagra."
    }
    if (pq && library == 'hnswlib' && params.input) {
        log.warn "CAGRA graphs usually outperform CPU hnswlib for PQ/OPQ and can be converted for CPU search with --cagra_to_hnsw true. Use --index hnswlib only when no GPU is available at build time. See docs/indexes.md."
    }
}

def validate_search_params() {
    if ((params.candidate_k as Integer) < (params.seed_k as Integer)) {
        error "--candidate_k must be >= --seed_k."
    }
    if ((params.cagra_graph_degree as Integer) > (params.cagra_intermediate_graph_degree as Integer)) {
        error "--cagra_graph_degree must be <= --cagra_intermediate_graph_degree."
    }
}

def collect_if_unused(unused, name, default_value, applies) {
    if (!applies && parameter_explicit(name, default_value)) {
        unused.add('--' + name)
    }
}

def warn_unused_index_params(library, kind) {
    def unused = []
    def faiss = library == 'faiss'
    def cagra = library == 'cagra'
    def ivf = library == 'ivf'
    def hnswlib = library == 'hnswlib'
    collect_if_unused(unused, 'faiss_index', 'flatip', faiss)
    collect_if_unused(unused, 'faiss_device', 'cpu', faiss)
    collect_if_unused(unused, 'faiss_gpu', false, faiss)
    collect_if_unused(unused, 'faiss_nlist', null, faiss && kind == 'IVFFlat')
    collect_if_unused(unused, 'faiss_nprobe', null, faiss && kind == 'IVFFlat')
    collect_if_unused(unused, 'faiss_hnsw_m', 32, faiss && kind == 'HNSW')
    collect_if_unused(unused, 'faiss_hnsw_ef_construction', 200, faiss && kind == 'HNSW')
    collect_if_unused(unused, 'faiss_hnsw_ef_search', 64, faiss && kind == 'HNSW')
    collect_if_unused(unused, 'cuvs_n_lists', null, ivf)
    collect_if_unused(unused, 'cuvs_n_probes', null, ivf)
    collect_if_unused(unused, 'cagra_intermediate_graph_degree', 128, cagra)
    collect_if_unused(unused, 'cagra_graph_degree', 64, cagra)
    collect_if_unused(unused, 'cagra_build_algo', 'nn_descent', cagra && params.quantize in ['none', 'sq'])
    collect_if_unused(unused, 'cagra_nndescent_iterations', 6, cagra && params.quantize in ['pq', 'opq'])
    collect_if_unused(unused, 'cagra_itopk_size', 64, cagra)
    collect_if_unused(unused, 'cagra_to_hnsw', false, cagra)
    collect_if_unused(unused, 'hnswlib_m', 32, hnswlib)
    collect_if_unused(unused, 'hnswlib_ef_construction', 200, hnswlib)
    collect_if_unused(unused, 'hnswlib_ef_search', 200, hnswlib)
    collect_if_unused(unused, 'hnswlib_random_seed', 1, hnswlib)
    collect_if_unused(unused, 'hnswlib_num_threads', 0, hnswlib)
    collect_if_unused(unused, 'pq_m', 16, params.quantize in ['pq', 'opq'])
    collect_if_unused(unused, 'pq_nbits', 4, params.quantize in ['pq', 'opq'])

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
    validate_faiss_device(library, kind)
    validate_cagra_device(library)
    validate_quantize_index(library, params.quantize)
    validate_search_params()
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
        library,
        params.quantize
    )

    graph_shards_ch = result.graphs
        .map { meta, tensors, sidecar ->
            [
                id               : meta.id,
                graphs           : tensors,
                metadata         : sidecar,
            ]
        }

    embedding_shards_ch = result.embeddings
        .map { meta, npz, manifest ->
            [
                id               : meta.id,
                embeddings       : npz,
                manifest         : manifest,
            ]
        }

    window_shards_ch = result.windows
        .map { meta, windows, windows_manifest ->
            [
                id               : meta.id,
                windows          : windows,
                windows_manifest : windows_manifest,
            ]
        }

    publish:
    graphs_shards    = BooleanParam.parse(params.save_graphs, false)
        ? graph_shards_ch
        : channel.empty()
    embeddings_shards = BooleanParam.parse(params.save_embeddings, false)
        ? embedding_shards_ch
        : channel.empty()
    windows_shards   = BooleanParam.parse(params.save_windows, false)
        ? window_shards_ch
        : channel.empty()
    database         = result.database
    quantized_windows = BooleanParam.parse(params.save_quantized_windows, false)
        ? result.quantized_windows
        : channel.empty()
    rerank_metrics   = result.rerank_metrics
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
    graphs_shards {
        path { sample ->
            sample.graphs >> 'graphs_shards/' + sample.id + '/'
            sample.metadata >> 'graphs_shards/' + sample.id + '/'
        }
    }
    embeddings_shards {
        path { sample ->
            sample.embeddings >> 'embeddings_shards/' + sample.id + '/'
            sample.manifest >> 'embeddings_shards/' + sample.id + '/'
        }
    }
    windows_shards {
        path { sample ->
            sample.windows >> 'windows_shards/' + sample.id + '/'
            sample.windows_manifest >> 'windows_shards/' + sample.id + '/'
        }
    }
    database {
        path { directory -> directory >> 'index' }
    }
    quantized_windows {
        path { directory -> directory >> 'windows_quantized' }
    }
    rerank_metrics {
        path { metrics -> metrics >> 'seeds/rerank_metrics.json' }
    }
    seeds {
        path { seed_file -> seed_file >> 'seeds/seeds.tsv' }
    }
    clusters {
        path { cluster_file -> cluster_file >> 'seeds/clusters.tsv' }
    }
    cluster_members {
        path { member_file -> member_file >> 'seeds/cluster_members.tsv' }
    }
    alignments {
        path '.'
    }
    alignment_text {
        path '.'
    }
    evd {
        path { evd_file ->
            evd_file >> 'index/evd.json'
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
