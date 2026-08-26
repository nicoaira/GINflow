include { BUILD_RNA_GRAPHS } from '../modules/build_rna_graphs/main'
include { EMBED_RNA_GRAPHS } from '../modules/embed_rna_graphs/main'
include { GENERATE_WINDOWS } from '../modules/generate_windows/main'
include { FIT_NODE_QUANTIZER; APPLY_NODE_QUANTIZER } from '../modules/quantize_nodes/main'
include { GENERATE_QUANTIZED_WINDOWS } from '../modules/generate_quantized_windows/main'

workflow PREPARE_WINDOWS {
    take:
    structures
    prefix
    shard_size
    quantize
    quantization_source

    main:
    ch_versions = channel.empty()
    ch_quantization = channel.empty()
    ch_quantized_windows = channel.empty()
    ch_quantized_npz = channel.empty()
    ch_quantized_manifests = channel.empty()
    ch_window_counts = channel.empty()
    ch_node_counts = channel.empty()
    ch_quantized_window_stats = channel.empty()

    ch_shards = channel
        .fromPath(structures, checkIfExists: true)
        .splitText(by: shard_size as Integer, file: true, keepHeader: true)
        .map { shard ->
            def id = prefix ? "${prefix}_${shard.baseName}" : shard.baseName
            tuple([id: id], shard)
        }

    BUILD_RNA_GRAPHS(ch_shards)
    ch_versions = ch_versions.mix(BUILD_RNA_GRAPHS.out.versions)

    ch_graphs_with_stats = BUILD_RNA_GRAPHS.out.graphs.map { meta, tensors, sidecar ->
        def stats = [node_count: 0L, edge_count: 0L, record_count: 0L, graph_bytes: tensors.size() as long]
        try {
            def payload = new groovy.json.JsonSlurperClassic().parseText(sidecar.text)
            stats.node_count = (payload.node_count ?: 0L) as long
            stats.edge_count = (payload.edge_count ?: 0L) as long
            stats.record_count = (payload.record_count ?: 0L) as long
        }
        catch (Exception _ignored) {
            // Stub metadata is empty.
        }
        tuple(meta + [resource_stats: stats], tensors, sidecar)
    }

    EMBED_RNA_GRAPHS(ch_graphs_with_stats)
    ch_versions = ch_versions.mix(EMBED_RNA_GRAPHS.out.versions)

    GENERATE_WINDOWS(EMBED_RNA_GRAPHS.out.embeddings)
    ch_versions = ch_versions.mix(GENERATE_WINDOWS.out.versions)

    ch_window_counts = GENERATE_WINDOWS.out.windows
        .map { meta, windows, manifest ->
            try {
                def payload = new groovy.json.JsonSlurperClassic().parseText(manifest.text)
                (payload.n_windows ?: 0) as long
            }
            catch (Exception _ignored) {
                0L
            }
        }
        .collect()
        .map { counts -> counts.sum() as long }

    ch_quantized_window_stats = EMBED_RNA_GRAPHS.out.embeddings
        .map { meta, npz, manifest ->
            def stats = [n_nodes: 0L, max_record_nodes: 0L, n_windows: 0L, max_record_windows: 0L]
            def window_size = Math.max(1L, params.window_size as long)
            def stride = Math.max(1L, params.window_stride as long)
            try {
                def payload = new groovy.json.JsonSlurperClassic().parseText(manifest.text)
                (payload.records ?: []).each { record ->
                    def shape = record.shape instanceof Collection ? record.shape : []
                    def length = Math.max(
                        0L,
                        (shape ? shape[0] : (record.core_length ?: record.node_count ?: record.length ?: 0L)) as long
                    )
                    def windows = length < window_size
                        ? 0L
                        : ((length - window_size) / stride) as long + 1L
                    stats.n_nodes += length
                    stats.max_record_nodes = Math.max(stats.max_record_nodes as long, length)
                    stats.n_windows += windows
                    stats.max_record_windows = Math.max(stats.max_record_windows as long, windows)
                }
            }
            catch (Exception _ignored) {
                // Stub manifests are empty.
            }
            stats
        }
        .collect()
        .map { shards ->
            [
                n_nodes: (shards.sum { shard_stats -> shard_stats.n_nodes as long } ?: 0L) as long,
                max_shard_nodes: (shards.collect { shard_stats -> shard_stats.n_nodes as long }.max() ?: 0L) as long,
                max_record_nodes: (shards.collect { shard_stats -> shard_stats.max_record_nodes as long }.max() ?: 0L) as long,
                n_windows: (shards.sum { shard_stats -> shard_stats.n_windows as long } ?: 0L) as long,
                max_shard_windows: (shards.collect { shard_stats -> shard_stats.n_windows as long }.max() ?: 0L) as long,
                max_record_windows: (shards.collect { shard_stats -> shard_stats.max_record_windows as long }.max() ?: 0L) as long,
            ]
        }

    ch_node_counts = ch_quantized_window_stats.map { stats -> stats.n_nodes as long }

    if (quantize) {
        ch_all_embeddings = EMBED_RNA_GRAPHS.out.embeddings
            .map { meta, npz, manifest -> npz }
            .collect()
        ch_all_manifests = EMBED_RNA_GRAPHS.out.embeddings
            .map { meta, npz, manifest -> manifest }
            .collect()

        if (quantization_source) {
            APPLY_NODE_QUANTIZER(ch_all_embeddings, ch_all_manifests, quantization_source)
            ch_quantization = APPLY_NODE_QUANTIZER.out.quantization
            ch_versions = ch_versions.mix(APPLY_NODE_QUANTIZER.out.versions)
        }
        else {
            FIT_NODE_QUANTIZER(ch_all_embeddings, ch_all_manifests, ch_node_counts)
            ch_quantization = FIT_NODE_QUANTIZER.out.quantization
            ch_versions = ch_versions.mix(FIT_NODE_QUANTIZER.out.versions)
        }

        GENERATE_QUANTIZED_WINDOWS(ch_quantization, ch_quantized_window_stats)
        ch_versions = ch_versions.mix(GENERATE_QUANTIZED_WINDOWS.out.versions)
        ch_window_counts = GENERATE_QUANTIZED_WINDOWS.out.summary.map { summary ->
            def payload = new groovy.json.JsonSlurperClassic().parseText(summary.text)
            (payload.n_windows ?: 0) as long
        }
        ch_quantized_windows = GENERATE_QUANTIZED_WINDOWS.out.windows.map { directory ->
            tuple([id: prefix ?: 'input'], directory)
        }
        ch_quantized_npz = GENERATE_QUANTIZED_WINDOWS.out.npz
        ch_quantized_manifests = GENERATE_QUANTIZED_WINDOWS.out.manifests
    }

    emit:
    graphs     = ch_graphs_with_stats
    embeddings = EMBED_RNA_GRAPHS.out.embeddings
    windows    = GENERATE_WINDOWS.out.windows
    quantization = ch_quantization
    quantized_windows = ch_quantized_windows
    quantized_npz = ch_quantized_npz
    quantized_manifests = ch_quantized_manifests
    window_counts = ch_window_counts
    node_counts = ch_node_counts
    versions   = ch_versions
}
