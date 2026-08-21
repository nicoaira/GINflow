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

    ch_shards = channel
        .fromPath(structures, checkIfExists: true)
        .splitText(by: shard_size as Integer, file: true, keepHeader: true)
        .map { shard ->
            def id = prefix ? "${prefix}_${shard.baseName}" : shard.baseName
            tuple([id: id], shard)
        }

    BUILD_RNA_GRAPHS(ch_shards)
    ch_versions = ch_versions.mix(BUILD_RNA_GRAPHS.out.versions)

    EMBED_RNA_GRAPHS(BUILD_RNA_GRAPHS.out.graphs)
    ch_versions = ch_versions.mix(EMBED_RNA_GRAPHS.out.versions)

    GENERATE_WINDOWS(EMBED_RNA_GRAPHS.out.embeddings)
    ch_versions = ch_versions.mix(GENERATE_WINDOWS.out.versions)

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
            FIT_NODE_QUANTIZER(ch_all_embeddings, ch_all_manifests)
            ch_quantization = FIT_NODE_QUANTIZER.out.quantization
            ch_versions = ch_versions.mix(FIT_NODE_QUANTIZER.out.versions)
        }

        GENERATE_QUANTIZED_WINDOWS(ch_quantization)
        ch_versions = ch_versions.mix(GENERATE_QUANTIZED_WINDOWS.out.versions)
        ch_quantized_windows = GENERATE_QUANTIZED_WINDOWS.out.windows.map { directory ->
            tuple([id: prefix ?: 'input'], directory)
        }
        ch_quantized_npz = GENERATE_QUANTIZED_WINDOWS.out.npz
        ch_quantized_manifests = GENERATE_QUANTIZED_WINDOWS.out.manifests
    }

    emit:
    graphs     = BUILD_RNA_GRAPHS.out.graphs
    embeddings = EMBED_RNA_GRAPHS.out.embeddings
    windows    = GENERATE_WINDOWS.out.windows
    quantization = ch_quantization
    quantized_windows = ch_quantized_windows
    quantized_npz = ch_quantized_npz
    quantized_manifests = ch_quantized_manifests
    versions   = ch_versions
}
