include { BUILD_RNA_GRAPHS } from '../modules/build_rna_graphs/main'
include { EMBED_RNA_GRAPHS } from '../modules/embed_rna_graphs/main'
include { GENERATE_WINDOWS } from '../modules/generate_windows/main'

workflow PREPARE_WINDOWS {
    take:
    structures
    prefix

    main:
    ch_versions = channel.empty()

    ch_shards = channel
        .fromPath(structures, checkIfExists: true)
        .splitText(by: params.shard_size, file: true, keepHeader: true)
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

    emit:
    graphs     = BUILD_RNA_GRAPHS.out.graphs
    embeddings = EMBED_RNA_GRAPHS.out.embeddings
    windows    = GENERATE_WINDOWS.out.windows
    versions   = ch_versions
}
