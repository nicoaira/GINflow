include { BUILD_RNA_GRAPHS } from '../modules/build_rna_graphs/main'
include { EMBED_RNA_GRAPHS } from '../modules/embed_rna_graphs/main'

workflow GINFLOW {

    take:
    structures

    main:

    ch_versions = channel.empty()

    ch_shards = channel
        .fromPath(structures, checkIfExists: true)
        .splitText(by: params.shard_size, file: true, keepHeader: true)
        .map { shard ->
            tuple([id: shard.baseName], shard)
        }

    BUILD_RNA_GRAPHS(ch_shards)
    ch_versions = ch_versions.mix(BUILD_RNA_GRAPHS.out.versions)

    EMBED_RNA_GRAPHS(BUILD_RNA_GRAPHS.out.graphs)
    ch_versions = ch_versions.mix(EMBED_RNA_GRAPHS.out.versions)

    emit:
    graphs     = BUILD_RNA_GRAPHS.out.graphs
    embeddings = EMBED_RNA_GRAPHS.out.embeddings
    versions   = ch_versions
}
