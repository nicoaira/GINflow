process RERANK_CANDIDATES {
    tag "${meta.id}"
    label 'process_medium'
    accelerator { params.exact_rerank_device == 'cuda' ? 1 : null }

    conda "${ task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml" }"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        (task.accelerator ? 'oras://community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:40ef78c1cb7c29cd' : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7920d351ee6307611f471013dcddff6d8a5c7ec7bf723d854e99a545376b63d8/data') :
        (task.accelerator ? 'community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:93cd6db656f6b1e4' : 'community.wave.seqera.io/library/python_numpy_faiss-cpu_mkl_libblas:078dd4f35c795d6e') }"

    input:
    tuple val(meta), path(candidates)
    path database
    path query_windows, stageAs: 'query_windows/*'
    path query_manifests, stageAs: 'query_manifests/*'

    output:
    tuple val(meta), path("*.seeds.tsv"), emit: seeds
    path "rerank.metrics.json", emit: metrics
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def workers = params.exact_rerank_workers > 0 ? params.exact_rerank_workers : task.cpus
    """
    rerank_candidates.py \\
        --candidates ${candidates} \\
        --database ${database} \\
        --query-windows query_windows/*.windows.npz \\
        --query-manifests query_manifests/*.windows.manifest.json \\
        --output ${prefix}.reranked.seeds.tsv \\
        --metrics rerank.metrics.json \\
        --output-k ${params.seed_k} \\
        --min-similarity ${params.seed_min_similarity} \\
        --batch-size ${params.exact_rerank_batch_size} \\
        --candidate-batch-size ${params.exact_rerank_candidate_batch_size} \\
        --workers ${workers} \\
        --device ${params.exact_rerank_device} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        ${task.accelerator ? 'cupy: \$(python3 -c "import cupy; print(cupy.__version__)")' : ''}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "query_id\\tquery_start\\tquery_end\\ttarget_id\\ttarget_start\\ttarget_end\\tscore\\trank" > ${prefix}.reranked.seeds.tsv
    echo '{"backend":"exact_original_window_rerank","device":"stub","n_seeds":0}' > rerank.metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        numpy: stub
    END_VERSIONS
    """
}
