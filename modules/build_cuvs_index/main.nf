process BUILD_CUVS_INDEX {
    tag "cuvs"
    label 'process_medium'
    accelerator 1

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:40ef78c1cb7c29cd' :
        'community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:93cd6db656f6b1e4' }"

    input:
    path windows, stageAs: 'windows/*'
    path manifests, stageAs: 'manifests/*'
    path embeddings, stageAs: 'embeddings/*'
    path metadata, stageAs: 'metadata/*'

    output:
    path "faiss", emit: database
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def n_lists  = params.cuvs_n_lists != null ? "--cuvs-n-lists ${params.cuvs_n_lists}" : ''
    def n_probes = params.cuvs_n_probes != null ? "--cuvs-n-probes ${params.cuvs_n_probes}" : ''
    """
    build_faiss.py \\
        --windows windows/*.windows.npz \\
        --manifests manifests/*.windows.manifest.json \\
        --embeddings embeddings/*.npz \\
        --graph-metadata metadata/*.json \\
        --outdir faiss \\
        --backend cuvs \\
        --cuvs-index-type ${params.cuvs_index} \\
        --cuvs-pq-bits ${params.cuvs_pq_bits} \\
        --cuvs-pq-dim ${params.cuvs_pq_dim} \\
        --cuvs-intermediate-graph-degree ${params.cuvs_intermediate_graph_degree} \\
        --cuvs-graph-degree ${params.cuvs_graph_degree} \\
        --cuvs-build-algo ${params.cuvs_build_algo} \\
        --cuvs-itopk-size ${params.cuvs_itopk_size} \\
        ${n_lists} \\
        ${n_probes} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cuvs: \$(python3 -c "from importlib.metadata import version; print(version('cuvs'))")
        cupy: \$(python3 -c "import cupy; print(cupy.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p faiss/cuvs
    touch faiss/cuvs/index.bin
    touch faiss/windows.tsv
    touch faiss/embeddings.npz
    touch faiss/records.tsv
    echo '{}' > faiss/meta.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cuvs: \$(python3 -c "from importlib.metadata import version; print(version('cuvs'))")
        cupy: \$(python3 -c "import cupy; print(cupy.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
