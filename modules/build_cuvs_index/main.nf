process BUILD_CUVS_INDEX {
    tag "cuvs"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:40ef78c1cb7c29cd' :
        'community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:93cd6db656f6b1e4' }"

    input:
    path windows, stageAs: 'windows/*'
    path manifests, stageAs: 'manifests/*'
    path embeddings, stageAs: 'embeddings/*'
    path metadata, stageAs: 'metadata/*'
    path quantization, stageAs: 'quantization'
    val n_windows
    val n_nodes

    output:
    path "index", emit: database
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def n_lists  = params.cuvs_n_lists != null ? "--cuvs-n-lists ${params.cuvs_n_lists}" : ''
    def n_probes = params.cuvs_n_probes != null ? "--cuvs-n-probes ${params.cuvs_n_probes}" : ''
    def quantization_arg = params.quantize == 'sq' ? '--quantization quantization' : ''
    def cagra_to_hnsw = params.cagra_to_hnsw
    """
    build_faiss.py \\
        --windows windows/*.windows.npz \\
        --manifests manifests/*.windows.manifest.json \\
        --embeddings embeddings/*.npz \\
        --graph-metadata metadata/*.json \\
        --outdir index \\
        --backend cuvs \\
        --cuvs-index-type ${params.index == 'ivf' ? 'ivf' : 'cagra'} \\
        --cuvs-intermediate-graph-degree ${params.cagra_intermediate_graph_degree} \\
        --cuvs-graph-degree ${params.cagra_graph_degree} \\
        --cuvs-build-algo ${params.cagra_build_algo} \\
        --cuvs-itopk-size ${params.cagra_itopk_size} \\
        ${cagra_to_hnsw ? '--cagra-to-hnsw' : ''} \\
        ${n_lists} \\
        ${n_probes} \\
        ${quantization_arg} \\
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
    mkdir -p index/cuvs
    touch index/cuvs/index.bin
    touch index/windows.tsv
    touch index/embeddings.npz
    touch index/records.tsv
    if [ "${params.quantize}" = "sq" ]; then
        mkdir -p index/quantization
        cp -a quantization/. index/quantization/
    fi
    echo '{"backend":"cuvs","index_type":"CAGRA"}' > index/meta.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cuvs: stub
        cupy: stub
        numpy: stub
    END_VERSIONS
    """
}
