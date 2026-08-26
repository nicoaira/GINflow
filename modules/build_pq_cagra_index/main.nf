process BUILD_PQ_CAGRA_INDEX {
    tag "pq-cagra"
    label 'process_high_memory'
    label 'process_gpu'
    accelerator 1

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_numpy_cudatoolkit_pq-cagra-adc:e0102895f642d35a' :
        'community.wave.seqera.io/library/python_numpy_cudatoolkit_pq-cagra-adc:ed61007abe908991' }"

    input:
    path quantized_windows, stageAs: 'quantized_windows'
    path quantization, stageAs: 'quantization'
    path embeddings, stageAs: 'embeddings/*'
    path metadata, stageAs: 'metadata/*'
    val n_windows
    val n_nodes

    output:
    path "faiss", emit: database
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    set -o pipefail
    mkdir -p faiss
    build_pq_cagra.py \\
        --windows-dir quantized_windows \\
        --quantization quantization \\
        --outdir faiss \\
        --graph-degree ${params.cagra_graph_degree} \\
        --intermediate-graph-degree ${params.cagra_intermediate_graph_degree} \\
        --nndescent-iterations ${params.cagra_nndescent_iterations} \\
        --itopk-size ${params.cagra_itopk_size} \\
        --embeddings embeddings/*.npz \\
        --graph-metadata metadata/*.json \\
        ${args} 2>&1 | tee progress.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p faiss/quantization
    touch faiss/cagra.index
    echo 'PQ_CAGRA_PROGRESS scope=task phase=complete percent=100'
    echo '{"backend":"cagra","index_type":"PQ_CAGRA"}' > faiss/meta.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        numpy: stub
    END_VERSIONS
    """
}
