process BUILD_HNSWLIB_INDEX {
    tag "hnswlib"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_numpy_hnswlib_cxx-compiler:ab5f45645254951d' :
        'community.wave.seqera.io/library/python_numpy_hnswlib_cxx-compiler:ab5f45645254951d' }"

    input:
    path quantized_windows, stageAs: 'quantized_windows'
    path quantization, stageAs: 'quantization'
    path embeddings, stageAs: 'embeddings/*'
    path metadata, stageAs: 'metadata/*'
    val n_windows
    path hnsw_bundle, stageAs: 'hnsw_bundle'

    output:
    path "faiss", emit: database
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def threads = params.hnswlib_num_threads != null ? "--num-threads ${params.hnswlib_num_threads}" : ''
    """
    build_hnswlib.py \\
        --windows-dir quantized_windows \\
        --quantization quantization \\
        --hnsw-bundle hnsw_bundle \\
        --outdir faiss \\
        --m ${params.hnswlib_m} \\
        --ef-construction ${params.hnswlib_ef_construction} \\
        --ef-search ${params.hnswlib_ef_search} \\
        --random-seed ${params.hnswlib_random_seed} \\
        --embeddings embeddings/*.npz \\
        --graph-metadata metadata/*.json \\
        ${threads} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hnswlib: 0.8.0-cpp
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        gxx: \$(g++ --version | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p faiss/quantization
    touch faiss/index.bin
    echo '{"backend":"hnswlib","index_type":"HNSWLIB_PQ"}' > faiss/meta.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hnswlib: stub
        numpy: stub
    END_VERSIONS
    """
}
