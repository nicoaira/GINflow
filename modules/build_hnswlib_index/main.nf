process BUILD_HNSWLIB_INDEX {
    tag "hnswlib"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e3/e39b53016c4ed7ede3a06391f841e0bf9d0711e35272ce1ab370ee149343a1fd/data' :
        'community.wave.seqera.io/library/python_numpy_libstdcxx-ng_pip_scann:b1bc94cdc1825d91' }"

    input:
    path quantized_windows, stageAs: 'quantized_windows'
    path quantization, stageAs: 'quantization'
    path embeddings, stageAs: 'embeddings/*'
    path metadata, stageAs: 'metadata/*'
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
    g++ -O3 -std=c++11 -fopenmp \\
        -I hnsw_bundle \\
        hnsw_bundle/compact_hnswlib.cpp \\
        -o compact_hnswlib

    build_compact_hnswlib.py \\
        --windows-dir quantized_windows \\
        --quantization quantization \\
        --executable ./compact_hnswlib \\
        --embeddings embeddings/*.npz \\
        --graph-metadata metadata/*.json \\
        --outdir faiss \\
        --m ${params.hnswlib_m} \\
        --ef-construction ${params.hnswlib_ef_construction} \\
        --ef-search ${params.hnswlib_ef_search} \\
        --random-seed ${params.hnswlib_random_seed} \\
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
    touch faiss/windows.tsv
    touch faiss/embeddings.npz
    touch faiss/records.tsv
    touch faiss/quantization/centroids.npy
    touch faiss/quantization/similarity.npy
    echo '{}' > faiss/quantization/quantization.json
    echo '{"backend":"hnswlib","index_type":"HNSWLIB_COMPACT","candidate_representation":"uint16_centroid_codes"}' > faiss/meta.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hnswlib: stub
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
