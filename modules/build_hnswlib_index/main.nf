process BUILD_HNSWLIB_INDEX {
    tag "hnswlib"
    label 'process_high_memory'

    conda "${ task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml" }"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        (task.accelerator ? 'oras://community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:40ef78c1cb7c29cd' : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e3/e39b53016c4ed7ede3a06391f841e0bf9d0711e35272ce1ab370ee149343a1fd/data') :
        (task.accelerator ? 'community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:93cd6db656f6b1e4' : 'community.wave.seqera.io/library/python_numpy_libstdcxx-ng_pip_scann:b1bc94cdc1825d91') }"

    input:
    path quantized_windows, stageAs: 'quantized_windows'
    path quantization, stageAs: 'quantization'
    path embeddings, stageAs: 'embeddings/*'
    path metadata, stageAs: 'metadata/*'
    path hnsw_bundle, stageAs: 'hnsw_bundle'
    path windows, stageAs: 'windows/*'
    path manifests, stageAs: 'manifests/*'

    output:
    path "faiss", emit: database
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def threads = params.hnswlib_num_threads != null ? "--num-threads ${params.hnswlib_num_threads}" : ''
    def gpu = task.accelerator as boolean
    def command = gpu ? """
    hnswlib_gpu.py build \\
        --windows-dir windows \\
        --manifests-dir manifests \\
        --quantization quantization \\
        --embeddings embeddings/*.npz \\
        --graph-metadata metadata/*.json \\
        --outdir faiss \\
        --intermediate-graph-degree ${params.hnswlib_gpu_intermediate_graph_degree} \\
        --graph-degree ${params.hnswlib_gpu_graph_degree} \\
        --build-algo ${params.hnswlib_gpu_build_algo} \\
        --itopk-size ${params.hnswlib_gpu_itopk_size} \\
        --search-batch-size ${params.hnswlib_gpu_search_batch_size} \\
        ${params.hnswlib_gpu_int8_scale != null ? "--int8-scale ${params.hnswlib_gpu_int8_scale}" : ''} \\
        ${args}
    """ : """
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
    """
    """
    ${command}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hnswlib: ${gpu ? 'cagra-gpu' : '0.8.0-cpp'}
        ${gpu ? "cuvs: \$(python3 -c \"from importlib.metadata import version; print(version('cuvs'))\")" : ''}
        ${gpu ? "cupy: \$(python3 -c \"import cupy; print(cupy.__version__)\")" : ''}
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        ${gpu ? '' : 'gxx: \$(g++ --version | head -1)'}
    END_VERSIONS
    """

    stub:
    def gpu = task.accelerator as boolean
    """
    mkdir -p faiss/quantization ${gpu ? 'faiss/cagra' : ''}
    ${gpu ? 'touch faiss/cagra/index.bin' : 'touch faiss/index.bin'}
    touch faiss/windows.tsv
    touch faiss/embeddings.npz
    touch faiss/records.tsv
    touch faiss/quantization/centroids.npy
    touch faiss/quantization/similarity.npy
    echo '{}' > faiss/quantization/quantization.json
    echo '${gpu ? '{"backend":"hnswlib","index_type":"HNSWLIB_GPU_CAGRA","candidate_representation":"original_normalized_window_int8"}' : '{"backend":"hnswlib","index_type":"HNSWLIB_COMPACT","candidate_representation":"uint16_centroid_codes"}'}' > faiss/meta.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hnswlib: ${gpu ? 'cagra-gpu-stub' : 'stub'}
        ${gpu ? 'cuvs: stub\n        cupy: stub' : ''}
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
