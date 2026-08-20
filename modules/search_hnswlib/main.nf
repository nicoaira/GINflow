process SEARCH_HNSWLIB {
    tag "${meta.id}"
    label 'process_medium'

    conda "${ task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml" }"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        (task.accelerator ? 'oras://community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:40ef78c1cb7c29cd' : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e3/e39b53016c4ed7ede3a06391f841e0bf9d0711e35272ce1ab370ee149343a1fd/data') :
        (task.accelerator ? 'community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:93cd6db656f6b1e4' : 'community.wave.seqera.io/library/python_numpy_libstdcxx-ng_pip_scann:b1bc94cdc1825d91') }"

    input:
    tuple val(meta), path(quantized_windows)
    path database
    path hnsw_bundle, stageAs: 'hnsw_bundle'
    path query_embeddings, stageAs: 'query_embeddings/*'
    path query_windows, stageAs: 'query_windows/*'
    path query_manifests, stageAs: 'query_manifests/*'

    output:
    tuple val(meta), path("*.seeds.tsv"), emit: seeds
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ef_search = params.hnswlib_ef_search != null ? "--ef-search ${params.hnswlib_ef_search}" : ''
    def threads = params.hnswlib_num_threads != null ? "--num-threads ${params.hnswlib_num_threads}" : ''
    def candidate_k = params.hnswlib_candidate_k != null ? "--candidate-k ${params.hnswlib_candidate_k}" : ''
    def rerank = params.hnswlib_rerank ? "--rerank --query-embeddings query_embeddings/*.npz" : ''
    def raw_min_score = (params.seed_min_similarity as BigDecimal) * (params.window_size as BigDecimal)
    def compact_min_score = params.hnswlib_rerank ? params.seed_min_similarity : raw_min_score
    def gpu = task.accelerator as boolean
    def command = gpu ? """
    hnswlib_gpu.py search \\
        --windows-dir query_windows \\
        --manifests-dir query_manifests \\
        --database ${database} \\
        --query-embeddings query_embeddings/*.npz \\
        --output ${prefix}.seeds.tsv \\
        --k ${params.seed_k} \\
        --candidate-k ${params.hnswlib_gpu_candidate_k} \\
        --min-similarity ${params.seed_min_similarity} \\
        --itopk-size ${params.hnswlib_gpu_itopk_size} \\
        --search-batch-size ${params.hnswlib_gpu_search_batch_size} \\
        ${args}
    """ : """
    g++ -O3 -std=c++11 -fopenmp \\
        -I hnsw_bundle \\
        hnsw_bundle/compact_hnswlib.cpp \\
        -o compact_hnswlib

    if python3 -c "import json; assert json.load(open('${database}/meta.json')).get('index_type') == 'HNSWLIB_COMPACT'"; then
        search_compact_hnswlib.py \\
            --executable ./compact_hnswlib \\
            --windows-dir quantized_windows \\
            --database ${database} \\
            --output ${prefix}.seeds.tsv \\
            --k ${params.seed_k} \\
            ${candidate_k} \\
            ${rerank} \\
            --min-similarity ${compact_min_score} \\
            ${ef_search} \\
            ${threads} \\
            ${args}
    else
        search_hnswlib.py \\
            --windows-dir quantized_windows \\
            --database ${database} \\
            --output ${prefix}.seeds.tsv \\
            --k ${params.seed_k} \\
            --min-similarity ${raw_min_score} \\
            ${ef_search} \\
            ${threads} \\
            ${args}
    fi
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
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gpu = task.accelerator as boolean
    """
    echo -e "query_id\tquery_start\tquery_end\ttarget_id\ttarget_start\ttarget_end\tscore\trank" > ${prefix}.seeds.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hnswlib: ${gpu ? 'cagra-gpu-stub' : 'stub'}
        ${gpu ? 'cuvs: stub\n        cupy: stub' : ''}
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
