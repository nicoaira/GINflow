process SEARCH_HNSWLIB {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e3/e39b53016c4ed7ede3a06391f841e0bf9d0711e35272ce1ab370ee149343a1fd/data' :
        'community.wave.seqera.io/library/python_numpy_libstdcxx-ng_pip_scann:b1bc94cdc1825d91' }"

    input:
    tuple val(meta), path(quantized_windows)
    path database
    path hnsw_bundle, stageAs: 'hnsw_bundle'
    path query_embeddings, stageAs: 'query_embeddings/*'

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
    """
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hnswlib: 0.8.0-cpp
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        gxx: \$(g++ --version | head -1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "query_id\\tquery_start\\tquery_end\\ttarget_id\\ttarget_start\\ttarget_end\\tscore\\trank" > ${prefix}.seeds.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hnswlib: stub
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
