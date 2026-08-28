process SEARCH_HNSWLIB {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_numpy_hnswlib_cxx-compiler:ab5f45645254951d' :
        'community.wave.seqera.io/library/python_numpy_hnswlib_cxx-compiler:ab5f45645254951d' }"

    input:
    tuple val(meta), path(windows), path(manifest)
    path database
    path hnsw_bundle, stageAs: 'hnsw_bundle'

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
    def separate_rerank = params.exact_rerank || params.hnswlib_rerank
    def search_k = separate_rerank ? params.candidate_k : params.seed_k
    def min_sim = (separate_rerank || params.quantize.toString().toLowerCase() in ['pq', 'opq']) ? '-inf' : String.valueOf(params.seed_min_similarity)
    """
    g++ -O3 -std=c++11 -fopenmp \\
        -I hnsw_bundle \\
        hnsw_bundle/pq_hnswlib.cpp \\
        -o pq_hnswlib

    mkdir -p query_windows
    cp ${windows} query_windows/
    cp ${manifest} query_windows/

    search_hnswlib.py \\
        --windows-dir query_windows \\
        --database ${database} \\
        --executable ./pq_hnswlib \\
        --output ${prefix}.seeds.tsv \\
        --k ${search_k} \\
        --min-similarity=${min_sim} \\
        ${ef_search} \\
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
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "query_id\\tquery_start\\tquery_end\\ttarget_id\\ttarget_start\\ttarget_end\\tscore\\trank" > ${prefix}.seeds.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hnswlib: stub
        numpy: stub
    END_VERSIONS
    """
}
