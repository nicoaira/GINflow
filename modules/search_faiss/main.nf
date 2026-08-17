def scann_requested() {
    return params.use_scann || ((params.index as String)?.trim()?.equalsIgnoreCase('scann'))
}

def ngt_requested() {
    return params.use_ngt || ((params.index as String)?.trim()?.equalsIgnoreCase('ngt'))
}

def cuvs_requested() {
    return params.use_cuvs || ((params.index as String)?.trim()?.equalsIgnoreCase('cuvs'))
}

def index_conda(moduleDir, task) {
    if (cuvs_requested()) {
        return "${moduleDir}/environment.cuvs.yml"
    }
    if (ngt_requested()) {
        return "${moduleDir}/environment.ngt.yml"
    }
    if (scann_requested()) {
        return "${moduleDir}/environment.scann.yml"
    }
    return task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml"
}

def index_container(task) {
    def sif = workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
    if (cuvs_requested()) {
        return sif
            ? 'oras://community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:40ef78c1cb7c29cd'
            : 'community.wave.seqera.io/library/python_numpy_cupy_cudatoolkit_pruned:93cd6db656f6b1e4'
    }
    if (ngt_requested()) {
        return sif
            ? 'oras://community.wave.seqera.io/library/python_numpy_ngt:ff987f0bb9f59553'
            : 'community.wave.seqera.io/library/python_numpy_ngt:9a0ca7a46e9c18b2'
    }
    if (scann_requested()) {
        return sif
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e3/e39b53016c4ed7ede3a06391f841e0bf9d0711e35272ce1ab370ee149343a1fd/data'
            : 'community.wave.seqera.io/library/python_numpy_libstdcxx-ng_pip_scann:b1bc94cdc1825d91'
    }
    if (task.accelerator) {
        return sif
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/40/405599a0a76073fa7c942f2773cad66490aabb99fb9825c79ba50ad1157e9a6a/data'
            : 'community.wave.seqera.io/library/python_numpy_faiss-gpu:7382ed4d7c6e6258'
    }
    return sif
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7920d351ee6307611f471013dcddff6d8a5c7ec7bf723d854e99a545376b63d8/data'
        : 'community.wave.seqera.io/library/python_numpy_faiss-cpu_mkl_libblas:078dd4f35c795d6e'
}

process SEARCH_FAISS {
    tag "${meta.id}"
    label 'process_low'

    conda { index_conda(moduleDir, task) }
    container { index_container(task) }

    input:
    tuple val(meta), path(windows), path(manifest)
    path database

    output:
    tuple val(meta), path("*.seeds.tsv"), emit: seeds
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def versions = cuvs_requested()
        ? 'cuvs: $(python3 -c "from importlib.metadata import version; print(version(\'cuvs\'))")'
        : ngt_requested()
        ? 'ngt: $(python3 -c "from importlib.metadata import version; print(version(\'ngt\'))")'
        : scann_requested()
        ? 'scann: $(python3 -c "from importlib.metadata import version; print(version(\'scann\'))")'
        : 'faiss: $(python3 -c "import faiss; print(faiss.__version__)")'
    if (cuvs_requested()) {
        def gpu_flag = task.accelerator ? '--gpu' : ''
        def n_probes = params.cuvs_n_probes != null ? "--cuvs-n-probes ${params.cuvs_n_probes}" : ''
        """
        search_faiss.py \\
            --windows ${windows} \\
            --manifest ${manifest} \\
            --database ${database} \\
            --output ${prefix}.seeds.tsv \\
            --k ${params.seed_k} \\
            --min-similarity ${params.seed_min_similarity} \\
            ${n_probes} \\
            ${gpu_flag} \\
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ${versions}
            cupy: \$(python3 -c "import cupy; print(cupy.__version__)")
            numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        END_VERSIONS
        """
    }
    else if (ngt_requested()) {
        """
        search_faiss.py \\
            --windows ${windows} \\
            --manifest ${manifest} \\
            --database ${database} \\
            --output ${prefix}.seeds.tsv \\
            --k ${params.seed_k} \\
            --min-similarity ${params.seed_min_similarity} \\
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ${versions}
            numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        END_VERSIONS
        """
    }
    else if (scann_requested()) {
        def scann_lts = params.scann_leaves_to_search != null ? "--scann-leaves-to-search ${params.scann_leaves_to_search}" : ''
        """
        search_faiss.py \\
            --windows ${windows} \\
            --manifest ${manifest} \\
            --database ${database} \\
            --output ${prefix}.seeds.tsv \\
            --k ${params.seed_k} \\
            --min-similarity ${params.seed_min_similarity} \\
            --scann-reorder ${params.scann_reorder} \\
            ${scann_lts} \\
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ${versions}
            numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        END_VERSIONS
        """
    }
    else {
        def gpu_flag  = task.accelerator ? '--gpu' : ''
        def nprobe    = params.faiss_nprobe != null ? "--nprobe ${params.faiss_nprobe}" : ''
        def ef_search = params.faiss_hnsw_ef_search != null ? "--hnsw-ef-search ${params.faiss_hnsw_ef_search}" : ''
        """
        search_faiss.py \\
            --windows ${windows} \\
            --manifest ${manifest} \\
            --database ${database} \\
            --output ${prefix}.seeds.tsv \\
            --k ${params.seed_k} \\
            --min-similarity ${params.seed_min_similarity} \\
            ${nprobe} \\
            ${ef_search} \\
            ${gpu_flag} \\
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ${versions}
            numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        END_VERSIONS
        """
    }

    stub:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def versions = cuvs_requested()
        ? 'cuvs: $(python3 -c "from importlib.metadata import version; print(version(\'cuvs\'))")'
        : ngt_requested()
        ? 'ngt: $(python3 -c "from importlib.metadata import version; print(version(\'ngt\'))")'
        : scann_requested()
        ? 'scann: $(python3 -c "from importlib.metadata import version; print(version(\'scann\'))")'
        : 'faiss: $(python3 -c "import faiss; print(faiss.__version__)")'
    """
    echo -e "query_id\\tquery_start\\tquery_end\\ttarget_id\\ttarget_start\\ttarget_end\\tscore\\trank" > ${prefix}.seeds.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${versions}
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
