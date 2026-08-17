def scann_requested() {
    return params.use_scann || ((params.index as String)?.trim()?.equalsIgnoreCase('scann'))
}

def ngt_requested() {
    return params.use_ngt || ((params.index as String)?.trim()?.equalsIgnoreCase('ngt'))
}

def index_conda(moduleDir, task) {
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

process BUILD_FAISS_INDEX {
    tag "faiss"
    label 'process_medium'

    conda { index_conda(moduleDir, task) }
    container { index_container(task) }

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
    def versions = ngt_requested()
        ? 'ngt: $(python3 -c "from importlib.metadata import version; print(version(\'ngt\'))")'
        : scann_requested()
        ? 'scann: $(python3 -c "from importlib.metadata import version; print(version(\'scann\'))")'
        : 'faiss: $(python3 -c "import faiss; print(faiss.__version__)")'
    if (ngt_requested()) {
        """
        build_faiss.py \\
            --windows windows/*.windows.npz \\
            --manifests manifests/*.windows.manifest.json \\
            --embeddings embeddings/*.npz \\
            --graph-metadata metadata/*.json \\
            --outdir faiss \\
            --backend ngt \\
            --ngt-index-type ${params.ngt_index} \\
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ${versions}
            numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        END_VERSIONS
        """
    }
    else if (scann_requested()) {
        def scann_soar   = params.scann_soar ? '--scann-soar' : ''
        def scann_leaves = params.scann_leaves != null ? "--scann-leaves ${params.scann_leaves}" : ''
        def scann_lts    = params.scann_leaves_to_search != null ? "--scann-leaves-to-search ${params.scann_leaves_to_search}" : ''
        """
        build_faiss.py \\
            --windows windows/*.windows.npz \\
            --manifests manifests/*.windows.manifest.json \\
            --embeddings embeddings/*.npz \\
            --graph-metadata metadata/*.json \\
            --outdir faiss \\
            --index-type ScaNN \\
            --scann-reorder ${params.scann_reorder} \\
            --scann-ah-dim ${params.scann_ah_dim} \\
            --scann-anisotropic ${params.scann_anisotropic} \\
            --num-neighbors ${params.seed_k} \\
            ${scann_soar} \\
            ${scann_leaves} \\
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
        def nlist     = params.faiss_nlist      != null ? "--nlist ${params.faiss_nlist}" : ''
        def nprobe    = params.faiss_nprobe     != null ? "--nprobe ${params.faiss_nprobe}" : ''
        def lsh_nbits = params.faiss_lsh_nbits  != null ? "--lsh-nbits ${params.faiss_lsh_nbits}" : ''
        """
        build_faiss.py \\
            --windows windows/*.windows.npz \\
            --manifests manifests/*.windows.manifest.json \\
            --embeddings embeddings/*.npz \\
            --graph-metadata metadata/*.json \\
            --outdir faiss \\
            --index-type ${params.faiss_index} \\
            --pq-m ${params.faiss_pq_m} \\
            --pq-nbits ${params.faiss_pq_nbits} \\
            --pq-m-refine ${params.faiss_pq_m_refine} \\
            --hnsw-m ${params.faiss_hnsw_m} \\
            --hnsw-ef-construction ${params.faiss_hnsw_ef_construction} \\
            --hnsw-ef-search ${params.faiss_hnsw_ef_search} \\
            --sq-type ${params.faiss_sq_type} \\
            ${nlist} \\
            ${nprobe} \\
            ${lsh_nbits} \\
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
    def versions = ngt_requested()
        ? 'ngt: $(python3 -c "from importlib.metadata import version; print(version(\'ngt\'))")'
        : scann_requested()
        ? 'scann: $(python3 -c "from importlib.metadata import version; print(version(\'scann\'))")'
        : 'faiss: $(python3 -c "import faiss; print(faiss.__version__)")'
    def stub_index_dir = ngt_requested() ? 'ngt' : 'scann'
    """
    mkdir -p faiss/${stub_index_dir}
    touch faiss/index.faiss
    touch faiss/windows.tsv
    touch faiss/embeddings.npz
    touch faiss/records.tsv
    echo '{}' > faiss/meta.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${versions}
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
