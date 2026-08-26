process ALIGN_CLUSTERS {
    tag "${clusters.baseName}"
    cpus { Math.max(1, params.align_cpus.toString().toInteger()) }
    time { 8.h * task.attempt }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c7db9e2e04ab9a21b4e3bf856e5c3a9e3eb1909cf3e251b366dca6e10ae8684a/data' :
        'community.wave.seqera.io/library/python_ginfinity-sw:151e772020737622' }"

    input:
    path clusters
    path query_embeddings, stageAs: 'query_emb/*'
    path query_metadata, stageAs: 'query_meta/*'
    path database
    path parameters
    path evd

    output:
    path "*.alignments.tsv",          emit: alignments
    path "*.alignments.txt",          emit: text
    path "*.alignment_stats.json",    emit: stats
    path "versions.yml",              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: clusters.baseName
    """
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export NUMBA_NUM_THREADS=${task.cpus}
    align_clusters.py \\
        --clusters ${clusters} \\
        --parameters ${parameters} \\
        --query-embeddings query_emb/*.npz \\
        --query-metadata query_meta/*.json \\
        --database ${database} \\
        --evd ${evd} \\
        --pad ${params.align_pad} \\
        --max-cells ${params.align_max_cells} \\
        --cpus ${task.cpus} \\
        --output ${prefix}.alignments.tsv \\
        --alignment-text ${prefix}.alignments.txt \\
        --stats-json ${prefix}.alignment_stats.json \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ginfinity-sw: \$(python3 -c "import ginfinity_sw; print(ginfinity_sw.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: clusters.baseName
    """
    echo -e "cluster_id\\tquery_id\\ttarget_id\\tscore\\ttotal_score\\tmax_score\\tbit_score\\tevalue\\tevalue_pair\\tlog_evalue\\talignment_count\\tquery_start\\tquery_end\\ttarget_start\\ttarget_end\\tquery_length\\ttarget_length\\tmatch_count\\taligned_columns\\tungapped_columns\\tbase_matches\\tstructure_matches\\tseed_count\\tmax_seed_score\\tquery_sequence\\tquery_structure\\ttarget_sequence\\ttarget_structure\\tquery_aligned\\ttarget_aligned" > ${prefix}.alignments.tsv
    touch ${prefix}.alignments.txt
    echo '{}' > ${prefix}.alignment_stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ginfinity-sw: \$(python3 -c "import ginfinity_sw; print(ginfinity_sw.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
