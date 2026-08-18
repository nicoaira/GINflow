#!/usr/bin/env nextflow

/*
 * Benchmark-only entry workflow. It intentionally stops after window creation:
 * index construction is performed independently by each benchmark backend.
 */

include { PREPARE_WINDOWS } from '../workflows/prepare_windows'

params.benchmark_prefix = null

workflow {
    main:
    if (!params.input) {
        error "--input is required for the benchmark window cache."
    }
    if (!params.benchmark_prefix) {
        error "--benchmark_prefix is required for stable cache artifact names."
    }

    def benchmark_shard_size = params.shard_size as int
    PREPARE_WINDOWS(params.input, params.benchmark_prefix, benchmark_shard_size)

    publish:
    graphs     = PREPARE_WINDOWS.out.graphs
    embeddings = PREPARE_WINDOWS.out.embeddings
    windows    = PREPARE_WINDOWS.out.windows
}

output {
    graphs {
        path 'graphs'
    }
    embeddings {
        path 'embeddings'
    }
    windows {
        path 'windows'
    }
}
