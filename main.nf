#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nicoaira/ginflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nicoaira/ginflow
----------------------------------------------------------------------------------------
*/

include { GINFLOW } from './workflows/ginflow'

workflow {
    main:
    result = GINFLOW(params.input)

    samples_ch = result.graphs
        .join(result.embeddings)
        .map { meta, tensors, sidecar, npz, manifest ->
            [
                id         : meta.id,
                graphs     : tensors,
                metadata   : sidecar,
                embeddings : npz,
                manifest   : manifest,
            ]
        }

    publish:
    samples = samples_ch
}

output {
    samples {
        path { sample ->
            sample.graphs >> "graphs/${sample.id}/"
            sample.metadata >> "graphs/${sample.id}/"
            sample.embeddings >> "embeddings/${sample.id}/"
            sample.manifest >> "embeddings/${sample.id}/"
        }
        index {
            path 'samples.csv'
            header true
        }
    }
}
