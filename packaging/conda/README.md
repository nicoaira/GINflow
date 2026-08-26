# Custom conda packages

GINflow publishes two packages to [`nicolas.aira`](https://anaconda.org/nicolas.aira) so
`-profile conda` and Wave-frozen `-profile docker` images need no extra local
installs.

| Package | Contents | Used by |
|---|---|---|
| `pq-cagra-adc` | GPU build + GPU/CPU ADC search | `BUILD_PQ_CAGRA_INDEX`, GPU `SEARCH_PQ_CAGRA` |
| `pq-cagra-adc-cpu` | CPU load + ADC walker only | CPU `SEARCH_PQ_CAGRA` (`--search_device cpu`) |

Rebuild (Python 3.12, Linux x86_64):

```bash
conda build packaging/conda/pq-cagra-adc-cpu
anaconda upload -u nicolas.aira --force "$(conda build packaging/conda/pq-cagra-adc-cpu --output)"

conda build packaging/conda/pq-cagra-adc
anaconda upload -u nicolas.aira --force "$(conda build packaging/conda/pq-cagra-adc --output)"
```

The GPU recipe uses CUDA 11.8 (`nvidia::cuda-nvcc=11.8.89` plus the matching
development headers and GCC 11.4 host compiler). After a new upload, freeze
Wave images from the module `environment*.yml` files:

```bash
wave --conda-file modules/build_pq_cagra_index/environment.yml \
    --conda-channels nicolas.aira,conda-forge --freeze --await 25m -o json
wave --conda-file modules/search_pq_cagra/environment.yml \
    --conda-channels nicolas.aira,conda-forge --freeze --await 25m -o json
wave --conda-file modules/search_pq_cagra/environment.yml \
    --conda-channels nicolas.aira,conda-forge --freeze --await 25m --singularity -o json
wave --conda-file modules/build_pq_cagra_index/environment.yml \
    --conda-channels nicolas.aira,conda-forge --freeze --await 25m --singularity -o json
```

Pin the resulting Docker and ORAS URLs in `modules/build_pq_cagra_index/main.nf`
and `modules/search_pq_cagra/main.nf`.

Current pins:

| Image | Tag |
|---|---|
| GPU Docker | `community.wave.seqera.io/library/python_numpy_cudatoolkit_pq-cagra-adc:ed61007abe908991` |
| GPU Singularity | `oras://community.wave.seqera.io/library/python_numpy_cudatoolkit_pq-cagra-adc:e0102895f642d35a` |
| CPU Docker | `community.wave.seqera.io/library/python_numpy_pq-cagra-adc-cpu:2ea6d576a29716e1` |
| CPU Singularity | `oras://community.wave.seqera.io/library/python_numpy_pq-cagra-adc-cpu:230079361cb63119` |
