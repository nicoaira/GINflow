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

The GPU recipe needs `nvcc` 11.8 (`nvidia::cuda-nvcc=11.8.89`). After a new
upload, freeze Wave images from the module `environment*.yml` files:

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
| GPU Docker | `community.wave.seqera.io/library/python_numpy_cudatoolkit_pq-cagra-adc:51e2d56cbda45e24` |
| GPU Singularity | `oras://community.wave.seqera.io/library/python_numpy_cudatoolkit_pq-cagra-adc:d24945e894e22c05` |
| CPU Docker | `community.wave.seqera.io/library/python_numpy_pq-cagra-adc-cpu:2ea6d576a29716e1` |
| CPU Singularity | `oras://community.wave.seqera.io/library/python_numpy_pq-cagra-adc-cpu:230079361cb63119` |
