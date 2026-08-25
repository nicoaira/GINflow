# Clustering and alignment

Window search only proposes **seeds**: short 11-nt hits. Clustering
turns those into HSP-like boxes. GINFINITY-SW then produces gapped
local alignments on a padded crop of each box. Pair-level collapse
makes one BLAST-style result per query–target pair.

```text
seeds.tsv
    → CLUSTER_SEEDS          (diagonal boxes)
    → ALIGN_CLUSTERS         (one HSP row per cluster; all queries, N CPU threads)
    → MERGE_ALIGNMENTS       (one row per query–target pair)
    → SPLIT_CLUSTERS         (only when plotting; one file per query)
```

## Seeds

`seeds.tsv` columns: `query_id`, `query_start`, `query_end`,
`target_id`, `target_start`, `target_end`, `score`, `rank`.

Coordinates are **0-based half-open** on the independent subject (the
full molecule, or the core of a slice). Self-hits are kept.

`--seed_min_similarity` (default 0.8) is a cosine threshold on the
**original** window after rerank. Without exact rerank, PQ/OPQ search
emits the top `--seed_k` ADC neighbours and does not apply this cutoff
(ADC is not cosine).

`--cluster_max_seed_rank` (default 10) drops worse-ranked seeds before
clustering. `0` keeps every seed.

## Diagonal clustering

`CLUSTER_SEEDS` (`bin/cluster_seeds.py`) groups seeds independently for
each `(query_id, target_id)`.

A cluster starts at the highest-scoring unused seed. It then repeatedly
adds the best remaining seed that satisfies all of:

| Constraint | Flag | Default | Meaning |
|---|---|---:|---|
| Box proximity | `--cluster_span` | 80 | Query and target intervals must overlap the current box expanded by this many nucleotides |
| Diagonal drift | `--cluster_diagonal_tolerance` | 12 | Seed diagonal `query_start - target_start` must fall within this many nt of the current `[diag_min, diag_max]` |
| Diagonal breadth | `--cluster_max_diagonal_span` | 96 | After adding the seed, `diag_max - diag_min` must not exceed this. `0` disables the cap |
| Minimum size | `--cluster_min_seeds` | 2 | Singletons are dropped |

The diagonal of an 11-nt seed is `query_start - target_start`. Seeds on
the same ungapped offset share a diagonal; a small tolerance allows
indels and embedding jitter without merging unrelated hits.

Published:

- `clusters.tsv` — one row per cluster (`cluster_id`, span, seed count,
  max score, diagonal stats)
- `cluster_members.tsv` — every seed that went into a cluster

These boxes are the SW crops, not the final alignments.

## Crop and Smith–Waterman

`ALIGN_CLUSTERS` (`bin/align_clusters.py`) loads original residue
embeddings for the query (from the query shards) and the target (from
the packed `faiss/` residue store, or a legacy per-id `embeddings.npz`).
Only identifiers present in the cluster table are loaded. It cuts a crop
around the cluster:

```text
query  [cluster.query_start - pad, cluster.query_end + pad)
target [cluster.target_start - pad, cluster.target_end + pad)
```

`--align_pad` defaults to 32. `--align_max_cells` (default 16 777 216)
is the maximum product of the two crop lengths; larger crops are
rejected so a runaway cluster cannot explode DP memory.

Scoring is **embedding-only cosine**, not a nucleotide matrix. Gap
open, gap extend, and the cosine → substitution map live in
`assets/alignment.json` and are passed to
[GINFINITY-SW](https://github.com/nicoaira/GINFINITY-SW). Independent
crops are aligned concurrently (`--align_cpus`).

Each cluster may yield several disjoint local HSPs inside the crop
(`align_multiple`). Filters for the SW call (and for EVD calibration):

| Flag | Default | Role |
|---|---:|---|
| `--align_max_alignments` | 16 | Max disjoint local HSPs kept from one crop |
| `--align_min_score` | 0.0 | Minimum HSP score |
| `--align_min_match_count` | 1 | Minimum matched columns |

`ALIGN_CLUSTERS` deliberately writes **one TSV row per seed cluster**.
It is not the pair-level result.

For sliced subjects, the sequence/structure written on the row is the
**core window**, with pairs that cross the cut rewritten as unpaired
(`.`) so the subject stays a legal, balanced structure.

## Pair collapse

`MERGE_ALIGNMENTS` (`bin/merge_alignments.py`) groups rows by
`(query_id, target_id)`:

- `total_score` — sum of disjoint HSP scores
- `max_score` — strongest HSP
- `score` — legacy alias of `total_score`
- `alignment_count` — number of HSPs
- `cluster_ids`, `hsp_scores`, `hsp_spans` — provenance
- `query_aligned` / `target_aligned` — gapped strings of the **first**
  HSP
- `evalue` — database E-value from `total_score`
- `evalue_pair` — same formula with target length instead of database
  residues
- `bit_score` — $(\lambda S_{\mathrm{total}} - \ln K) / \ln 2$

Rows are ranked by **ascending database E-value**.

`alignments.txt` is a human-readable rendering: one pair header, then
the six-line GINFINITY-SW block for each HSP.

How E-values are fit: [E-values](statistics.md).
