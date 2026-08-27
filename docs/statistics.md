# E-values

GINflow reports BLAST-style database E-values for each query–target
pair, using the **sum of disjoint local HSP scores** as the raw score.

$$
E = K\, m\, N\, \exp(-\lambda S_{\mathrm{total}})
$$

| Symbol | Meaning |
|---|---|
| $S_{\mathrm{total}}$ | Sum of disjoint GINFINITY-SW HSP scores for the pair |
| $m$ | Query length (nucleotides) |
| $N$ | Searchable database residues (`evd.json` → `database_residues`) |
| $\lambda$, $K$ | Karlin–Altschul parameters fit at database-build time |

`evalue_pair` uses the target length instead of $N$. Bit score:

$$
S_{\mathrm{bit}} = \frac{\lambda S_{\mathrm{total}} - \ln K}{\ln 2}
$$

The legacy conversion $K = e^{-\lambda\mu}$ is **not** used. $K$
comes from a length-aware Gumbel MLE so that
$\mu = \ln(Kmn)/\lambda$.

## Why a custom null

GINFINITY-SW scores cosine similarity of residue embeddings, not BLOSUM
or a nucleotide matrix. Published Karlin–Altschul parameters do not
apply. GINflow estimates $(\lambda, K)$ from a **reverse-sequence
null** that keeps local embedding correlation (the reversed query is
still a real RNA embedding) but destroys homology.

## How `ESTIMATE_EVD` samples

`bin/estimate_evd.py`, process `ESTIMATE_EVD`:

1. Load `index/embeddings.npz` (original residue vectors).
2. Draw `--evd_samples` (default 1000) random pairs of distinct
   database molecules.
3. Reverse the query embedding along the sequence axis; leave the
   target as-is.
4. Crop both to at most `--evd_max_length` (default 400) so each DP
   stays $O(L^2)$. `--align_max_cells` still caps the product.
5. Run the same multi-HSP SW used in search (`align_multiple` with
   `--align_max_alignments`, `--align_min_score`,
   `--align_min_match_count`).
6. Record `total_score` and the two crop lengths.

If too few samples succeed (or too few scores are positive), the
process errors: the scoring system is then too harsh for a Gumbel tail,
or the database is too small.

`--evd_seed` (default 1) makes the draw reproducible.
`-profile smoke_test` lowers `--evd_samples` to 200.

## Fitting λ and K

Positive null scores are treated as observations from

$$
S - \frac{\ln(mn)}{\lambda} \sim \mathrm{Gumbel}\left(\frac{\ln K}{\lambda},\; \frac{1}{\lambda}\right)
$$

The fitter (`fit_karlin_altschul`):

1. Moment estimate of a two-parameter Gumbel on the raw scores.
2. A few iterations that peel off $\ln(mn)/\lambda$ and re-fit.
3. Nelder–Mead on the Gumbel negative log-likelihood in
   $(\ln\lambda, \ln K)$.
4. A Kolmogorov–Smirnov statistic against the standard Gumbel CDF
   $F(z)=\exp(-e^{-z})$, stored as `ks_statistic`.

`index/evd.json` records `lambda`, `K`, `database_residues`, the
Gumbel loc/scale, sample counts, and that KS value.

On a query-only run, GINflow reuses `index/evd.json` from `--database`
when present. Otherwise it fits EVD from the bundled embeddings before
alignment.

## How to read the numbers

- **Database E (`evalue`)** — expected number of pairs in a database
  of size $N$ with score ≥ this pair. This is the rank key.
- **Pair E (`evalue_pair`)** — same score against a database that is
  only this one target. Useful when you care about “is this pair
  better than chance on these two lengths?” rather than “would this
  show up in a large search?”
- **Bits** — log-scaled score; higher is better, independent of
  database size.
- **`total_score` vs `max_score`** — several weak HSPs can outrank one
  strong HSP on database E because $S_{\mathrm{total}}$ is the sum.
  Look at `alignment_count` and `hsp_spans` before treating a pair as
  a single motif.

There is no pipeline-wide E-value cutoff. Filter in
`report.html` or downstream from `alignments.tsv`.
