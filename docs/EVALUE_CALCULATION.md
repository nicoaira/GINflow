# E-value Calculation in GINflow

## Overview

GINflow implements a **BLAST-like E-value calculation** to assess the statistical significance of sequence alignments based on protein language model embeddings. The E-value (Expectation value) represents the expected number of alignments with a score at least as good as the observed score that would occur by chance in a database search.

**Key Concept**: A lower E-value indicates a more statistically significant alignment. For example:
- E-value = 0.001: Expect 1 such alignment per 1,000 random searches (highly significant)
- E-value = 1.0: Expect 1 such alignment per search (marginally significant)
- E-value = 10.0: Expect 10 such alignments per search (likely noise)

---

## Two-Phase Calculation Strategy

### Phase 1: EVD Parameter Estimation (Pre-computation)
**Module**: `calculate_evd_params` ([bin/calculate_evd_params.py](../bin/calculate_evd_params.py))

Estimates database-wide Extreme Value Distribution parameters (λ and K) that characterize the null distribution of alignment scores. These parameters are computed **once per database** and reused across all queries.

### Phase 2: E-value Computation (Per-alignment)
**Module**: `align_candidates` ([bin/align_candidates.py](../bin/align_candidates.py))

Uses the pre-computed parameters to calculate E-values for each alignment during the alignment phase.

---

## Phase 1: EVD Parameter Estimation

### 1.1 Background Statistics Sampling

**Purpose**: Establish a baseline similarity distribution for unrelated sequence pairs.

**Method** (function: `sample_background`, [calculate_evd_params.py:451-538](../bin/calculate_evd_params.py#L451-L538)):

1. Randomly sample `--background-samples` (default: 10,000) pairs of nodes from the embedding database
2. For each pair:
   - Select random node `i` from random sequence A
   - Select random node `j` from random sequence B
   - Compute dot product: $\text{sim} = \text{embedding}_A[i] \cdot \text{embedding}_B[j]$
3. Estimate population parameters from the sample:
   $$\mu_0 = \text{mean}(\text{all\_similarities})$$
   $$\sigma_0 = \text{std}(\text{all\_similarities})$$

**Output**: Background mean $\mu_0$ and standard deviation $\sigma_0$

**Implementation Notes**:
- Uses multiprocessing for parallel sampling across workers
- Supports shared memory mode to reduce memory overhead when using multiple workers
- Each worker uses a deterministic random seed (`base_seed + worker_id`) for reproducibility

---

### 1.2 Z-score Transformation

**Purpose**: Convert raw dot products to normalized scores comparable across different sequence pairs.

**Formula** (function: `compute_score`, [calculate_evd_params.py:541-546](../bin/calculate_evd_params.py#L541-L546)):

$$z\text{-score} = \gamma \times \frac{\text{dot} - \mu_0}{\sigma_0}$$

$$\text{score} = \text{clamp}(z\text{-score}, \text{score}_{\min}, \text{score}_{\max})$$

Where:
- $\gamma$ (gamma) = scaling factor (default: 1.5)
- $\text{dot}$ = raw dot product between embedding vectors
- $\mu_0, \sigma_0$ = background statistics from step 1.1
- $\text{score}_{\min}, \text{score}_{\max}$ = clamping bounds (default: -4, 8)

**Rationale**:
- Normalization accounts for embedding space properties (dimensions, typical magnitudes)
- Gamma scaling emphasizes differences (higher $\gamma$ → larger score differences)
- Clamping prevents extreme outliers from dominating the alignment

---

### 1.3 Null Distribution Generation

**Purpose**: Generate a distribution of scores from random (unrelated) alignments.

**Method** (function: `estimate_evd_parameters`, [calculate_evd_params.py:852-1041](../bin/calculate_evd_params.py#L852-L1041)):

1. Sample `--evd-samples` (default: 1,000) random sequence pairs
2. For each pair:
   - **Shuffle** the embedding vectors to destroy sequential information while preserving composition
     - Example: $[v_1, v_2, v_3, v_4] \rightarrow [v_3, v_1, v_4, v_2]$
     - Function: `shuffle_embeddings` ([calculate_evd_params.py:549-563](../bin/calculate_evd_params.py#L549-L563))
   - Perform **banded Smith-Waterman alignment** with the same parameters used for real alignments:
     - Affine gap penalties: gap_open (default: 12), gap_extend (default: 2)
     - Band width (default: 96)
     - X-drop termination (default: 50)
     - Scoring uses z-scores from step 1.2
   - Record the alignment score (if > 0)

3. Result: Collection of null scores $S_{\text{null}} = \{s_1, s_2, \ldots, s_n\}$

**Why Shuffling?**
- Preserves the embedding composition (same vectors, same overall statistics)
- Destroys biological sequential relationships
- Creates realistic "random" sequences that match database characteristics
- Avoids biases from truly random synthetic sequences

**Smith-Waterman Implementation**:
- Function: `smith_waterman` ([calculate_evd_params.py:566-730](../bin/calculate_evd_params.py#L566-L730))
- Banded algorithm for efficiency (only computes DP cells within a diagonal band)
- Three-state model for affine gaps (M=match, X=gap in query, Y=gap in target)
- Local alignment (scores can start/end anywhere)
- X-drop early termination when score drops too far from the maximum

---

### 1.4 Gumbel Distribution Fitting

**Purpose**: Model the null score distribution with an Extreme Value Distribution (EVD).

**Theory**: The maximum score from independent local alignments follows a Gumbel (Type I EVD) distribution. This is the theoretical foundation of BLAST statistics (Karlin-Altschul statistics).

**Fitting Process** ([calculate_evd_params.py:1025-1031](../bin/calculate_evd_params.py#L1025-L1031)):

1. Fit Gumbel distribution to null scores using Maximum Likelihood Estimation (MLE):
   ```python
   from scipy.stats import gumbel_r
   μ_gumbel, β_gumbel = gumbel_r.fit(S_null)
   ```

   The Gumbel PDF is:
   $$f(x) = \frac{1}{\beta} \exp\left(-\frac{x - \mu}{\beta}\right) \exp\left(-\exp\left(-\frac{x - \mu}{\beta}\right)\right)$$

2. Convert Gumbel parameters to **Karlin-Altschul parameters**:
   $$\lambda = \frac{1}{\beta_{\text{gumbel}}}$$
   $$K = \exp(-\lambda \times \mu_{\text{gumbel}})$$

**Parameter Interpretation**:
- $\lambda$ **(lambda)**: Score decay rate; higher $\lambda$ means scores decrease faster with randomness
- $K$: Search space correction factor; accounts for finite sequence lengths and edge effects

**Output** ([calculate_evd_params.py:1102-1126](../bin/calculate_evd_params.py#L1102-L1126)):

Saved to `evd_params.json`:
```json
{
  "evd_lambda": 0.005560,
  "evd_K": 0.5144,
  "background_mu": 0.0580,
  "background_sigma": 0.1202,
  "gamma": 1.5,
  "band_width": 96,
  "gap_open": 12.0,
  "gap_extend": 2.0,
  "xdrop": 50.0,
  "score_min": -4.0,
  "score_max": 8.0,
  "evd_samples": 1000,
  "background_samples": 5000,
  "random_seed": 42,
  "num_sequences": 119
}
```

**Important**: All alignment parameters (gap penalties, band width, etc.) are saved because the EVD parameters are **only valid** for alignments performed with identical settings.

---

## Phase 2: E-value Computation During Alignment

### 2.1 Real Alignment Scoring

**Module**: `align_candidates.py`

1. Load pre-computed EVD parameters from `evd_params.json` ([align_candidates.py:945-993](../bin/align_candidates.py#L945-L993))
2. Sample background statistics from the current embedding set (same as Phase 1.1)
3. For each cluster of seeds:
   - Extract query and target embedding regions
   - Perform banded Smith-Waterman alignment
   - Obtain raw alignment score $S$

The alignment uses **identical parameters** to Phase 1.3 to ensure statistical validity.

---

### 2.2 Bit Score Calculation

**Purpose**: Normalize scores to a database-independent scale.

**Formula** (function: `calculate_evalue`, [align_candidates.py:377-415](../bin/align_candidates.py#L377-L415)):

$$S' = \frac{\lambda \times S - \ln(K)}{\ln(2)}$$

Where:
- $S$ = raw alignment score
- $\lambda$ = lambda parameter from Phase 1
- $K$ = K parameter from Phase 1
- $\ln$ = natural logarithm

**Interpretation**:
- Bit score is database-independent (unlike raw scores)
- Represents the score in units of "doublings" (powers of 2)
- Higher bit score = more significant alignment
- Typical range: 0-1000+ (no theoretical upper bound)

---

### 2.3 E-value Calculation

**Formula** ([align_candidates.py:410-413](../bin/align_candidates.py#L410-L413)):

$$E = m \times n \times 2^{-S'}$$

Where:
- $m$ = effective query length (database size in nucleotides)
- $n$ = effective target length (query length in nucleotides)
- $S'$ = bit score from step 2.2
- $m \times n$ = total search space

**Detailed Calculation** ([align_candidates.py:377-415](../bin/align_candidates.py#L377-L415)):

```python
def calculate_evalue(
    score: float,
    query_len: int,
    target_len: int,
    lambda_param: float,
    K_param: float,
    m: int,  # total database length
    n: int   # total query length
) -> Tuple[float, float]:
    # Effective lengths (account for edge effects)
    # Use full lengths as effective lengths
    m_eff = max(1, m)
    n_eff = max(1, n)

    # Search space
    search_space = m_eff * n_eff

    # Bit score
    bit_score = (lambda_param * score - np.log(K_param)) / np.log(2.0)

    # E-value
    evalue = search_space * np.power(2.0, -bit_score)

    return evalue, bit_score
```

**Search Space Interpretation**:
- $m \times n$ represents the total number of possible alignment start positions
- For a database with $N$ sequences of average length $L_{\text{db}}$ and a query of length $L_q$:
  - $m = N \times L_{\text{db}}$ (total database size)
  - $n = L_q$ (query size)
  - Search space $\approx N \times L_{\text{db}} \times L_q$

**Example Calculation**:

Given:
- Raw score $S = 150$
- $\lambda = 0.00556$, $K = 0.5144$ (from evd_params.json)
- Query length = 500 nt
- Database size = 100,000 nt (e.g., 500 sequences × 200 nt average)

**Step 1 - Bit score:**
$$S' = \frac{0.00556 \times 150 - \ln(0.5144)}{\ln(2)}$$
$$S' = \frac{0.834 - (-0.665)}{0.693}$$
$$S' = \frac{1.499}{0.693} = 2.16$$

**Step 2 - Search space:**
$$m \times n = 100{,}000 \times 500 = 50{,}000{,}000$$

**Step 3 - E-value:**
$$E = 50{,}000{,}000 \times 2^{-2.16}$$
$$E = 50{,}000{,}000 \times 0.224$$
$$E \approx 11{,}200{,}000$$

This high E-value indicates the alignment is not statistically significant.

---

## Mathematical Foundations

### Karlin-Altschul Statistics

The E-value calculation is based on **Karlin-Altschul statistics** developed for BLAST:

**Key Theorem**: For local alignments of random sequences, the number of alignments with score $\geq S$ follows a Poisson distribution with mean:

$$E(S) = K \times m \times n \times e^{-\lambda S}$$

For large databases, this approximates the probability that at least one alignment has score $\geq S$:

$$P(\text{max score} \geq S) \approx 1 - e^{-E} \approx E \quad \text{(when } E \ll 1\text{)}$$

**Conversion to Bit Scores**: To make E-values database-independent, we use bit scores:

$$S' = \frac{\lambda S - \ln K}{\ln 2}$$
$$E = m \times n \times 2^{-S'}$$

This formulation has the property that each bit of score doubles the significance (halves the E-value).

---

### Gumbel Distribution Properties

The **Gumbel distribution** (Type I Extreme Value Distribution) models the maximum of a sample of random variables:

**PDF**:
$$f(x; \mu, \beta) = \frac{1}{\beta} \exp(-(z + e^{-z}))$$
where $z = \frac{x - \mu}{\beta}$

**CDF**:
$$F(x; \mu, \beta) = \exp\left(-\exp\left(-\frac{x - \mu}{\beta}\right)\right)$$

**Parameters**:
- $\mu$ (location): mode of the distribution
- $\beta$ (scale): controls spread; larger $\beta$ = more spread out

**Relation to Karlin-Altschul**:
$$\lambda = \frac{1}{\beta}$$
$$K = \exp(-\lambda\mu)$$

Therefore:
$$P(\text{score} \geq S) \approx 1 - \exp\left(-\exp(-\lambda(S - \mu))\right)$$
$$\approx \exp(-\lambda(S - \mu)) \quad \text{(for large } S\text{)}$$
$$= \exp(-\lambda S) \times \exp(\lambda\mu)$$
$$= K^{-1} \times \exp(-\lambda S)$$

When multiplied by search space $mn$, this gives the E-value formula.

---

## Implementation Details

### Shared Memory Optimization

For large databases, EVD parameter calculation can be memory-intensive when using multiple workers. GINflow implements a **shared memory strategy** ([calculate_evd_params.py:148-219](../bin/calculate_evd_params.py#L148-L219)):

**Problem**: Each worker process normally receives a copy of the entire embedding database
- With 8 workers: memory usage = $9\times$ (1 main process + 8 workers)

**Solution**: Store embeddings in shared memory accessible to all workers
- Memory usage $\approx 1\times$ (only one copy in shared memory)
- Workers attach to shared memory without copying data

**Implementation**:
1. Convert embeddings dict to flat numpy array in `/dev/shm`
2. Store metadata (sequence IDs, offsets, lengths) separately
3. Workers create `SharedEmbeddings` accessor using the shared memory name
4. Cleanup after processing completes

**Requirements**:
- Sufficient `/dev/shm` space (checked automatically)
- Docker containers need `--shm-size` flag (default: 6GB in GINflow)

---

### Sequence Sampling Strategy

**Problem**: EVD calculation can be slow for very large databases (millions of sequences).

**Solution**: Random sampling ([calculate_evd_params.py:129-145](../bin/calculate_evd_params.py#L129-L145))

**Method**:
1. Collect all sequence IDs from embeddings file (lightweight pass)
2. If `--sampled-sequences` > 0 and less than total:
   - Randomly sample specified number of sequences
   - Load only embeddings for sampled sequences
3. Perform EVD calculation on the sample

**Rationale**:
- EVD parameters characterize the scoring scheme, not individual sequences
- A representative sample (e.g., 50,000 sequences) provides accurate parameter estimates
- Dramatically reduces computation time for large databases

**Trade-offs**:
- Faster computation
- Reduced memory usage
- Slightly increased variance in parameter estimates (negligible with large samples)

---

### Numerical Stability Considerations

**Small σ protection** ([calculate_evd_params.py:536-538](../bin/calculate_evd_params.py#L536-L538)):
```python
if sigma < 1e-6:
    sigma = 1.0
```
Prevents division by zero when embeddings are nearly identical.

**Score clamping**:
```python
score = np.clip(z_score, score_min, score_max)
```
Prevents extreme z-scores from causing numerical overflow in exponential calculations.

**Minimum score threshold** ([calculate_evd_params.py:1017-1023](../bin/calculate_evd_params.py#L1017-L1023)):
```python
if len(null_scores) < 10:
    return 0.1, 0.01  # Default fallback parameters
```
Ensures sufficient data for robust parameter estimation.

---

## Parameter Configuration

### EVD Parameter Calculation (`calculate_evd_params`)

Configured in [modules/calculate_evd_params/main.nf](../modules/calculate_evd_params/main.nf):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--background-samples` | 10,000 | Number of random node pairs for background statistics |
| `--evd-samples` | 1,000 | Number of random alignments for null distribution |
| `--sampled-sequences` | 50,000 | Max sequences to sample (0 = use all) |
| `--gamma` | 1.5 | Z-score scaling factor |
| `--band-width` | 96 | Alignment band width (nt) |
| `--gap-open` | 12 | Gap opening penalty |
| `--gap-extend` | 2 | Gap extension penalty |
| `--xdrop` | 50 | X-drop termination threshold |
| `--score-min` | -4 | Minimum score clamp |
| `--score-max` | 8 | Maximum score clamp |
| `--random-seed` | 42 | Random seed for reproducibility |
| `--workers` | 6 | Number of parallel workers |

### E-value Calculation (`align_candidates`)

Configured in [modules/align_candidates/main.nf](../modules/align_candidates/main.nf):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--calculate-evalue` | `true` | Enable E-value calculation |
| `--evd-lambda` | (from JSON) | Pre-computed λ parameter |
| `--evd-K` | (from JSON) | Pre-computed K parameter |

**Note**: E-value calculation requires `evd_params.json` from the `calculate_evd_params` step.

---

## Output Format

### EVD Parameters File (`evd_params.json`)

```json
{
  "evd_lambda": 0.005560146185312788,
  "evd_K": 0.5144443054339665,
  "background_mu": 0.058035311607411134,
  "background_sigma": 0.12017227046781641,
  "gamma": 1.5,
  "band_width": 96,
  "gap_open": 12.0,
  "gap_extend": 2.0,
  "xdrop": 50.0,
  "score_min": -4.0,
  "score_max": 8.0,
  "evd_samples": 1000,
  "background_samples": 5000,
  "random_seed": 42,
  "num_sequences": 119,
  "sampled_sequences": 50000,
  "total_sequences_in_database": 500,
  "sequences_sampled": true
}
```

### Alignment Output (`alignments.tsv`)

E-values are added as a column in the alignment results:

```
query_id  target_id  score  evalue      bit_score  alignment_length  ...
seq_001   seq_042    156.2  0.000234    45.3       180              ...
seq_001   seq_103    142.8  0.00156     42.1       165              ...
seq_002   seq_089    138.5  0.00421     40.8       172              ...
```

**Column Descriptions**:
- `score`: Raw Smith-Waterman alignment score (z-scored)
- `evalue`: Expected number of alignments with score ≥ observed score by chance
- `bit_score`: Database-independent normalized score (optional, if calculated)

---

## Interpretation Guidelines

### E-value Thresholds

**Standard BLAST Guidelines** (adapted for GINflow):

| E-value Range | Interpretation | Action |
|---------------|----------------|--------|
| < 0.001 | Highly significant | Strong homology evidence |
| 0.001 - 0.01 | Significant | Probable homology |
| 0.01 - 0.1 | Marginally significant | Possible homology, verify manually |
| 0.1 - 1.0 | Weak signal | Likely noise, low confidence |
| > 1.0 | Not significant | Random similarity |

**Context-Specific Adjustments**:
- **Small databases** (< 1,000 sequences): Use stricter threshold (e.g., 0.0001)
- **Large databases** (> 100,000 sequences): Can use relaxed threshold (e.g., 0.01)
- **Hypothesis generation**: Consider E < 1.0 for exploratory analysis
- **High-confidence filtering**: Use E < 0.001 for publication-quality results

---

### Bit Score Interpretation

**Advantage**: Database-independent comparison across different searches.

**Rough Guidelines**:
- **Bit score > 50**: Strong similarity (equivalent to E < 0.001 in typical databases)
- **Bit score 30-50**: Moderate similarity
- **Bit score < 30**: Weak similarity

**Property**: Each additional bit halves the E-value:
$$\begin{align}
\text{Bit score} = N \quad &\rightarrow \quad E \approx X \\
\text{Bit score} = N+1 \quad &\rightarrow \quad E \approx X/2 \\
\text{Bit score} = N+10 \quad &\rightarrow \quad E \approx X/1024
\end{align}$$

---

## Assumptions and Limitations

### Assumptions

1. **Independence**: Alignments are independent (approximation; real alignments may overlap)
2. **Gumbel Distribution**: Maximum scores follow Gumbel statistics (valid for local alignments)
3. **Sequence Composition**: Database sequences have similar composition to sampled sequences
4. **Parameter Consistency**: Alignment parameters match between EVD estimation and real searches
5. **Large Database**: Formula most accurate when E << m × n

### Limitations

1. **Edge Effects**:
   - Current implementation uses full sequence lengths
   - Could be improved with effective length corrections

2. **Compositional Bias**:
   - No composition-based statistics (unlike BLAST's `-comp_based_stats`)
   - Could affect E-values for sequences with unusual embedding distributions

3. **Multiple Testing**:
   - E-values don't account for multiple hypothesis testing across many queries
   - Consider Bonferroni correction: $E_{\text{corrected}} = E \times \text{num\_queries}$

4. **Small Database Inaccuracy**:
   - Poisson approximation breaks down when E >> 1
   - For small databases, use more conservative thresholds

5. **Parameter Mismatch**:
   - E-values invalid if alignment parameters differ from EVD estimation
   - Always use `evd_params.json` from matching parameter set

---

## Troubleshooting

### Problem: E-values are all very large (> 100)

**Possible Causes**:
- Database too small (small m × n search space)
- Alignment parameters too strict (low scores)
- EVD parameters estimated from different settings

**Solutions**:
1. Verify `evd_params.json` matches alignment parameters
2. Check database size (should have hundreds to thousands of sequences)
3. Review alignment parameter settings (gap penalties, gamma, etc.)

---

### Problem: E-values are all very small (< 0.0001)

**Possible Causes**:
- Database too large (inflated m × n)
- Alignment parameters too permissive (inflated scores)
- Bug in length calculations

**Solutions**:
1. Verify sequence length calculations are correct
2. Check for parameter mismatches
3. Inspect raw alignment scores for sanity

---

### Problem: EVD parameter estimation fails

**Error**: "Not enough sequences for EVD estimation"

**Solutions**:
1. Ensure database has at least 20 sequences with length ≥ 20
2. Check embedding loading is successful
3. Verify TSV format matches expected columns

**Error**: "Insufficient /dev/shm space"

**Solutions**:
1. Increase Docker `--shm-size` (e.g., `--shm-size=8g`)
2. Reduce `--sampled-sequences` to use fewer sequences
3. Disable shared memory by setting `workers=1`

---

### Problem: E-values don't match expectations

**Check**:
1. **Reproducibility**: Run with same `--random-seed` - results should be identical
2. **Null Distribution**: Inspect `evd_params.json` - typical λ ≈ 0.001-0.01, K ≈ 0.1-1.0
3. **Bit Scores**: Should be positive for significant alignments
4. **Score Distribution**: Plot alignment scores - should show tail above background

**Validation**:
```python
import numpy as np

# Check if E-value formula is working
score = 150
lambda_param = 0.00556
K_param = 0.5144
m = 100000
n = 500

bit_score = (lambda_param * score - np.log(K_param)) / np.log(2)
evalue = m * n * np.power(2, -bit_score)

print(f"Bit score: {bit_score:.2f}")
print(f"E-value: {evalue:.2e}")
```

---

## References

### Primary Literature

1. **Karlin, S. & Altschul, S.F. (1990)**
   "Methods for assessing the statistical significance of molecular sequence features by using general scoring schemes"
   *Proc. Natl. Acad. Sci. USA* 87:2264-2268
   [doi:10.1073/pnas.87.6.2264](https://doi.org/10.1073/pnas.87.6.2264)

2. **Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990)**
   "Basic local alignment search tool"
   *J. Mol. Biol.* 215:403-410
   [doi:10.1016/S0022-2836(05)80360-2](https://doi.org/10.1016/S0022-2836(05)80360-2)

3. **Altschul, S.F., Madden, T.L., Schäffer, A.A., Zhang, J., Zhang, Z., Miller, W. & Lipman, D.J. (1997)**
   "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs"
   *Nucleic Acids Res.* 25:3389-3402
   [doi:10.1093/nar/25.17.3389](https://doi.org/10.1093/nar/25.17.3389)

### Statistical Theory

4. **Gumbel, E.J. (1958)**
   *Statistics of Extremes*
   Columbia University Press, New York

5. **Dembo, A., Karlin, S. & Zeitouni, O. (1994)**
   "Limit distribution of maximal non-aligned two-sequence segmental score"
   *Ann. Probab.* 22:2022-2039

### Implementation Guides

6. **NCBI BLAST Documentation**
   "The Statistics of Sequence Similarity Scores"
   [https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html](https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html)

7. **Eddy, S.R. (2004)**
   "What is a hidden Markov model?"
   *Nature Biotechnology* 22:1315-1316
   (Contains accessible explanation of E-value interpretation)

---

## Related Documentation

- [Alignment Phase Details](ALIGNMENT_README.md) - Full description of the alignment algorithm
- [EVALUE_IMPLEMENTATION_GUIDE.md](../EVALUE_IMPLEMENTATION_GUIDE.md) - Development guide for E-value integration
- [Nextflow Configuration](../nextflow.config) - Parameter settings and defaults

---

## Summary

GINflow's E-value calculation follows the proven BLAST statistical framework:

1. **Estimate null distribution**: Align shuffled sequences to get score distribution
2. **Fit Gumbel model**: Extract λ and K parameters via maximum likelihood
3. **Calculate bit scores**: Normalize raw scores to database-independent scale
4. **Compute E-values**: Multiply survival probability by search space

This approach provides statistically rigorous significance assessment for embedding-based alignments, enabling confident identification of true homologs while filtering random similarities.

**Key Insight**: E-values answer the question: *"How many alignments this good would I expect to see by chance?"* - providing an intuitive, database-aware measure of alignment significance.
