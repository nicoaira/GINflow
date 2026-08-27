# Plotting and report

Every query run writes `report.html`. Structure and SW-matrix plots are
optional and are inlined in that page when requested.

## `report.html`

`WRITE_REPORT` (`bin/write_report.py`) builds a standalone HTML file:
no extra assets at view time. Open it in a browser.

- Hits grouped by query, ranked by database E-value
- Filterable by E-value
- Pagination: 10 / 25 / 50 / 100 / 150 per page (default 10)
- Theme: `--report_theme light` (default) or `dark`. A masthead toggle
  switches themes while viewing
- Each card: database E, pair E, bits, total score, max HSP score, HSP
  count, aligned span on both molecules
- Floating **Back to top** control after scrolling
- Masthead logo: `docs/images/ginflow_logo.svg`; favicon:
  `docs/images/ginflow_icon.svg`

RNArtistCore and SW plots use a two-column Query | Target (or cosine |
SW scores) grid. R4RNA is one full-width alignment arc plot per HSP.

## Structure plots

`--plot_backend`:

| Value | What is drawn |
|---|---|
| `none` (default) | No structure plots |
| `rnartistcore` | RNArtistCore 2D for query and target, per HSP |
| `r4rna` | One R4RNA alignment-coordinate SVG per HSP |
| `both` | Both of the above |

**RNArtistCore** colours the aligned span with
`--plot_highlight_colour` (default `#24B064`); the rest of the molecule
is gray. Conda package `nicolas.aira::rnartistcore=0.4.6` (OpenJDK 17–21).

**R4RNA** draws both structures on the alignment coordinates: query
arcs above, target arcs flipped below, identity ribbon between the
backbones. Shared pairs and matching bases use the highlight colour.
Conda package `nicolas.aira::r-r4rna=2.0.9`.

## Smith–Waterman matrices

`--plot_sw` (default `false`) writes two SVGs per plotted HSP:

1. Crop cosine matrix (query residues × target residues)
2. Substitution-score matrix with the SW traceback overlaid

Uses the GINFINITY-SW container (no extra image). The traceback uses
`--plot_highlight_colour`.

## How many plots

`--plot_max_pairs` (default 25) counts **unique query–target pairs per
query**, not HSP rows and not SVG files. Every HSP belonging to a
selected pair is rendered, so the SVG count can be larger than 25.

Each query gets its own draw task. Draw processes use `task.cpus`
workers (6 with the default `process_medium` label).

Published under `plots/`:

```
plots/rnartistcore/*.svg
plots/r4rna/*.svg
plots/sw/*.svg
```

Filenames use the alignment `baseName` so two slices of the same
accession do not overwrite each other.

## Typical command

```bash
nextflow run nicoaira/ginflow \
    -profile docker \
    --query queries.tsv \
    --database results/index \
    --plot_backend both \
    --plot_sw \
    --plot_max_pairs 25 \
    --report_theme light \
    --outdir search_results
```

Plots are the expensive optional tail of a query run. Keep
`--plot_backend none` while iterating on index settings; turn plots on
for a publication figure or a small query set.
