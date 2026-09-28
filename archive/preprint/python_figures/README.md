# Python figure suite

Python alternative to `scripts/Cx_pipiens_pangenome_figures.R`. The two
suites are kept side-by-side intentionally: the R script remains the
reference for the figures in the original manuscript, while this Python
suite generates the publication-quality versions for the revised paper
using the same data files.

## Layout

```
scripts/python_figures/
├── env.yaml          # conda env (matplotlib, seaborn, upsetplot, toytree, …)
├── style.py          # shared rcParams, colour palette, helpers
├── data_loaders.py   # read partitions, tree, CAFE5, ANI, TE, coords
├── figure1.py        # Fig 1: UpSet + accumulation + composition + form-specific
├── figure2.py        # Fig 2: ML tree + CAFE5 divergent bars
├── figure3.py        # Fig 3: ANI heatmap + TE composition + protein counts
├── figure4.py        # Fig 4: 2x3 synteny dot-plot grid
├── figure5.py        # Fig 5: pangenome flower diagram
├── run_all.py        # run all five in order
└── README.md         # this file
```

## Quick start

```bash
# Build the conda env (one time)
mamba env create -f scripts/python_figures/env.yaml
conda activate cxpan_figures

# Render every figure
python scripts/python_figures/run_all.py
```

Outputs land in `results/figures/`:

```
results/figures/figure1_pangenome.{pdf,png}
results/figures/figure2_tree_cafe.{pdf,png}
results/figures/figure3_ani_te_proteins.{pdf,png}
results/figures/figure4_synteny.{pdf,png}
results/figures/figure5_flower.{pdf,png}
```

Each figure script is independently runnable if you only want one:

```bash
python scripts/python_figures/figure2.py \
    --repo-root . \
    --out results/figures/figure2_tree_cafe
```

## Design conventions

The suite follows the conventions current in the comparative-genomics /
pangenome literature ca. 2022–2025:

- **Colours.** Okabe-Ito palette, colour-blind safe. Each genome gets a
  stable colour throughout (`style.GENOME_COLORS`); compartment colours
  (core / shell / cloud) are muted blue / orange / grey.
- **Tip ordering.** Every panel that lists genomes (CAFE5 bars, ANI
  heatmap rows, TE bars, protein counts) uses the same tip order
  matching the species-tree topology recovered from the corrected
  supermatrix: (Cx. pallens, Cx. quinquefasciatus, Cx. molestus,
  Cx. pipiens).
- **Format.** PDF (vector) for publication submission plus a 300-dpi PNG
  preview. To produce a 600-dpi LZW TIFF, edit the `save_publication`
  call in any figure script and add `tiff=True`.
- **Fonts.** Arial 9 pt body, TrueType-embedded (`pdf.fonttype = 42`) so
  the MDPI typesetter accepts it.
- **Page widths.** Figure widths are set in millimetres matching MDPI
  Insects column widths (85 mm single, 180 mm double).

## What each figure plots

### Figure 1 — Pangenome architecture
1. UpSet plot of orthogroup intersections across the four genomes, with
   the all-4 (core) intersection highlighted.
2. Pangenome accumulation curve in Tettelin style: power-law fit for
   the pan-genome, exponential decay for the core.
3. Per-genome stacked bar of core / shell / cloud orthogroup membership.
4. Bar chart of form-specific (cloud) orthogroup counts per genome.

### Figure 2 — Phylogenomics + gene-family evolution
- Left: maximum-likelihood species tree with ultrafast-bootstrap support
  labels on internal branches.
- Right: divergent horizontal bar chart of per-branch CAFE5
  expansions (green) and contractions (red). The y-axis order matches
  the tree's tip order.

### Figure 3 — Genome-level comparisons
- ANI heatmap (square, values annotated, `rocket_r` colourmap zoomed to
  the actual data range so off-diagonal differences are visible).
- TE composition stacked bars per genome.
- Predicted protein-coding gene counts per genome.

### Figure 4 — Synteny dot plots
2×3 grid of pairwise nucmer alignments. Collinear alignments are dark
grey; inverted alignments are red. Filters: ≥5 kb alignment length and
≥90% identity (configurable in `figure4.py`).

### Figure 5 — Flower diagram
Central circle shows the core orthogroup count; petals around it show
each genome's form-specific (cloud) count.

## When the upstream pipeline hasn't finished

Each loader in `data_loaders.py` falls back to manually-curated numbers
from the manuscript when the corresponding result file doesn't exist
yet, so figures render with paper-baseline data even before the real
results land. Once the pipeline completes, just rerun `run_all.py` and
the figures regenerate against the new data.
