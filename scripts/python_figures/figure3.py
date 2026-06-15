"""
Figure 3 — Genome-level comparisons (composite 3 panels).

  (a) Pairwise ANI heatmap (4x4, values annotated)
  (b) TE composition stacked bars (interspersed / simple / low-complexity / unclassified)
  (c) Predicted protein-coding gene counts per genome
"""
from __future__ import annotations
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import seaborn as sns

from style import (apply_style, fig_size, save_publication, panel_letter,
                   GENOME_COLORS, GENOME_PRETTY, INGROUP)
from data_loaders import compute_ani_from_coords, load_te_summary, load_protein_counts


# ----------------------------------------------------------------------------
# Panel A — ANI heatmap
# ----------------------------------------------------------------------------
def draw_ani(ax, ani_df: pd.DataFrame):
    """Annotated ANI heatmap with viridis-like palette, scaled to the
    actual data range (avoids wasting dynamic range on the 100% diagonal)."""
    # Mask the diagonal so it doesn't dominate vmin/vmax
    mat = ani_df.reindex(index=INGROUP, columns=INGROUP).copy()
    off_diag = mat.where(~np.eye(len(mat), dtype=bool))
    vmin = float(off_diag.min().min()) - 0.05
    vmax = float(off_diag.max().max()) + 0.05

    sns.heatmap(mat, ax=ax, annot=True, fmt=".1f",
                cmap="rocket_r", vmin=vmin, vmax=vmax,
                annot_kws={"fontsize": 7, "color": "white"}, square=True,
                cbar_kws={"label": "ANI (%)", "shrink": 0.7, "pad": 0.04},
                linewidths=0.4, linecolor="white")
    # Use shorter abbreviations so labels fit
    short_labels = {
        "Cx_quinquefasciatus": "Cx. quinq.",
        "Cx_pallens":          "Cx. p. pallens",
        "Cx_molestus":         "Cx. p. molestus",
        "Cx_pipiens":          "Cx. p. pipiens",
    }
    ax.set_xticklabels([short_labels[g] for g in INGROUP],
                       rotation=40, ha="right", style="italic", fontsize=7)
    ax.set_yticklabels([short_labels[g] for g in INGROUP],
                       rotation=0, style="italic", fontsize=7)


# ----------------------------------------------------------------------------
# Panel B — TE stacked composition
# ----------------------------------------------------------------------------
def draw_te(ax, te_df: pd.DataFrame):
    """Stacked horizontal bars, one per genome."""
    te = te_df.reindex(INGROUP).fillna(0)
    classes = [c for c in ("Retroelements", "DNA_transposons", "Rolling_circles",
                           "Unclassified", "Interspersed", "Simple_repeats",
                           "Low_complexity") if c in te.columns]
    palette = {
        "Retroelements":    "#1F77B4",
        "DNA_transposons":  "#FF7F0E",
        "Rolling_circles":  "#2CA02C",
        "Unclassified":     "#8C564B",
        "Interspersed":     "#1F77B4",  # used as fallback if class breakdown unavailable
        "Simple_repeats":   "#9467BD",
        "Low_complexity":   "#7F7F7F",
    }

    ys = np.arange(len(te))
    left = np.zeros(len(te))
    for cls in classes:
        ax.barh(ys, te[cls], left=left, color=palette.get(cls, "#888"),
                edgecolor="white", linewidth=0.5, label=cls.replace("_", " "))
        left += te[cls].values

    ax.set_yticks(ys)
    short_labels = {
        "Cx_quinquefasciatus": "Cx. quinq.",
        "Cx_pallens":          "Cx. p. pallens",
        "Cx_molestus":         "Cx. p. molestus",
        "Cx_pipiens":          "Cx. p. pipiens",
    }
    ax.set_yticklabels([short_labels[g] for g in te.index],
                       fontsize=7, style="italic")
    ax.invert_yaxis()
    ax.set_xlabel("Percent of genome (%)")
    # Move legend outside the plot to the right so it doesn't overlap bars
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0),
              fontsize=7, ncol=1, frameon=False, borderaxespad=0.3,
              title=None)


# ----------------------------------------------------------------------------
# Panel C — protein counts
# ----------------------------------------------------------------------------
def draw_proteins(ax, counts: pd.Series):
    """Bar chart of predicted protein counts per genome."""
    counts = counts.reindex(INGROUP)
    bars = ax.bar(range(len(counts)), counts.values,
                  color=[GENOME_COLORS[g] for g in counts.index],
                  edgecolor="none")
    for i, v in enumerate(counts.values):
        ax.text(i, v + max(counts.values) * 0.01, f"{int(v):,}",
                ha="center", va="bottom", fontsize=8)
    short_labels = {
        "Cx_quinquefasciatus": "Cx. quinq.",
        "Cx_pallens":          "Cx. p. pallens",
        "Cx_molestus":         "Cx. p. molestus",
        "Cx_pipiens":          "Cx. p. pipiens",
    }
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels([short_labels[g] for g in counts.index],
                       rotation=40, ha="right", style="italic", fontsize=7)
    ax.set_ylabel("Predicted proteins")
    ax.set_ylim(0, max(counts) * 1.12)


# ----------------------------------------------------------------------------
# Main composite
# ----------------------------------------------------------------------------
def build_figure(repo_root: Path, out_basename: str):
    apply_style()
    fig = plt.figure(figsize=fig_size(180, 110))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.1, 1.3, 1.0], wspace=0.65,
                           left=0.10, right=0.97, top=0.92, bottom=0.32)
    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1])
    ax_c = fig.add_subplot(gs[2])

    draw_ani(ax_a, compute_ani_from_coords(repo_root))
    panel_letter(ax_a, "A", x=-0.20, y=1.05)

    draw_te(ax_b, load_te_summary(repo_root))
    panel_letter(ax_b, "B", x=-0.30, y=1.05)

    draw_proteins(ax_c, load_protein_counts(repo_root))
    panel_letter(ax_c, "C", x=-0.25, y=1.05)

    save_publication(fig, out_basename)
    plt.close(fig)
    print(f"Wrote {out_basename}.pdf and {out_basename}.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path(".").resolve())
    ap.add_argument("--out", type=str, default="results/figures/figure3_ani_te_proteins")
    args = ap.parse_args()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    build_figure(args.repo_root, args.out)


if __name__ == "__main__":
    main()
