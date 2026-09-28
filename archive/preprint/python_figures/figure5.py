"""
Figure 5 — Pangenome flower diagram.

Central circle shows the core orthogroup count shared by all 4 genomes.
Each petal represents one genome's form-specific (cloud) orthogroups.
A central caption labels the shell orthogroup count.

Many recent pangenome papers are dropping flower diagrams in favour of the
information already in UpSet plots, but the flower is still a recognisable
"glance" figure for editors and is included here for completeness.
"""
from __future__ import annotations
import argparse
from pathlib import Path
from math import cos, sin, pi

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

from style import (apply_style, fig_size, save_publication, panel_letter,
                   GENOME_COLORS, GENOME_PRETTY, INGROUP)
from data_loaders import load_partitions, load_pangenome_summary


def compute_flower_data(repo_root: Path):
    """Return (core_count, shell_count, form_specific_dict, total_pan)."""
    summary = load_pangenome_summary(repo_root)
    core  = int(summary.loc[summary["compartment"] == "core",  "n_orthogroups"].iloc[0])
    shell = int(summary.loc[summary["compartment"] == "shell", "n_orthogroups"].iloc[0])
    cloud = int(summary.loc[summary["compartment"] == "cloud", "n_orthogroups"].iloc[0])
    total = core + shell + cloud

    # Form-specific per genome
    fs = {}
    for g in INGROUP:
        key = f"cloud_{g}_specific"
        row = summary[summary["compartment"] == key]
        fs[g] = int(row["n_orthogroups"].iloc[0]) if len(row) else 0
    return core, shell, cloud, total, fs


def draw_flower(ax, core, shell, cloud, total, fs):
    """Render the flower."""
    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-2.5, 2.5)
    ax.set_aspect("equal")
    ax.axis("off")

    # Petals (ellipses) rotated around centre
    petal_a = 1.8       # major radius
    petal_b = 0.7       # minor radius
    petal_offset = 1.2  # how far from centre the centroid sits
    for i, g in enumerate(INGROUP):
        angle = 90 + i * (360 / len(INGROUP))  # start at top, go clockwise
        rad = np.radians(angle)
        cx = petal_offset * cos(rad)
        cy = petal_offset * sin(rad)
        # Convert to matplotlib angle convention
        ell = mpatches.Ellipse((cx, cy), width=petal_a * 2, height=petal_b * 2,
                               angle=angle - 90, facecolor=GENOME_COLORS[g],
                               alpha=0.55, edgecolor="black", linewidth=0.6,
                               zorder=2)
        ax.add_patch(ell)
        # Label at petal tip
        tip_x = (petal_offset + petal_a * 0.85) * cos(rad)
        tip_y = (petal_offset + petal_a * 0.85) * sin(rad)
        ax.text(tip_x, tip_y, f"{GENOME_PRETTY[g]}\n{fs[g]} form-specific",
                ha="center", va="center", fontsize=8, style="italic",
                zorder=5)

    # Central circle with core count
    core_circle = mpatches.Circle((0, 0), 0.85, facecolor="#FAFAFA",
                                  edgecolor="black", linewidth=1.0, zorder=4)
    ax.add_patch(core_circle)
    ax.text(0, 0.15, f"{core:,}", ha="center", va="center",
            fontsize=14, fontweight="bold", zorder=5)
    ax.text(0, -0.20, "core\northogroups", ha="center", va="center",
            fontsize=8, zorder=5)

    # Shell count below
    ax.text(0, -2.3, f"shell: {shell:,}  |  pan-genome total: {total:,}",
            ha="center", va="center", fontsize=8, color="#555")


def build_figure(repo_root: Path, out_basename: str):
    apply_style()
    core, shell, cloud, total, fs = compute_flower_data(repo_root)
    fig, ax = plt.subplots(figsize=fig_size(120, 120))
    draw_flower(ax, core, shell, cloud, total, fs)
    save_publication(fig, out_basename)
    plt.close(fig)
    print(f"Wrote {out_basename}.pdf and {out_basename}.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path(".").resolve())
    ap.add_argument("--out", type=str, default="results/figures/figure5_flower")
    args = ap.parse_args()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    build_figure(args.repo_root, args.out)


if __name__ == "__main__":
    main()
