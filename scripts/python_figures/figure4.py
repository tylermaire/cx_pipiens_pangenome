"""
Figure 4 — Pairwise synteny dot plots in a 2x3 grid.

Each panel shows a nucmer alignment for one pair of ingroup genomes, drawn from
the show-coords -THrd output. Forward-strand alignments are plotted in dark
grey, inverted alignments in red.
"""
from __future__ import annotations
import argparse
from pathlib import Path
from itertools import combinations

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd

from style import (apply_style, fig_size, save_publication, panel_letter,
                   GENOME_PRETTY, INGROUP, FORWARD_COLOR, INVERTED_COLOR)
from data_loaders import load_coords


def draw_dot_panel(ax, coords: pd.DataFrame, ref_label: str, qry_label: str,
                   min_len_kb: float = 5.0, min_idy: float = 90.0):
    """Plot one pairwise dot plot.

    `coords` columns from show-coords -THrd:
        S1 E1 S2 E2 LEN1 LEN2 IDY R Q
    Inverted alignments have S2 > E2 (query end < query start).
    """
    if coords is None or len(coords) == 0:
        ax.text(0.5, 0.5, "No alignments", transform=ax.transAxes,
                ha="center", va="center", fontsize=9, color="grey")
        ax.set_title(f"{ref_label} (x) vs {qry_label} (y)", fontsize=8)
        return

    # Filter short / low-identity alignments
    df = coords.copy()
    df["len_kb"] = df["LEN1"] / 1000
    df = df[(df["len_kb"] >= min_len_kb) & (df["IDY"] >= min_idy)]
    if len(df) == 0:
        ax.text(0.5, 0.5, "No alignments after filter",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=8, color="grey")
        return

    df["inverted"] = (df["E2"] - df["S2"]) < 0

    # Forward
    fwd = df[~df["inverted"]]
    for _, r in fwd.iterrows():
        ax.plot([r["S1"]/1e6, r["E1"]/1e6], [r["S2"]/1e6, r["E2"]/1e6],
                color=FORWARD_COLOR, lw=0.4, alpha=0.7)
    # Inverted
    inv = df[df["inverted"]]
    for _, r in inv.iterrows():
        ax.plot([r["S1"]/1e6, r["E1"]/1e6], [r["S2"]/1e6, r["E2"]/1e6],
                color=INVERTED_COLOR, lw=0.5, alpha=0.85)

    ax.set_title(f"{ref_label} (x) vs {qry_label} (y)",
                 fontsize=8, style="italic", pad=2)
    ax.set_xlabel(f"{ref_label} (Mb)", fontsize=7)
    ax.set_ylabel(f"{qry_label} (Mb)", fontsize=7)
    ax.tick_params(axis="both", which="major", labelsize=6)
    ax.set_aspect("equal", adjustable="datalim")
    # Light grid
    ax.grid(True, color="#DDDDDD", lw=0.3, alpha=0.7)


def build_figure(repo_root: Path, out_basename: str):
    apply_style()
    pairs = list(combinations(INGROUP, 2))[:6]   # 6 pairs from 4 genomes

    fig = plt.figure(figsize=fig_size(180, 130))
    gs = gridspec.GridSpec(2, 3, hspace=0.45, wspace=0.4,
                           left=0.07, right=0.97, top=0.93, bottom=0.10)

    letters = ["A", "B", "C", "D", "E", "F"]
    for i, (ref, qry) in enumerate(pairs):
        row, col = divmod(i, 3)
        ax = fig.add_subplot(gs[row, col])
        coords = load_coords(repo_root, ref, qry)
        if coords is None:
            # Try reverse direction (nucmer is sensitive to order)
            coords = load_coords(repo_root, qry, ref)
        draw_dot_panel(ax, coords, GENOME_PRETTY[ref], GENOME_PRETTY[qry])
        panel_letter(ax, letters[i], x=-0.20, y=1.10, fontsize=10)

    # Bottom legend
    fig.text(0.5, 0.04,
             "Collinear alignments (dark grey)  |  inverted alignments (red)  "
             "|  ≥5 kb and ≥90% identity",
             ha="center", fontsize=8, color="#555")

    save_publication(fig, out_basename)
    plt.close(fig)
    print(f"Wrote {out_basename}.pdf and {out_basename}.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path(".").resolve())
    ap.add_argument("--out", type=str, default="results/figures/figure4_synteny")
    args = ap.parse_args()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    build_figure(args.repo_root, args.out)


if __name__ == "__main__":
    main()
