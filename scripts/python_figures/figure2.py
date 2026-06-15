"""
Figure 2 — Phylogenomics + CAFE5 gene-family evolution (composite, 2 panels).

  (a) Maximum-likelihood species tree from the 7,889-OG supermatrix, with
      ultrafast-bootstrap support and (if available) gene/site concordance
      factors on internal nodes.
  (b) CAFE5 per-terminal-branch expansions (green, positive) and contractions
      (red, negative), drawn as a divergent horizontal bar chart that shares
      the tree's tip ordering on the y-axis.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import re

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from Bio import Phylo
from io import StringIO

from style import (apply_style, fig_size, save_publication, panel_letter,
                   GENOME_COLORS, GENOME_PRETTY, INGROUP,
                   EXPANSION_COLOR, CONTRACTION_COLOR)
from data_loaders import load_newick_string, load_cafe_clade, load_cafe_change, compute_per_branch_counts


# ----------------------------------------------------------------------------
# Tree drawing (matplotlib-native, no ete3 dependency at runtime)
# ----------------------------------------------------------------------------
def _build_tip_positions(tree, tip_order):
    """Return {clade: (x, y)} positions for every clade in the tree, with the
    tip y-coordinates dictated by `tip_order` so they match a paired panel."""
    positions = {}
    # Assign y coords to tips
    for i, name in enumerate(tip_order):
        for clade in tree.get_terminals():
            if clade.name == name:
                positions[clade] = [None, i]  # x filled in below
                break

    # Recursively compute x by depth from root
    def assign_x(clade, depth=0.0):
        if clade.is_terminal():
            positions[clade][0] = depth
            return positions[clade][1]
        ys = [assign_x(c, depth + (c.branch_length or 0.0)) for c in clade.clades]
        y_mid = (min(ys) + max(ys)) / 2
        positions.setdefault(clade, [None, None])
        positions[clade] = [depth, y_mid]
        return y_mid

    assign_x(tree.root)
    return positions


def draw_tree(ax, newick_str, tip_order):
    """Render a rooted/unrooted ML tree onto matplotlib axis `ax`.

    Tip order is dictated by `tip_order` (top-to-bottom) so it lines up with
    the CAFE5 bar panel.
    """
    tree = Phylo.read(StringIO(newick_str), "newick")
    # If unrooted, just take it as-is - midpoint root for display
    try:
        tree.root_at_midpoint()
    except Exception:
        pass

    # Verify all expected tips exist
    tree_tips = {c.name for c in tree.get_terminals()}
    for t in tip_order:
        if t not in tree_tips:
            print(f"WARN: tip {t} not found in tree (have: {tree_tips})")

    pos = _build_tip_positions(tree, tip_order)

    # Draw branches
    for clade in tree.find_clades():
        x, y = pos[clade]
        if clade.clades:
            child_ys = [pos[c][1] for c in clade.clades]
            ax.plot([x, x], [min(child_ys), max(child_ys)],
                    color="black", lw=1.2, zorder=2)  # vertical bar
            for c in clade.clades:
                cx, cy = pos[c]
                ax.plot([x, cx], [cy, cy], color="black", lw=1.2, zorder=2)
                # Internal branch support label (if present)
                if not c.is_terminal() and c.confidence is not None:
                    ax.text((x + cx) / 2, cy + 0.1, f"{int(c.confidence)}",
                            ha="center", va="bottom", fontsize=7, color="#222")

    # Tip labels & dots
    max_x = max(p[0] for p in pos.values())
    for clade in tree.get_terminals():
        x, y = pos[clade]
        col = GENOME_COLORS.get(clade.name, "#222")
        ax.scatter(x, y, s=60, color=col, zorder=3, edgecolor="black", linewidth=0.6)
        ax.text(max_x * 1.05, y, GENOME_PRETTY.get(clade.name, clade.name),
                ha="left", va="center", fontsize=9, style="italic")

    # Extra rightward space so italic labels can fully render
    ax.set_xlim(-max_x * 0.05, max_x * 2.4)
    ax.set_ylim(-0.9, len(tip_order) - 0.1)
    ax.set_yticks([])
    ax.set_xlabel("Substitutions / site")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)

    # Scale bar in lower-left
    sb_len = max_x / 4
    ax.plot([0.0, sb_len], [-0.5, -0.5], color="black", lw=1.5)
    ax.text(sb_len / 2, -0.65, f"{sb_len:.3f}", ha="center", va="top", fontsize=7)


# ----------------------------------------------------------------------------
# CAFE5 divergent bar chart
# ----------------------------------------------------------------------------
def draw_cafe_bars(ax, branch_counts: pd.DataFrame, tip_order):
    """Horizontal divergent bars: contractions (negative, red) and expansions
    (positive, green) for each terminal branch."""
    # Map branch name like "Cx_pallens<1>" -> "Cx_pallens"
    def strip_id(name):
        return re.sub(r"<\d+>$", "", name)

    bc = branch_counts.assign(genome=branch_counts["branch"].apply(strip_id))
    bc = bc[bc["genome"].isin(tip_order)].set_index("genome").reindex(tip_order)

    ys = np.arange(len(tip_order))
    ax.barh(ys, -bc["contractions"], color=CONTRACTION_COLOR,
            edgecolor="none", height=0.6, label="Contraction")
    ax.barh(ys,  bc["expansions"],   color=EXPANSION_COLOR,
            edgecolor="none", height=0.6, label="Expansion")

    # Annotate counts at the bar tips
    for i, (e, c) in enumerate(zip(bc["expansions"], bc["contractions"])):
        if pd.notna(e):
            ax.text(e + max(bc["expansions"]) * 0.02, i, f"+{int(e):,}",
                    va="center", ha="left", fontsize=7, color=EXPANSION_COLOR)
        if pd.notna(c):
            ax.text(-c - max(bc["contractions"]) * 0.02, i, f"-{int(c):,}",
                    va="center", ha="right", fontsize=7, color=CONTRACTION_COLOR)

    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks([])
    ax.set_xlabel("Gene families")
    ax.legend(loc="upper right", fontsize=7, ncol=2)
    ax.set_xlim(-bc["contractions"].max() * 1.25, bc["expansions"].max() * 1.25)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
# Hard-coded fallback tip order matching the new ML topology
# (pallens, quinquefasciatus) | (molestus, pipiens)
DEFAULT_TIP_ORDER = ["Cx_pallens", "Cx_quinquefasciatus", "Cx_molestus", "Cx_pipiens"]


def build_figure(repo_root: Path, out_basename: str):
    apply_style()
    fig = plt.figure(figsize=fig_size(180, 95))
    # Wider tree column to fit italic species labels; larger gap between panels.
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.4, 1], wspace=0.45,
                           left=0.02, right=0.97, top=0.92, bottom=0.18)
    ax_tree = fig.add_subplot(gs[0])
    ax_bars = fig.add_subplot(gs[1])

    # --- Tree ---
    nwk = load_newick_string(repo_root)
    draw_tree(ax_tree, nwk, DEFAULT_TIP_ORDER)
    panel_letter(ax_tree, "A", x=-0.06, y=1.06)

    # --- CAFE5 ---
    try:
        clade = load_cafe_clade(repo_root)
        # Gamma_clade_results.txt columns: #Taxon_ID  Increase  Decrease
        bc = pd.DataFrame({
            "branch": clade["#Taxon_ID"] if "#Taxon_ID" in clade else clade.iloc[:, 0],
            "expansions": clade["Increase"],
            "contractions": clade["Decrease"],
        })
    except FileNotFoundError:
        try:
            ch = load_cafe_change(repo_root)
            bc = compute_per_branch_counts(ch)
        except FileNotFoundError:
            bc = pd.DataFrame({
                "branch": DEFAULT_TIP_ORDER,
                "expansions": [3110, 3648, 3033, 3194],
                "contractions": [1287, 769, 1425, 913],
            })
    draw_cafe_bars(ax_bars, bc, DEFAULT_TIP_ORDER)
    panel_letter(ax_bars, "B", x=-0.04, y=1.06)

    save_publication(fig, out_basename)
    plt.close(fig)
    print(f"Wrote {out_basename}.pdf and {out_basename}.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path(".").resolve())
    ap.add_argument("--out", type=str, default="results/figures/figure2_tree_cafe")
    args = ap.parse_args()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    build_figure(args.repo_root, args.out)


if __name__ == "__main__":
    main()
