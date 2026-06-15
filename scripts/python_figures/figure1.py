"""
Figure 1 — Pangenome architecture (composite 2x2 panel).

  (a) UpSet plot of orthogroup intersections across the 4 ingroup genomes
  (b) Pangenome accumulation curve (Tettelin-style; pan/core fits)
  (c) Per-genome stacked bar of core/shell/cloud composition
  (d) Form-specific (cloud) orthogroup counts per genome
"""
from __future__ import annotations
import argparse
from pathlib import Path
from itertools import combinations

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from style import (apply_style, fig_size, save_publication, panel_letter,
                   GENOME_COLORS, GENOME_PRETTY, COMPARTMENT_COLORS, INGROUP)
from data_loaders import load_partitions, load_pangenome_summary, build_membership_matrix


# ----------------------------------------------------------------------------
# Panel A — UpSet plot
# ----------------------------------------------------------------------------
def draw_upset(ax_upper, ax_lower, ax_matrix, membership: pd.DataFrame):
    """Draw a minimal hand-rolled UpSet plot.

    Composed of three sub-axes:
      ax_upper   : intersection-size bar chart
      ax_lower   : set-size bar chart on the left
      ax_matrix  : dot matrix on bottom-right
    """
    # Compute non-empty intersections, sorted by degree then size
    intersections = []
    for r in range(1, len(INGROUP) + 1):
        for combo in combinations(INGROUP, r):
            mask = np.ones(len(membership), dtype=bool)
            for g in INGROUP:
                mask &= membership[g].values if g in combo else ~membership[g].values
            n = int(mask.sum())
            if n > 0:
                intersections.append((combo, n))
    # Sort by degree desc (so all-4 first), then size desc within each degree
    intersections.sort(key=lambda kv: (-len(kv[0]), -kv[1]))
    n_inter = len(intersections)

    # Top: bar chart of intersection sizes
    xs = np.arange(n_inter)
    sizes = [n for _, n in intersections]
    colors = ["#1F77B4" if len(c) == len(INGROUP) else "#888888"
              for c, _ in intersections]
    bars = ax_upper.bar(xs, sizes, color=colors, width=0.85, edgecolor="none")
    ax_upper.set_ylabel("Orthogroup\nintersection size")
    ax_upper.set_xticks([])
    ax_upper.set_xlim(-0.5, n_inter - 0.5)
    # Annotate count above bars
    for x, s in zip(xs, sizes):
        if s > 0:
            ax_upper.text(x, s + max(sizes) * 0.02, f"{s:,}", ha="center",
                          va="bottom", fontsize=6, rotation=0)

    # Left: set-size bars (per-genome totals)
    set_sizes = {g: int(membership[g].sum()) for g in INGROUP}
    ys = np.arange(len(INGROUP))[::-1]
    ax_lower.barh(ys, [set_sizes[g] for g in INGROUP],
                  color=[GENOME_COLORS[g] for g in INGROUP], edgecolor="none")
    ax_lower.invert_xaxis()
    ax_lower.set_yticks(ys)
    ax_lower.set_yticklabels([GENOME_PRETTY[g] for g in INGROUP],
                             fontsize=8, style="italic")
    ax_lower.set_xlabel("Set size")

    # Matrix: dots showing genome membership of each intersection
    for i, (combo, _) in enumerate(intersections):
        for j, g in enumerate(INGROUP[::-1]):
            if g in combo:
                ax_matrix.scatter(i, j, s=42, color=GENOME_COLORS[g], zorder=3)
            else:
                ax_matrix.scatter(i, j, s=42, color="#DDDDDD", zorder=2)
        # Connect dots within a single intersection
        in_combo_y = [len(INGROUP) - 1 - INGROUP.index(g) for g in combo]
        if len(in_combo_y) > 1:
            ax_matrix.plot([i, i], [min(in_combo_y), max(in_combo_y)],
                           color="#444444", lw=1.2, zorder=1)
    ax_matrix.set_xlim(-0.5, n_inter - 0.5)
    ax_matrix.set_ylim(-0.5, len(INGROUP) - 0.5)
    ax_matrix.set_xticks([])
    ax_matrix.set_yticks(range(len(INGROUP)))
    ax_matrix.set_yticklabels([])
    ax_matrix.tick_params(left=False, bottom=False)
    for s in ("top", "right", "bottom", "left"):
        ax_matrix.spines[s].set_visible(False)


# ----------------------------------------------------------------------------
# Panel B — accumulation curve
# ----------------------------------------------------------------------------
def _pan_power(x, a, b, c):    return a * (x ** b) + c          # Heaps' law
def _core_decay(x, a, b, c):   return a * np.exp(-x / b) + c    # exponential decay

def draw_accumulation(ax, membership: pd.DataFrame, rng=None):
    """Compute and plot Tettelin-style pan + core accumulation curves."""
    rng = rng or np.random.default_rng(42)
    n_total = len(INGROUP)
    pan_pts, core_pts = [], []
    for k in range(1, n_total + 1):
        for combo in combinations(INGROUP, k):
            sub = membership[list(combo)]
            # An orthogroup is in the pan-genome of `combo` if any member has it
            pan = int(sub.any(axis=1).sum())
            # Core within `combo` = present in every selected genome
            core = int(sub.all(axis=1).sum())
            pan_pts.append((k, pan))
            core_pts.append((k, core))

    pan_df = pd.DataFrame(pan_pts, columns=["k", "n"])
    core_df = pd.DataFrame(core_pts, columns=["k", "n"])

    # Scatter all sampled subset values
    ax.scatter(pan_df["k"] + 0.07, pan_df["n"], s=14, color="#444",
               alpha=0.35, zorder=2, label=None)
    ax.scatter(core_df["k"] - 0.07, core_df["n"], s=14, color="#444",
               alpha=0.35, zorder=2, label=None)

    # Mean per k
    pan_mean = pan_df.groupby("k").mean()
    core_mean = core_df.groupby("k").mean()
    ks = np.array(pan_mean.index)

    # Fit (only meaningful if N>=3)
    if len(ks) >= 3:
        try:
            popt_pan, _ = curve_fit(_pan_power, ks, pan_mean["n"].values,
                                    p0=[1000, 0.3, 11000], maxfev=20000)
            xfine = np.linspace(1, n_total, 100)
            ax.plot(xfine, _pan_power(xfine, *popt_pan), color="#1F77B4",
                    lw=2, label=f"Pan-genome\n  $y = {popt_pan[0]:.0f}\\,x^{{{popt_pan[1]:.2f}}} + {popt_pan[2]:.0f}$")
        except Exception:
            ax.plot(ks, pan_mean["n"], "o-", color="#1F77B4", lw=2, label="Pan-genome")
        try:
            popt_core, _ = curve_fit(_core_decay, ks, core_mean["n"].values,
                                     p0=[5000, 2.0, 11000], maxfev=20000)
            ax.plot(xfine, _core_decay(xfine, *popt_core), color="#D62728",
                    lw=2, label=f"Core\n  $y = {popt_core[0]:.0f}\\,e^{{-x/{popt_core[1]:.2f}}} + {popt_core[2]:.0f}$")
        except Exception:
            ax.plot(ks, core_mean["n"], "o-", color="#D62728", lw=2, label="Core")
    else:
        ax.plot(ks, pan_mean["n"], "o-", color="#1F77B4", lw=2, label="Pan-genome")
        ax.plot(ks, core_mean["n"], "o-", color="#D62728", lw=2, label="Core")

    ax.set_xlabel("Number of genomes sampled")
    ax.set_ylabel("Orthogroup count")
    ax.set_xticks(range(1, n_total + 1))
    ax.legend(loc="center right", fontsize=7, frameon=False)


# ----------------------------------------------------------------------------
# Panel C — per-genome stacked composition
# ----------------------------------------------------------------------------
def draw_composition(ax, partitions: pd.DataFrame):
    """Stacked bars: each genome's orthogroup count partitioned core/shell/cloud."""
    rows = []
    for g in INGROUP:
        sub = partitions[partitions[g] > 0]
        for comp in ("core", "shell", "cloud"):
            n = int((sub["compartment"] == comp).sum())
            rows.append({"genome": g, "compartment": comp, "count": n})
    df = pd.DataFrame(rows).pivot(index="genome", columns="compartment", values="count")
    df = df.reindex(INGROUP)[["core", "shell", "cloud"]]

    bottom = np.zeros(len(df))
    for comp in ("core", "shell", "cloud"):
        ax.bar(range(len(df)), df[comp], bottom=bottom,
               color=COMPARTMENT_COLORS[comp], label=comp.capitalize(),
               edgecolor="white", linewidth=0.5)
        # Inline count text inside each segment
        for i, v in enumerate(df[comp]):
            if v > 200:
                ax.text(i, bottom[i] + v / 2, f"{v:,}", ha="center",
                        va="center", fontsize=7, color="white" if comp == "core" else "black")
        bottom += df[comp].values
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels([GENOME_PRETTY[g] for g in INGROUP],
                       rotation=30, ha="right", style="italic", fontsize=8)
    ax.set_ylabel("Number of orthogroups")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0),
              fontsize=8, title="Compartment", title_fontsize=8)


# ----------------------------------------------------------------------------
# Panel D — form-specific cloud counts
# ----------------------------------------------------------------------------
def draw_form_specific(ax, partitions: pd.DataFrame):
    """Bar chart of form-specific (cloud) orthogroups per genome."""
    cloud = partitions[partitions["compartment"] == "cloud"]
    fs_counts = []
    for g in INGROUP:
        others = [o for o in INGROUP if o != g]
        mask = (cloud[g] > 0) & (cloud[others].sum(axis=1) == 0)
        fs_counts.append((g, int(mask.sum())))

    df = pd.DataFrame(fs_counts, columns=["genome", "count"]).set_index("genome")
    df = df.reindex(INGROUP)
    bars = ax.bar(range(len(df)), df["count"],
                  color=[GENOME_COLORS[g] for g in df.index],
                  edgecolor="none")
    for i, v in enumerate(df["count"]):
        ax.text(i, v + max(df["count"]) * 0.02, f"{v}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels([GENOME_PRETTY[g] for g in df.index],
                       rotation=30, ha="right", style="italic", fontsize=8)
    ax.set_ylabel("Form-specific orthogroups")


# ----------------------------------------------------------------------------
# Main composite layout
# ----------------------------------------------------------------------------
def build_figure(repo_root: Path, out_basename: str):
    apply_style()
    partitions = load_partitions(repo_root)
    membership = build_membership_matrix(partitions)

    fig = plt.figure(figsize=fig_size(180, 165))
    outer = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1.4, 1],
                              hspace=0.45, wspace=0.35,
                              left=0.08, right=0.97, top=0.95, bottom=0.10)

    # --- Panel A: UpSet (top-left, occupies whole top-left quadrant) ---
    upset_gs = outer[0, 0].subgridspec(2, 2, height_ratios=[2, 1],
                                       width_ratios=[1, 3.5],
                                       wspace=0.05, hspace=0.05)
    ax_upper  = fig.add_subplot(upset_gs[0, 1])
    ax_lower  = fig.add_subplot(upset_gs[1, 0])
    ax_matrix = fig.add_subplot(upset_gs[1, 1])
    draw_upset(ax_upper, ax_lower, ax_matrix, membership)
    panel_letter(ax_upper, "A", x=-0.10, y=1.10)

    # --- Panel B: accumulation curve (top-right) ---
    ax_b = fig.add_subplot(outer[0, 1])
    draw_accumulation(ax_b, membership)
    panel_letter(ax_b, "B", x=-0.18, y=1.04)

    # --- Panel C: stacked composition (bottom-left) ---
    ax_c = fig.add_subplot(outer[1, 0])
    draw_composition(ax_c, partitions)
    panel_letter(ax_c, "C", x=-0.12, y=1.04)

    # --- Panel D: form-specific cloud (bottom-right) ---
    ax_d = fig.add_subplot(outer[1, 1])
    draw_form_specific(ax_d, partitions)
    panel_letter(ax_d, "D", x=-0.18, y=1.04)

    save_publication(fig, out_basename)
    plt.close(fig)
    print(f"Wrote {out_basename}.pdf and {out_basename}.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path(".").resolve())
    ap.add_argument("--out", type=str, default="results/figures/figure1_pangenome")
    args = ap.parse_args()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    build_figure(args.repo_root, args.out)


if __name__ == "__main__":
    main()
