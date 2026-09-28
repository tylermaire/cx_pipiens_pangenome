#!/usr/bin/env python3
"""
make_figure9_anvio.py — Build Figure 9 (anvi'o pangenome composition) from the
anvi-summarize gene-clusters TSV.

Reads anvio_data/Culex_pipiens_complex_gene_clusters_summary.txt and writes
figures/Figure_9_anvio_pangenome.{png,pdf}.

The figure has three panels:
  (a) UpSet-style intersection plot — the 15 non-empty subsets of {4 genomes}
      sorted by intersection size, with intersection counts above each bar and
      the genome-membership matrix below.
  (b) Per-genome cluster count (one bar per genome).
  (c) Cluster size distribution (log-scaled).
"""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
IN = ROOT / "anvio_data" / "Culex_pipiens_complex_gene_clusters_summary.txt"
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)

SPECIES_ORDER = ["Cx_molestus", "Cx_pallens", "Cx_pipiens", "Cx_quinquefasciatus"]
SPECIES_LABELS = {
    "Cx_molestus": "Cx. molestus",
    "Cx_pallens": "Cx. pallens",
    "Cx_pipiens": "Cx. pipiens",
    "Cx_quinquefasciatus": "Cx. quinquefasciatus",
}
SPECIES_COLORS = {
    "Cx_molestus": "#E69F00",
    "Cx_pallens": "#56B4E9",
    "Cx_pipiens": "#009E73",
    "Cx_quinquefasciatus": "#CC79A7",
}

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.dpi": 300, "savefig.bbox": "tight",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def italic(s):
    return r"$\it{" + s.replace(" ", r"\ ") + "}$"


def main():
    print("Loading anvi'o gene clusters summary...")
    df = pd.read_csv(IN, sep="\t", usecols=["gene_cluster_id", "genome_name"])
    print(f"  {len(df):,} rows ({df['gene_cluster_id'].nunique():,} unique gene clusters)")

    # Presence/absence matrix: clusters x genomes
    pa = (df.groupby(["gene_cluster_id", "genome_name"]).size().unstack(fill_value=0) > 0).astype(int)
    pa = pa.reindex(columns=SPECIES_ORDER, fill_value=0)
    total_clusters = len(pa)
    print(f"  {total_clusters:,} total gene clusters across {pa.shape[1]} genomes")

    # ---- All 2^4 - 1 = 15 non-empty intersection patterns ----
    pattern_keys = []
    for cluster in pa.index:
        key = tuple(pa.loc[cluster].values)  # e.g. (1,1,1,1)
        pattern_keys.append(key)
    pa["pattern"] = pattern_keys
    intersections = pa["pattern"].value_counts().sort_values(ascending=False)
    print(f"  {len(intersections)} distinct intersection patterns:")
    for pat, n in intersections.items():
        members = [SPECIES_ORDER[i] for i, v in enumerate(pat) if v == 1]
        print(f"    {pat}  n={n:>6,}  ({', '.join(members)})")

    # ---- Cluster-size distribution (number of genes in each cluster) ----
    sizes = df.groupby("gene_cluster_id").size()

    # ---- Per-genome cluster count ----
    per_genome = pa[SPECIES_ORDER].sum(axis=0)

    # =====================================================================
    # Figure
    # =====================================================================
    fig = plt.figure(figsize=(14, 8.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.3, 1.0], height_ratios=[3.5, 1.6],
                          hspace=0.55, wspace=0.40)

    # Panel A (top-left): UpSet intersection bars
    ax_bars = fig.add_subplot(gs[0, 0])
    n_patterns = len(intersections)
    x = np.arange(n_patterns)
    counts = intersections.values

    bars = ax_bars.bar(x, counts, color="#444444", edgecolor="white", linewidth=0.5)
    for i, c in enumerate(counts):
        ax_bars.text(i, c + total_clusters * 0.012, f"{c:,}",
                     ha="center", va="bottom", fontsize=9, rotation=0)
    ax_bars.set_ylim(0, counts.max() * 1.15)
    ax_bars.set_ylabel("Gene clusters in intersection", fontsize=11)
    ax_bars.set_xticks([])
    ax_bars.set_title("(a) Cluster intersections", fontsize=11, fontweight="bold",
                      loc="left", pad=8)
    ax_bars.spines["bottom"].set_visible(False)

    # Panel A (matrix below): genome-membership dots/lines
    ax_mat = fig.add_subplot(gs[1, 0], sharex=ax_bars)
    n_genomes = len(SPECIES_ORDER)
    y = np.arange(n_genomes)
    ax_mat.set_xlim(-0.5, n_patterns - 0.5)
    ax_mat.set_ylim(-0.6, n_genomes - 0.4)
    ax_mat.invert_yaxis()

    for i, (pat, _) in enumerate(intersections.items()):
        members = [j for j, v in enumerate(pat) if v == 1]
        # Vertical line linking included genomes
        if len(members) >= 2:
            ax_mat.plot([i, i], [min(members), max(members)], color="#444444", lw=1.5, zorder=1)
        # Filled dots for included, light dots for excluded
        for j in range(n_genomes):
            if j in members:
                ax_mat.scatter(i, j, s=85, color=SPECIES_COLORS[SPECIES_ORDER[j]],
                               edgecolor="#222", lw=0.6, zorder=2)
            else:
                ax_mat.scatter(i, j, s=85, color="#dddddd",
                               edgecolor="#888", lw=0.4, zorder=2)

    ax_mat.set_yticks(range(n_genomes))
    ax_mat.set_yticklabels([italic(SPECIES_LABELS[s]) for s in SPECIES_ORDER], fontsize=10)
    ax_mat.set_xticks([])
    ax_mat.set_xlabel("Intersection (sorted by size, left = largest)", fontsize=10.5)
    for spine_name in ["top", "right", "bottom", "left"]:
        ax_mat.spines[spine_name].set_visible(False)
    ax_mat.tick_params(left=False)

    # Panel B (top-right): per-genome cluster counts
    ax_pg = fig.add_subplot(gs[0, 1])
    bars = ax_pg.barh(
        range(n_genomes),
        per_genome.values,
        color=[SPECIES_COLORS[s] for s in per_genome.index],
        edgecolor="white", lw=0.6,
    )
    for i, v in enumerate(per_genome.values):
        ax_pg.text(v + total_clusters * 0.005, i, f"{v:,}",
                   va="center", fontsize=9.5)
    ax_pg.set_yticks(range(n_genomes))
    ax_pg.set_yticklabels([italic(SPECIES_LABELS[s]) for s in per_genome.index])
    ax_pg.invert_yaxis()
    ax_pg.set_xlim(0, per_genome.max() * 1.18)
    ax_pg.set_xlabel("Gene clusters per genome", fontsize=10.5)
    ax_pg.set_title("(b) Per-genome cluster count", fontsize=11, fontweight="bold",
                    loc="left", pad=8)

    # Panel C (bottom-right): cluster-size distribution
    ax_sz = fig.add_subplot(gs[1, 1])
    bins = np.logspace(0, np.log10(sizes.max() + 1), 30)
    ax_sz.hist(sizes.values, bins=bins, color="#666666", edgecolor="white", lw=0.4)
    ax_sz.set_xscale("log")
    ax_sz.set_xlabel("Genes per cluster", fontsize=10.5)
    ax_sz.set_ylabel("# clusters", fontsize=10)
    ax_sz.set_title("(c) Cluster-size distribution", fontsize=11, fontweight="bold",
                    loc="left", pad=8)
    ax_sz.text(0.97, 0.94,
               f"median {int(np.median(sizes))}  |  max {int(sizes.max())}",
               transform=ax_sz.transAxes, ha="right", va="top", fontsize=9, color="#444")

    fig.savefig(OUT / "Figure_9_anvio_pangenome.png")
    fig.savefig(OUT / "Figure_9_anvio_pangenome.pdf")
    print("  wrote Figure_9_anvio_pangenome.{png,pdf}")
    plt.close(fig)


if __name__ == "__main__":
    main()
