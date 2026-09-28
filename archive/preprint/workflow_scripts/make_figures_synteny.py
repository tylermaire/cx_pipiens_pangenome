#!/usr/bin/env python3
"""
make_figures_synteny.py — Figure 7 (ANI heatmap) and Figure 8 (synteny composite).

Reads from results/synteny/ and writes Figure_7_ani.{png,pdf} and
Figure_8_synteny.{png,pdf} to figures/.
"""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results"
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)

SPECIES = ["Cx_molestus", "Cx_pallens", "Cx_pipiens", "Cx_quinquefasciatus"]
SPECIES_LABEL = {
    "Cx_molestus": "Cx. molestus",
    "Cx_pallens": "Cx. pallens",
    "Cx_pipiens": "Cx. pipiens",
    "Cx_quinquefasciatus": "Cx. quinquefasciatus",
    "Cx_tarsalis": "Cx. tarsalis",
}

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.dpi": 300, "savefig.bbox": "tight",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def italic(s):
    return r"$\it{" + s.replace(" ", r"\ ") + "}$"


def fig7_ani_heatmap():
    print("Figure 7: ANI heatmap")
    df = pd.read_csv(RES / "synteny" / "ani_matrix.tsv", sep="\t", index_col=0)
    # Restrict to the four ingroup species (tarsalis row/col is NA)
    df = df.reindex(SPECIES)[SPECIES]
    df = df.astype(float)

    fig, ax = plt.subplots(figsize=(6.5, 5.2))
    mask = np.eye(len(df), dtype=bool)
    sns.heatmap(
        df, annot=True, fmt=".2f", cmap="YlGnBu_r",
        vmin=92.5, vmax=95.0, cbar_kws={"label": "skani ANI (%)"},
        linewidths=0.6, linecolor="white",
        annot_kws={"size": 11, "weight": "bold"},
        mask=mask, ax=ax,
    )
    # Diagonal cells filled with neutral gray and 100 labels
    for i in range(len(df)):
        ax.add_patch(plt.Rectangle((i, i), 1, 1, fc="#e0e0e0", ec="white", lw=0.6))
        ax.text(i + 0.5, i + 0.5, "100", ha="center", va="center",
                fontsize=11, fontweight="bold", color="#444444")

    labels = [italic(SPECIES_LABEL[s]) for s in df.index]
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_yticklabels(labels, rotation=0)
    ax.set_xlabel("")
    ax.set_ylabel("")

    fig.savefig(OUT / "Figure_7_ani.png")
    fig.savefig(OUT / "Figure_7_ani.pdf")
    print("  wrote Figure_7_ani.png and Figure_7_ani.pdf")
    plt.close(fig)


def fig8_synteny_composite():
    print("Figure 8: Synteny composite")
    plots_dir = RES / "synteny" / "synteny_dotplots"
    pairs = [
        ("Cx_pallens", "Cx_pipiens", "(a) pallens vs pipiens"),
        ("Cx_pallens", "Cx_quinquefasciatus", "(b) pallens vs quinquefasciatus"),
        ("Cx_molestus", "Cx_pallens", "(c) molestus vs pallens"),
        ("Cx_molestus", "Cx_pipiens", "(d) molestus vs pipiens"),
        ("Cx_pipiens", "Cx_quinquefasciatus", "(e) pipiens vs quinquefasciatus"),
        ("Cx_molestus", "Cx_quinquefasciatus", "(f) molestus vs quinquefasciatus"),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 9))
    for ax, (s1, s2, label) in zip(axes.flat, pairs):
        # Try both orderings of the filename
        candidates = [
            plots_dir / f"{s1}_vs_{s2}.png",
            plots_dir / f"{s2}_vs_{s1}.png",
        ]
        path = next((p for p in candidates if p.exists()), None)
        if path is None:
            ax.text(0.5, 0.5, f"missing:\n{s1}_vs_{s2}", ha="center", va="center")
            ax.set_xticks([]); ax.set_yticks([])
            continue
        img = mpimg.imread(path)
        ax.imshow(img)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(label, fontsize=10.5, pad=6)
        for spine in ax.spines.values():
            spine.set_visible(False)

    fig.tight_layout()
    fig.savefig(OUT / "Figure_8_synteny.png")
    fig.savefig(OUT / "Figure_8_synteny.pdf")
    print("  wrote Figure_8_synteny.png and Figure_8_synteny.pdf")
    plt.close(fig)


if __name__ == "__main__":
    sns.set_style("white")
    fig7_ani_heatmap()
    fig8_synteny_composite()
