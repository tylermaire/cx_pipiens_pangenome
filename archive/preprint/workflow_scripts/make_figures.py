#!/usr/bin/env python3
"""
make_figures.py — Generate publication figures from the clean pipeline run.

Outputs PNG (300 dpi) and PDF for each figure into the project's figures/
directory. All figures intentionally use journal-style typography: italic
species names, color-blind safe palettes, large readable fonts.

Run from the project root:
    python3 workflow/scripts/make_figures.py
"""
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results"
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)

# ---- Style ---------------------------------------------------------------
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# Color palette (Okabe–Ito, color-blind safe) for the four ingroup species
SPECIES_ORDER = ["Cx_molestus", "Cx_pallens", "Cx_pipiens", "Cx_quinquefasciatus"]
SPECIES_LABELS = {
    "Cx_molestus": "Cx. molestus",
    "Cx_pallens": "Cx. pallens",
    "Cx_pipiens": "Cx. pipiens",
    "Cx_quinquefasciatus": "Cx. quinquefasciatus",
    "Cx_tarsalis": "Cx. tarsalis",
}
SPECIES_COLORS = {
    "Cx_molestus": "#E69F00",
    "Cx_pallens": "#56B4E9",
    "Cx_pipiens": "#009E73",
    "Cx_quinquefasciatus": "#CC79A7",
    "Cx_tarsalis": "#999999",
}
COMP_COLORS = {"core": "#2C7FB8", "shell": "#7FCDBB", "cloud": "#EDF8B1"}


def save(fig, name):
    """Save both PNG and PDF copies of a figure."""
    fig.savefig(OUT / f"{name}.png")
    fig.savefig(OUT / f"{name}.pdf")
    print(f"  wrote {name}.png and {name}.pdf")


def italic(label):
    """Return matplotlib-italic-formatted species name."""
    spaced = label.replace(" ", r"\ ")
    return r"$\it{" + spaced + "}$"


# ==========================================================================
# Figure 1: Assembly quality (BUSCO + QUAST)
# ==========================================================================
def fig1_assembly_quality():
    print("Figure 1: Assembly quality")
    # BUSCO protein-mode summaries
    busco = []
    for s in SPECIES_ORDER + ["Cx_tarsalis"]:
        f = RES / "busco_proteins" / s / f"short_summary.specific.diptera_odb10.{s}.json"
        if not f.exists():
            continue
        with open(f) as fh:
            j = json.load(fh)
        r = j["results"]
        busco.append({
            "sample": s,
            "Single copy": r["Single copy percentage"],
            "Multi copy": r["Multi copy percentage"],
            "Fragmented": r["Fragmented percentage"],
            "Missing": r["Missing percentage"],
        })
    busco_df = pd.DataFrame(busco).set_index("sample")

    # QUAST
    quast = pd.read_csv(RES / "quast" / "transposed_report.tsv", sep="\t")
    quast = quast.set_index("Assembly")

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.6),
                              gridspec_kw={"width_ratios": [1.15, 1.05]})

    # --- Panel A: BUSCO stacked horizontal bar ----------------------------
    ax = axes[0]
    cats = ["Single copy", "Multi copy", "Fragmented", "Missing"]
    cat_colors = ["#1f78b4", "#a6cee3", "#ff7f00", "#e31a1c"]
    samples_with_tarsalis = SPECIES_ORDER + ["Cx_tarsalis"]
    busco_df = busco_df.reindex(samples_with_tarsalis)
    left = np.zeros(len(busco_df))
    for cat, col in zip(cats, cat_colors):
        vals = busco_df[cat].values
        ax.barh(range(len(busco_df)), vals, left=left, color=col, label=cat,
                edgecolor="white", linewidth=0.4)
        for i, v in enumerate(vals):
            if v >= 4:
                ax.text(left[i] + v/2, i, f"{v:.1f}", ha="center", va="center",
                        fontsize=8, color="white" if cat == "Single copy" else "black")
        left += vals
    ax.set_yticks(range(len(busco_df)))
    ax.set_yticklabels([italic(SPECIES_LABELS[s]) for s in busco_df.index])
    ax.invert_yaxis()
    ax.set_xlabel("BUSCO category (% of 3,285 diptera_odb10 markers)")
    ax.set_xlim(0, 100)
    ax.set_title("(a) Protein-mode BUSCO completeness")
    # Legend below the panel so it doesn't overlap data
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18),
              ncol=4, frameon=False, fontsize=9)

    # --- Panel B: Assembly statistics summary -----------------------------
    ax = axes[1]
    ax.axis("off")
    stats_rows = []
    for s in samples_with_tarsalis:
        if s not in quast.index:
            continue
        row = quast.loc[s]
        stats_rows.append([
            SPECIES_LABELS[s],
            f"{int(row['Total length']) / 1e6:.1f}",
            f"{int(row['# contigs'])}",
            f"{int(row['N50']) / 1e6:.1f}",
            f"{float(row['GC (%)']):.2f}",
        ])
    cols = ["Species", "Size (Mb)", "Contigs", "N50 (Mb)", "GC (%)"]
    col_widths = [0.42, 0.16, 0.14, 0.16, 0.12]
    tbl = ax.table(cellText=stats_rows, colLabels=cols, loc="center",
                   cellLoc="center", colColours=["#dddddd"] * len(cols),
                   colWidths=col_widths)
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9.5)
    tbl.scale(1, 1.6)
    for j, col in enumerate(cols):
        tbl[(0, j)].set_text_props(weight="bold")
    # Italicize species column
    for i in range(1, len(stats_rows) + 1):
        cell = tbl[(i, 0)]
        cell.set_text_props(style="italic")
    ax.set_title("(b) Assembly statistics")
    save(fig, "Figure_1_assembly_quality")
    plt.close(fig)


# ==========================================================================
# Figure 2: Pangenome composition
# ==========================================================================
def fig2_pangenome():
    print("Figure 2: Pangenome composition")
    summary = pd.read_csv(RES / "pangenome" / "pangenome_summary.tsv", sep="\t")
    parts = pd.read_csv(RES / "pangenome" / "partitioned_orthogroups.tsv", sep="\t")
    parts["comp_main"] = parts["compartment"].str.replace(
        r"_Cx_.*_specific$", "", regex=True)

    fig = plt.figure(figsize=(13, 4.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.15, 1.05], wspace=0.35)

    # --- Panel A: Donut of pangenome partition ----------------------------
    ax = fig.add_subplot(gs[0])
    main = summary.iloc[:3]  # core, shell, cloud
    sizes = main["n_orthogroups"].values
    pct = sizes / sizes.sum() * 100
    colors = [COMP_COLORS[c] for c in main["compartment"]]
    # Labels for the small slices go outside; large slice gets its label inside.
    def make_label(c, n, p):
        return f"{c.title()}\n{n:,}\n({p:.1f}%)"
    slice_labels = [make_label(c, n, p)
                    for c, n, p in zip(main["compartment"], sizes, pct)]
    wedges, texts = ax.pie(
        sizes, colors=colors, startangle=90,
        labels=slice_labels,
        labeldistance=1.18,
        wedgeprops={"edgecolor": "white", "linewidth": 2},
        textprops={"fontsize": 9.5, "ha": "center"})
    centre = plt.Circle((0, 0), 0.55, color="white")
    ax.add_artist(centre)
    ax.text(0, 0.08, f"{sizes.sum():,}", ha="center", va="center",
            fontsize=17, fontweight="bold")
    ax.text(0, -0.14, "orthogroups", ha="center", va="center", fontsize=10)
    ax.set_title("(a) Pangenome partition", pad=14)

    # --- Panel B: Per-species gene counts in each compartment -------------
    ax = fig.add_subplot(gs[1])
    rows = []
    for s in SPECIES_ORDER:
        for comp in ["core", "shell", "cloud"]:
            sub = parts[parts["comp_main"] == comp]
            rows.append({"species": s, "compartment": comp, "genes": int(sub[s].sum())})
    bar = pd.DataFrame(rows).pivot(index="species", columns="compartment",
                                    values="genes")
    bar = bar.reindex(SPECIES_ORDER)[["core", "shell", "cloud"]]
    bar.plot(kind="bar", stacked=True, ax=ax,
             color=[COMP_COLORS["core"], COMP_COLORS["shell"], COMP_COLORS["cloud"]],
             edgecolor="white", linewidth=0.5, legend=False)
    ax.set_xticklabels([italic(SPECIES_LABELS[s]) for s in bar.index],
                       rotation=20, ha="right")
    ax.set_ylabel("Genes per pangenome compartment")
    ax.set_xlabel("")
    # Annotate totals
    for i, s in enumerate(bar.index):
        total = bar.loc[s].sum()
        ax.text(i, total + 250, f"{total:,}", ha="center", fontsize=9)
    ax.legend(["Core", "Shell", "Cloud"], title="", loc="upper center",
              bbox_to_anchor=(0.5, -0.22), ncol=3, frameon=False, fontsize=9)
    ax.set_title("(b) Genes per species by compartment")
    ax.set_ylim(0, bar.sum(axis=1).max() * 1.18)

    # --- Panel C: Species-specific cloud orthogroups ----------------------
    ax = fig.add_subplot(gs[2])
    cloud_rows = summary[summary["compartment"].str.contains("_specific", na=False)]
    cloud_rows = cloud_rows.copy()
    cloud_rows["species"] = cloud_rows["compartment"].str.replace(
        r"cloud_(Cx_.*?)_specific", r"\1", regex=True)
    cloud_rows = cloud_rows.set_index("species").reindex(SPECIES_ORDER)
    colors = [SPECIES_COLORS[s] for s in cloud_rows.index]
    bars = ax.bar(range(len(cloud_rows)), cloud_rows["n_orthogroups"],
                  color=colors, edgecolor="white", linewidth=0.6)
    for bar_, n in zip(bars, cloud_rows["n_orthogroups"]):
        ax.text(bar_.get_x() + bar_.get_width() / 2, bar_.get_height() + 4,
                f"{int(n)}", ha="center", fontsize=10)
    ax.set_xticks(range(len(cloud_rows)))
    ax.set_xticklabels([italic(SPECIES_LABELS[s]) for s in cloud_rows.index],
                       rotation=20, ha="right")
    ax.set_ylabel("Species-specific cloud orthogroups")
    ax.set_ylim(0, cloud_rows["n_orthogroups"].max() * 1.25)
    ax.set_title("(c) Form-specific cloud orthogroups")

    save(fig, "Figure_2_pangenome")
    plt.close(fig)


# ==========================================================================
# Figure 3: Species tree with concordance factors
# ==========================================================================
def fig3_species_tree():
    print("Figure 3: Species tree")
    # Hand-render the 4-taxon tree from the Newick. With only 4 taxa we know
    # the layout exactly, so we draw it cleanly with matplotlib rather than
    # pulling in ete3 / dendropy.
    nwk = (RES / "phylo" / "concord.cf.tree").read_text().strip()
    # Newick: (Cx_molestus:0.0325,(Cx_pallens:0.0235,Cx_quinquefasciatus:0.0689)
    #          100/59.3/45.4:0.0141, Cx_pipiens:0.0339);
    m = re.search(r"\(Cx_molestus:([\d\.]+),"
                  r"\(Cx_pallens:([\d\.]+),Cx_quinquefasciatus:([\d\.]+)\)"
                  r"([\d\./]+):([\d\.]+),"
                  r"Cx_pipiens:([\d\.]+)\);", nwk)
    if not m:
        raise RuntimeError(f"Could not parse Newick: {nwk}")
    bl_mol, bl_pal, bl_qui, support, bl_inner, bl_pip = m.groups()
    bl_mol, bl_pal, bl_qui, bl_inner, bl_pip = (
        float(bl_mol), float(bl_pal), float(bl_qui),
        float(bl_inner), float(bl_pip))

    fig, ax = plt.subplots(figsize=(9.5, 6))
    # Layout: unrooted as 3-way star with an internal node for (pallens, quinque)
    # Place taxa around a center.
    # We'll draw it as a rooted-looking tree with Cx_molestus at top, the
    # (pallens, quinque) clade at right, Cx_pipiens at bottom — but make
    # clear via the figure caption that the rooting is arbitrary (no usable
    # outgroup in current run).
    center = (0, 0)
    # Tip vertical positions.
    y_mol = 1.4
    y_pal = 0.55
    y_qui = -0.55
    y_pip = -1.4
    # Scale branch lengths
    scale = 30
    x_inner = bl_inner * scale  # internal node x for (pallens, quinque)

    # Internal node between root and (pallens, quinque)
    inner_xy = (x_inner, (y_pal + y_qui) / 2)

    # Branches from root
    root_x = 0
    # molestus
    ax.plot([root_x, bl_mol * scale], [y_mol, y_mol], "k-", linewidth=2.5)
    ax.plot([root_x, root_x], [y_mol, inner_xy[1]], "k-", linewidth=2.5)
    # (pallens, quinque) inner branch
    ax.plot([root_x, inner_xy[0]], [inner_xy[1], inner_xy[1]],
            "k-", linewidth=2.5)
    # pipiens
    ax.plot([root_x, root_x], [inner_xy[1], y_pip], "k-", linewidth=2.5)
    ax.plot([root_x, bl_pip * scale], [y_pip, y_pip], "k-", linewidth=2.5)
    # From inner node to pallens / quinque
    ax.plot([inner_xy[0], inner_xy[0]], [y_pal, y_qui], "k-", linewidth=2.5)
    ax.plot([inner_xy[0], inner_xy[0] + bl_pal * scale], [y_pal, y_pal],
            "k-", linewidth=2.5)
    ax.plot([inner_xy[0], inner_xy[0] + bl_qui * scale], [y_qui, y_qui],
            "k-", linewidth=2.5)

    # Tip labels
    tip_kw = dict(va="center", fontsize=12, fontstyle="italic")
    ax.text(bl_mol * scale + 0.1, y_mol, "Cx. molestus", **tip_kw)
    ax.text(inner_xy[0] + bl_pal * scale + 0.1, y_pal, "Cx. pallens", **tip_kw)
    ax.text(inner_xy[0] + bl_qui * scale + 0.1, y_qui, "Cx. quinquefasciatus",
            **tip_kw)
    ax.text(bl_pip * scale + 0.1, y_pip, "Cx. pipiens", **tip_kw)

    # Support label as plain bold text floating above the inner branch.
    # No box — the box obscured the branch joints. This is the standard
    # convention used in published phylogenies.
    ax.text(x_inner / 2, inner_xy[1] + 0.12, support,
            ha="center", va="bottom", fontsize=10.5,
            fontweight="bold", color="#222222")

    # Scale bar (pushed below the lowest tip with clearance)
    sb_x0, sb_x1, sb_y = -0.05, 0.05 * scale, -1.80
    ax.plot([sb_x0, sb_x0 + 0.05 * scale], [sb_y, sb_y], "k-", linewidth=2)
    ax.text((sb_x0 + sb_x0 + 0.05 * scale) / 2, sb_y - 0.10, "0.05",
            ha="center", va="top", fontsize=10)

    # Footnote about rooting
    ax.text(0, -2.15,
            "Branch labels: ultrafast bootstrap / gene concordance / site concordance.\n"
            "Tree rooting is illustrative; Cx. tarsalis outgroup not available in "
            "the current run (see paper_audit.md).",
            ha="center", va="top", fontsize=8.5, color="#444444", style="italic")

    ax.set_xlim(-0.5, 3.8)
    ax.set_ylim(-2.4, 1.7)
    ax.axis("off")
    save(fig, "Figure_3_species_tree")
    plt.close(fig)


# ==========================================================================
# Figure 4: CAFE gene family expansions and contractions per lineage
# ==========================================================================
def fig4_cafe():
    print("Figure 4: CAFE gene family evolution")
    clade = (RES / "cafe" / "output" / "Gamma_clade_results.txt").read_text()
    # Parse the simple tab-separated text
    rows = []
    for line in clade.strip().splitlines()[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        name = re.sub(r"<\d+>$", "", parts[0])
        rows.append({"sample": name, "expansion": int(parts[1]),
                     "contraction": int(parts[2])})
    df = pd.DataFrame(rows).set_index("sample").reindex(SPECIES_ORDER)

    # Lambda / alpha
    gamma = (RES / "cafe" / "output" / "Gamma_results.txt").read_text()
    lam = float(re.search(r"Lambda:\s*([\d\.]+)", gamma).group(1))
    alpha_m = re.search(r"Alpha:\s*([\d\.]+)", gamma)
    alpha = float(alpha_m.group(1)) if alpha_m else float("nan")
    n_sig = sum(1 for _ in open(RES / "cafe" / "significant_families.tsv")) - 1

    fig, ax = plt.subplots(figsize=(9, 4.8))
    x = np.arange(len(df))
    w = 0.38
    bars_e = ax.bar(x - w / 2, df["expansion"], w, color="#2ca25f",
                    edgecolor="white", linewidth=0.5, label="Expansion")
    bars_c = ax.bar(x + w / 2, df["contraction"], w, color="#cb181d",
                    edgecolor="white", linewidth=0.5, label="Contraction")
    for bar in bars_e:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 50,
                f"{int(bar.get_height()):,}", ha="center", fontsize=10)
    for bar in bars_c:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 50,
                f"{int(bar.get_height()):,}", ha="center", fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels([italic(SPECIES_LABELS[s]) for s in df.index],
                       rotation=15, ha="right")
    ax.set_ylabel("Gene family expansions / contractions")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18),
              ncol=2, frameon=False, fontsize=10)
    ax.set_ylim(0, df.values.max() * 1.18)
    save(fig, "Figure_4_cafe")
    plt.close(fig)


# ==========================================================================
# Figure 5: Key gene families heatmap
# ==========================================================================
def fig5_key_families():
    print("Figure 5: Key gene families")
    df = pd.read_csv(RES / "functional" / "key_families_wide.tsv", sep="\t")
    df = df.set_index("family")[SPECIES_ORDER]
    # Sort families by mean copy number descending
    df = df.loc[df.mean(axis=1).sort_values(ascending=False).index]

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    sns.heatmap(df, annot=True, fmt="d", cmap="YlGnBu", cbar_kws={"label": "Copy number"},
                linewidths=0.4, linecolor="white", ax=ax,
                xticklabels=[italic(SPECIES_LABELS[s]) for s in df.columns],
                annot_kws={"fontsize": 11})
    ax.set_xlabel("")
    ax.set_ylabel("Gene family")
    plt.xticks(rotation=20, ha="right")
    plt.yticks(rotation=0)
    save(fig, "Figure_5_key_families")
    plt.close(fig)


# ==========================================================================
# Figure 6: TE composition per genome
# ==========================================================================
def fig6_te_composition():
    print("Figure 6: TE composition")
    rows = []
    for s in SPECIES_ORDER:
        tbl = RES / "repeats" / s / f"{s}.fasta.tbl"
        if not tbl.exists():
            continue
        text = tbl.read_text()
        total = float(re.search(r"total length:\s+(\d+)\s+bp", text).group(1))
        interspersed = float(re.search(
            r"Total interspersed repeats:\s+(\d+)\s+bp", text).group(1))
        simple = float(re.search(r"Simple repeats:\s+\d+\s+(\d+)\s+bp", text).group(1))
        low = float(re.search(r"Low complexity:\s+\d+\s+(\d+)\s+bp", text).group(1))
        unmasked = total - interspersed - simple - low
        rows.append({
            "sample": s,
            "Interspersed TE": interspersed / total * 100,
            "Simple repeats": simple / total * 100,
            "Low complexity": low / total * 100,
            "Unmasked": unmasked / total * 100,
        })
    df = pd.DataFrame(rows).set_index("sample").reindex(SPECIES_ORDER)

    fig, ax = plt.subplots(figsize=(9, 4.6))
    cats = ["Interspersed TE", "Simple repeats", "Low complexity", "Unmasked"]
    colors = ["#6a3d9a", "#fdbf6f", "#fb9a99", "#cccccc"]
    left = np.zeros(len(df))
    for cat, col in zip(cats, colors):
        ax.barh(range(len(df)), df[cat], left=left, color=col, label=cat,
                edgecolor="white", linewidth=0.5)
        for i, v in enumerate(df[cat]):
            if v > 4:
                ax.text(left[i] + v / 2, i, f"{v:.1f}%",
                        ha="center", va="center", fontsize=9,
                        color="white" if cat == "Interspersed TE" else "black")
        left += df[cat].values
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels([italic(SPECIES_LABELS[s]) for s in df.index])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Genome fraction (%)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18),
              ncol=4, frameon=False, fontsize=9)
    save(fig, "Figure_6_te_composition")
    plt.close(fig)


if __name__ == "__main__":
    sns.set_style("white")
    fig1_assembly_quality()
    fig2_pangenome()
    fig3_species_tree()
    fig4_cafe()
    fig5_key_families()
    fig6_te_composition()
    print("All figures written to:", OUT)
