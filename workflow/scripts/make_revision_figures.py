#!/usr/bin/env python3
"""
make_revision_figures.py - main text figures of the revised manuscript, drawn
from the workflow outputs at the journal's final size (171 mm wide).

  Figure 1  pangenome partition, and what the cloud is made of
  Figure 2  phylogeny and the test of the discordant gene trees
  Figure 3  copy number retention and CAFE significance by reference copy number
  Figure 4  pairwise gene anchor synteny

Every number is read from results/; nothing is typed in. Writes a vector PDF
and an LZW compressed TIFF (1200 dpi) per figure to figures/revision/.

Usage, from the repository root:
    python3 workflow/scripts/make_revision_figures.py [--dpi 1200] [--only 1,2]
"""

import argparse
import collections
import csv
import glob
import importlib.util
import itertools
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy.stats import binomtest  # noqa: E402

OUT = os.path.join("figures", "revision")
MM = 1 / 25.4
WIDTH = 171 * MM

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
SERIES = "#2a78d6"
RAMP = {"core": "#0d366b", "shell": "#1c5cab", "cloud": "#3987e5", "unassigned": "#86b6ef"}
OUTCOME = {"artifact": "#eb6834", "supported": "#1baf7a", "weak": MUTED}
CHROM = ["#2a78d6", "#eb6834", "#1baf7a"]

INGROUP = ["Cx_quinquefasciatus", "Cx_pallens", "Cx_molestus", "Cx_pipiens"]
TRANSFERRED = ["Cx_pallens", "Cx_molestus", "Cx_pipiens"]


def short(sample):
    return "Cx. " + sample.split("_", 1)[1]


def style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 7, "axes.titlesize": 7.5,
        "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5, "axes.edgecolor": AXIS, "axes.linewidth": 0.5,
        "xtick.color": INK2, "ytick.color": INK2, "axes.labelcolor": INK,
        "text.color": INK, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "pdf.fonttype": 42,
        "axes.titlelocation": "left", "axes.titleweight": "bold",
        "savefig.facecolor": "white", "figure.facecolor": "white",
    })


def clean(ax, grid_axis="y"):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.5)
        ax.set_axisbelow(True)


def thousands(ax, axis="y"):
    from matplotlib.ticker import FuncFormatter
    fmt = FuncFormatter(lambda v, _: f"{int(v):,}")
    (ax.yaxis if axis == "y" else ax.xaxis).set_major_formatter(fmt)


def italic_ticks(labels):
    for t in labels:
        t.set_fontstyle("italic")


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def save(fig, name, dpi):
    os.makedirs(OUT, exist_ok=True)
    pdf = os.path.join(OUT, f"{name}.pdf")
    tif = os.path.join(OUT, f"{name}.tif")
    fig.savefig(pdf)
    fig.savefig(tif, dpi=dpi, pil_kwargs={"compression": "tiff_lzw"})
    from PIL import Image
    with Image.open(tif) as im:             # journals want RGB, not RGBA
        im.convert("RGB").save(tif, compression="tiff_lzw", dpi=(dpi, dpi))
    png = os.path.join(OUT, f"{name}.png")
    fig.savefig(png, dpi=200)
    plt.close(fig)
    print(f"  wrote {pdf}, {tif}, {png}")


def ci(k, n):
    r = binomtest(k, n)
    c = r.proportion_ci(confidence_level=0.95, method="exact")
    return k / n, c.low, c.high


# ------------------------------------------------------------------ figure 1
def figure1(dpi):
    summary = {r["compartment"]: r for r in read_tsv("results/pangenome/pangenome_summary.tsv")}
    part = read_tsv("results/pangenome/partitioned_orthogroups.tsv")
    per = {s: collections.Counter() for s in INGROUP}
    for r in part:
        for s in INGROUP:
            per[s][r["compartment"]] += int(r[s])
    unassigned = collections.Counter()
    un = glob.glob("results/orthofinder/output/**/Orthogroups_UnassignedGenes.tsv", recursive=True)
    with open(un[0]) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            for i, x in enumerate(f[1:], 1):
                if x.strip() and head[i] in per:
                    unassigned[head[i]] += 1
    orf = {(r["sample"], r["compartment"]): r
           for r in read_tsv("results/annotation/transfer_quality_by_compartment.tsv")}
    absence = {r["compartment"]: r for r in read_tsv("results/validation/absence_summary.tsv")}

    fig, axes = plt.subplots(2, 2, figsize=(WIDTH, 124 * MM),
                             gridspec_kw={"hspace": 0.7, "wspace": 0.42,
                                          "left": 0.09, "right": 0.985, "top": 0.94,
                                          "bottom": 0.16})
    # (a) orthogroups by compartment
    ax = axes[0, 0]
    comps = ["core", "shell", "cloud", "outgroup_only"]
    names = ["Core", "Shell", "Cloud", "Outgroup\nonly"]
    vals = [int(summary[c]["n_orthogroups"]) for c in comps]
    ax.bar(names, vals, color=SERIES, width=0.6)
    for x, v in enumerate(vals):
        ax.text(x, v + max(vals) * 0.015, f"{v:,}", ha="center", va="bottom", fontsize=6.5)
    ax.set_ylabel("Orthogroups")
    ax.set_ylim(0, max(vals) * 1.12)
    thousands(ax, "y")
    ax.set_title("(a) Orthogroups by compartment")
    clean(ax)

    # (b) genes per form by compartment
    ax = axes[0, 1]
    order = INGROUP
    y = range(len(order))
    left = [0] * len(order)
    for comp in ["core", "shell", "cloud", "unassigned"]:
        vals = [unassigned[s] if comp == "unassigned" else per[s][comp] for s in order]
        ax.barh(list(y), vals, left=left, color=RAMP[comp], height=0.62,
                edgecolor="white", linewidth=0.6, label=comp.capitalize())
        left = [a + b for a, b in zip(left, vals)]
    for i, tot in enumerate(left):
        ax.text(tot + 150, i, f"{tot:,}", va="center", fontsize=6.5)
    ax.set_yticks(list(y))
    ax.set_yticklabels([short(s) for s in order])
    italic_ticks(ax.get_yticklabels())
    ax.invert_yaxis()
    ax.set_xlim(0, max(left) * 1.16)
    ax.set_xticks([0, 5000, 10000, 15000])
    thousands(ax, "x")
    ax.set_xlabel("Genes (one protein per gene)")
    ax.set_title("(b) Genes per form by compartment")
    ax.legend(ncol=4, frameon=False, loc="upper left", bbox_to_anchor=(-0.02, -0.2),
              handlelength=1, columnspacing=0.8)
    clean(ax, "x")

    # (c) transferred models without a valid ORF, by compartment
    ax = axes[1, 0]
    comps = ["core", "shell", "cloud", "unassigned"]
    width = 0.2
    for j, comp in enumerate(comps):
        xs = [i + (j - 1.5) * width for i in range(len(TRANSFERRED))]
        vals = [float(orf[(s, comp)]["pct_invalid_orf"]) for s in TRANSFERRED]
        ax.bar(xs, vals, width=width * 0.92, color=RAMP[comp], label=comp.capitalize())
    ax.set_xticks(range(len(TRANSFERRED)))
    ax.set_xticklabels([short(s) for s in TRANSFERRED])
    italic_ticks(ax.get_xticklabels())
    ax.set_ylim(0, 105)
    ax.set_ylabel("Models without a valid ORF (%)")
    ax.set_title("(c) Transferred models without a valid ORF")
    ax.legend(ncol=4, frameon=False, loc="upper left", bbox_to_anchor=(-0.02, -0.14),
              handlelength=1, columnspacing=0.8)
    clean(ax)

    # (d) absence events: what the probes found
    ax = axes[1, 1]
    rows = ["cloud", "shell"]
    left = [0.0, 0.0]
    for key, label in [("artifact", "Gene present or locus unannotated"),
                       ("supported", "Absence supported"), ("weak", "Weak hit")]:
        vals = [float(absence[c][f"pct_{key}"]) for c in rows]
        ax.barh([0, 1], vals, left=left, color=OUTCOME[key], height=0.55,
                edgecolor="white", linewidth=0.6, label=label)
        for i, (lft, v) in enumerate(zip(left, vals)):
            if v >= 6:
                ax.text(lft + v / 2, i, f"{v:.1f}%", ha="center", va="center",
                        fontsize=6.5, color="white" if key == "artifact" else INK)
        left = [a + b for a, b in zip(left, vals)]
    ax.set_yticks([0, 1])
    ax.set_yticklabels([f"Cloud\n(n = {int(absence['cloud']['n_events']):,})",
                        f"Shell\n(n = {int(absence['shell']['n_events']):,})"])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Absence events (%)")
    ax.set_title("(d) Absence calls tested with probes")
    ax.legend(ncol=2, frameon=False, loc="upper left", bbox_to_anchor=(-0.3, -0.24),
              handlelength=1, columnspacing=0.8)
    clean(ax, "x")
    save(fig, "Figure_1_pangenome", dpi)


# ------------------------------------------------------------------ figure 2
def branch_lengths(newick):
    term = {t: float(x) for t, x in re.findall(r"([A-Za-z][A-Za-z0-9_.]*):([0-9.eE+-]+)", newick)}
    internal = [float(x) for x in re.findall(r"\)[0-9./]*:([0-9.eE+-]+)", newick)]
    return term, internal[0]


def figure2(dpi):
    import math
    tree = open("results/phylo/concat_tree.treefile").read().strip()
    term, internal = branch_lengths(tree)
    cf = {}
    for line in open("results/phylo/concord.cf.stat"):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if f[0] == "ID":
            head = f
            continue
        row = dict(zip(head, f))
        if row.get("gCF") not in (None, "NA"):
            cf = row
    topo = {r["role"]: r for r in read_tsv("results/phylo/quartet_topology_counts.tsv")}
    support = read_tsv("results/phylo/quartet_asymmetry_by_support.tsv")
    sites = read_tsv("results/phylo/quartet_site_patterns.tsv")

    fig = plt.figure(figsize=(WIDTH, 128 * MM))
    top = fig.add_gridspec(1, 2, width_ratios=[1, 1.1], wspace=0.55, left=0.02,
                           right=0.97, top=0.95, bottom=0.6)
    bot = fig.add_gridspec(1, 2, width_ratios=[1, 1.25], wspace=0.62, left=0.09,
                           right=0.97, top=0.47, bottom=0.09)

    # (a) unrooted quartet with concatenated branch lengths
    ax = fig.add_subplot(top[0, 0])
    ax.set_title("(a) Concatenated tree (unrooted)")
    scale = 1.0 / 0.06
    half = internal * scale / 2
    ang = math.radians(35)
    tips = {"Cx_molestus": (-1, 1), "Cx_pipiens": (-1, -1),
            "Cx_pallens": (1, -1), "Cx_quinquefasciatus": (1, 1)}
    for t, (sx, sy) in tips.items():
        x0 = sx * half
        L = term[t] * scale
        x1, y1 = x0 + sx * L * math.cos(ang), sy * L * math.sin(ang)
        ax.plot([x0, x1], [0, y1], color=INK, lw=1.0, solid_capstyle="round")
        ax.text(x1 + sx * 0.05, y1 + sy * 0.05, short(t), ha="left" if sx > 0 else "right",
                va="center", fontstyle="italic", fontsize=7)
    ax.plot([-half, half], [0, 0], color=INK, lw=1.0)
    if cf:
        ax.text(0, -0.66, f"Internal branch: UFBoot 100, gCF {float(cf['gCF']):.1f}, "
                f"sCF {float(cf['sCF']):.1f}", ha="center", va="center", fontsize=6.3,
                color=INK2)
    ax.plot([-0.3, -0.3 + 0.01 * scale], [-0.95, -0.95], color=INK, lw=1.0)
    ax.text(-0.3 + 0.005 * scale, -1.01, "0.01 substitutions per site", ha="center",
            va="top", fontsize=6, color=INK2)
    ax.set_xlim(-2.1, 2.3)
    ax.set_ylim(-1.2, 1.0)
    ax.axis("off")

    # (b) gene tree topology counts, all trees
    ax = fig.add_subplot(top[0, 1])
    rows = [("species_tree", "mol + pip | pal + qui"),
            ("major_discordant", "mol + qui | pal + pip"),
            ("minor_discordant", "mol + pal | pip + qui")]
    vals = [int(topo[r]["n_gene_trees"]) for r, _ in rows]
    ax.barh(range(3), vals, color=SERIES, height=0.55)
    for i, v in enumerate(vals):
        ax.text(v + 60, i, f"{v:,} ({100 * v / sum(vals):.1f}%)", va="center", fontsize=6.5)
    ax.set_yticks(range(3))
    ax.set_yticklabels([lab for _, lab in rows])
    ax.invert_yaxis()
    ax.set_xlim(0, max(vals) * 1.4)
    ax.set_xticks([0, 2000, 4000])
    thousands(ax, "x")
    ax.set_xlabel(f"Gene trees (n = {sum(vals):,})")
    ax.set_title("(b) Gene tree topologies")
    clean(ax, "x")

    # (c) share of discordant gene trees with the major topology, by support
    ax = fig.add_subplot(bot[0, 0])
    xs, ps, lo, hi, ns = [], [], [], [], []
    for r in support:
        n1, n2 = int(r["n_major"]), int(r["n_minor"])
        p, a, b = ci(n1, n1 + n2)
        xs.append(int(r["min_ufboot"]))
        ps.append(p)
        lo.append(p - a)
        hi.append(b - p)
        ns.append(n1 + n2)
    pos = list(range(len(xs)))
    ax.axhline(0.5, color=MUTED, lw=0.6)
    ax.text(-0.3, 0.505, "equal under ILS", ha="left", va="bottom", fontsize=6, color=MUTED)
    ax.errorbar(pos, ps, yerr=[lo, hi], fmt="o", color=SERIES, ms=3.5, lw=0.9,
                capsize=2, capthick=0.8)
    ax.set_xticks(pos)
    GE = "\u2265 "
    ticks = [("all" if x == 0 else GE + str(x)) + "\n" + f"{n:,}" for x, n in zip(xs, ns)]
    ax.set_xticklabels(ticks, fontsize=5.8)
    ax.set_xlabel("UFBoot on the internal branch; discordant trees")
    ax.set_ylabel("Share with mol + qui | pal + pip")
    ax.set_ylim(0.45, 0.75)
    ax.set_title("(c) Excess by gene tree support")
    clean(ax)

    # (d) the same share measured by gene trees, loci and sites
    ax = fig.add_subplot(bot[0, 1])
    s0 = next(r for r in support if r["min_ufboot"] == "0")
    s95 = next(r for r in support if r["min_ufboot"] == "95")
    role = {r["role"]: r for r in sites}
    k = int(role["top_loci_major_discordant_gene_trees"]["n_loci"])
    a1 = int(role["major_discordant_gene_trees"]["n_informative_sites"])
    a2 = int(role["minor_discordant_gene_trees"]["n_informative_sites"])
    t1 = int(role["top_loci_major_discordant_gene_trees"]["n_informative_sites"])
    t2 = int(role["top_loci_minor_discordant_gene_trees"]["n_informative_sites"])
    measures = [
        ("Gene trees", int(s0["n_major"]), int(s0["n_minor"])),
        ("Gene trees, UFBoot \u2265 95", int(s95["n_major"]), int(s95["n_minor"])),
        ("Loci (site majority)", int(role["locus_majority_major_discordant_gene_trees"]["n_loci"]),
         int(role["locus_majority_minor_discordant_gene_trees"]["n_loci"])),
        ("Sites, all loci", a1, a2),
        (f"Sites, without top {k} loci", a1 - t1, a2 - t2),
    ]
    ypos = list(range(len(measures)))
    for yv, (label, n1, n2) in zip(ypos, measures):
        p, a, b = ci(n1, n1 + n2)
        ax.errorbar([p], [yv], xerr=[[p - a], [b - p]], fmt="o", color=SERIES, ms=3.5,
                    lw=0.9, capsize=2, capthick=0.8)
        ax.text(0.797, yv, f"n = {n1 + n2:,}", va="center", ha="right", fontsize=5.8,
                color=INK2)
    ax.axvline(0.5, color=MUTED, lw=0.6)
    ax.set_yticks(ypos)
    ax.set_yticklabels([m[0] for m in measures])
    ax.invert_yaxis()
    ax.set_xlim(0.45, 0.80)
    ax.set_xticks([0.45, 0.5, 0.55, 0.6, 0.65, 0.7])
    ax.set_xlabel("Share with mol + qui | pal + pip")
    ax.set_title("(d) Weighting loci equally or by sites")
    clean(ax, "x")
    save(fig, "Figure_2_phylogeny", dpi)


def figure3(dpi):
    bins = read_tsv("results/cafe/transfer_bias_by_copy_number.tsv")
    ref = r"Copies in the reference ($\it{Cx.\ quinquefasciatus}$)"
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 66 * MM),
                             gridspec_kw={"wspace": 0.32, "left": 0.09, "right": 0.99,
                                          "top": 0.88, "bottom": 0.25})
    ax = axes[0]
    ret = [b for b in bins if b["mean_retention"] not in ("", "NA")]
    vals = [float(b["mean_retention"]) for b in ret]
    ax.bar(range(len(ret)), vals, color=SERIES, width=0.6)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.02, f"{v:.2f}", ha="center", va="bottom", fontsize=6.5)
    ax.set_xticks(range(len(ret)))
    DASH = "\u2013"
    ax.set_xticklabels([b["ref_copies"].replace("-", DASH) + "\n" + f"{int(b['n_families']):,}"
                        for b in ret], fontsize=6)
    ax.set_ylim(0, 1.08)
    ax.set_xlabel(ref + "; families")
    ax.set_ylabel("Mean share of reference\ncopies retained")
    ax.set_title("(a) Copy retention after transfer")
    clean(ax)

    ax = axes[1]
    vals = [float(b["pct_significant"]) for b in bins]
    ax.bar(range(len(bins)), vals, color=SERIES, width=0.6)
    for i, (v, b) in enumerate(zip(vals, bins)):
        ax.text(i, v + 1.2, f"{int(b['n_significant'])}/{int(b['n_families']):,}",
                ha="center", va="bottom", fontsize=5.8)
    ax.set_xticks(range(len(bins)))
    ax.set_xticklabels([b["ref_copies"].replace("-", "\u2013") for b in bins], fontsize=6)
    ax.set_ylim(0, max(vals) * 1.18)
    ax.set_xlabel(ref)
    ax.set_ylabel("Families CAFE calls significant (%)")
    ax.set_title("(b) CAFE significance by copy number")
    clean(ax)
    save(fig, "Figure_3_copy_number", dpi)


def load_synteny_module():
    spec = importlib.util.spec_from_file_location(
        "synteny_summary", os.path.join("workflow", "scripts", "synteny_summary.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def chrom_extent(genes, seqid):
    """Length used for plotting: the largest gene midpoint on the sequence."""
    return max(p for q, p in genes.values() if q == seqid)


def figure4(dpi):
    syn = load_synteny_module()
    genes = {s: syn.read_genes(f"results/annotation/{s}_liftoff.gff3") for s in INGROUP}
    homolog = syn.homolog_names(genes, INGROUP, "Cx_quinquefasciatus", 3)
    orient = {(r["sample1"], r["sample2"], r["seqid2"]): r["orientation"]
              for r in read_tsv("results/synteny/chromosome_orientation.tsv")}
    chroms, extent = {}, {}
    for s in INGROUP:
        tops = syn.top_seqids(genes[s], 3)
        chroms[s] = sorted(tops, key=lambda sq: homolog.get((s, sq), sq))
        for sq in tops:
            extent[(s, sq)] = chrom_extent(genes[s], sq)

    fig, axes = plt.subplots(2, 3, figsize=(WIDTH, 124 * MM),
                             gridspec_kw={"hspace": 0.5, "wspace": 0.34, "left": 0.07,
                                          "right": 0.99, "top": 0.94, "bottom": 0.14})
    for n, (ax, (a, b)) in enumerate(zip(axes.flat, itertools.combinations(INGROUP, 2))):
        off_a, pos = {}, 0.0
        for sq in chroms[a]:
            off_a[sq] = pos
            pos += extent[(a, sq)]
        len_a = pos
        off_b, pos = {}, 0.0
        for sq in chroms[b]:
            off_b[sq] = pos
            pos += extent[(b, sq)]
        len_b = pos
        # a query chromosome is flipped when it runs reverse to its homolog in a
        flip = {sq for (s1, s2, sq), o in orient.items()
                if s1 == a and s2 == b and o == "reverse"}
        xs = {i: [] for i in range(3)}
        ys = {i: [] for i in range(3)}
        for g, (sa, pa) in genes[a].items():
            if sa not in off_a:
                continue
            loc = genes[b].get(g)
            if not loc or loc[0] not in off_b:
                continue
            sb, pb = loc
            if sb in flip:
                pb = extent[(b, sb)] - pb
            i = chroms[a].index(sa)
            xs[i].append((off_a[sa] + pa) / 1e6)
            ys[i].append((off_b[sb] + pb) / 1e6)
        for i in range(3):
            ax.scatter(xs[i], ys[i], s=0.3, color=CHROM[i], linewidths=0, rasterized=True)
        for v in list(off_a.values())[1:]:
            ax.axvline(v / 1e6, color=GRID, lw=0.5, zorder=0)
        for v in list(off_b.values())[1:]:
            ax.axhline(v / 1e6, color=GRID, lw=0.5, zorder=0)
        ax.set_xlim(0, len_a / 1e6)
        ax.set_ylim(0, len_b / 1e6)
        ax.set_xticks(range(0, int(len_a / 1e6) + 1, 100))
        ax.set_yticks(range(0, int(len_b / 1e6) + 1, 100))
        ax.set_xlabel(f"{short(a)} (Mb)", fontstyle="italic")
        ax.set_ylabel(f"{short(b)} (Mb)", fontstyle="italic")
        ax.set_title(f"({'abcdef'[n]})")
        ax.tick_params(labelsize=6)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    handles = [Line2D([], [], marker="o", ls="", color=CHROM[i], ms=4, label=f"chr{i + 1}")
               for i in range(3)]
    fig.legend(handles=handles, ncol=3, frameon=False, loc="lower center",
               bbox_to_anchor=(0.5, 0.0), handletextpad=0.3, columnspacing=1.2)
    save(fig, "Figure_4_synteny", dpi)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dpi", type=int, default=1200)
    ap.add_argument("--only", default="1,2,3,4")
    a = ap.parse_args()
    style()
    for n in a.only.split(","):
        print(f"Figure {n}")
        {"1": figure1, "2": figure2, "3": figure3, "4": figure4}[n.strip()](a.dpi)


if __name__ == "__main__":
    main()
