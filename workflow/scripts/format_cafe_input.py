#!/usr/bin/env python3
"""format_cafe_input.py - gene counts and an ultrametric tree for CAFE 5.

Gene counts: the OrthoFinder counts of the ingroup forms only (the outgroup
column is dropped), families with more than 100 genes in a form removed, and
families with no ingroup gene removed. The orthogroup id goes in Desc, since
CAFE names families by Desc.

Tree: CAFE needs a rooted ultrametric tree, and neither divergence times nor
calibrated branch lengths are known, so branch lengths are assumed. The
topology comes from the rooted five taxon tree (input.rooted, rooted with the
outgroup, which is then dropped). When the root falls between the two pairs
of forms, the tree is the one used since V3, ((A:1,B:1):0.5,(C:1,D:1):0.5);
when it falls elsewhere, node heights count the levels below each node and
branch lengths follow from them, which keeps the tree ultrametric. Without a
rooted tree (older runs), the unrooted four taxon ML tree is rooted between
its two pairs, as before.
"""
import glob
import os
import re

import pandas as pd

of_dir = snakemake.input.counts
tree_file = snakemake.input.tree
rooted_file = getattr(snakemake.input, "rooted", None)
ingroup = list(getattr(snakemake.params, "ingroup", []))
outgroup = list(getattr(snakemake.params, "outgroup", []))

# --- gene counts ---------------------------------------------------------
gc_files = glob.glob(os.path.join(of_dir, "**", "Orthogroups.GeneCount.tsv"), recursive=True)
gene_counts = pd.read_csv(gc_files[0], sep="\t", index_col=0)
if "Total" in gene_counts.columns:
    gene_counts = gene_counts.drop(columns=["Total"])
if ingroup:
    missing = [s for s in ingroup if s not in gene_counts.columns]
    if missing:
        raise SystemExit(f"ingroup samples missing from the OrthoFinder counts: {missing}")
    gene_counts = gene_counts[ingroup]
max_per_species = gene_counts.max(axis=1)
filtered = gene_counts[max_per_species <= 100].copy()
filtered = filtered[filtered.sum(axis=1) > 0]
print(f"Total orthogroups: {len(gene_counts)}")
print(f"With ingroup genes and at most 100 per form: {len(filtered)}")
cafe_table = filtered.copy()
cafe_table.insert(0, "Desc", cafe_table.index)
cafe_table.index.name = "Family ID"
cafe_table.to_csv(snakemake.output.counts, sep="\t")


# --- rooted topology -------------------------------------------------------
def unrooted_pairs(path):
    """The two pairs of the unrooted four taxon ML tree."""
    nwk = open(path).read().strip().rstrip(";")
    topo = re.sub(r":[0-9.eE+-]+", "", nwk)
    topo = re.sub(r"\)[0-9./]+", ")", topo)
    taxa = re.findall(r"[A-Za-z_][A-Za-z0-9_]*", topo)
    m = re.search(r"\(([^()]+)\)", topo)
    inner = [t.strip() for t in m.group(1).split(",")] if m else taxa[:2]
    outer = [t for t in taxa if t not in inner]
    return inner, outer


def rooted_ingroup(path, out):
    """Nested tuples of the ingroup topology, rooted with the outgroup."""
    import sys
    sys.path.insert(0, os.path.join(os.getcwd(), "workflow", "scripts"))
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from rooting import parse, ingroup_clades, nested
    if len(out) != 1:
        raise SystemExit(f"expected one outgroup, got {out}")
    tree = parse(open(path).read())
    clades = [c for c, _, _ in ingroup_clades(tree, out[0])]
    ingroup = sorted(set(tree.tips()) - {out[0]})
    return nested(clades, ingroup)


def height(node):
    return 0 if isinstance(node, str) else 1 + max(height(c) for c in node)


def ultrametric(node, parent_height=None):
    """Newick with branch lengths from node heights (levels below a node)."""
    h = height(node)
    if isinstance(node, str):
        text = node
    else:
        text = "(" + ",".join(ultrametric(c, h) for c in node) + ")"
    if parent_height is None:
        return text + ";"
    return f"{text}:{float(parent_height - h)}"


if rooted_file and os.path.exists(rooted_file) and outgroup:
    topo = rooted_ingroup(rooted_file, outgroup)
    print(f"Rooted ingroup topology: {topo}")
    balanced = (not isinstance(topo, str) and len(topo) == 2
                and all(not isinstance(c, str) and len(c) == 2
                        and all(isinstance(t, str) for t in c) for c in topo))
    if balanced:
        (a, b), (c, d) = topo
        tree = f"(({a}:1.0,{b}:1.0):0.5,({c}:1.0,{d}:1.0):0.5);"
        print("Root between the two pairs: tips 1, the two root branches 0.5 (assumed)")
    else:
        tree = ultrametric(topo)
        print("Root not between the two pairs: branch lengths from node levels (assumed)")
else:
    inner, outer = unrooted_pairs(tree_file)
    print(f"No rooted tree; rooting the ML tree between {inner} and {outer}")
    tree = (f"(({inner[0]}:1.0,{inner[1]}:1.0):0.5,"
            f"({outer[0]}:1.0,{outer[1]}:1.0):0.5);")
with open(snakemake.output.tree, "w") as f:
    f.write(tree + "\n")
print(f"Ultrametric tree: {tree}")
