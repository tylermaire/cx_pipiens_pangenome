import pandas as pd
import glob, os, re

of_dir    = snakemake.input.counts
tree_file = snakemake.input.tree

# --- Build CAFE5 gene-count input ---
gc_files = glob.glob(os.path.join(of_dir, "**", "Orthogroups.GeneCount.tsv"), recursive=True)
gene_counts = pd.read_csv(gc_files[0], sep="\t", index_col=0)
if "Total" in gene_counts.columns:
    gene_counts = gene_counts.drop(columns=["Total"])
max_per_species = gene_counts.max(axis=1)
filtered = gene_counts[max_per_species <= 100].copy()
filtered = filtered[filtered.sum(axis=1) > 0]
print(f"Total orthogroups: {len(gene_counts)}")
print(f"After filtering: {len(filtered)}")
cafe_table = filtered.copy()
cafe_table.insert(0, "Desc", "n/a")
cafe_table.index.name = "Family ID"
cafe_table.to_csv(snakemake.output.counts, sep="\t")

# --- Build CAFE5 ultrametric tree using the ML topology ---
# Read the IQ-TREE ML tree; extract topology only (ignore branch lengths).
with open(tree_file) as f:
    ml_nwk = f.read().strip().rstrip(";")

# Strip branch lengths and support values, leave only taxon names + topology
topo = re.sub(r":[0-9.eE+-]+", "", ml_nwk)
topo = re.sub(r"\)\d+", ")", topo)
print(f"ML topology (no branch lengths): {topo}")

# Convert unrooted 4-taxon Newick to a rooted ultrametric tree.
# Identify the 4 taxa and the bipartition.
import io
taxa = re.findall(r"[A-Za-z_][A-Za-z0-9_]*", topo)
print(f"Taxa: {taxa}")

# Parse the bipartition: the internal pair is whichever two are inside ()
import re as _re
m = _re.search(r"\(([^()]+)\)", topo)
if m:
    inner_pair = [t.strip() for t in m.group(1).split(",")]
    outer_pair = [t for t in taxa if t not in inner_pair]
else:
    inner_pair = taxa[:2]
    outer_pair = taxa[2:]
print(f"Bipartition: {inner_pair} | {outer_pair}")

# Build a rooted ultrametric tree where both clades are equal depth (1.5):
# (((A,B):0.5):1.0,((C,D):0.5):1.0); with leaf branches of 1.0 each
ultrametric = (
    f"(({inner_pair[0]}:1.0,{inner_pair[1]}:1.0):0.5,"
    f"({outer_pair[0]}:1.0,{outer_pair[1]}:1.0):0.5);"
)
with open(snakemake.output.tree, "w") as f:
    f.write(ultrametric + "\n")
print(f"Ultrametric tree: {ultrametric}")
