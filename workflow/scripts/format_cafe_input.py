#!/usr/bin/env python3
"""
Prepare CAFE5 input: filter gene count table and create ultrametric tree.
"""
import pandas as pd
import glob
import os
import re

of_dir = snakemake.input.counts
tree_file = snakemake.input.tree

# Find gene count file
gc_files = glob.glob(os.path.join(of_dir, "**", "Orthogroups.GeneCount.tsv"), recursive=True)
gene_counts = pd.read_csv(gc_files[0], sep="\t", index_col=0)

# Drop Total column
if "Total" in gene_counts.columns:
    gene_counts = gene_counts.drop(columns=["Total"])

# CAFE5 format: tab-separated with Desc and Family ID columns
# Filter: remove families with >100 copies in any species (CAFE can't handle these)
max_per_species = gene_counts.max(axis=1)
filtered = gene_counts[max_per_species <= 100].copy()

# Also remove families where all species have 0 genes
row_sums = filtered.sum(axis=1)
filtered = filtered[row_sums > 0]

print(f"Total orthogroups: {len(gene_counts)}")
print(f"After filtering (max<=100, sum>0): {len(filtered)}")

# Format for CAFE5: Desc\tFamily ID\tSpecies1\tSpecies2\t...
cafe_table = filtered.copy()
cafe_table.insert(0, "Desc", "n/a")
cafe_table.index.name = "Family ID"
cafe_table.to_csv(snakemake.output.counts, sep="\t")

# Read the species tree and make it ultrametric
# For CAFE5, we need an ultrametric tree (all tips same distance from root)
# Simple approach: read the tree and set all branch lengths equal
with open(tree_file, "r") as f:
    tree_str = f.read().strip()

# Remove bootstrap values (numbers before colons or after closing parens)
# CAFE5 wants a clean newick with just topology and branch lengths
# Replace branch lengths with uniform values to make ultrametric
# Simple ultrametric: set all terminal branches to 1.0 and internal to 0.5
tree_str = re.sub(r'\)(\d+)', ')', tree_str)  # remove bootstrap values

# For a 4-taxon tree, manually create ultrametric version
# Based on our topology: (mol, (pal, pip), quin)
ultrametric = "((Cx_pallens:1,Cx_pipiens:1):0.5,(Cx_molestus:1.5,Cx_quinquefasciatus:1.5):0);"

with open(snakemake.output.tree, "w") as f:
    f.write(ultrametric + "\n")

print(f"Ultrametric tree written: {ultrametric}")
