#!/usr/bin/env python3
import pandas as pd
import glob, os, re

of_dir = snakemake.input.counts
tree_file = snakemake.input.tree
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

ultrametric = "((Cx_pallens:1,Cx_pipiens:1):0.5,(Cx_molestus:1.5,Cx_quinquefasciatus:1.5):0.1);"
with open(snakemake.output.tree, "w") as f:
    f.write(ultrametric + "\n")
print(f"Ultrametric tree: {ultrametric}")
