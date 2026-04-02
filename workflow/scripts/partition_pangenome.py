#!/usr/bin/env python3
import pandas as pd
import os
import glob

of_dir = snakemake.input[0]
gc_files = glob.glob(os.path.join(of_dir, "**", "Orthogroups.GeneCount.tsv"), recursive=True)
if not gc_files:
    raise FileNotFoundError(f"Could not find Orthogroups.GeneCount.tsv in {of_dir}")

gene_counts = pd.read_csv(gc_files[0], sep="\t", index_col=0)
if "Total" in gene_counts.columns:
    gene_counts = gene_counts.drop(columns=["Total"])

all_samples = list(gene_counts.columns)
outgroup = [s for s in all_samples if "tarsalis" in s.lower()]
ingroup = [s for s in all_samples if s not in outgroup]
n_ingroup = len(ingroup)

print(f"Ingroup samples: {ingroup}")
print(f"Outgroup samples: {outgroup}")
print(f"Total orthogroups: {len(gene_counts)}")

ingroup_presence = (gene_counts[ingroup] > 0).sum(axis=1)
gene_counts["ingroup_count"] = ingroup_presence
gene_counts["compartment"] = "unassigned"
gene_counts.loc[ingroup_presence == n_ingroup, "compartment"] = "core"
gene_counts.loc[(ingroup_presence >= 2) & (ingroup_presence < n_ingroup), "compartment"] = "shell"
gene_counts.loc[ingroup_presence == 1, "compartment"] = "cloud"
gene_counts.loc[ingroup_presence == 0, "compartment"] = "outgroup_only"

gene_counts.to_csv(snakemake.output.table, sep="\t")

summary_rows = []
for comp in ["core", "shell", "cloud", "outgroup_only"]:
    subset = gene_counts[gene_counts["compartment"] == comp]
    summary_rows.append({"compartment": comp, "n_orthogroups": len(subset), "total_genes_ingroup": int(subset[ingroup].sum().sum())})

cloud = gene_counts[gene_counts["compartment"] == "cloud"]
for sp in ingroup:
    sp_specific = cloud[(cloud[sp] > 0) & (cloud[ingroup].drop(columns=[sp]).sum(axis=1) == 0)]
    summary_rows.append({"compartment": f"cloud_{sp}_specific", "n_orthogroups": len(sp_specific), "total_genes_ingroup": int(sp_specific[sp].sum())})

summary = pd.DataFrame(summary_rows)
summary.to_csv(snakemake.output.summary, sep="\t", index=False)

print("\n=== PANGENOME SUMMARY ===")
for _, row in summary.iterrows():
    print(f"  {row['compartment']}: {row['n_orthogroups']} orthogroups, {row['total_genes_ingroup']} genes")
print("=========================\n")
