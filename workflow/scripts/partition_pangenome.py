#!/usr/bin/env python3
"""
Partition OrthoFinder orthogroups into core, shell, and cloud pangenome compartments.
"""
import pandas as pd
import os
import glob

# Find the gene count file within OrthoFinder output
# OrthoFinder creates a dated subdirectory, so we need to search for it
of_dir = snakemake.input[0]
gene_count_files = glob.glob(
    os.path.join(of_dir, "**", "Orthogroups.GeneCount.tsv"),
    recursive=True
)

if not gene_count_files:
    raise FileNotFoundError(
        f"Could not find Orthogroups.GeneCount.tsv in {of_dir}"
    )

gene_counts = pd.read_csv(gene_count_files[0], sep="\t", index_col=0)

# Drop the 'Total' column if present
if "Total" in gene_counts.columns:
    gene_counts = gene_counts.drop(columns=["Total"])

# Identify ingroup samples (exclude outgroup - Cx_tarsalis)
all_samples = list(gene_counts.columns)
outgroup = [s for s in all_samples if "tarsalis" in s.lower()]
ingroup = [s for s in all_samples if s not in outgroup]

print(f"Ingroup samples: {ingroup}")
print(f"Outgroup samples: {outgroup}")
print(f"Total orthogroups: {len(gene_counts)}")

# Count how many ingroup species have at least 1 gene per orthogroup
ingroup_presence = (gene_counts[ingroup] > 0).sum(axis=1)
n_ingroup = len(ingroup)

# Classify each orthogroup
gene_counts["ingroup_count"] = ingroup_presence
gene_counts["compartment"] = "unassigned"
gene_counts.loc[ingroup_presence == n_ingroup, "compartment"] = "core"
gene_counts.loc[
    (ingroup_presence >= 2) & (ingroup_presence < n_ingroup), "compartment"
] = "shell"
gene_counts.loc[ingroup_presence == 1, "compartment"] = "cloud"
gene_counts.loc[ingroup_presence == 0, "compartment"] = "outgroup_only"

# Save full partitioned table
gene_counts.to_csv(snakemake.output.table, sep="\t")

# Generate summary
summary_rows = []

# Overall counts
for comp in ["core", "shell", "cloud", "outgroup_only"]:
    subset = gene_counts[gene_counts["compartment"] == comp]
    n_orthogroups = len(subset)
    total_genes = subset[ingroup].sum().sum()
    summary_rows.append({
        "compartment": comp,
        "n_orthogroups": n_orthogroups,
        "total_genes_ingroup": int(total_genes),
    })

summary = pd.DataFrame(summary_rows)

# Add per-species counts for each compartment
for comp in ["core", "shell", "cloud"]:
    subset = gene_counts[gene_counts["compartment"] == comp]
    for sp in ingroup:
        col_name = f"{comp}_{sp}"
        summary.loc[summary["compartment"] == comp, col_name] = int(
            (subset[sp] > 0).sum()
        )

# Add species-specific cloud counts
cloud = gene_counts[gene_counts["compartment"] == "cloud"]
for sp in ingroup:
    # Orthogroups where ONLY this species has genes
    sp_specific = cloud[
        (cloud[sp] > 0) & (cloud[ingroup].drop(columns=[sp]).sum(axis=1) == 0)
    ]
    summary_rows.append({
        "compartment": f"cloud_{sp}_specific",
        "n_orthogroups": len(sp_specific),
        "total_genes_ingroup": int(sp_specific[sp].sum()),
    })

summary = pd.DataFrame(summary_rows)
summary.to_csv(snakemake.output.summary, sep="\t", index=False)

# Print summary to log
print("\n=== PANGENOME SUMMARY ===")
for _, row in summary.iterrows():
    print(f"  {row['compartment']}: {row['n_orthogroups']} orthogroups, "
          f"{row['total_genes_ingroup']} genes")
print("=========================\n")
