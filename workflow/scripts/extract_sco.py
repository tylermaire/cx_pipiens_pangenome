#!/usr/bin/env python3
import os, glob
import pandas as pd
from Bio import SeqIO

of_dir = snakemake.input[0]
outdir = snakemake.output[0]
os.makedirs(outdir, exist_ok=True)

gc_files = glob.glob(os.path.join(of_dir, "**", "Orthogroups.GeneCount.tsv"), recursive=True)
og_files = glob.glob(os.path.join(of_dir, "**", "Orthogroups.tsv"), recursive=True)

gene_counts = pd.read_csv(gc_files[0], sep="\t", index_col=0)
if "Total" in gene_counts.columns:
    gene_counts = gene_counts.drop(columns=["Total"])

ingroup = [c for c in gene_counts.columns if "tarsalis" not in c.lower()]
sco_mask = True
for sp in ingroup:
    sco_mask = sco_mask & (gene_counts[sp] == 1)
sco_ogs = gene_counts[sco_mask].index.tolist()
print(f"Found {len(sco_ogs)} single-copy orthogroups across ingroup")

og_membership = pd.read_csv(og_files[0], sep="\t", index_col=0)

all_seqs = {}
for fa_file in glob.glob("results/proteins/*.fa"):
    for rec in SeqIO.parse(fa_file, "fasta"):
        all_seqs[rec.id] = rec
seq_dir_list = glob.glob(os.path.join(of_dir, "**", "WorkingDirectory"), recursive=True)
if seq_dir_list:
    for fa_file in glob.glob(os.path.join(seq_dir_list[0], "Species*.fa")):
        for rec in SeqIO.parse(fa_file, "fasta"):
            all_seqs[rec.id] = rec

written = 0
for og in sco_ogs:
    if og not in og_membership.index:
        continue
    outfile = os.path.join(outdir, f"{og}.fa")
    seqs_written = 0
    with open(outfile, "w") as f:
        for sp in ingroup:
            if sp not in og_membership.columns:
                continue
            genes = og_membership.loc[og, sp]
            if pd.isna(genes):
                continue
            gene_list = [g.strip() for g in str(genes).split(",")]
            for gene_id in gene_list:
                if gene_id in all_seqs:
                    f.write(f">{sp}\n{str(all_seqs[gene_id].seq)}\n")
                    seqs_written += 1
    if seqs_written == len(ingroup):
        written += 1
    else:
        os.remove(outfile)
print(f"Wrote {written} complete SCO FASTA files to {outdir}")
