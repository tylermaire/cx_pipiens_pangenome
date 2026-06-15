import glob, os, pandas as pd
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

# Per-species sequence dictionary to avoid ID collisions
per_species_seqs = {}
for fa_file in glob.glob("results/proteins/*.fa"):
    sp = os.path.basename(fa_file).replace(".fa", "")
    per_species_seqs[sp] = {
        rec.id: str(rec.seq) for rec in SeqIO.parse(fa_file, "fasta")
    }
print("Loaded proteomes:", {sp: len(d) for sp, d in per_species_seqs.items()})

written = 0
for og in sco_ogs:
    if og not in og_membership.index:
        continue
    seqs_for_this_og = {}
    for sp in ingroup:
        if sp not in og_membership.columns:
            continue
        genes = og_membership.loc[og, sp]
        if pd.isna(genes):
            continue
        gene_id = str(genes).split(",")[0].strip()
        if sp in per_species_seqs and gene_id in per_species_seqs[sp]:
            seqs_for_this_og[sp] = per_species_seqs[sp][gene_id]
    if len(seqs_for_this_og) == len(ingroup):
        outfile = os.path.join(outdir, f"{og}.fa")
        with open(outfile, "w") as f:
            for sp in ingroup:
                f.write(f">{sp}\n{seqs_for_this_og[sp]}\n")
        written += 1

print(f"Wrote {written} complete SCO FASTA files to {outdir}")
