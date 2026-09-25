#!/usr/bin/env python3
"""
cloud_composition.py - what the cloud orthogroups are made of.

A cloud orthogroup holds genes from exactly one ingroup genome. In this study
that reading has two problems, and this script measures both.

1. Outgroup genes. OrthoFinder orthogroups contain at least two genes, and
   unassigned singletons are not partitioned. A cloud orthogroup with a single
   ingroup gene therefore has to contain an outgroup gene (Cx. perexiguus from
   V5, Cx. tarsalis before), so "form specific" orthogroups can be shared with
   the outgroup.

2. Same gene elsewhere. Every model in the transferred genomes was placed by
   Liftoff from a Cx. quinquefasciatus gene and keeps that gene's ID. If a
   cloud gene's ID also carries a protein coding model in another ingroup
   genome, the gene is present there and the orthogroup is a clustering split,
   not a presence/absence difference.

Outputs:
  cloud_composition.tsv          one row per cloud orthogroup
  cloud_composition_summary.tsv  one row per form, plus a total

Snakemake provides:
    input.table, input.of, input.gffs, input.proteins
    params.ingroup, params.reference
    output.per_og, output.summary
"""

import glob
import os
import re

import pandas as pd

ATTR_ID = re.compile(r"(?:^|;)ID=([^;]+)")
ATTR_PARENT = re.compile(r"(?:^|;)Parent=([^;,]+)")
TRANSCRIPT_TYPES = {"mRNA", "transcript"}


def fasta_ids(path):
    with open(path) as fh:
        return {line[1:].split()[0] for line in fh if line.startswith(">")}


def transcript_to_gene(path):
    """{transcript: gene}, periods removed to match the protein FASTA headers."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in TRANSCRIPT_TYPES:
                continue
            m = ATTR_ID.search(f[8])
            if not m:
                continue
            tid = m.group(1).replace(".", "")
            p = ATTR_PARENT.search(f[8])
            out[tid] = p.group(1).replace(".", "") if p else tid
    return out


def load_members(of_dir):
    files = glob.glob(os.path.join(of_dir, "**", "Orthogroups.tsv"), recursive=True)
    if not files:
        raise SystemExit(f"No Orthogroups.tsv found under {of_dir}")
    return pd.read_csv(files[0], sep="\t", index_col=0).fillna("")


def split_genes(cell):
    return [g.strip() for g in str(cell).split(",") if g.strip()]


def compose(table, members, tx_to_gene, coding, ingroup, reference):
    """Per cloud orthogroup composition.

    table       partitioned_orthogroups.tsv as a DataFrame (index = orthogroup)
    members     Orthogroups.tsv as a DataFrame (index = orthogroup)
    tx_to_gene  {form: {transcript: gene}}
    coding      {form: set of transcripts with a protein}
    """
    outgroup = [c for c in members.columns if c not in ingroup]
    gene_model = {}           # (form, gene) -> transcript with a protein
    tx_to_og = {}             # (form, transcript) -> orthogroup
    for form in ingroup:
        for tid in coding[form]:
            gene_model.setdefault((form, tx_to_gene[form].get(tid, tid)), tid)
    for og, row in members.iterrows():
        for form in ingroup:
            if form in members.columns:
                for g in split_genes(row[form]):
                    tx_to_og[(form, g)] = og

    compartment = table["compartment"].to_dict()
    rows = []
    cloud = table[table["compartment"] == "cloud"]
    for og, r in cloud.iterrows():
        form = next(s for s in ingroup if r[s] > 0)
        genes = split_genes(members.at[og, form]) if og in members.index else []
        gene_ids = sorted({tx_to_gene[form].get(g, g) for g in genes})
        n_out = int(sum(r[c] for c in outgroup if c in r.index))
        where = []            # (other form, transcript, orthogroup)
        for other in ingroup:
            if other == form:
                continue
            for gid in gene_ids:
                t = gene_model.get((other, gid))
                if t:
                    where.append((other, t, tx_to_og.get((other, t), "unassigned")))
        other_ogs = sorted({w[2] for w in where})
        rows.append({
            "orthogroup": og,
            "form": form,
            "is_reference_form": form == reference,
            "n_ingroup_genes": int(r[form]),
            "n_outgroup_genes": n_out,
            "genes": ",".join(genes),
            "gene_ids": ",".join(gene_ids),
            "same_gene_forms": ",".join(sorted({w[0] for w in where})),
            "same_gene_orthogroups": ",".join(other_ogs),
            "same_gene_compartments": ",".join(sorted({compartment.get(o, o)
                                                       for o in other_ogs})),
            "in_reference_elsewhere": any(w[0] == reference for w in where),
            "status": "same_gene_elsewhere" if where else "no_counterpart",
        })
    return pd.DataFrame(rows)


def summarise(per_og, ingroup):
    rows = []
    for form in ingroup + ["total"]:
        sub = per_og if form == "total" else per_og[per_og["form"] == form]
        rows.append({
            "form": form,
            "n_cloud_orthogroups": len(sub),
            "with_outgroup_gene": int((sub["n_outgroup_genes"] > 0).sum()),
            "single_ingroup_gene": int((sub["n_ingroup_genes"] == 1).sum()),
            "same_gene_elsewhere": int((sub["status"] == "same_gene_elsewhere").sum()),
            "same_gene_in_reference": int(sub["in_reference_elsewhere"].sum()),
            "no_counterpart": int((sub["status"] == "no_counterpart").sum()),
        })
    out = pd.DataFrame(rows)
    out["pct_same_gene_elsewhere"] = (100.0 * out["same_gene_elsewhere"]
                                      / out["n_cloud_orthogroups"].clip(lower=1)).round(1)
    return out


def main():
    ingroup = list(snakemake.params.ingroup)
    reference = snakemake.params.reference
    table = pd.read_csv(snakemake.input.table, sep="\t", index_col=0)
    members = load_members(snakemake.input.of)
    gffs = {os.path.basename(p).replace("_liftoff.gff3", ""): p
            for p in snakemake.input.gffs}
    prots = {os.path.basename(p)[:-3]: p for p in snakemake.input.proteins}
    tx_to_gene = {f: transcript_to_gene(gffs[f]) for f in ingroup}
    coding = {f: fasta_ids(prots[f]) for f in ingroup}

    per_og = compose(table, members, tx_to_gene, coding, ingroup, reference)
    per_og.to_csv(snakemake.output.per_og, sep="\t", index=False)
    summary = summarise(per_og, ingroup)
    summary.to_csv(snakemake.output.summary, sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
