#!/usr/bin/env python3
"""
extract_absence_queries.py - build the query sets needed to test every absence
call in the pangenome.

An "absence event" is a (orthogroup, genome) pair where the orthogroup table
says the orthogroup has no gene in that genome. Cloud orthogroups generate
three absence events each, shell orthogroups one or two. Every one of those is
a claim that a gene is missing, and every one is testable.

For each orthogroup this writes:
  queries.faa   the longest protein in the orthogroup (protein-level probe)
  queries.fna   that gene's genomic locus, introns included (DNA-level probe)
  manifest.tsv  orthogroup, compartment, representative, present and absent forms

The two probes fail differently. Protein alignment tolerates divergence but can
land on a paralog; DNA alignment is positionally stricter but degrades faster
with sequence distance. Agreement between them is what makes a call solid.

Snakemake provides:
    input.table, input.of, params.ingroup, params.flank, output.dir
"""

import glob
import os

import pandas as pd

table = snakemake.input.table
of_dir = snakemake.input.of
ingroup = list(snakemake.params.ingroup)
flank = int(getattr(snakemake.params, "flank", 0))
out_dir = snakemake.output.dir

os.makedirs(out_dir, exist_ok=True)


def read_fasta(path):
    seqs, rid, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if rid is not None:
                    seqs[rid] = "".join(chunks)
                rid = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if rid is not None:
        seqs[rid] = "".join(chunks)
    return seqs


def load_fai(fasta):
    """Return {seqid: (offset, linebases, linewidth)} from the .fai index."""
    idx = {}
    with open(fasta + ".fai") as fh:
        for line in fh:
            name, length, offset, lb, lw = line.split("\t")[:5]
            idx[name] = (int(offset), int(lb), int(lw), int(length))
    return idx


def fetch(fh, idx, seqid, start, end):
    """1-based inclusive slice of a FASTA sequence using its .fai index."""
    if seqid not in idx:
        return None
    offset, lb, lw, length = idx[seqid]
    start = max(1, start)
    end = min(length, end)
    if end < start:
        return None
    out = []
    # walk line by line so we never load a whole chromosome
    first_line = (start - 1) // lb
    pos_in_line = (start - 1) % lb
    fh.seek(offset + first_line * lw + pos_in_line)
    need = end - start + 1
    while need > 0:
        take = min(need, lb - pos_in_line)
        chunk = fh.read(take)
        out.append(chunk)
        need -= take
        fh.read(lw - lb)  # skip the newline(s)
        pos_in_line = 0
    return "".join(out)


def mrna_spans(gff):
    """{stripped mRNA id: (seqid, start, end)} — periods removed to match FASTA."""
    spans = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("mRNA", "transcript"):
                continue
            attr = f[8]
            i = attr.find("ID=")
            if i < 0:
                continue
            tid = attr[i + 3:].split(";")[0].replace(".", "")
            spans[tid] = (f[0], int(f[3]), int(f[4]))
    return spans


part = pd.read_csv(table, sep="\t")
part = part[part["compartment"].isin(["cloud", "shell"])]

og_files = glob.glob(os.path.join(of_dir, "**", "Orthogroups.tsv"), recursive=True)
if not og_files:
    raise SystemExit(f"No Orthogroups.tsv found under {of_dir}")
members = pd.read_csv(og_files[0], sep="\t", index_col=0).fillna("")

proteomes = {s: read_fasta(f"results/proteins/{s}.fa") for s in ingroup}
spans = {s: mrna_spans(f"results/annotation/{s}_liftoff.gff3") for s in ingroup}
fais = {s: load_fai(f"resources/genomes/{s}.fasta") for s in ingroup}
handles = {s: open(f"resources/genomes/{s}.fasta") for s in ingroup}

faa = open(os.path.join(out_dir, "queries.faa"), "w")
fna = open(os.path.join(out_dir, "queries.fna"), "w")
man = open(os.path.join(out_dir, "manifest.tsv"), "w")
man.write("orthogroup\tcompartment\trep_gene\trep_form\tpresent_forms\tabsent_forms\n")

n_dna = n_events = 0
for _, row in part.iterrows():
    og = row["Orthogroup"]
    if og not in members.index:
        continue
    present = [s for s in ingroup if row[s] > 0]
    absent = [s for s in ingroup if row[s] == 0]
    if not present or not absent:
        continue

    # representative = longest protein anywhere in the orthogroup
    best = (0, None, None)
    for form in present:
        for gene in [g.strip() for g in str(members.at[og, form]).split(",") if g.strip()]:
            seq = proteomes[form].get(gene)
            if seq and len(seq) > best[0]:
                best = (len(seq), gene, form)
    _, gene, form = best
    if gene is None:
        continue

    faa.write(f">{og} gene={gene} form={form}\n")
    seq = proteomes[form][gene]
    for i in range(0, len(seq), 60):
        faa.write(seq[i:i + 60] + "\n")

    loc = spans[form].get(gene)
    if loc:
        seqid, start, end = loc
        dna = fetch(handles[form], fais[form], seqid, start - flank, end + flank)
        if dna:
            fna.write(f">{og} gene={gene} form={form} loc={seqid}:{start}-{end}\n")
            for i in range(0, len(dna), 60):
                fna.write(dna[i:i + 60] + "\n")
            n_dna += 1

    man.write(f"{og}\t{row['compartment']}\t{gene}\t{form}\t"
              f"{','.join(present)}\t{','.join(absent)}\n")
    n_events += len(absent)

for h in handles.values():
    h.close()
faa.close(); fna.close(); man.close()

print(f"orthogroups queried: {len(part)}")
print(f"genomic loci recovered: {n_dna}")
print(f"absence events to test: {n_events}")
