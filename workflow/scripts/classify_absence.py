#!/usr/bin/env python3
"""
classify_absence.py - decide, for every absence call in the pangenome, whether
the gene is really absent or merely unannotated.

Each absence event is an (orthogroup, genome) pair the orthogroup table says is
empty. Evidence is combined from four sources:

  protein   miniprot alignment of the orthogroup's representative protein
  dna       minimap2 alignment of that gene's genomic locus
  context   every gene the target genome annotates under the hit interval
  identity  pairwise protein identity between the query and EACH of those genes

The identity step is what separates a clustering split from a paralog. A hit
that lands on annotated ground proves nothing by itself: a spliced alignment can
span a long interval and overlap several unrelated genes. Only if one of those
genes is near identical to the query is it the same gene.

Calls:

  absent              no probe finds the locus. The absence is real
  paralog_only        a hit exists, but every gene under it is only distantly
                      similar. The hit is a paralog; the absence is real
  unannotated_locus   sequence is there, nothing annotated on it. Liftoff placed
                      no model, so the absence is an annotation gap
  clustered_elsewhere sequence is there and annotated, and that gene is near
                      identical to the query. The gene is present; orthology
                      clustering split it
  weak                a hit below threshold; reported separately, never binned

'absent' and 'paralog_only' both support a real presence/absence difference.

Snakemake provides:
    input.manifest, input.prot_pafs, input.dna_pafs, input.of
    params.min_id, params.min_cov, params.dna_min_cov, params.same_gene_id,
    params.ingroup, params.max_overlap
    output.calls, output.summary
"""

import collections
import glob
import os

import pandas as pd

try:
    from Bio import Align
    from Bio.Align import substitution_matrices
except ImportError:
    raise SystemExit("Biopython required in the phylo environment")

min_id = float(getattr(snakemake.params, "min_id", 0.80))
min_cov = float(getattr(snakemake.params, "min_cov", 0.70))
dna_min_cov = float(getattr(snakemake.params, "dna_min_cov", 0.50))
same_gene_id = float(getattr(snakemake.params, "same_gene_id", 0.90))
max_overlap = int(getattr(snakemake.params, "max_overlap", 10))
ingroup = list(snakemake.params.ingroup)


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


def parse_paf(path):
    """{(query, target_genome): (identity, coverage, seqid, start, end)}."""
    base = os.path.basename(path).replace(".paf", "")
    target = base.split("_vs_")[-1]
    best = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            q, qlen, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
            tname, ts, te = f[5], int(f[7]), int(f[8])
            nmatch, blocklen = int(f[9]), int(f[10])
            ident = nmatch / max(blocklen, 1)
            cov = (qe - qs) / max(qlen, 1)
            key = (q, target)
            if key not in best or ident * cov > best[key][0] * best[key][1]:
                best[key] = (ident, cov, tname, ts, te)
    return best


aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -11
aligner.extend_gap_score = -1
aligner.mode = "global"

_cache = {}


def protein_identity(a, b, key):
    """Identical residues over the shorter sequence, cached by gene pair."""
    if key in _cache:
        return _cache[key]
    if not a or not b or min(len(a), len(b)) / max(len(a), len(b)) < 0.30:
        _cache[key] = 0.0
        return 0.0
    aln = aligner.align(a, b)[0]
    same = 0
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        same += sum(1 for x, y in zip(a[s1:e1], b[s2:e2]) if x == y)
    val = same / min(len(a), len(b))
    _cache[key] = val
    return val


# --- annotation intervals and orthogroup membership per genome ---
og_files = glob.glob(os.path.join(snakemake.input.of, "**", "Orthogroups.tsv"),
                     recursive=True)
members = pd.read_csv(og_files[0], sep="\t", index_col=0).fillna("")
gene_to_og = {}
for og, row in members.iterrows():
    for form in ingroup:
        if form in members.columns:
            for g in [x.strip() for x in str(row[form]).split(",") if x.strip()]:
                gene_to_og[(form, g)] = og

intervals = {}
for form in ingroup:
    by_seq = collections.defaultdict(list)
    with open(f"results/annotation/{form}_liftoff.gff3") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("mRNA", "transcript"):
                continue
            i = f[8].find("ID=")
            if i < 0:
                continue
            tid = f[8][i + 3:].split(";")[0].replace(".", "")
            by_seq[f[0]].append((int(f[3]), int(f[4]), tid))
    for k in by_seq:
        by_seq[k].sort()
    intervals[form] = by_seq

proteomes = {f: read_fasta(f"results/proteins/{f}.fa") for f in ingroup}


def overlapping(form, seqid, start, end):
    """All annotated transcripts overlapping the interval, longest overlap first."""
    hits = []
    for s, e, tid in intervals[form].get(seqid, []):
        if s > end:
            break
        if e >= start:
            hits.append((min(e, end) - max(s, start), tid))
    hits.sort(reverse=True)
    return [tid for _, tid in hits[:max_overlap]]


prot, dna = {}, {}
for p in snakemake.input.prot_pafs:
    prot.update(parse_paf(p))
for p in snakemake.input.dna_pafs:
    dna.update(parse_paf(p))

man = pd.read_csv(snakemake.input.manifest, sep="\t")

rows = []
for n, r in enumerate(man.itertuples(), 1):
    og, comp = r.orthogroup, r.compartment
    qseq = proteomes.get(r.rep_form, {}).get(r.rep_gene)
    for target in str(r.absent_forms).split(","):
        target = target.strip()
        if not target:
            continue
        pi, pc, pseq, ps, pe = prot.get((og, target), (0.0, 0.0, "", 0, 0))
        di, dc, dseq, ds, de = dna.get((og, target), (0.0, 0.0, "", 0, 0))

        prot_ok = pi >= min_id and pc >= min_cov
        dna_ok = di >= min_id and dc >= dna_min_cov

        if not (prot_ok or dna_ok):
            call = "weak" if (pi > 0 or di > 0) else "absent"
            rows.append((og, comp, target, round(pi, 3), round(pc, 3),
                         round(di, 3), round(dc, 3), call, "", 0.0))
            continue

        seqid, s, e = (pseq, ps, pe) if prot_ok else (dseq, ds, de)
        genes = overlapping(target, seqid, s, e)
        # compare against EVERY overlapping gene, keep the best match
        best_id, best_gene = 0.0, ""
        for g in genes:
            val = protein_identity(qseq, proteomes[target].get(g),
                                   (r.rep_gene, target, g))
            if val > best_id:
                best_id, best_gene = val, g

        if not genes:
            call = "unannotated_locus"
        elif best_id >= same_gene_id:
            call = "clustered_elsewhere"
        else:
            call = "paralog_only"

        rows.append((og, comp, target, round(pi, 3), round(pc, 3),
                     round(di, 3), round(dc, 3), call, best_gene, round(best_id, 4)))
    if n % 500 == 0:
        print(f"  ...{n}/{len(man)} orthogroups", flush=True)

cols = ["orthogroup", "compartment", "absent_from", "prot_identity",
        "prot_coverage", "dna_identity", "dna_coverage", "call",
        "best_overlapping_gene", "best_protein_identity"]
out = pd.DataFrame(rows, columns=cols)
out.to_csv(snakemake.output.calls, sep="\t", index=False)

summary = (out.groupby(["compartment", "call"]).size()
           .unstack(fill_value=0).reset_index())
summary.to_csv(snakemake.output.summary, sep="\t", index=False)

real = {"absent", "paralog_only"}
print(f"\n=== absence validation (hit id>={min_id} cov>={min_cov}; "
      f"same-gene protein identity >={same_gene_id}) ===")
for comp in ["cloud", "shell"]:
    sub = out[out["compartment"] == comp]
    if sub.empty:
        continue
    n = len(sub)
    print(f"\n{comp}: {n} absence events")
    for call, k in sub["call"].value_counts().items():
        print(f"   {call:22s} {k:6d}  ({k / n:6.1%})")
    good = sub["call"].isin(real).sum()
    print(f"   -> {good}/{n} ({good / n:.1%}) of {comp} absences are supported")

print("\n'absent' and 'paralog_only' support a real presence/absence difference;")
print("'clustered_elsewhere' and 'unannotated_locus' are artifacts.")
