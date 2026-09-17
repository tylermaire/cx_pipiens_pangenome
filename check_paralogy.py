#!/usr/bin/env python3
"""
check_paralogy.py - decide whether 'clustered_elsewhere' calls are the same gene
split across two orthogroups, or a paralog standing in for a truly absent gene.

For each such call we have a query protein (the orthogroup's representative) and
the annotated gene its alignment landed on in the target genome. If those two
proteins are near identical, they are the same gene and the orthogroup was split
by clustering, so the recorded absence is an artifact. If they are only distantly
similar, the hit is a paralog and the absence is real.

At 93 to 95% ingroup ANI, true orthologs between these forms sit well above 90%
protein identity, and paralogs are typically far below it. The script reports the
full distribution so the threshold is a reading of the data, not an assumption.

Run from the repository root:
    python3 check_paralogy.py
"""

import collections
import sys

import pandas as pd

try:
    from Bio import Align
    from Bio.Align import substitution_matrices
except ImportError:
    sys.exit("Biopython required: run inside the phylo conda env")

CALLS = "results/validation/absence_calls.tsv"
MANIFEST = "results/validation/queries/manifest.tsv"


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


calls = pd.read_csv(CALLS, sep="\t")
man = pd.read_csv(MANIFEST, sep="\t").set_index("orthogroup")
sub = calls[calls["call"] == "clustered_elsewhere"].copy()
print(f"clustered_elsewhere events: {len(sub)}")

forms = sorted(set(man["rep_form"]) | set(sub["absent_from"]))
prot = {f: read_fasta(f"results/proteins/{f}.fa") for f in forms}

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -11
aligner.extend_gap_score = -1
aligner.mode = "global"


def identity(a, b):
    """Fraction of identical residues over the shorter sequence."""
    if not a or not b:
        return None
    # cheap length guard: wildly different lengths are not the same gene
    if min(len(a), len(b)) / max(len(a), len(b)) < 0.30:
        return 0.0
    aln = aligner.align(a, b)[0]
    ia, ib = aln.aligned
    same = 0
    for (s1, e1), (s2, e2) in zip(ia, ib):
        same += sum(1 for x, y in zip(a[s1:e1], b[s2:e2]) if x == y)
    return same / min(len(a), len(b))


rows = []
missing = 0
for i, r in enumerate(sub.itertuples(), 1):
    og, target, other = r.orthogroup, r.absent_from, r.overlapping_gene
    if og not in man.index or not isinstance(other, str):
        missing += 1
        continue
    rep_gene, rep_form = man.at[og, "rep_gene"], man.at[og, "rep_form"]
    qa = prot.get(rep_form, {}).get(rep_gene)
    qb = prot.get(target, {}).get(other)
    if qa is None or qb is None:
        missing += 1
        continue
    pid = identity(qa, qb)
    rows.append((og, r.compartment, target, rep_gene, other,
                 len(qa), len(qb), round(pid, 4)))
    if i % 500 == 0:
        print(f"  ...{i}/{len(sub)}", flush=True)

out = pd.DataFrame(rows, columns=["orthogroup", "compartment", "absent_from",
                                  "query_gene", "overlapping_gene",
                                  "query_len", "other_len", "identity"])
out.to_csv("results/validation/paralogy_check.tsv", sep="\t", index=False)
print(f"\npairs compared: {len(out)}   unresolvable: {missing}")

print("\nprotein identity distribution")
bins = [(0.95, 1.01, ">=95%  same gene"),
        (0.90, 0.95, "90-95%  same gene"),
        (0.70, 0.90, "70-90%  ambiguous"),
        (0.40, 0.70, "40-70%  likely paralog"),
        (0.00, 0.40, "<40%    paralog or spurious")]
for lo, hi, lab in bins:
    n = ((out["identity"] >= lo) & (out["identity"] < hi)).sum()
    print(f"  {lab:28s} {n:6d}  ({n/max(len(out),1):6.1%})")

for comp in ["cloud", "shell"]:
    c = out[out["compartment"] == comp]
    if c.empty:
        continue
    split = (c["identity"] >= 0.90).sum()
    print(f"\n{comp}: {split}/{len(c)} ({split/len(c):.1%}) of clustered_elsewhere "
          f"calls are the SAME gene (>=90% identity), i.e. a clustering split")
