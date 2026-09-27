#!/usr/bin/env python3
"""locus_quality.py - gene model quality of each five taxon locus, the lowest
identity among its ingroup sequences, and the loci the rooted tree uses.

Why: where the transferred models of two or three ingroup forms share an
error relative to the reference (most often a reading frame shifted where
RefSeq corrected the reference model for an error in the reference genome),
the misaligned codons after it look like derived alleles shared by those
forms. In the first V5 run such loci carried most of the sites that placed
the root on the reference's branch (rooting_site_patterns.tsv), so the rooted
tree uses loci whose four ingroup models are intact (config rooted: loci).

Outputs
-------
  quality    one row per locus with a codon alignment or a trimmed protein
             alignment: reason ('intact' or why not, as
             quartet_asymmetry.locus_reasons), intact (as d_statistics,
             quartet_asymmetry.intact_loci), min_ingroup_identity (lowest
             pairwise identity of the ingroup sequences over codon alignment
             columns where both carry A, C, G or T), codon_columns,
             tree_alignment (a trimmed protein alignment of at least 50
             columns exists) and tree (the locus is in the rooted tree)
  tree_loci  the loci of the rooted tree, one per line

Snakemake provides input.aln (codon/ and trimmed/ inside), input.of,
input.gffs; params.ingroup, params.outgroup, params.reference, params.loci
('intact' or 'all'); output.quality, output.tree_loci.
"""
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.join(os.getcwd(), "workflow", "scripts"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from quartet_asymmetry import intact_loci, locus_reasons  # noqa: E402

BASES = frozenset("ACGT")


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def min_identity(seqs, taxa):
    """Lowest pairwise identity among taxa over columns where both carry a
    base; None when a pair shares no such column."""
    low = None
    for i, a in enumerate(taxa):
        for b in taxa[i + 1:]:
            pairs = [(x, y) for x, y in zip(seqs[a], seqs[b]) if x in BASES and y in BASES]
            if not pairs:
                return None
            ident = sum(x == y for x, y in pairs) / len(pairs)
            low = ident if low is None else min(low, ident)
    return low


def run(aln_dir, of_dir, gffs, ingroup, outgroup, reference, loci, out_quality, out_tree):
    if loci not in ("intact", "all"):
        raise SystemExit(f"config rooted: loci must be 'intact' or 'all', not {loci!r}")
    intact = intact_loci(of_dir, {s: gffs[s] for s in ingroup}, reference)
    reasons = locus_reasons(of_dir, {s: gffs[s] for s in ingroup}, reference)
    codon = {os.path.basename(p)[:-4]: p
             for p in glob.glob(os.path.join(aln_dir, "codon", "*.fna"))}
    trimmed = {os.path.basename(p)[:-5]
               for p in glob.glob(os.path.join(aln_dir, "trimmed", "*.trim"))}
    rows = []
    for og in sorted(set(codon) | trimmed):
        ident, cols = None, 0
        if og in codon:
            seqs = read_fasta(codon[og])
            if set(ingroup) | {outgroup} <= set(seqs):
                ident = min_identity(seqs, list(ingroup))
                cols = len(seqs[outgroup])
        ok = intact.get(og)
        reason = reasons.get(og, "not a single copy orthogroup of the ingroup")
        if ok is not None and ok != (reason == "intact"):
            reason = "not intact"            # never expected: the two checks share their rules
        in_tree = og in trimmed and (loci == "all" or ok is True)
        rows.append({"orthogroup": og, "reason": reason,
                     "intact": "NA" if ok is None else ok,
                     "min_ingroup_identity": "NA" if ident is None else round(ident, 5),
                     "codon_columns": cols, "tree_alignment": og in trimmed, "tree": in_tree})
    with open(out_quality, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["orthogroup"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with open(out_tree, "w") as fh:
        for r in rows:
            if r["tree"]:
                fh.write(r["orthogroup"] + "\n")
    n_tree = sum(r["tree"] for r in rows)
    n_aln = sum(r["tree_alignment"] for r in rows)
    print(f"{len(rows)} five taxon loci; {sum(r['intact'] is True for r in rows)} intact; "
          f"rooted tree ({loci}): {n_tree} of {n_aln} loci with a trimmed alignment")


def main():
    if "snakemake" not in globals():
        sys.exit("Run via Snakemake (see the docstring).")
    sm = globals()["snakemake"]
    gffs = {os.path.basename(g).split("_liftoff")[0]: g for g in sm.input.gffs}
    run(sm.input.aln, sm.input.of, gffs, list(sm.params.ingroup), sm.params.outgroup,
        sm.params.reference, sm.params.loci, sm.output.quality, sm.output.tree_loci)


if __name__ == "__main__":
    main()
