#!/usr/bin/env python3
"""check_seqids.py - stop the run when a GFF names sequences its FASTA lacks.

    python check_seqids.py genome.fasta genes.gff3

A gene set and a genome taken from different sources can name the same
sequence differently (chromosome names against INSDC accessions). gffread
then silently extracts nothing for the unmatched sequences, so the mismatch
is checked here, where it is cheap, instead of surfacing as a small proteome.
Exit status 1 when fewer than 99% of the GFF's features lie on sequences in
the FASTA, or when the GFF has no coding sequence.
"""
import sys


def fasta_names(path):
    with open(path) as fh:
        return {line[1:].split()[0] for line in fh if line.startswith(">")}


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    fasta, gff = sys.argv[1:]
    names = fasta_names(fasta)
    n_feat = n_on = n_cds = 0
    missing = set()
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 9:
                continue
            n_feat += 1
            if f[2] == "CDS":
                n_cds += 1
            if f[0] in names:
                n_on += 1
            else:
                missing.add(f[0])
    share = n_on / n_feat if n_feat else 0.0
    print(f"{fasta}: {len(names)} sequences; {gff}: {n_feat} features, {n_cds} CDS, "
          f"{100 * share:.2f}% on sequences in the FASTA")
    if missing:
        print(f"  sequences named in the GFF but absent from the FASTA ({len(missing)}): "
              f"{', '.join(sorted(missing)[:10])}")
    if not n_cds:
        sys.exit("the GFF has no CDS features")
    if share < 0.99:
        sys.exit("the GFF and the FASTA name their sequences differently; take both from "
                 "the same source")


if __name__ == "__main__":
    main()
