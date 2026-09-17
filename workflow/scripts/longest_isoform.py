#!/usr/bin/env python3
"""
longest_isoform.py - reduce a gffread protein FASTA to one sequence per gene.

Background
----------
`gffread -y` emits one protein per *transcript*. Where the GFF3 carries
multiple mRNAs for a gene, every isoform enters the FASTA, so downstream
orthology and gene-family analyses treat isoforms as separate genes. This
script keeps only the longest protein per gene.

Gene membership is taken from the GFF3 itself: each mRNA's ID is mapped to
its Parent gene. Sequences whose transcript ID is not found in the GFF are
kept as-is rather than silently dropped.

Usage
-----
    python longest_isoform.py --gff ann.gff3 --fasta proteins_all.fa \
        --out proteins.fa [--report counts.tsv] [--strip-periods]

--strip-periods reproduces the `sed 's/\\.//g'` step in the current
extract_proteins rule. Apply it here instead of afterwards so that FASTA
headers and GFF IDs are normalised the same way.
"""

import argparse
import gzip
import re
import sys
from collections import defaultdict

ATTR_ID = re.compile(r"(?:^|;)ID=([^;]+)")
ATTR_PARENT = re.compile(r"(?:^|;)Parent=([^;,]+)")


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def transcript_to_gene(gff_path, strip_periods):
    """Map every mRNA/transcript ID to its parent gene ID."""
    mapping = {}
    feature_types = {"mRNA", "transcript", "V_gene_segment", "C_gene_segment"}
    with opener(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] not in feature_types:
                continue
            tid = ATTR_ID.search(parts[8])
            gid = ATTR_PARENT.search(parts[8])
            if not tid:
                continue
            t = tid.group(1)
            g = gid.group(1) if gid else t
            if strip_periods:
                t, g = t.replace(".", ""), g.replace(".", "")
            mapping[t] = g
    return mapping


def read_fasta(path, strip_periods):
    """Yield (header_line, sequence) pairs."""
    header, chunks = None, []
    with opener(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line.replace(".", "") if strip_periods else line
                chunks = []
            else:
                chunks.append(line.replace(".", "") if strip_periods else line)
    if header is not None:
        yield header, "".join(chunks)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True, help="GFF3 used to extract the proteins")
    ap.add_argument("--fasta", required=True, help="gffread -y output (all isoforms)")
    ap.add_argument("--out", required=True, help="filtered FASTA, one protein per gene")
    ap.add_argument("--report", help="optional TSV of input/output counts")
    ap.add_argument("--strip-periods", action="store_true",
                    help="remove '.' from IDs and sequences, as the current rule does")
    args = ap.parse_args()

    t2g = transcript_to_gene(args.gff, args.strip_periods)
    if not t2g:
        sys.exit(f"No mRNA features parsed from {args.gff}; check the GFF3 feature types.")

    best = {}          # gene -> (length, header, sequence)
    unmapped = []      # sequences whose transcript ID was not in the GFF
    n_in = 0

    for header, seq in read_fasta(args.fasta, args.strip_periods):
        n_in += 1
        tid = header[1:].split()[0]
        gene = t2g.get(tid)
        if gene is None:
            unmapped.append((header, seq))
            continue
        length = len(seq.replace("*", ""))
        # tie-break on transcript ID so the choice is deterministic across runs
        if gene not in best or (length, tid) > (best[gene][0], best[gene][3]):
            best[gene] = (length, header, seq, tid)

    with open(args.out, "w") as out:
        for gene in sorted(best):
            _, header, seq, _ = best[gene]
            out.write(header + "\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")
        for header, seq in unmapped:
            out.write(header + "\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")

    n_out = len(best) + len(unmapped)
    msg = (f"{args.fasta}: {n_in} proteins in, {n_out} out "
           f"({len(best)} genes, {len(unmapped)} unmapped kept), "
           f"{n_in - n_out} isoforms dropped")
    print(msg, file=sys.stderr)

    if args.report:
        with open(args.report, "w") as rep:
            rep.write("file\tproteins_in\tgenes\tunmapped\tproteins_out\tisoforms_dropped\n")
            rep.write(f"{args.fasta}\t{n_in}\t{len(best)}\t{len(unmapped)}\t"
                      f"{n_out}\t{n_in - n_out}\n")


if __name__ == "__main__":
    main()
