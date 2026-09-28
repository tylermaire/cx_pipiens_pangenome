#!/usr/bin/env python3
"""codon_align.py - codon alignments from protein alignments, trimmed alike.

    python codon_align.py --sco DIR --alignments DIR --out DIR --summary TSV

For every locus, DIR/<og>.aln is the MAFFT alignment of <og>.faa and
DIR/<og>.cols the column map trimAl printed (-colnumbering) when it trimmed
that alignment. <og>.fna in --sco holds each sequence's codons in frame,
codon i coding residue i of <og>.faa (extract_sco5.py). Each aligned residue
is replaced by its codon and each gap by '---', and only the protein columns
trimAl kept are written, so the codon alignment and the trimmed protein
alignment hold the same columns.

A locus is skipped, and counted, when a sequence in the alignment is not the
translation it should be (the alignment and the CDS disagree), when trimAl
kept no column, or when files are missing.
"""
import argparse
import glob
import os
import sys


def read_fasta(path):
    seqs, order, name = {}, [], None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                order.append(name)
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}, order


def read_columns(path):
    """Kept column indices (0 based) from trimAl's '#ColumnsMap' line."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#ColumnsMap"):
                text = line.split("\t", 1)[1] if "\t" in line else line[len("#ColumnsMap"):]
                return [int(x) for x in text.replace(",", " ").split()]
    return []


def back_translate(aligned, codons):
    """Codon string for an aligned protein, or None when they disagree."""
    residues = aligned.replace("-", "")
    if len(codons) != 3 * len(residues):
        return None
    out, i = [], 0
    for aa in aligned:
        if aa == "-":
            out.append("---")
        else:
            out.append(codons[3 * i:3 * i + 3])
            i += 1
    return out


def codon_alignment(aln_path, cols_path, fna_path, faa_path=None):
    """({name: trimmed codon string}, order) or (None, reason)."""
    aln, order = read_fasta(aln_path)
    nt, _ = read_fasta(fna_path)
    faa = read_fasta(faa_path)[0] if faa_path and os.path.exists(faa_path) else None
    keep = read_columns(cols_path)
    if not keep:
        return None, "no_columns_kept"
    width = len(next(iter(aln.values()))) if aln else 0
    if not aln or any(len(s) != width for s in aln.values()):
        return None, "ragged_alignment"
    if max(keep) >= width:
        return None, "column_map_out_of_range"
    out = {}
    for name in order:
        aligned = aln[name].upper()
        if name not in nt:
            return None, "missing_cds"
        if faa is not None and aligned.replace("-", "") != faa.get(name, "").upper():
            return None, "alignment_not_translation"
        codons = back_translate(aligned, nt[name].upper())
        if codons is None:
            return None, "cds_length_mismatch"
        out[name] = "".join(codons[c] for c in keep)
    return out, order


def run(sco_dir, aln_dir, out_dir, summary):
    os.makedirs(out_dir, exist_ok=True)
    loci = sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(aln_dir, "*.aln")))
    tally = {"alignments": len(loci), "written": 0}
    for og in loci:
        cols = os.path.join(aln_dir, f"{og}.cols")
        fna = os.path.join(sco_dir, f"{og}.fna")
        if not os.path.exists(cols) or not os.path.exists(fna):
            tally["missing_files"] = tally.get("missing_files", 0) + 1
            continue
        seqs, order = codon_alignment(os.path.join(aln_dir, f"{og}.aln"), cols, fna,
                                      os.path.join(sco_dir, f"{og}.faa"))
        if seqs is None:
            tally[order] = tally.get(order, 0) + 1
            continue
        with open(os.path.join(out_dir, f"{og}.fna"), "w") as fh:
            for name in order:
                fh.write(f">{name}\n{seqs[name]}\n")
        tally["written"] += 1
    with open(summary, "w") as fh:
        fh.write("item\tvalue\n")
        for k, v in tally.items():
            fh.write(f"{k}\t{v}\n")
    print(f"codon alignments: {tally}")
    if loci and tally["written"] == 0:
        sys.exit("no codon alignment was written")
    return tally


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sco", required=True, help="extract_sco5.py output (<og>.faa, <og>.fna)")
    ap.add_argument("--alignments", required=True, help="<og>.aln and <og>.cols")
    ap.add_argument("--out", required=True, help="directory for <og>.fna codon alignments")
    ap.add_argument("--summary", required=True)
    a = ap.parse_args()
    run(a.sco, a.alignments, a.out, a.summary)


if __name__ == "__main__":
    main()
