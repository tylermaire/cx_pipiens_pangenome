#!/usr/bin/env python3
"""
anvio_inputs.py - anvi'o contigs and external gene calls from the workflow's
gene sets, so that anvi'o clusters exactly the proteins OrthoFinder used.

For each ingroup form, every gene (one per gene, the longest isoform, as in
results/proteins/<form>.fa) becomes one contig holding its coding sequence
(results/cds/<form>.fa), with one external gene call covering the contig and
carrying the protein itself in the aa_sequence column. anvi'o then uses these
proteins as given instead of translating: transferred models whose coding
sequence is not a multiple of three (a frameshift in transfer) keep the
protein gffread made from them, and are flagged partial. A few transferred
models are only a few nucleotides long; anvi'o refuses contigs shorter than
its k-mer size (4), so those contigs are padded with N after the coding
sequence (the gene call still covers the coding sequence only). Contig names
are replaced by simple names (anvi'o accepts letters, digits and
underscores); the map back to the transcript identifiers is written beside
them.

Outputs, in --out:
    <form>.contigs.fa        one contig per gene
    <form>.gene_calls.tsv    external gene calls with aa_sequence
    <form>.gene_map.tsv      contig, gene_callers_id, transcript, partial
    external_genomes.txt     name and contigs database path for each form

Usage (repository root, results/ holding the run):
    python workflow/scripts/anvio_inputs.py --out results/anvio/inputs \\
        --db-dir results/anvio
"""
import argparse
import os
import sys

FORMS = ["Cx_quinquefasciatus", "Cx_pallens", "Cx_molestus", "Cx_pipiens"]
ABBR = {"Cx_quinquefasciatus": "qui", "Cx_pallens": "pal", "Cx_molestus": "mol",
        "Cx_pipiens": "pip"}
MIN_CONTIG = 4          # anvi-gen-contigs-database default k-mer size


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def write_form(form, proteins, cds, out):
    missing = [t for t in proteins if t not in cds]
    if missing:
        sys.exit(f"{form}: {len(missing)} proteins without a coding sequence, e.g. {missing[:3]}")
    fa = open(os.path.join(out, f"{form}.contigs.fa"), "w")
    gc = open(os.path.join(out, f"{form}.gene_calls.tsv"), "w")
    gm = open(os.path.join(out, f"{form}.gene_map.tsv"), "w")
    gc.write("gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tcall_type\tsource\t"
             "version\taa_sequence\n")
    gm.write("contig\tgene_callers_id\ttranscript\tpartial\n")
    n_partial = n_padded = 0
    for i, (tx, aa) in enumerate(proteins.items()):
        contig = f"{ABBR[form]}_{i + 1:06d}"
        nt = cds[tx].upper()
        aa = aa.replace("*", "").replace(".", "")
        if not aa:
            sys.exit(f"{form}: empty protein for {tx}")
        partial = int(len(nt) % 3 != 0)
        n_partial += partial
        seq = nt + "N" * max(0, MIN_CONTIG - len(nt))
        n_padded += len(seq) > len(nt)
        fa.write(f">{contig}\n")
        for j in range(0, len(seq), 80):
            fa.write(seq[j:j + 80] + "\n")
        gc.write(f"{i}\t{contig}\t0\t{len(nt)}\tf\t{partial}\t1\tworkflow\tv5\t{aa}\n")
        gm.write(f"{contig}\t{i}\t{tx}\t{partial}\n")
    for fh in (fa, gc, gm):
        fh.close()
    return len(proteins), n_partial, n_padded


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--proteins", default="results/proteins")
    ap.add_argument("--cds", default="results/cds")
    ap.add_argument("--out", default="results/anvio/inputs")
    ap.add_argument("--db-dir", default="results/anvio",
                    help="where the contigs databases will be written")
    ap.add_argument("--forms", default=",".join(FORMS))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    forms = [f.strip() for f in args.forms.split(",") if f.strip()]
    total = 0
    with open(os.path.join(args.out, "external_genomes.txt"), "w") as eg:
        eg.write("name\tcontigs_db_path\n")
        for form in forms:
            proteins = read_fasta(os.path.join(args.proteins, f"{form}.fa"))
            cds = read_fasta(os.path.join(args.cds, f"{form}.fa"))
            n, n_partial, n_padded = write_form(form, proteins, cds, args.out)
            total += n
            eg.write(f"{form}\t{os.path.abspath(os.path.join(args.db_dir, form + '.db'))}\n")
            print(f"{form}: {n:,} genes, {n_partial:,} with a coding sequence not a multiple of three, "
                  f"{n_padded:,} contigs padded to {MIN_CONTIG} nt")
    print(f"total {total:,} genes")


if __name__ == "__main__":
    main()
