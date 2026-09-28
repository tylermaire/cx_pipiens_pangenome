#!/usr/bin/env python3
"""extract_sco5.py - single copy orthologs of all five taxa, protein and CDS.

The rooted analyses (the rooted species tree and the D statistics) need loci
with one gene in each ingroup form and one in the outgroup. For each such
orthogroup this writes

    <og>.faa   the five proteins, named by sample
    <og>.fna   the five coding sequences, in frame, with every stop codon
               removed, so that each protein is the exact translation of its
               CDS (codon i codes residue i)

The protein written is the translation of the CDS in the frame that best
reproduces the protein the pangenome used (results/proteins). The frame is
not always 0: a partial model can start mid codon, and gffread -x ignores the
phase. Stop codons are removed because longest_isoform.py removes stop
symbols (including in frame stops of broken transferred models) from the
proteins. A sequence whose translation matches its protein at less than 95%
identity in every frame removes the locus, and is counted.

Snakemake provides:
    input.of        OrthoFinder output directory
    input.proteins  results/proteins/<sample>.fa for every sample
    input.cds       results/cds/<sample>.fa for every sample
    params.samples  every sample, ingroup first, outgroup last
    params.reference  the Liftoff reference (its gene id locates each locus)
    output.fastas   directory for <og>.faa and <og>.fna
    output.loci     one row per locus written: og and the gene of each sample
    output.accounting  counts of candidate, written and removed loci
"""
import csv
import difflib
import glob
import itertools
import os
import sys

BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: AMINO[16 * i + 4 * j + k]
               for i, a in enumerate(BASES) for j, b in enumerate(BASES)
               for k, c in enumerate(BASES)}
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T", "R": "AG", "Y": "CT",
         "S": "CG", "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
         "H": "ACT", "V": "ACG", "N": "ACGT"}
MIN_IDENTITY = 0.95


def translate_codon(codon, _cache={}):
    """Standard code; an ambiguous codon translates when every reading of it
    gives one amino acid (GCN is A), otherwise X."""
    codon = codon.upper()
    if codon in CODON_TABLE:
        return CODON_TABLE[codon]
    if codon in _cache:
        return _cache[codon]
    options = [IUPAC.get(b) for b in codon]
    if len(codon) != 3 or any(o is None for o in options):
        aa = "X"
    else:
        aas = {CODON_TABLE["".join(p)] for p in itertools.product(*options)}
        aa = aas.pop() if len(aas) == 1 else "X"
    _cache[codon] = aa
    return aa


def in_frame(cds, frame):
    """(protein, codons) for one frame, stop codons and a trailing partial
    codon removed."""
    codons = [cds[i:i + 3] for i in range(frame, len(cds) - 2, 3)]
    keep = [(translate_codon(c), c) for c in codons]
    keep = [(a, c) for a, c in keep if a != "*"]
    return "".join(a for a, _ in keep), [c for _, c in keep]


def identity(a, b):
    if a == b:
        return 1.0
    if not a or not b:
        return 0.0
    sm = difflib.SequenceMatcher(None, a, b, autojunk=False)
    same = sum(block.size for block in sm.get_matching_blocks())
    return same / max(len(a), len(b))


def best_frame(cds, protein):
    """(identity, frame, protein, codons) of the frame closest to protein."""
    best = None
    for frame in (0, 1, 2):
        prot, codons = in_frame(cds, frame)
        ident = identity(prot, protein)
        if best is None or ident > best[0]:
            best = (ident, frame, prot, codons)
        if ident == 1.0:
            break
    return best


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def newest(paths):
    return max(paths, key=os.path.getmtime) if paths else None


def orthogroup_tables(of_dir):
    """Orthogroups.tsv and Orthogroups.GeneCount.tsv of the newest result."""
    counts = newest(glob.glob(os.path.join(of_dir, "**", "Orthogroups.GeneCount.tsv"),
                              recursive=True))
    if counts is None:
        sys.exit(f"no Orthogroups.GeneCount.tsv under {of_dir}")
    members = os.path.join(os.path.dirname(counts), "Orthogroups.tsv")
    if not os.path.exists(members):
        sys.exit(f"no Orthogroups.tsv beside {counts}")
    return members, counts


def wrap(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def run(of_dir, samples, protein_files, cds_files, outdir, loci_out, accounting_out):
    os.makedirs(outdir, exist_ok=True)
    members, counts = orthogroup_tables(of_dir)
    single = []
    with open(counts) as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        missing = [s for s in samples if s not in rows.fieldnames]
        if missing:
            sys.exit(f"samples missing from {counts}: {missing}")
        for r in rows:
            if all(r[s] == "1" for s in samples):
                single.append(r["Orthogroup"])
    single_set = set(single)
    genes = {}
    with open(members) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["Orthogroup"] in single_set:
                genes[r["Orthogroup"]] = {s: r[s].strip() for s in samples}
    proteins = {s: read_fasta(protein_files[s]) for s in samples}
    cds = {s: read_fasta(cds_files[s]) for s in samples}
    tally = {"single_copy_in_all_taxa": len(single), "written": 0,
             "removed_no_protein": 0, "removed_no_cds": 0, "removed_cds_mismatch": 0}
    per_sample = {s: {"exact": 0, "frame_not_0": 0, "below_identity": 0} for s in samples}
    loci_rows = []
    for og in single:
        record, reason = {}, None
        for s in samples:
            gene = genes.get(og, {}).get(s, "")
            prot = proteins[s].get(gene)
            if prot is None:
                reason = "removed_no_protein"
                break
            nt = cds[s].get(gene)
            if nt is None:
                reason = "removed_no_cds"
                break
            prot = prot.upper().replace("*", "").replace(".", "")
            ident, frame, trans, codons = best_frame(nt.upper(), prot)
            if ident < MIN_IDENTITY:
                per_sample[s]["below_identity"] += 1
                reason = "removed_cds_mismatch"
                break
            per_sample[s]["exact"] += ident == 1.0
            per_sample[s]["frame_not_0"] += frame != 0
            record[s] = (gene, trans, "".join(codons))
        if reason:
            tally[reason] += 1
            continue
        with open(os.path.join(outdir, f"{og}.faa"), "w") as fa, \
                open(os.path.join(outdir, f"{og}.fna"), "w") as fn:
            for s in samples:
                gene, trans, nt = record[s]
                fa.write(f">{s}\n{wrap(trans)}\n")
                fn.write(f">{s}\n{wrap(nt)}\n")
        loci_rows.append({"orthogroup": og, **{s: record[s][0] for s in samples}})
        tally["written"] += 1
    with open(loci_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["orthogroup"] + samples, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(loci_rows)
    with open(accounting_out, "w") as fh:
        fh.write("item\tsample\tvalue\n")
        for k, v in tally.items():
            fh.write(f"{k}\t\t{v}\n")
        for s in samples:
            for k, v in per_sample[s].items():
                fh.write(f"cds_translation_{k}\t{s}\t{v}\n")
    print(f"{len(single)} orthogroups single copy in all {len(samples)} taxa; "
          f"{tally['written']} written to {outdir}")
    for k in ("removed_no_protein", "removed_no_cds", "removed_cds_mismatch"):
        print(f"  {k}: {tally[k]}")
    return tally


def main():
    if "snakemake" not in globals():
        sys.exit("Run via Snakemake (see the docstring).")
    sm = globals()["snakemake"]
    samples = list(sm.params.samples)

    def by_sample(paths):
        out = {os.path.basename(p).rsplit(".", 1)[0]: p for p in paths}
        missing = [s for s in samples if s not in out]
        if missing:
            sys.exit(f"no file for {missing} among {list(paths)}")
        return out

    run(sm.input.of, samples, by_sample(sm.input.proteins), by_sample(sm.input.cds),
        sm.output.fastas, sm.output.loci, sm.output.accounting)


if __name__ == "__main__":
    main()
