#!/usr/bin/env python3
"""
transfer_quality.py - how many of the protein models kept for analysis are
broken.

extract_proteins translates every transcript with gffread and then strips the
'.' stop symbols from the sequences, so a transferred model whose reading frame
shifted (an indel in the target that Liftoff did not repair; the workflow runs
Liftoff without -polish) enters OrthoFinder as a chimera of two frames with no
sign that anything is wrong. This script re-reads the unstripped translation
for the models that were kept and counts:

  internal stop     an in frame stop codon before the last residue
  no start Met      the model does not begin with methionine
  no terminal stop  the model does not end in a stop codon
  complete          none of the three

Usage:
    transfer_quality.py --raw gffread.faa --kept proteins.fa --sample NAME --out out.tsv
"""

import argparse


def read_fasta(path, strip_ids):
    seqs, rid, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if rid is not None:
                    seqs[rid] = "".join(chunks)
                rid = line[1:].split()[0]
                if strip_ids:
                    rid = rid.replace(".", "")
                chunks = []
            else:
                chunks.append(line.strip())
    if rid is not None:
        seqs[rid] = "".join(chunks)
    return seqs


def assess(raw, kept_ids):
    n = internal = no_met = no_stop = complete = 0
    for tid in kept_ids:
        seq = raw.get(tid)
        if seq is None:
            continue
        n += 1
        has_stop = seq.endswith(".") or seq.endswith("*")
        body = seq[:-1] if has_stop else seq
        bad_internal = ("." in body) or ("*" in body)
        starts_met = seq[:1].upper() == "M"
        internal += bad_internal
        no_met += not starts_met
        no_stop += not has_stop
        complete += (not bad_internal) and starts_met and has_stop
    pct = lambda k: round(100.0 * k / n, 3) if n else 0.0
    return {"n_models": n,
            "n_internal_stop": internal, "pct_internal_stop": pct(internal),
            "n_no_start_met": no_met, "pct_no_start_met": pct(no_met),
            "n_no_terminal_stop": no_stop, "pct_no_terminal_stop": pct(no_stop),
            "n_complete": complete, "pct_complete": pct(complete)}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw", required=True, help="gffread -y output, stops kept")
    ap.add_argument("--kept", required=True, help="results/proteins/{sample}.fa")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    raw = read_fasta(a.raw, strip_ids=True)
    kept = read_fasta(a.kept, strip_ids=False)
    row = {"sample": a.sample, **assess(raw, kept.keys())}
    row["n_kept_not_found"] = len(kept) - row["n_models"]
    with open(a.out, "w") as fh:
        fh.write("\t".join(row) + "\n")
        fh.write("\t".join(str(v) for v in row.values()) + "\n")
    print("\t".join(f"{k}={v}" for k, v in row.items()))


if __name__ == "__main__":
    main()
