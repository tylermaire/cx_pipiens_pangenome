#!/usr/bin/env python3
"""
skani_matrix.py - pairwise whole-genome ANI across the five-taxon panel.

The ANI matrix is cited in the Results but was previously produced outside the
workflow, with two versions in the repository disagreeing by 0.2-0.3% per cell
and no record of which parameters produced which. This rule replaces both.

skani is run once per pair with explicit parameters (recorded in the output
header) rather than relying on defaults that may change between versions.
Pairs for which skani reports no value are written as NA rather than dropped:
an outgroup can be divergent enough that genome-wide identity falls below the
threshold at which skani will report (Cx. tarsalis, the outgroup up to V4,
was), and that is a result worth showing.

Snakemake provides:
    input.genomes  one FASTA per sample
    params.samples, params.extra (skani flags)
    output.matrix  square TSV, 100.00 on the diagonal, NA where unreported
    output.long    one row per pair, with the raw skani fields retained
"""

import itertools
import os
import subprocess
import sys


def run_skani(a_path, b_path, extra):
    """Return (ani, af_query, af_ref) for one pair, or (None, None, None)."""
    cmd = ["skani", "dist", a_path, b_path] + list(extra)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print(f"    skani failed: {proc.stderr.strip()[:200]}", file=sys.stderr)
        return None, None, None
    lines = [l for l in proc.stdout.splitlines() if l.strip()]
    if len(lines) < 2:            # header only: no alignment above threshold
        return None, None, None
    header = lines[0].split("\t")
    fields = lines[1].split("\t")
    row = dict(zip(header, fields))

    def num(*names):
        for n in names:
            if n in row:
                try:
                    return float(row[n])
                except ValueError:
                    return None
        return None

    return num("ANI"), num("Align_fraction_query"), num("Align_fraction_ref")


def main():
    if "snakemake" not in globals():
        sys.exit("Run via Snakemake.")

    samples = list(snakemake.params.samples)
    extra = list(getattr(snakemake.params, "extra", []))
    paths = {}
    for p in snakemake.input.genomes:
        s = os.path.basename(p).replace(".fasta", "")
        paths[s] = p
    missing = [s for s in samples if s not in paths]
    if missing:
        sys.exit(f"No genome FASTA for: {missing}")

    ani = {}
    rows = []
    for a, b in itertools.combinations(samples, 2):
        print(f"  skani {a} vs {b}", flush=True)
        val, afq, afr = run_skani(paths[a], paths[b], extra)
        ani[(a, b)] = ani[(b, a)] = val
        rows.append({"sample1": a, "sample2": b,
                     "ani": "NA" if val is None else f"{val:.2f}",
                     "align_fraction_query": "" if afq is None else f"{afq:.4f}",
                     "align_fraction_ref": "" if afr is None else f"{afr:.4f}"})
        print(f"    ANI = {'NA (below reporting threshold)' if val is None else f'{val:.4f}'}")

    with open(snakemake.output.matrix, "w") as out:
        out.write("sample\t" + "\t".join(samples) + "\n")
        for a in samples:
            cells = []
            for b in samples:
                if a == b:
                    cells.append("100.00")
                else:
                    v = ani.get((a, b))
                    cells.append("NA" if v is None else f"{v:.4f}")
            out.write(a + "\t" + "\t".join(cells) + "\n")

    with open(snakemake.output.long, "w") as out:
        out.write(f"# skani dist {' '.join(extra)}\n")
        out.write("sample1\tsample2\tani\talign_fraction_query\talign_fraction_ref\n")
        for r in rows:
            out.write("\t".join([r["sample1"], r["sample2"], r["ani"],
                                 r["align_fraction_query"],
                                 r["align_fraction_ref"]]) + "\n")

    n_na = sum(1 for r in rows if r["ani"] == "NA")
    print(f"\n{len(rows)} pairs, {n_na} below skani's reporting threshold")
    if n_na:
        print("NA pairs are reported as such rather than omitted; genome-wide")
        print("identity below roughly 80% is outside the range skani reports.")


if __name__ == "__main__":
    main()
