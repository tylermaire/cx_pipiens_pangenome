#!/usr/bin/env python3
"""filter_cds.py - keep the coding sequences of the proteins used downstream.

    python filter_cds.py --cds all_cds.fa --proteins proteins.fa --out cds.fa \
        [--strip-periods]

gffread -x writes the spliced coding sequence of every transcript. The
protein files keep one isoform per gene (longest_isoform.py) and, with
--strip-periods, drop the periods from identifiers. This applies the same
identifier change to the CDS headers and keeps exactly the transcripts that
the protein file holds, so a gene's protein and its CDS share one name.
"""
import argparse


def read_fasta(path, strip_periods=False):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                if strip_periods:
                    name = name.replace(".", "")
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cds", required=True)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--strip-periods", action="store_true")
    a = ap.parse_args()
    keep = {name for name, _ in read_fasta(a.proteins)}
    n = 0
    with open(a.out, "w") as out:
        for name, seq in read_fasta(a.cds, a.strip_periods):
            if name in keep:
                out.write(f">{name}\n")
                seq = seq.upper()
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i + 60] + "\n")
                n += 1
    print(f"{a.out}: {n} of {len(keep)} proteins have a coding sequence")


if __name__ == "__main__":
    main()
