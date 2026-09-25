#!/usr/bin/env python3
"""normalize_gff.py - make a native (Ensembl) GFF3 look like the others.

    python normalize_gff.py in.gff3 out.gff3

Ensembl writes typed identifiers (ID=transcript:ENSX0001;Parent=gene:ENSX0002).
gffread copies the ID into the protein FASTA header, and a colon in a sequence
name breaks Newick gene trees downstream. This removes the type prefix from
every ID and Parent value (gene:, transcript:, CDS:, exon: and the UTR types)
and replaces any colon still left with an underscore. Nothing else changes:
coordinates, phases and all other attributes are copied as they are.
"""
import re
import sys

PREFIX = re.compile(r"^(?:gene|transcript|CDS|exon|five_prime_UTR|three_prime_UTR):")


def clean(value):
    return PREFIX.sub("", value).replace(":", "_")


def normalize_attributes(attrs):
    out = []
    for kv in attrs.split(";"):
        if "=" not in kv:
            out.append(kv)
            continue
        key, value = kv.split("=", 1)
        if key in ("ID", "Parent"):
            value = ",".join(clean(v) for v in value.split(","))
        out.append(f"{key}={value}")
    return ";".join(out)


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, dst = sys.argv[1:]
    n = 0
    with open(src) as fh, open(dst, "w") as out:
        for line in fh:
            if line.startswith("#") or not line.strip():
                out.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                out.write(line)
                continue
            f[8] = normalize_attributes(f[8])
            out.write("\t".join(f) + "\n")
            n += 1
    print(f"{src}: {n} features written to {dst}")


if __name__ == "__main__":
    main()
