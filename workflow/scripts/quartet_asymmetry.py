#!/usr/bin/env python3
"""
quartet_asymmetry.py - test whether the two discordant gene-tree topologies
occur at equal frequency.

Why this test
-------------
Under incomplete lineage sorting alone, the two topologies that disagree with
the species tree are expected in EQUAL frequency. A significant excess of one
over the other is the signal that introgression, not ILS alone, is producing
the discordance. This is the quartet-count analogue of the ABBA/BABA test.

Runs on files already in the repository; nothing needs recomputing:
    results/phylo/all_gene_trees.nwk

Usage
-----
    python quartet_asymmetry.py \
        --trees results/phylo/all_gene_trees.nwk \
        --sister Cx_pallens,Cx_quinquefasciatus
"""

import argparse
import collections
import re
import sys

try:
    from scipy.stats import binomtest
except ImportError:
    binomtest = None


def parse_cherry(newick):
    """Return the frozenset of the first two-taxon clade in a 4-taxon tree."""
    t = re.sub(r":[0-9eE.+-]+", "", newick)      # branch lengths
    t = re.sub(r"\)[0-9.]+", ")", t)             # support values
    m = re.findall(r"\(([A-Za-z0-9_.-]+),([A-Za-z0-9_.-]+)\)", t)
    if not m:
        return None
    pair = frozenset(m[0])
    return pair if len(pair) == 2 else None


def canonical(pair, taxa):
    """A split is the same whichever side you name it from; pick one name."""
    other = frozenset(taxa - pair)
    return min([pair, other], key=lambda s: sorted(s)[0])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trees", required=True, help="newick gene trees, one per line")
    ap.add_argument("--sister", required=True,
                    help="comma-separated sister pair in the species tree")
    args = ap.parse_args()

    sister = frozenset(args.sister.split(","))
    if len(sister) != 2:
        sys.exit("--sister needs exactly two comma-separated names")

    counts, unparsed = collections.Counter(), 0
    taxa = set()
    for line in open(args.trees):
        line = line.strip()
        if not line:
            continue
        pair = parse_cherry(line)
        if pair is None:
            unparsed += 1
            continue
        taxa |= set(re.findall(r"[A-Za-z][A-Za-z0-9_.-]*", re.sub(r":[0-9eE.+-]+", "", line)))
        counts[pair] += 1

    if len(taxa) != 4:
        sys.exit(f"expected 4 taxa, found {len(taxa)}: {sorted(taxa)}")

    splits = collections.Counter()
    for pair, n in counts.items():
        splits[canonical(pair, taxa)] += n
    total = sum(splits.values())

    def label(key):
        a = "+".join(sorted(key))
        b = "+".join(sorted(taxa - key))
        return f"({a}) | ({b})"

    print(f"gene trees: {total} parsed, {unparsed} unparsed\n")
    print("topology frequencies")
    for key, n in splits.most_common():
        tag = "  <- species tree" if key == canonical(sister, taxa) else ""
        print(f"  {label(key):55s} {n:6d}  {100*n/total:6.2f}%{tag}")

    conc = canonical(sister, taxa)
    disc = sorted(((k, v) for k, v in splits.items() if k != conc),
                  key=lambda kv: -kv[1])
    (k1, n1), (k2, n2) = disc
    d = (n1 - n2) / (n1 + n2)

    print("\nsymmetry test on the two discordant topologies")
    print(f"  major discordant  {label(k1):55s} {n1}")
    print(f"  minor discordant  {label(k2):55s} {n2}")
    print(f"  asymmetry (n1-n2)/(n1+n2) = {d:.4f}")
    if binomtest is not None:
        r = binomtest(n1, n1 + n2, 0.5)
        ci = r.proportion_ci(0.95)
        print(f"  exact binomial: n = {n1+n2}, p = {r.pvalue:.3e}, "
              f"95% CI {ci.low:.4f} to {ci.high:.4f}")
        print("\n  Equal frequencies are the ILS expectation. Rejecting equality is "
              "evidence\n  for introgression, but an unrooted quartet cannot say which "
              "pair exchanged\n  genes: the excess topology groups two non-sister pairs "
              "at once.")
    else:
        print("  install scipy for the exact binomial test")


if __name__ == "__main__":
    main()
