#!/usr/bin/env python3
"""
quartet_asymmetry.py - which discordant gene tree topology is in excess, and
does the excess survive filtering on gene tree support and a site level count.

Why this test
-------------
With four taxa there is one internal branch and three possible unrooted
topologies. Under incomplete lineage sorting alone the two topologies that
disagree with the species tree are expected in EQUAL frequency. A significant
excess of one is evidence that something beyond ILS (introgression, or a
systematic error) produces part of the discordance. It is the quartet count
analogue of the ABBA/BABA test.

What an unrooted quartet cannot tell you: the excess topology pairs two taxa on
each side, so gene flow between either pair produces it. It detects departure
from ILS; it does not say which lineages exchanged genes.

Outputs
-------
  topology_counts   the three splits with counts, their role (species tree,
                    major and minor discordant) and, when concord.cf.stat is
                    given, which of IQ-TREE's gDF1 and gDF2 each one is
  by_support        the binomial test repeated on gene trees whose internal
                    branch has ultrafast bootstrap support at or above
                    0, 50, 70, 80, 90 and 95. Estimation error is concentrated
                    in poorly supported trees, so an excess that holds or grows
                    with support is not a product of noisy gene trees
  site_patterns     parsimony informative sites (two states, two taxa each) in
                    the trimmed alignments, counted per split, with the same
                    binomial test on the two discordant splits. Summed site
                    counts weight each locus by its number of informative
                    sites, so a few loci with hundreds of such sites (often
                    misaligned or misannotated) can outweigh thousands of
                    typical loci. The table therefore also gives each locus
                    one vote (the split with the most informative sites in
                    that locus, ties and loci without such sites set aside),
                    with the binomial test on the two discordant splits, and
                    reports how much of the site total the top 1% of loci hold

Usage
-----
Inside the workflow this runs as a Snakemake script. Standalone:

    python quartet_asymmetry.py --trees results/phylo/all_gene_trees.nwk \\
        --species-tree results/phylo/concat_tree.treefile \\
        [--cf-stat results/phylo/concord.cf.stat] \\
        [--alignments results/phylo/trimmed] [--outdir results/phylo]
"""

import argparse
import collections
import glob
import os
import re
import sys

try:
    from scipy.stats import binomtest
except ImportError:
    sys.exit("scipy required for the exact binomial test")

SUPPORT_THRESHOLDS = [0, 50, 70, 80, 90, 95]
GAP_CHARS = set("-X?*.BZJUO")
CHERRY = re.compile(r"\(([A-Za-z][A-Za-z0-9_.]*)(?::[0-9.eE+-]+)?,"
                    r"([A-Za-z][A-Za-z0-9_.]*)(?::[0-9.eE+-]+)?\)"
                    r"([0-9.]+)?(?:/[^:),;]*)?")
TAXON = re.compile(r"[(,]([A-Za-z][A-Za-z0-9_.]*)")


def parse_tree(newick):
    """(taxa, cherry pair, internal branch support) for a four taxon tree."""
    taxa = frozenset(TAXON.findall(newick))
    m = CHERRY.search(newick)
    if not m:
        return taxa, None, None
    pair = frozenset(m.group(1, 2))
    support = float(m.group(3)) if m.group(3) else None
    return taxa, (pair if len(pair) == 2 else None), support


def canonical(pair, taxa):
    """A split is the same whichever side names it; use the side holding the
    alphabetically first taxon."""
    other = frozenset(taxa - pair)
    return pair if min(taxa) in pair else other


def label(split, taxa):
    a = "+".join(sorted(split))
    b = "+".join(sorted(taxa - split))
    return f"{a} | {b}"


def binom(n1, n2):
    if n1 + n2 == 0:
        return float("nan"), float("nan"), float("nan")
    r = binomtest(n1, n1 + n2, 0.5)
    ci = r.proportion_ci(0.95)
    return r.pvalue, round(ci.low, 4), round(ci.high, 4)


def read_gene_trees(path):
    trees, taxa = [], set()
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        t, pair, sup = parse_tree(line)
        if pair is None:
            continue
        taxa |= t
        trees.append((pair, sup))
    return trees, frozenset(taxa)


def read_cf_stat(path):
    """{'gDF1_N': int, 'gDF2_N': int, ...} from IQ-TREE's concord.cf.stat."""
    header, row = None, None
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if header is None:
            header = f
        else:
            row = f
            break
    if not header or not row:
        return {}
    return dict(zip(header, row))


def site_patterns(aln_dir, taxa, per_locus=None):
    """Count parsimony informative site patterns per split over *.trim files.
    If per_locus is a list, one Counter per alignment is appended to it."""
    order = sorted(taxa)
    counts = collections.Counter()
    n_files = n_sites = 0
    for path in sorted(glob.glob(os.path.join(aln_dir, "*.trim"))):
        seqs, name = {}, None
        for line in open(path):
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
        seqs = {k: "".join(v).upper() for k, v in seqs.items()}
        if set(seqs) != set(order):
            continue
        n_files += 1
        local = collections.Counter()
        cols = zip(*(seqs[t] for t in order))
        for col in cols:
            n_sites += 1
            if any(c in GAP_CHARS for c in col):
                continue
            states = collections.Counter(col)
            if len(states) != 2 or sorted(states.values()) != [2, 2]:
                continue
            first = col[0]
            pair = frozenset(t for t, c in zip(order, col) if c == first)
            local[canonical(pair, taxa)] += 1
        counts.update(local)
        if per_locus is not None:
            per_locus.append(local)
    return counts, n_files, n_sites


def per_locus_rows(per_locus, splits, roles, taxa, top_share=0.01):
    """Rows giving each locus one vote, and the site share of the top loci."""
    conc, major, minor = splits
    votes = collections.Counter()
    for loc in per_locus:
        best = max(loc[s] for s in splits)
        if best == 0:
            votes["none"] += 1
            continue
        winners = [s for s in splits if loc[s] == best]
        votes[winners[0] if len(winners) == 1 else "tie"] += 1
    decided = sum(votes[s] for s in splits)
    rows = []
    for s, role in zip(splits, roles):
        rows.append({"split": label(s, taxa), "role": f"locus_majority_{role}",
                     "n_loci": votes[s],
                     "pct": round(100.0 * votes[s] / decided, 2) if decided else 0.0})
    m1, m2 = votes[major], votes[minor]
    p, lo, hi = binom(m1, m2)
    rows.append({"split": "locus majority: gene tree major vs minor", "role": "binomial_test",
                 "n_loci": m1 + m2,
                 "pct": round(100.0 * m1 / (m1 + m2), 2) if m1 + m2 else float("nan"),
                 "binomial_p": p, "ci95_low": lo, "ci95_high": hi})
    rows.append({"split": "loci without a majority split", "role": "locus_majority_undecided",
                 "n_loci": votes["tie"] + votes["none"],
                 "n_loci_no_informative_sites": votes["none"]})

    totals = sorted((sum(loc[s] for s in splits) for loc in per_locus), reverse=True)
    all_sites = sum(totals)
    k = max(1, round(top_share * len(per_locus)))
    ranked = sorted(per_locus, key=lambda loc: -sum(loc[s] for s in splits))[:k]
    top = collections.Counter()
    for loc in ranked:
        top.update({s: loc[s] for s in splits})
    top_sum = sum(top.values())
    median = (totals[len(totals) // 2] if len(totals) % 2 else
              (totals[len(totals) // 2 - 1] + totals[len(totals) // 2]) / 2) if totals else 0
    rows.append({"split": "informative sites per locus", "role": "median",
                 "n_loci": len(per_locus), "n_informative_sites": median})
    rows.append({"split": f"top {round(100 * top_share)}% of loci by informative sites",
                 "role": "concentration", "n_loci": k, "n_informative_sites": top_sum,
                 "pct": round(100.0 * top_sum / all_sites, 2) if all_sites else 0.0})
    for s, role in zip(splits, roles):
        rows.append({"split": label(s, taxa), "role": f"top_loci_{role}",
                     "n_loci": k, "n_informative_sites": top[s],
                     "pct": round(100.0 * top[s] / top_sum, 2) if top_sum else 0.0})
    return rows


def analyse(trees, taxa, sister, cf=None, site_counts=None, per_locus=None):
    """Return (topology rows, support rows, site rows) as lists of dicts."""
    if len(taxa) != 4:
        raise SystemExit(f"expected 4 taxa, found {len(taxa)}: {sorted(taxa)}")
    conc = canonical(frozenset(sister), taxa)
    splits = [conc] + sorted({canonical(frozenset(p), taxa) for p in
                              [frozenset(x) for x in _pairs(taxa)]} - {conc},
                             key=lambda s: label(s, taxa))

    counts = collections.Counter(canonical(p, taxa) for p, _ in trees)
    disc = sorted(splits[1:], key=lambda s: -counts[s])
    major, minor = disc
    total = sum(counts.values())

    cf_label = {}
    if cf:
        g1, g2 = int(float(cf.get("gDF1_N", -1))), int(float(cf.get("gDF2_N", -1)))
        for s in disc:
            if counts[s] == g1 and counts[s] != g2:
                cf_label[s] = "gDF1"
            elif counts[s] == g2 and counts[s] != g1:
                cf_label[s] = "gDF2"

    topo = []
    for s, role in [(conc, "species_tree"), (major, "major_discordant"),
                    (minor, "minor_discordant")]:
        topo.append({"split": label(s, taxa), "role": role, "n_gene_trees": counts[s],
                     "pct": round(100.0 * counts[s] / total, 2) if total else 0.0,
                     "iqtree_label": cf_label.get(s, "gCF" if s == conc else "")})

    support_rows = []
    for thr in SUPPORT_THRESHOLDS:
        sub = collections.Counter(canonical(p, taxa) for p, sup in trees
                                  if thr == 0 or (sup is not None and sup >= thr))
        n = sum(sub.values())
        n1, n2 = sub[major], sub[minor]
        p, lo, hi = binom(n1, n2)
        support_rows.append({
            "min_ufboot": thr, "n_gene_trees": n, "n_concordant": sub[conc],
            "pct_concordant": round(100.0 * sub[conc] / n, 2) if n else float("nan"),
            "major_split": label(major, taxa), "n_major": n1,
            "minor_split": label(minor, taxa), "n_minor": n2,
            "major_share_of_discordant": round(n1 / (n1 + n2), 4) if n1 + n2 else float("nan"),
            "binomial_p": p, "ci95_low": lo, "ci95_high": hi})

    site_rows = []
    if site_counts is not None:
        n_all = sum(site_counts.values())
        s1, s2 = site_counts[major], site_counts[minor]
        p, lo, hi = binom(s1, s2)
        for s, role in [(conc, "species_tree"), (major, "major_discordant_gene_trees"),
                        (minor, "minor_discordant_gene_trees")]:
            site_rows.append({"split": label(s, taxa), "role": role,
                              "n_informative_sites": site_counts[s],
                              "pct": round(100.0 * site_counts[s] / n_all, 2) if n_all else 0.0})
        site_rows.append({"split": "discordant sites: gene tree major vs minor",
                          "role": "binomial_test", "n_informative_sites": s1 + s2,
                          "pct": round(100.0 * s1 / (s1 + s2), 2) if s1 + s2 else float("nan"),
                          "binomial_p": p, "ci95_low": lo, "ci95_high": hi})
        if per_locus:
            site_rows.extend(per_locus_rows(
                per_locus, [conc, major, minor],
                ["species_tree", "major_discordant_gene_trees", "minor_discordant_gene_trees"],
                taxa))
    return topo, support_rows, site_rows


def _pairs(taxa):
    t = sorted(taxa)
    return [(t[0], t[1]), (t[0], t[2]), (t[0], t[3])]


def sister_from_species_tree(path):
    _, pair, _ = parse_tree(open(path).read().strip())
    if pair is None:
        raise SystemExit(f"could not read a cherry from {path}")
    return sorted(pair)


def write_tsv(rows, path):
    import csv
    if not rows:
        open(path, "w").close()
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", extrasaction="ignore",
                           lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.3e}" if k == "binomial_p" and isinstance(v, float)
                            else v) for k, v in r.items()})


def run(gene_trees, species_tree, cf_stat, aln_dir, out_topo, out_support, out_sites):
    trees, taxa = read_gene_trees(gene_trees)
    sister = sister_from_species_tree(species_tree)
    cf = read_cf_stat(cf_stat) if cf_stat and os.path.exists(cf_stat) else None
    site_counts, per_locus = None, []
    if aln_dir and os.path.isdir(aln_dir):
        site_counts, n_files, n_sites = site_patterns(aln_dir, taxa, per_locus)
        print(f"site patterns: {sum(site_counts.values())} informative of {n_sites} "
              f"columns in {n_files} alignments")
    topo, support, sites = analyse(trees, taxa, sister, cf, site_counts, per_locus)
    write_tsv(topo, out_topo)
    write_tsv(support, out_support)
    if out_sites:
        write_tsv(sites, out_sites)

    print(f"gene trees: {len(trees)}; species tree split: {label(frozenset(sister), taxa)}\n")
    for r in topo:
        print(f"  {r['split']:55s} {r['n_gene_trees']:6d} {r['pct']:6.2f}%  "
              f"{r['role']} {r['iqtree_label']}")
    print("\nby internal branch support")
    for r in support:
        print(f"  UFBoot>={r['min_ufboot']:3d}  n={r['n_gene_trees']:5d}  "
              f"major {r['n_major']:5d}  minor {r['n_minor']:5d}  "
              f"share {r['major_share_of_discordant']}  p={r['binomial_p']:.2e}")
    for r in sites:
        print(f"  sites  {r['split']:55s} {r['role']:40s} loci {r.get('n_loci', '')!s:>6} "
              f"sites {r.get('n_informative_sites', '')!s:>7} pct {r.get('pct', '')}")


def main():
    if "snakemake" in globals():
        sm = globals()["snakemake"]
        run(sm.params.gene_trees, sm.input.tree, sm.params.cf_stat,
            sm.params.trimmed_dir, sm.output.topology, sm.output.support,
            sm.output.sites)
        return
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trees", required=True, help="newick gene trees, one per line")
    ap.add_argument("--species-tree", required=True, help="concatenated ML tree")
    ap.add_argument("--cf-stat", help="IQ-TREE concord.cf.stat, to label gDF1 and gDF2")
    ap.add_argument("--alignments", help="directory of trimmed *.trim alignments")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    run(a.trees, a.species_tree, a.cf_stat, a.alignments,
        os.path.join(a.outdir, "quartet_topology_counts.tsv"),
        os.path.join(a.outdir, "quartet_asymmetry_by_support.tsv"),
        os.path.join(a.outdir, "quartet_site_patterns.tsv"))


if __name__ == "__main__":
    main()
