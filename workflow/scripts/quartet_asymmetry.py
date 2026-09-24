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
                    reports how much of the site total the top 1% of loci hold,
                    and how many of those loci, against the rest, contain a
                    transferred model without a valid ORF
  robustness        the gene tree test on subsets that remove the two obvious
                    sources of spurious discordance: gene trees whose internal
                    branch sits at IQ-TREE's minimum length (no substitution
                    supports any resolution, so the topology is arbitrary),
                    and loci with a broken gene model (a transferred model
                    without a valid ORF, or any model whose reference source is
                    partial or carries a RefSeq sequence exception)

Usage
-----
Inside the workflow this runs as a Snakemake script. Standalone:

    python quartet_asymmetry.py --trees results/phylo/all_gene_trees.nwk \\
        --species-tree results/phylo/concat_tree.treefile \\
        [--cf-stat results/phylo/concord.cf.stat] \\
        [--alignments results/phylo/trimmed] [--outdir results/phylo] \\
        [--tree-ids results/phylo/gene_tree_ids.txt \\
         --orthogroups results/orthofinder/output \\
         --gffs results/annotation/<form>_liftoff.gff3 ... --reference <form>]

The robustness table needs the locus of every gene tree. all_gene_trees.nwk is
the concatenation of results/phylo/gene_trees/*.treefile in sorted order, so
the loci are read from those file names, or from gene_tree_ids.txt (the same
listing, shipped with patch runs instead of the IQ-TREE files).
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
INTERNAL_LENGTH = re.compile(r":([0-9.eE+-]+)")
MIN_BRANCH = 1.1e-6            # IQ-TREE floors branch lengths at 1e-6
ATTR_ID = re.compile(r"(?:^|;)ID=([^;]+)")
TRANSCRIPT_TYPES = {"mRNA", "transcript"}


def parse_tree(newick):
    """(taxa, cherry pair, internal branch support) for a four taxon tree."""
    taxa = frozenset(TAXON.findall(newick))
    m = CHERRY.search(newick)
    if not m:
        return taxa, None, None
    pair = frozenset(m.group(1, 2))
    support = float(m.group(3)) if m.group(3) else None
    return taxa, (pair if len(pair) == 2 else None), support


def internal_length(newick):
    """Length of the internal branch of a four taxon tree, or None. IQ-TREE
    writes unrooted trees as (A,(B,C)support:length,D), so the length follows
    the cherry."""
    m = CHERRY.search(newick)
    if not m:
        return None
    ml = INTERNAL_LENGTH.match(newick, m.end())
    return float(ml.group(1)) if ml else None


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


def tree_loci(gene_trees):
    """Locus names in the order of all_gene_trees.nwk: the sorted
    gene_trees/*.treefile names, or the shipped gene_tree_ids.txt listing."""
    phylo_dir = os.path.dirname(gene_trees)
    files = sorted(glob.glob(os.path.join(phylo_dir, "gene_trees", "*.treefile")))
    names = [os.path.basename(p) for p in files]
    listing = os.path.join(phylo_dir, "gene_tree_ids.txt")
    if not names and os.path.exists(listing):
        names = [line.strip() for line in open(listing) if line.strip()]
    return [n.split(".", 1)[0] for n in names]


def read_gene_trees_by_locus(path, loci):
    """[(locus, pair, support, internal length)] for trees with a cherry."""
    lines = [line.strip() for line in open(path) if line.strip()]
    if len(lines) != len(loci):
        raise SystemExit(f"{len(lines)} gene trees but {len(loci)} locus names")
    out = []
    for locus, line in zip(loci, lines):
        _, pair, sup = parse_tree(line)
        if pair is not None:
            out.append((locus, pair, sup, internal_length(line)))
    return out


def gff_transcripts(gff):
    """{transcript (periods removed, as in the protein files): attributes}."""
    out = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in TRANSCRIPT_TYPES:
                continue
            m = ATTR_ID.search(f[8])
            if m:
                out[m.group(1).replace(".", "")] = dict(
                    kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
    return out


def intact_loci(of_dir, gffs, reference):
    """{orthogroup: True if every ingroup model is intact} for single copy
    orthogroups. A model is intact when its reference source model is neither
    partial nor carries a RefSeq exception and, for a transferred model, when
    Liftoff found a valid ORF. gffs maps form to its annotation GFF."""
    tables = glob.glob(os.path.join(of_dir, "**", "Orthogroups.tsv"), recursive=True)
    tables = [t for t in tables if os.sep + "Orthogroups" + os.sep in t] or tables
    if not tables:
        raise SystemExit(f"no Orthogroups.tsv under {of_dir}")
    attrs = {form: gff_transcripts(path) for form, path in gffs.items()}
    ref = attrs[reference]

    def clean_source(t):
        a = ref.get(t)
        return a is not None and a.get("partial") != "true" and "exception" not in a

    def intact(form, t):
        if not clean_source(t):
            return False
        return form == reference or attrs[form].get(t, {}).get("valid_ORF") == "True"

    out = {}
    with open(tables[0]) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        cols = {form: head.index(form) for form in gffs}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            genes = {form: [g.strip() for g in (f[i] if i < len(f) else "").split(",")
                            if g.strip()] for form, i in cols.items()}
            if all(len(g) == 1 for g in genes.values()):
                out[f[0]] = all(intact(form, g[0]) for form, g in genes.items())
    return out


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


def site_patterns(aln_dir, taxa, per_locus=None, names=None):
    """Count parsimony informative site patterns per split over *.trim files.
    If per_locus is a list, one Counter per alignment is appended to it, and
    if names is a list, the alignment's locus name alongside."""
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
        if names is not None:
            names.append(os.path.basename(path).split(".", 1)[0])
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
    # the summed site test again, without the top loci
    rest = collections.Counter()
    for loc in per_locus:
        rest.update({s: loc[s] for s in splits})
    rest.subtract(top)
    r1, r2 = rest[major], rest[minor]
    p, lo, hi = binom(r1, r2)
    rows.append({"split": f"discordant sites outside the top {round(100 * top_share)}% of loci: "
                          "gene tree major vs minor",
                 "role": "binomial_test", "n_loci": len(per_locus) - k,
                 "n_informative_sites": r1 + r2,
                 "pct": round(100.0 * r1 / (r1 + r2), 2) if r1 + r2 else float("nan"),
                 "binomial_p": p, "ci95_low": lo, "ci95_high": hi})
    return rows


def main_splits(trees, taxa, sister):
    """(species tree split, major discordant, minor discordant), the major
    being the discordant split with more gene trees."""
    if len(taxa) != 4:
        raise SystemExit(f"expected 4 taxa, found {len(taxa)}: {sorted(taxa)}")
    conc = canonical(frozenset(sister), taxa)
    others = sorted({canonical(frozenset(p), taxa) for p in
                     [frozenset(x) for x in _pairs(taxa)]} - {conc},
                    key=lambda s: label(s, taxa))
    counts = collections.Counter(canonical(p, taxa) for p, _ in trees)
    major, minor = sorted(others, key=lambda s: -counts[s])
    return conc, major, minor


def analyse(trees, taxa, sister, cf=None, site_counts=None, per_locus=None):
    """Return (topology rows, support rows, site rows) as lists of dicts."""
    conc, major, minor = main_splits(trees, taxa, sister)
    disc = [major, minor]
    counts = collections.Counter(canonical(p, taxa) for p, _ in trees)
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


def subset_row(name, pairs, taxa, splits):
    """Gene tree counts and the binomial test for one subset of gene trees."""
    conc, major, minor = splits
    sub = collections.Counter(canonical(p, taxa) for p in pairs)
    n = sum(sub.values())
    n1, n2 = sub[major], sub[minor]
    p, lo, hi = binom(n1, n2)
    return {"subset": name, "n_gene_trees": n, "n_concordant": sub[conc],
            "pct_concordant": round(100.0 * sub[conc] / n, 2) if n else float("nan"),
            "n_major": n1, "n_minor": n2,
            "major_share_of_discordant": round(n1 / (n1 + n2), 4) if n1 + n2 else float("nan"),
            "binomial_p": p, "ci95_low": lo, "ci95_high": hi}


def robustness_rows(by_locus, taxa, splits, intact=None):
    """The gene tree test after removing unresolved gene trees (internal
    branch at the minimum length) and loci with a broken model."""
    def resolved(t):
        return t[3] is not None and t[3] > MIN_BRANCH

    rows = [subset_row("all gene trees", [t[1] for t in by_locus], taxa, splits),
            subset_row("internal branch above the minimum length",
                       [t[1] for t in by_locus if resolved(t)], taxa, splits),
            subset_row("internal branch at the minimum length",
                       [t[1] for t in by_locus if t[3] is not None and not resolved(t)],
                       taxa, splits)]
    if intact:
        good = [t for t in by_locus if intact.get(t[0]) is True]
        bad = [t for t in by_locus if intact.get(t[0]) is False]
        rows += [
            subset_row("loci with four intact models", [t[1] for t in good], taxa, splits),
            subset_row("loci with four intact models, internal branch above the minimum length",
                       [t[1] for t in good if resolved(t)], taxa, splits),
            subset_row("loci with four intact models, UFBoot >= 95",
                       [t[1] for t in good if t[2] is not None and t[2] >= 95], taxa, splits),
            subset_row("loci with a model that is not intact", [t[1] for t in bad], taxa, splits),
        ]
    return rows


def locus_quality_rows(per_locus, names, splits, intact, top_share=0.01):
    """Share of loci with a model that is not intact, among the loci holding
    the most informative sites and among the rest. The ranking is the one
    per_locus_rows uses, so the top loci are the same loci."""
    order = sorted(range(len(per_locus)),
                   key=lambda i: -sum(per_locus[i][s] for s in splits))
    k = max(1, round(top_share * len(per_locus)))
    top, rest = set(order[:k]), set(order[k:])
    rows = []
    for label_, idx in ((f"top {round(100 * top_share)}% of loci by informative sites", top),
                        (f"loci outside the top {round(100 * top_share)}%", rest)):
        known = [i for i in idx if names[i] in intact]
        broken = sum(1 for i in known if not intact[names[i]])
        rows.append({"split": label_, "role": "loci_with_a_model_not_intact",
                     "n_loci": len(known), "n_loci_not_intact": broken,
                     "pct": round(100.0 * broken / len(known), 2) if known else float("nan")})
    return rows


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


def run(gene_trees, species_tree, cf_stat, aln_dir, out_topo, out_support, out_sites,
        out_robust=None, of_dir=None, gffs=None, reference=None):
    trees, taxa = read_gene_trees(gene_trees)
    sister = sister_from_species_tree(species_tree)
    cf = read_cf_stat(cf_stat) if cf_stat and os.path.exists(cf_stat) else None
    site_counts, per_locus, names = None, [], []
    if aln_dir and os.path.isdir(aln_dir):
        site_counts, n_files, n_sites = site_patterns(aln_dir, taxa, per_locus, names)
        print(f"site patterns: {sum(site_counts.values())} informative of {n_sites} "
              f"columns in {n_files} alignments")
    topo, support, sites = analyse(trees, taxa, sister, cf, site_counts, per_locus)
    splits = main_splits(trees, taxa, sister)
    intact = intact_loci(of_dir, gffs, reference) if of_dir and gffs and reference else None
    if intact is not None:
        print(f"single copy loci with four intact models: "
              f"{sum(intact.values())} of {len(intact)}")
        if per_locus:
            sites += locus_quality_rows(per_locus, names, splits, intact)
    write_tsv(topo, out_topo)
    write_tsv(support, out_support)
    if out_sites:
        write_tsv(sites, out_sites)
    robust = []
    if out_robust:
        by_locus = read_gene_trees_by_locus(gene_trees, tree_loci(gene_trees))
        robust = robustness_rows(by_locus, taxa, splits, intact)
        write_tsv(robust, out_robust)

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
    for r in robust:
        print(f"  {r['subset']:70s} n={r['n_gene_trees']:5d} concordant {r['pct_concordant']}% "
              f"major {r['n_major']} minor {r['n_minor']} share {r['major_share_of_discordant']} "
              f"p={r['binomial_p']:.2e}")


def main():
    if "snakemake" in globals():
        sm = globals()["snakemake"]
        gffs = getattr(sm.input, "gffs", None)
        run(sm.params.gene_trees, sm.input.tree, sm.params.cf_stat,
            sm.params.trimmed_dir, sm.output.topology, sm.output.support,
            sm.output.sites, getattr(sm.output, "robustness", None),
            getattr(sm.input, "of", None) or getattr(sm.params, "orthogroups", None),
            {os.path.basename(g).split("_liftoff")[0]: g for g in gffs} if gffs else None,
            getattr(sm.params, "reference", None))
        return
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trees", required=True, help="newick gene trees, one per line")
    ap.add_argument("--species-tree", required=True, help="concatenated ML tree")
    ap.add_argument("--cf-stat", help="IQ-TREE concord.cf.stat, to label gDF1 and gDF2")
    ap.add_argument("--alignments", help="directory of trimmed *.trim alignments")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--orthogroups", help="OrthoFinder output directory")
    ap.add_argument("--gffs", nargs="*", help="<form>_liftoff.gff3 of the four ingroup forms")
    ap.add_argument("--reference", help="the form whose annotation was transferred")
    a = ap.parse_args()
    run(a.trees, a.species_tree, a.cf_stat, a.alignments,
        os.path.join(a.outdir, "quartet_topology_counts.tsv"),
        os.path.join(a.outdir, "quartet_asymmetry_by_support.tsv"),
        os.path.join(a.outdir, "quartet_site_patterns.tsv"),
        os.path.join(a.outdir, "quartet_robustness.tsv"), a.orthogroups,
        {os.path.basename(g).split("_liftoff")[0]: g for g in a.gffs} if a.gffs else None,
        a.reference)


if __name__ == "__main__":
    main()
