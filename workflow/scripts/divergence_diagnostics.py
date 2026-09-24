#!/usr/bin/env python3
"""
divergence_diagnostics.py - how much of the phylogenomic divergence is real.

The concatenated tree gives ingroup branch lengths of several percent amino
acid divergence, far more than forms that interbreed should show. Branch
lengths from a concatenated partition model are pulled up by a minority of
loci with very long branches, which is what broken gene models (a frameshifted
or truncated transfer, or mismatched isoforms) produce. This script separates
the typical locus from the outliers.

Outputs:
  branch_lengths   per taxon: branch in the concatenated tree, and across the
                   per locus gene trees the median, mean, 90th percentile,
                   share above 0.1 and share at IQ-TREE's minimum (no
                   substitution inferred on the branch); the internal branch
                   likewise
  pairwise         per pair of taxa: protein identity over the trimmed SCO
                   alignments (columns where both have a residue), median,
                   mean, 5th percentile and share of loci below 95%
  loci             SCO FASTAs, trimmed alignments kept, gene trees built, and
                   how many alignments without a tree are invariant or have
                   fewer than four distinct sequences (IQ-TREE needs four)

Snakemake provides:
    input.tree, input.concord
    params.gene_trees, params.trimmed_dir, params.sco_dir
    output.branch_lengths, output.pairwise, output.loci
"""

import csv
import glob
import itertools
import os
import re
import statistics
import sys

MIN_BRANCH = 1.1e-6          # IQ-TREE floors branch lengths at 1e-6
LONG_BRANCH = 0.1
TERMINAL = re.compile(r"([A-Za-z][A-Za-z0-9_.]*):([0-9.eE+-]+)")
INTERNAL = re.compile(r"\)[0-9./]*:([0-9.eE+-]+)")
GAP = set("-X?*.")


def quantile(values, q):
    if not values:
        return float("nan")
    v = sorted(values)
    k = (len(v) - 1) * q
    lo, hi = int(k), min(int(k) + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo]) * (k - lo)


def tree_branches(newick):
    """({taxon: terminal length}, [internal lengths]) for one newick string."""
    term = {t: float(x) for t, x in TERMINAL.findall(newick)}
    internal = [float(x) for x in INTERNAL.findall(newick)]
    return term, internal


def branch_table(concat_newick, gene_tree_lines):
    concat_term, concat_int = tree_branches(concat_newick)
    per_taxon = {t: [] for t in concat_term}
    internal = []
    n = 0
    for line in gene_tree_lines:
        line = line.strip()
        if not line:
            continue
        n += 1
        term, ints = tree_branches(line)
        for t, x in term.items():
            per_taxon.setdefault(t, []).append(x)
        internal.extend(ints)
    rows = []
    for name, vals, concat in ([(t, per_taxon[t], concat_term.get(t)) for t in sorted(per_taxon)]
                               + [("internal", internal, concat_int[0] if concat_int else None)]):
        rows.append({
            "branch": name,
            "concat_tree": round(concat, 5) if concat is not None else "",
            "n_gene_trees": len(vals),
            "gene_tree_median": round(statistics.median(vals), 5) if vals else "",
            "gene_tree_mean": round(statistics.mean(vals), 5) if vals else "",
            "gene_tree_p90": round(quantile(vals, 0.9), 5) if vals else "",
            "share_above_0.1": round(sum(x > LONG_BRANCH for x in vals) / len(vals), 4) if vals else "",
            "share_minimum": round(sum(x <= MIN_BRANCH for x in vals) / len(vals), 4) if vals else "",
        })
    return rows, n


def read_alignment(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def pair_identity(a, b):
    same = both = 0
    for x, y in zip(a, b):
        if x in GAP or y in GAP:
            continue
        both += 1
        same += x == y
    return same / both if both else None


def pairwise_table(aln_paths):
    per_pair = {}
    invariant, few_distinct = set(), set()
    for path in aln_paths:
        seqs = read_alignment(path)
        if len(seqs) < 2:
            continue
        locus = os.path.basename(path).rsplit(".", 1)[0]
        n_distinct = len(set(seqs.values()))
        if n_distinct == 1:
            invariant.add(locus)
        if n_distinct < 4:
            few_distinct.add(locus)
        for a, b in itertools.combinations(sorted(seqs), 2):
            ident = pair_identity(seqs[a], seqs[b])
            if ident is not None:
                per_pair.setdefault((a, b), []).append(ident)
    rows = []
    for (a, b), vals in sorted(per_pair.items()):
        rows.append({
            "taxon_a": a, "taxon_b": b, "n_loci": len(vals),
            "median_identity": round(statistics.median(vals), 5),
            "mean_identity": round(statistics.mean(vals), 5),
            "p05_identity": round(quantile(vals, 0.05), 5),
            "share_below_0.95": round(sum(v < 0.95 for v in vals) / len(vals), 4),
        })
    return rows, invariant, few_distinct


def write(rows, path):
    if not rows:
        open(path, "w").close()
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def run(concat_tree, gene_trees, trimmed_dir, sco_dir, out_branch, out_pair, out_loci):
    gene_lines = open(gene_trees).read().splitlines() if os.path.exists(gene_trees) else []
    branch_rows, n_trees = branch_table(open(concat_tree).read().strip(), gene_lines)
    write(branch_rows, out_branch)

    trim = sorted(glob.glob(os.path.join(trimmed_dir, "*.trim")))
    pair_rows, invariant, few_distinct = pairwise_table(trim)
    write(pair_rows, out_pair)

    phylo_dir = os.path.dirname(gene_trees)
    tree_ids = {os.path.basename(p).rsplit(".", 1)[0]
                for p in glob.glob(os.path.join(phylo_dir, "gene_trees", "*.treefile"))}
    listing = os.path.join(phylo_dir, "gene_tree_ids.txt")
    if not tree_ids and os.path.exists(listing):
        # patch runs ship a listing instead of the 70,000 IQ-TREE files
        tree_ids = {line.strip().rsplit(".", 1)[0] for line in open(listing) if line.strip()}
    trimmed_ids = {os.path.basename(p).rsplit(".", 1)[0] for p in trim}
    no_tree = trimmed_ids - tree_ids if tree_ids else set()
    n_sco = len(glob.glob(os.path.join(sco_dir, "*.fa"))) if sco_dir else 0
    sco_listing = os.path.join(phylo_dir, "sco_fasta_ids.txt")
    if not n_sco and os.path.exists(sco_listing):
        n_sco = sum(1 for line in open(sco_listing) if line.strip().endswith(".fa"))
    loci = [
        {"measure": "sco_fastas", "value": n_sco},
        {"measure": "trimmed_alignments_kept", "value": len(trim)},
        {"measure": "gene_trees_in_all_gene_trees", "value": n_trees},
        {"measure": "trimmed_without_gene_tree", "value": len(no_tree) if tree_ids else ""},
        {"measure": "trimmed_without_gene_tree_invariant",
         "value": len(no_tree & invariant) if tree_ids else ""},
        {"measure": "trimmed_without_gene_tree_fewer_than_4_distinct",
         "value": len(no_tree & few_distinct) if tree_ids else ""},
        {"measure": "trimmed_invariant_total", "value": len(invariant)},
        {"measure": "trimmed_fewer_than_4_distinct_total", "value": len(few_distinct)},
    ]
    write(loci, out_loci)

    for r in branch_rows:
        print(f"  {r['branch']:22s} concat {r['concat_tree']!s:>9}  median {r['gene_tree_median']!s:>9}"
              f"  mean {r['gene_tree_mean']!s:>9}  >0.1 {r['share_above_0.1']}")
    for r in pair_rows:
        print(f"  {r['taxon_a']} vs {r['taxon_b']}: median identity {r['median_identity']}"
              f"  share < 95% {r['share_below_0.95']}")
    for r in loci:
        print(f"  {r['measure']}: {r['value']}")


def main():
    if "snakemake" in globals():
        sm = globals()["snakemake"]
        run(sm.input.tree, sm.params.gene_trees, sm.params.trimmed_dir,
            sm.params.sco_dir, sm.output.branch_lengths, sm.output.pairwise,
            sm.output.loci)
        return
    if len(sys.argv) != 5:
        sys.exit("usage: divergence_diagnostics.py CONCAT_TREE GENE_TREES TRIMMED_DIR OUTDIR")
    concat, genes, trimmed, outdir = sys.argv[1:]
    run(concat, genes, trimmed, None,
        os.path.join(outdir, "branch_length_summary.tsv"),
        os.path.join(outdir, "sco_pairwise_identity.tsv"),
        os.path.join(outdir, "locus_accounting.tsv"))


if __name__ == "__main__":
    main()
