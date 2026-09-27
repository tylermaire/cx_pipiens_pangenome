#!/usr/bin/env python3
"""
cafe_exact_pvalues.py - exact CAFE5 family P values on a four taxon tree.

CAFE5 gives each gene family the probability, under the fitted birth and
death model, of a family whose likelihood is at most the observed family's,
given the root size; the family's P value is the largest of these over root
sizes 1 to rint(1.25 x its largest count). CAFE5 estimates each probability
from 1,000 families it simulates per root size (compute_pvalues(..., 1000) in
execute.cpp) with an unseeded random number generator. The estimates
therefore change from run to run, by about 0.007 near P = 0.05, and every
family with the same counts gets the same estimate, so a count pattern shared
by hundreds of families moves all of them across the threshold together.
In the V5 run CAFE5 estimated 0.042 for the 572 families with one copy in one
form of each pair and none in the other two; the exact value is 0.0555, and
in V4 CAFE5 had estimated 0.054 for the same patterns.

With four ingroup forms on a tree ((A,B),(C,D)), the two pairs below the root
evolve independently once the root size is fixed, so the probability of a
configuration of counts is the product of the two pair probabilities and the
distribution can be summed exactly instead of sampled. This script does that
with CAFE5's definition of the P value and CAFE5's model:

  - transition probabilities of the equal rate birth and death model
    (the_probability_of_going_from_parent_fam_size_to_c in probability.cpp),
    with the base lambda, which CAFE5 also uses for P values under the gamma
    model (execute.cpp);
  - family sizes truncated as CAFE5 truncates them (user_data.cpp): the
    largest count in the table plus max(50, that count / 5);
  - configurations with the same likelihood as the observed family counted as
    at least as extreme (CAFE5 counts simulated likelihoods <= the observed);
  - root sizes 1 to rint(1.25 x the family's largest count), largest P kept.

Only families CAFE5 tested (those in Gamma_family_results.txt) are scored.

Snakemake provides:
    input.counts    results/cafe/gene_counts_filtered.tsv
    input.tree      results/cafe/ultrametric_tree.nwk
    input.cafe_dir  results/cafe/output
    output.table    results/cafe/family_pvalues.tsv
"""

import csv
import os
import re
import sys
from math import lgamma

import numpy as np

TIE = 1e-9   # relative tolerance for likelihoods equal up to rounding


# ---------------------------------------------------------------- model
def transition(lam, t, n):
    """P[s, c], s and c in 0..n: CAFE5's birth and death probability of going
    from s to c copies along a branch of length t with rate lam,
        sum over j of C(s, j) C(s + c - j - 1, s - 1) a^(s + c - 2j) (1 - 2a)^j
    with a = lam t / (1 + lam t)."""
    a = lam * t / (1 + lam * t)
    coeff = 1 - 2 * a
    if not (0 < a and coeff > 0):
        raise ValueError(f"branch saturated: lambda {lam}, length {t}")
    la, lc = np.log(a), np.log(coeff)
    lf = np.array([lgamma(k + 1) for k in range(2 * n + 2)])   # log k!
    P = np.zeros((n + 1, n + 1))
    P[0, 0] = 1.0
    c = np.arange(n + 1)[:, None]
    j = np.arange(n + 1)[None, :]
    for s in range(1, n + 1):
        ok = j <= np.minimum(s, c)
        jj = np.where(ok, j, 0)
        log_terms = (lf[s] - lf[jj] - lf[s - jj]                      # C(s, j)
                     + lf[s + c - jj - 1] - lf[s - 1]                # C(s+c-j-1, s-1)
                     - lf[np.where(ok, c - jj, 0)]
                     + (s + c - 2 * jj) * la + jj * lc)
        P[s] = np.where(ok, np.exp(np.where(ok, log_terms, 0.0)), 0.0).sum(axis=1)
    return P


# ---------------------------------------------------------------- tree
def parse_newick(text):
    """Nested tuples: a leaf is (name, length); an inner node is
    (children, length)."""
    s = text.strip().rstrip(";")
    pos = 0

    def length():
        nonlocal pos
        if pos < len(s) and s[pos] == ":":
            m = re.compile(r":([0-9.eE+-]+)").match(s, pos)
            pos = m.end()
            return float(m.group(1))
        return 0.0

    def node():
        nonlocal pos
        if s[pos] == "(":
            pos += 1
            kids = [node()]
            while s[pos] == ",":
                pos += 1
                kids.append(node())
            if s[pos] != ")":
                raise ValueError(f"bad Newick at {pos}: {text}")
            pos += 1
            m = re.compile(r"[^:,();]*").match(s, pos)   # optional label
            pos = m.end()
            return (kids, length())
        m = re.compile(r"[^:,();]+").match(s, pos)
        pos = m.end()
        return (m.group(0).strip(), length())

    return node()


def two_pairs(newick):
    """[((A, tA), (B, tB), t_AB), ((C, tC), (D, tD), t_CD)] for
    ((A,B),(C,D)); anything else raises ValueError."""
    root, _ = parse_newick(newick)
    if not isinstance(root, list) or len(root) != 2:
        raise ValueError("the root must have two children")
    pairs = []
    for kids, t in root:
        if not isinstance(kids, list) or len(kids) != 2 or \
                any(isinstance(k[0], list) for k in kids):
            raise ValueError("each child of the root must be a pair of leaves")
        pairs.append((kids[0], kids[1], t))
    return pairs


# ---------------------------------------------------------------- P values
class PairDistribution:
    """Probabilities of the counts of one pair given the root size."""

    def __init__(self, lam, t_stem, t_a, t_b, n):
        self.T0 = transition(lam, t_stem, n)
        self.Ta = transition(lam, t_a, n)
        self.Tb = self.Ta if t_b == t_a else transition(lam, t_b, n)
        self.cache = {}

    def matrix(self, s):
        """F[x, y] = P(first leaf x, second leaf y | root size s)."""
        if s not in self.cache:
            F = self.Ta.T @ (self.T0[s][:, None] * self.Tb)
            order = np.sort(F.ravel())
            self.cache[s] = (F, order, np.cumsum(order))
        return self.cache[s]


def p_value_at_root(d1, d2, counts, s):
    """P(a configuration at most as likely as `counts` | root size s)."""
    F1, v1, _ = d1.matrix(s)
    F2, v2, c2 = d2.matrix(s)
    obs = F1[counts[0], counts[1]] * F2[counts[2], counts[3]]
    with np.errstate(divide="ignore"):
        limit = np.where(v1 > 0, obs * (1 + TIE) / v1, np.inf)
    idx = np.searchsorted(v2, limit, side="right")
    below = np.where(idx > 0, c2[np.maximum(idx - 1, 0)], 0.0)
    return float(np.sum(v1 * below))


def family_p_value(d1, d2, counts, largest, max_root, cache):
    """CAFE5's P: the largest P over root sizes 1..rint(1.25 x largest)."""
    top = min(int(np.rint(largest * 1.25)), max_root)
    best = 0.0
    for s in range(1, top + 1):
        key = (counts, s)
        if key not in cache:
            cache[key] = p_value_at_root(d1, d2, counts, s)
        best = max(best, cache[key])
    return min(best, 1.0)


def read_lambda(cafe_dir):
    for name in ("Gamma_results.txt", "Base_results.txt"):
        path = os.path.join(cafe_dir, name)
        if os.path.exists(path):
            for line in open(path):
                if line.startswith("Lambda:"):
                    values = line.split(":", 1)[1].split()
                    if len(values) != 1:
                        sys.exit(f"ERROR: {path} has several lambdas; exact P values "
                                 "are implemented for one lambda only")
                    return float(values[0]), name.split("_")[0]
    sys.exit(f"ERROR: no Lambda line in {cafe_dir}/Gamma_results.txt or Base_results.txt")


def read_cafe_pvalues(cafe_dir, model):
    path = os.path.join(cafe_dir, f"{model}_family_results.txt")
    out = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        out[f[0].strip()] = float(f[1])
    return out


def run(counts_path, tree_path, cafe_dir, out_path):
    try:
        pairs = two_pairs(open(tree_path).read())
    except ValueError as err:
        sys.exit(f"ERROR: exact P values need a four taxon tree with two pairs below "
                 f"the root ({err}); set cafe: pvalues: cafe in config/config.yaml "
                 f"to use CAFE5's simulated estimates")
    lam, model = read_lambda(cafe_dir)
    cafe_p = read_cafe_pvalues(cafe_dir, model)

    with open(counts_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        species_cols = [c for c in reader.fieldnames if c not in ("Desc", "Family ID")]
        rows = {r["Family ID"].strip(): r for r in reader}
    # CAFE5 takes the largest count over every column of the table, tree or not
    largest_all = max(int(r[c]) for r in rows.values() for c in species_cols)
    n = largest_all + max(50, largest_all // 5)
    max_root = max(30, int(np.rint(largest_all * 1.25)))

    (a, ta), (b, tb), t1 = pairs[0]
    (c, tc), (d, td), t2 = pairs[1]
    leaves = [a, b, c, d]
    missing = [x for x in leaves if x not in species_cols]
    if missing:
        sys.exit(f"ERROR: tree leaves {missing} are not columns of {counts_path}")
    d1 = PairDistribution(lam, t1, ta, tb, n)
    d2 = d1 if (t1, ta, tb) == (t2, tc, td) else PairDistribution(lam, t2, tc, td, n)

    absent = [f for f in cafe_p if f not in rows]
    if absent:
        sys.exit(f"ERROR: {len(absent)} families in the CAFE5 results are not in "
                 f"{counts_path}, e.g. {absent[:3]}")

    cache, out = {}, []
    for fam, pc in cafe_p.items():
        r = rows[fam]
        counts = tuple(int(r[x]) for x in leaves)
        largest = max(int(r[x]) for x in species_cols)
        out.append((fam, counts, pc, family_p_value(d1, d2, counts, largest, max_root, cache)))

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write("orthogroup\t" + "\t".join(leaves) + "\tcafe_pvalue\texact_pvalue\n")
        for fam, counts, pc, pe in out:
            fh.write(f"{fam}\t" + "\t".join(map(str, counts)) + f"\t{pc}\t{pe:.6g}\n")

    diff = np.array([pc - pe for _, _, pc, pe in out])
    print(f"lambda {lam} ({model} model); {len(out):,} tested families; "
          f"sizes truncated at {n}, root sizes up to {max_root}")
    print(f"CAFE5 estimate minus exact P: mean {diff.mean():+.4f}, "
          f"sd {diff.std():.4f}, largest |difference| {np.abs(diff).max():.4f}")
    return leaves, out


if __name__ == "__main__":
    if "snakemake" in globals():
        run(snakemake.input.counts, snakemake.input.tree,
            snakemake.input.cafe_dir, snakemake.output.table)
    elif len(sys.argv) == 5:
        run(*sys.argv[1:5])
    else:
        sys.exit("usage: cafe_exact_pvalues.py COUNTS TREE CAFE_DIR OUT")
