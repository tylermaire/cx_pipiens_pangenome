#!/usr/bin/env python3
"""d_statistics.py - ABBA BABA tests with the outgroup, from codon alignments.

For a test D(P1, P2; P3, O), with P1 and P2 sisters relative to P3 and O the
outgroup, a site where O carries allele A and exactly one other allele B is
present among P1, P2 and P3 is

    ABBA  P1 = A, P2 = B, P3 = B
    BABA  P1 = B, P2 = A, P3 = B
    BBAA  P1 = B, P2 = B, P3 = A   (agrees with the species tree)

D = (ABBA - BABA) / (ABBA + BABA). Incomplete lineage sorting alone makes
ABBA and BABA equally frequent (D = 0). D > 0 means P2 and P3 share more
derived alleles than P1 and P3 do; D < 0 the reverse.

Sites are codon alignment columns (codon_align.py) where all five taxa carry
A, C, G or T. Each test is run on all such sites and on third codon positions
only, and on all loci and on loci whose four ingroup gene models are intact
(quartet_asymmetry.intact_loci: the reference model neither partial nor
carrying a RefSeq exception, and each transferred model with a valid ORF).

Standard errors come from a weighted block jackknife (Busing et al. 1999)
over blocks of the genome: loci are placed by their reference gene in
windows of block_size bp of the reference assembly, each window a block,
weighted by its number of sites. Z = D / SE; p is two sided from the normal
distribution, and the Bonferroni value multiplies it by the number of tests
in the same locus set and site class.

The tests: params.tests lists planned tests as [P1, P2, P3]. Each is checked
against the rooted species tree (input.tree, rooted with the outgroup); a
planned test whose P1 and P2 are not sisters relative to P3 is still run but
flagged, and the test the tree implies for those three taxa is added. With
no planned tests, one test per trio of ingroup taxa is taken from the tree,
P1 and P2 in alphabetical order.

Outputs: output.table (one row per test, locus set and site class) and
output.per_locus (site pattern counts per locus and test, all sites).
"""
import collections
import csv
import glob
import math
import os
import sys

sys.path.insert(0, os.path.join(os.getcwd(), "workflow", "scripts"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rooting import parse, ingroup_clades, sister_pair  # noqa: E402

NUCLEOTIDES = frozenset("ACGT")
SITE_CLASSES = ("all_sites", "third_positions")


# ---------------------------------------------------------------- tests
def plan_tests(tree_path, outgroup, ingroup, planned):
    """[(test id, P1, P2, P3, consistent with the rooted tree, origin)]."""
    clades = [c for c, _, _ in ingroup_clades(parse(open(tree_path).read()), outgroup)]
    tests, seen = [], set()

    def implied(trio):
        pair = sister_pair(clades, trio)
        if pair is None:
            return None
        p1, p2 = sorted(pair)
        (p3,) = set(trio) - pair
        return p1, p2, p3

    for p1, p2, p3 in planned or []:
        for name in (p1, p2, p3):
            if name not in ingroup:
                raise SystemExit(f"planned test taxon {name} is not an ingroup sample")
        pair = sister_pair(clades, (p1, p2, p3))
        ok = pair == frozenset((p1, p2))
        tests.append([p1, p2, p3, ok, "planned"])
        seen.add(frozenset((p1, p2, p3)))
        if not ok:
            imp = implied((p1, p2, p3))
            if imp and list(imp) != [p1, p2, p3]:
                tests.append([*imp, True, "implied by the rooted tree"])
    if not planned:
        import itertools
        for trio in itertools.combinations(sorted(ingroup), 3):
            if frozenset(trio) in seen:
                continue
            imp = implied(trio)
            if imp:
                tests.append([*imp, True, "implied by the rooted tree"])
    return [(f"T{i + 1}", *t) for i, t in enumerate(tests)]


# ---------------------------------------------------------------- loci
def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}


def reference_positions(gff):
    """{transcript id with periods removed: (seqid, start)} for mRNAs."""
    out = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("mRNA", "transcript"):
                continue
            for kv in f[8].split(";"):
                if kv.startswith("ID="):
                    out[kv[3:].replace(".", "")] = (f[0], int(f[3]))
                    break
    return out


def locus_blocks(loci_table, reference, ref_gff, block_size):
    """{orthogroup: (block id, seqid, start)}; a locus whose reference gene
    has no position is its own block."""
    pos = reference_positions(ref_gff) if ref_gff and os.path.exists(ref_gff) else {}
    out = {}
    with open(loci_table) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            og = r["orthogroup"]
            p = pos.get(r.get(reference, "").replace(".", ""))
            if p is None:
                out[og] = (f"unplaced:{og}", "", "")
            else:
                out[og] = (f"{p[0]}:{p[1] // block_size}", p[0], p[1])
    return out


def count_locus(seqs, taxa_order, tests, outgroup):
    """{(test id, site class): Counter(ABBA, BABA, BBAA)} and site counts."""
    counts = {(t[0], c): collections.Counter() for t in tests for c in SITE_CLASSES}
    n_sites = collections.Counter()
    cols = zip(*(seqs[t] for t in taxa_order))
    idx = {t: i for i, t in enumerate(taxa_order)}
    o = idx[outgroup]
    tix = [(t[0], idx[t[1]], idx[t[2]], idx[t[3]]) for t in tests]
    for i, col in enumerate(cols):
        if not NUCLEOTIDES.issuperset(col):
            continue
        classes = ("all_sites", "third_positions") if i % 3 == 2 else ("all_sites",)
        for c in classes:
            n_sites[c] += 1
        a = col[o]
        for tid, i1, i2, i3 in tix:
            x1, x2, x3 = col[i1], col[i2], col[i3]
            derived = {x for x in (x1, x2, x3) if x != a}
            if len(derived) != 1:
                continue
            if x1 == a and x2 == x3:
                pattern = "ABBA"
            elif x2 == a and x1 == x3:
                pattern = "BABA"
            elif x3 == a and x1 == x2:
                pattern = "BBAA"
            else:
                continue
            for c in classes:
                counts[(tid, c)][pattern] += 1
    return counts, n_sites


# ---------------------------------------------------------------- statistics
def d_value(abba, baba):
    return (abba - baba) / (abba + baba) if abba + baba else float("nan")


def block_jackknife(blocks):
    """Weighted delete one block jackknife (Busing et al. 1999) of D.
    blocks: [(abba, baba, n_sites)]. Returns (D, SE, number of blocks)."""
    blocks = [b for b in blocks if b[2] > 0]
    A = sum(b[0] for b in blocks)
    B = sum(b[1] for b in blocks)
    n = sum(b[2] for b in blocks)
    d = d_value(A, B)
    g = len(blocks)
    if g < 2 or d != d:
        return d, float("nan"), g
    loo, weights = [], []
    for a_j, b_j, m_j in blocks:
        rest = (A - a_j) + (B - b_j)
        if rest == 0:
            return d, float("nan"), g
        loo.append(((A - a_j) - (B - b_j)) / rest)
        weights.append(m_j)
    d_jack = g * d - sum((1 - m / n) * dj for m, dj in zip(weights, loo))
    var = 0.0
    for m, dj in zip(weights, loo):
        h = n / m
        if h <= 1:
            return d, float("nan"), g
        tau = h * d - (h - 1) * dj
        var += (tau - d_jack) ** 2 / (h - 1)
    return d, math.sqrt(var / g), g


def two_sided_p(z):
    return math.erfc(abs(z) / math.sqrt(2)) if z == z else float("nan")


# ---------------------------------------------------------------- run
def run(codon_dir, tree, loci_table, ref_gff, outgroup, ingroup, reference, planned,
        block_size, intact, out_table, out_per_locus):
    tests = plan_tests(tree, outgroup, ingroup, planned)
    for t in tests:
        print(f"{t[0]}: D({t[1]}, {t[2]}; {t[3]}, {outgroup})  "
              f"{'consistent with' if t[4] else 'NOT consistent with'} the rooted tree ({t[5]})")
    blocks = locus_blocks(loci_table, reference, ref_gff, block_size)
    taxa_order = list(ingroup) + [outgroup]
    per_locus, skipped = [], 0
    for path in sorted(glob.glob(os.path.join(codon_dir, "*.fna"))):
        og = os.path.basename(path)[:-4]
        seqs = read_fasta(path)
        if set(seqs) != set(taxa_order) or len({len(s) for s in seqs.values()}) != 1:
            skipped += 1
            continue
        counts, n_sites = count_locus(seqs, taxa_order, tests, outgroup)
        block, seqid, start = blocks.get(og, (f"unplaced:{og}", "", ""))
        per_locus.append({"orthogroup": og, "block": block, "seqid": seqid, "start": start,
                          "intact": intact.get(og, "NA") if intact is not None else "NA",
                          "counts": counts, "n_sites": n_sites})
    print(f"{len(per_locus)} codon alignments counted, {skipped} skipped")

    rows = []
    locus_sets = [("all_loci", lambda r: True)]
    if intact is not None:
        locus_sets.append(("intact_loci", lambda r: r["intact"] is True))
    for set_name, keep in locus_sets:
        loci = [r for r in per_locus if keep(r)]
        for c in SITE_CLASSES:
            n_tests = len(tests)
            for tid, p1, p2, p3, ok, origin in tests:
                by_block = collections.defaultdict(lambda: [0, 0, 0])
                tot = collections.Counter()
                for r in loci:
                    k = r["counts"][(tid, c)]
                    b = by_block[r["block"]]
                    b[0] += k["ABBA"]
                    b[1] += k["BABA"]
                    b[2] += r["n_sites"][c]
                    tot.update(k)
                    tot["sites"] += r["n_sites"][c]
                d, se, g = block_jackknife([tuple(v) for v in by_block.values()])
                z = d / se if se == se and se > 0 else float("nan")
                p = two_sided_p(z)
                if d != d:
                    excess = "NA"
                elif d > 0:
                    excess = f"{p2}+{p3}"
                elif d < 0:
                    excess = f"{p1}+{p3}"
                else:
                    excess = "none"
                rows.append({
                    "analysis": f"{tid} {set_name} {c}", "test": tid,
                    "statistic": f"D({p1},{p2};{p3},{outgroup})",
                    "locus_set": set_name, "site_class": c,
                    "P1": p1, "P2": p2, "P3": p3, "O": outgroup,
                    "consistent_with_rooted_tree": ok, "origin": origin,
                    "n_loci": len(loci), "n_blocks": g, "n_sites": tot["sites"],
                    "BBAA": tot["BBAA"], "ABBA": tot["ABBA"], "BABA": tot["BABA"],
                    "D": round(d, 5) if d == d else "NA",
                    "se": round(se, 5) if se == se else "NA",
                    "Z": round(z, 3) if z == z else "NA",
                    "p": f"{p:.3e}" if p == p else "NA",
                    "p_bonferroni": f"{min(1.0, p * n_tests):.3e}" if p == p else "NA",
                    "abs_Z_at_least_3": (abs(z) >= 3) if z == z else "NA",
                    "excess_derived_sharing": excess})
    with open(out_table, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with open(out_per_locus, "w", newline="") as fh:
        head = ["orthogroup", "block", "seqid", "start", "intact", "n_sites", "n_third"]
        for t in tests:
            head += [f"{t[0]}_ABBA", f"{t[0]}_BABA", f"{t[0]}_BBAA"]
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(head)
        for r in per_locus:
            row = [r["orthogroup"], r["block"], r["seqid"], r["start"], r["intact"],
                   r["n_sites"]["all_sites"], r["n_sites"]["third_positions"]]
            for t in tests:
                k = r["counts"][(t[0], "all_sites")]
                row += [k["ABBA"], k["BABA"], k["BBAA"]]
            w.writerow(row)
    for r in rows:
        print(f"  {r['analysis']:32s} {r['statistic']:62s} ABBA {r['ABBA']:6d} "
              f"BABA {r['BABA']:6d} BBAA {r['BBAA']:7d} D {r['D']} Z {r['Z']} "
              f"blocks {r['n_blocks']}")
    return rows


def main():
    if "snakemake" not in globals():
        sys.exit("Run via Snakemake (see the docstring).")
    sm = globals()["snakemake"]
    ingroup = list(sm.params.ingroup)
    reference = sm.params.reference
    gffs = {os.path.basename(g).split("_liftoff")[0]: g for g in sm.input.gffs}
    intact = None
    if getattr(sm.input, "of", None):
        from quartet_asymmetry import intact_loci
        intact = intact_loci(sm.input.of, {s: gffs[s] for s in ingroup}, reference)
        print(f"loci with four intact ingroup models: {sum(intact.values())} of {len(intact)}")
    run(sm.params.codon, sm.input.tree, sm.input.loci, gffs.get(reference),
        sm.params.outgroup, ingroup, reference, [list(t) for t in sm.params.tests],
        int(sm.params.block_size), intact, sm.output.table, sm.output.per_locus)


if __name__ == "__main__":
    main()
