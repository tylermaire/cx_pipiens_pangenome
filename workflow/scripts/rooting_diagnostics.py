#!/usr/bin/env python3
"""rooting_diagnostics.py - which loci carry the signal that roots the
ingroup, and do the root and the D statistics hold on cleaner locus sets.

Sites: codon alignment columns where all five taxa carry A, C, G or T and
exactly two alleles occur. The outgroup allele is taken as ancestral and the
site is labelled by the ingroup taxa carrying the other (derived) allele;
sites with one or all ingroup taxa derived are not counted. A site derived in
two taxa supports the clade of those two, one derived in all but one taxon
supports a root on the branch to that taxon.

Where transferred models of several forms share an error relative to the
reference (a reading frame shifted where RefSeq corrected the reference
model, say), every misaligned codon after it counts as a derived allele
shared by those forms. The tables show this by gene model quality.

Outputs
-------
  patterns   locus_class, site_class, n_loci, derived_in, n_sites,
             per_100_loci: every derived allele pattern in every locus class
             (all loci; intact and not; each reason a locus is not intact;
             and loci above ingroup identity thresholds)
  gene_trees locus_set, rank, rooted_topology, is_species_tree, n_resolved,
             pct_resolved, n_resolved_total: the most frequent rooted
             topologies of the resolved gene trees (both internal branches
             above IQ-TREE's minimum length) in each locus set
  dstat      locus_set, test, statistic, n_loci, n_blocks, ABBA, BABA, D, se,
             Z: each test of d_statistics on each locus set, with the same
             weighted block jackknife, from d_statistics' per locus counts

Snakemake provides input.quality, input.aln, input.gene_trees, input.ids,
input.tree, input.dstat, input.per_locus; params.ingroup, params.outgroup,
params.identity (thresholds); output.patterns, output.gene_trees,
output.dstat.
"""
import collections
import csv
import itertools
import os
import sys

sys.path.insert(0, os.path.join(os.getcwd(), "workflow", "scripts"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rooting import parse, ingroup_clades, rooted_topology  # noqa: E402
from d_statistics import block_jackknife  # noqa: E402

BASES = frozenset("ACGT")
MIN_BRANCH = 1.1e-6            # IQ-TREE floors branch lengths at 1e-6
TOP_TOPOLOGIES = 5


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def patterns(seqs, ingroup, outgroup, start=0, step=1):
    """{frozenset of derived ingroup taxa: sites}."""
    c = collections.Counter()
    out = seqs[outgroup]
    cols = [seqs[t] for t in ingroup]
    for i in range(start, len(out), step):
        o = out[i]
        col = [s[i] for s in cols]
        if o not in BASES or any(x not in BASES for x in col) or len(set(col) - {o}) != 1:
            continue
        derived = frozenset(t for t, x in zip(ingroup, col) if x != o)
        if 1 < len(derived) < len(ingroup):
            c[derived] += 1
    return c


def fmt(taxa):
    return "+".join(sorted(taxa))


def locus_sets(quality, thresholds):
    """[(name, set of loci)]: the classes of the patterns table."""
    have = [r for r in quality if r["codon"]]
    classes = [("all", {r["orthogroup"] for r in have}),
               ("intact", {r["orthogroup"] for r in have if r["intact"] == "True"}),
               ("not intact", {r["orthogroup"] for r in have if r["intact"] == "False"})]
    by_reason = collections.defaultdict(set)
    for r in have:
        if r["intact"] == "False":
            by_reason[r["reason"]].add(r["orthogroup"])
    classes += [(f"not intact: {k}", v)
                for k, v in sorted(by_reason.items(), key=lambda kv: (-len(kv[1]), kv[0]))]
    return classes, identity_sets(have, thresholds)


def identity_sets(rows, thresholds):
    out = []
    ident = {r["orthogroup"]: float(r["min_ingroup_identity"]) for r in rows
             if r["min_ingroup_identity"] not in ("", "NA")}
    for t in thresholds:
        pct = f"{100 * t:g}%"
        out.append((f"identity >= {pct}", {o for o, v in ident.items() if v >= t}))
    for t in thresholds:
        pct = f"{100 * t:g}%"
        out.append((f"intact, identity >= {pct}",
                    {r["orthogroup"] for r in rows if r["intact"] == "True"
                     and r["orthogroup"] in ident and ident[r["orthogroup"]] >= t}))
    return out


def run(quality_path, aln_dir, gene_trees, ids_path, tree_path, dstat_path, per_locus_path,
        ingroup, outgroup, thresholds, out_patterns, out_gene_trees, out_dstat):
    quality = list(csv.DictReader(open(quality_path), delimiter="\t"))
    counts = {}
    for r in quality:
        path = os.path.join(aln_dir, "codon", f"{r['orthogroup']}.fna")
        r["codon"] = False
        if os.path.exists(path):
            seqs = read_fasta(path)
            if set(ingroup) | {outgroup} <= set(seqs):
                counts[r["orthogroup"]] = (patterns(seqs, ingroup, outgroup),
                                           patterns(seqs, ingroup, outgroup, 2, 3))
                r["codon"] = True
    classes, id_sets = locus_sets(quality, thresholds)
    all_patterns = [frozenset(c) for k in range(2, len(ingroup))
                    for c in itertools.combinations(sorted(ingroup), k)]
    rows = []
    for name, loci in classes + id_sets:
        for site_class, idx in (("all_sites", 0), ("third_positions", 1)):
            tot = collections.Counter()
            for og in loci:
                tot.update(counts[og][idx])
            for p in all_patterns:
                rows.append({"locus_class": name, "site_class": site_class, "n_loci": len(loci),
                             "derived_in": fmt(p), "n_sites": tot[p],
                             "per_100_loci": round(100.0 * tot[p] / len(loci), 1) if loci else 0.0})
    write(rows, out_patterns)

    sets = [classes[0], classes[1], classes[2]] + id_sets
    species = species_topology(tree_path, ingroup, outgroup)
    topo = {}
    trees = [l.strip() for l in open(gene_trees) if l.strip()]
    ids = [l.strip() for l in open(ids_path) if l.strip()]
    for og, text in zip(ids, trees):
        clades = ingroup_clades(parse(text), outgroup)
        if clades and all(l is not None and l > MIN_BRANCH for _, _, l in clades):
            topo[og] = rooted_topology([c for c, _, _ in clades], ingroup)
    rows = []
    for name, loci in sets:
        tt = collections.Counter(topo[o] for o in loci if o in topo)
        n = sum(tt.values())
        ranked = tt.most_common(TOP_TOPOLOGIES)
        if species not in [t for t, _ in ranked] and tt[species]:
            ranked.append((species, tt[species]))
        for rank, (t, k) in enumerate(ranked, 1):
            rows.append({"locus_set": name, "rank": rank, "rooted_topology": t,
                         "is_species_tree": t == species, "n_resolved": k,
                         "pct_resolved": round(100.0 * k / n, 2) if n else 0.0,
                         "n_resolved_total": n})
    write(rows, out_gene_trees)

    tests = [(r["test"], r["statistic"]) for r in csv.DictReader(open(dstat_path), delimiter="\t")
             if r["locus_set"] == "all_loci" and r["site_class"] == "all_sites"]
    per = {r["orthogroup"]: r for r in csv.DictReader(open(per_locus_path), delimiter="\t")}
    rows = []
    for name, loci in sets:
        use = [per[o] for o in sorted(loci) if o in per]
        for tid, statistic in tests:
            blocks = collections.defaultdict(lambda: [0, 0, 0])
            for r in use:
                b = blocks[r["block"]]
                b[0] += int(r[f"{tid}_ABBA"])
                b[1] += int(r[f"{tid}_BABA"])
                b[2] += int(r["n_sites"])
            d, se, g = block_jackknife([tuple(v) for v in blocks.values()])
            z = d / se if se == se and se > 0 else float("nan")
            rows.append({"locus_set": name, "test": tid, "statistic": statistic,
                         "n_loci": len(use), "n_blocks": g,
                         "ABBA": sum(b[0] for b in blocks.values()),
                         "BABA": sum(b[1] for b in blocks.values()),
                         "D": round(d, 5) if d == d else "NA",
                         "se": round(se, 5) if se == se else "NA",
                         "Z": round(z, 3) if z == z else "NA"})
    write(rows, out_dstat)
    report(out_patterns, out_gene_trees, out_dstat)


def species_topology(tree_path, ingroup, outgroup):
    clades = ingroup_clades(parse(open(tree_path).read().strip()), outgroup)
    return rooted_topology([c for c, _, _ in clades], ingroup)


def write(rows, path):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["empty"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def report(patterns_path, trees_path, dstat_path):
    rows = [r for r in csv.DictReader(open(patterns_path), delimiter="\t")
            if r["site_class"] == "all_sites"]
    pats = list(dict.fromkeys(r["derived_in"] for r in rows))
    print("derived allele patterns per 100 loci, all sites")
    for name in dict.fromkeys(r["locus_class"] for r in rows):
        sel = {r["derived_in"]: r for r in rows if r["locus_class"] == name}
        print(f"  {name[:40]:40s} {sel[pats[0]]['n_loci']:>6} " +
              " ".join(f"{sel[p]['per_100_loci']:>7}" for p in pats))
    for r in csv.DictReader(open(trees_path), delimiter="\t"):
        if r["rank"] in ("1", "2"):
            print(f"  {r['locus_set']:28s} {r['rooted_topology']:70s} {r['pct_resolved']}%")
    for r in csv.DictReader(open(dstat_path), delimiter="\t"):
        print(f"  {r['locus_set']:28s} {r['test']:4s} {r['statistic']:62s} D {r['D']} Z {r['Z']}")


def main():
    if "snakemake" not in globals():
        sys.exit("Run via Snakemake (see the docstring).")
    sm = globals()["snakemake"]
    run(sm.input.quality, sm.input.aln, sm.input.gene_trees, sm.input.ids, sm.input.tree,
        sm.input.dstat, sm.input.per_locus, list(sm.params.ingroup), sm.params.outgroup,
        [float(x) for x in sm.params.identity], sm.output.patterns, sm.output.gene_trees,
        sm.output.dstat)


if __name__ == "__main__":
    main()
