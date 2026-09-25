#!/usr/bin/env python3
"""rooted_summary.py - the species tree rooted with the outgroup, its support,
and the rooted topologies of the five taxon gene trees.

Outputs
-------
  summary      item/value rows: loci and sites in the concatenated alignment,
               the rooted ingroup topology, where the root falls, and for each
               internal branch (named by its ingroup clade) the ultrafast
               bootstrap, gene and site concordance factors and length
  topologies   one row per rooted ingroup topology seen among the gene trees:
               its quartet split, whether it is the species tree, and counts
               among all gene trees, gene trees with both internal branches
               above IQ-TREE's minimum length (resolved), and gene trees with
               both internal branches at UFBoot 95 or more

Under incomplete lineage sorting alone the two quartet splits that disagree
with the species tree are equally frequent; the rooted topology adds where
the root falls within each split.

Snakemake provides input.tree, input.gene_trees, input.ids, input.accounting,
input.alignments; params.outgroup, params.ingroup, params.trimmed (the
trimmed protein alignments of the concatenated tree), params.cf_stat,
params.cf_branch; output.summary, output.topologies.
"""
import collections
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.join(os.getcwd(), "workflow", "scripts"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rooting import parse, ingroup_clades, rooted_topology, quartet_split  # noqa: E402

MIN_BRANCH = 1.1e-6            # IQ-TREE floors branch lengths at 1e-6


def fasta_width(path):
    """Columns of the first sequence of an alignment."""
    n, seen = 0, False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if seen:
                    break
                seen = True
            elif seen:
                n += len(line.strip())
    return n


def read_cf_stat(path):
    """{branch id: row} from IQ-TREE's .cf.stat."""
    if not path or not os.path.exists(path):
        return {}
    rows = [l.rstrip("\n").split("\t") for l in open(path)
            if l.strip() and not l.startswith("#")]
    if not rows:
        return {}
    head = rows[0]
    return {r[0]: dict(zip(head, r)) for r in rows[1:]}


def support(label):
    """UFBoot from an internal label ('100' or '100/85.2/70.1')."""
    if not label:
        return None
    try:
        return float(label.split("/")[0])
    except ValueError:
        return None


def clade_name(clade):
    return "(" + ",".join(sorted(clade)) + ")"


def root_position(clades, ingroup):
    """Describe the root: the two ingroup clades that meet at it."""
    ingroup = frozenset(ingroup)
    top = sorted((c for c in clades if not any(c < d for d in clades)), key=sorted)
    if len(top) == 2 and top[0] | top[1] == ingroup:
        return f"between {clade_name(top[0])} and {clade_name(top[1])}"
    if len(top) == 1:
        rest = sorted(ingroup - top[0])
        return f"between {clade_name(top[0])} and {', '.join(rest)}"
    return "unresolved"


def species_tree_rows(tree_path, cf_stat, cf_branch, outgroup, ingroup):
    tree = parse(open(tree_path).read())
    clades = ingroup_clades(tree, outgroup)
    topo = rooted_topology([c for c, _, _ in clades], ingroup)
    rows = [("rooted ingroup topology", topo),
            ("root position", root_position([c for c, _, _ in clades], ingroup))]
    cf = read_cf_stat(cf_stat)
    ids = {}
    if cf_branch and os.path.exists(cf_branch):
        for c, label, _ in ingroup_clades(parse(open(cf_branch).read()), outgroup):
            ids[c] = label
    for c, label, length in sorted(clades, key=lambda x: (len(x[0]), sorted(x[0]))):
        name = clade_name(c)
        rows.append((f"{name} UFBoot", support(label)))
        rows.append((f"{name} branch length", length))
        stat = cf.get(ids.get(c, ""), {})
        for k in ("gCF", "gCF_N", "gDF1", "gDF1_N", "gDF2", "gDF2_N", "gDFP", "gDFP_N",
                  "gN", "sCF", "sDF1", "sDF2", "sN"):
            rows.append((f"{name} {k}", stat.get(k, "NA") if stat else "NA"))
    return rows, topo, frozenset(c for c, _, _ in clades)


def gene_tree_rows(trees_path, ids_path, outgroup, ingroup, species_topo):
    lines = [l.strip() for l in open(trees_path) if l.strip()]
    ids = [l.strip() for l in open(ids_path) if l.strip()] if ids_path else []
    if ids and len(ids) != len(lines):
        raise SystemExit(f"{len(lines)} gene trees but {len(ids)} ids")
    counts = {k: collections.Counter() for k in ("all", "resolved", "ufboot95")}
    splits = {}
    skipped = 0
    for line in lines:
        tree = parse(line)
        if set(tree.tips()) != set(ingroup) | {outgroup}:
            skipped += 1
            continue
        clades = ingroup_clades(tree, outgroup)
        topo = rooted_topology([c for c, _, _ in clades], ingroup)
        splits[topo] = quartet_split([c for c, _, _ in clades], ingroup)
        counts["all"][topo] += 1
        if clades and all(l is not None and l > MIN_BRANCH for _, _, l in clades):
            counts["resolved"][topo] += 1
        sups = [support(lab) for _, lab, _ in clades]
        if clades and all(s is not None and s >= 95 for s in sups):
            counts["ufboot95"][topo] += 1
    species_split = splits.get(species_topo)
    rows = []
    totals = {k: sum(v.values()) for k, v in counts.items()}
    for topo, n in counts["all"].most_common():
        row = {"rooted_topology": topo, "quartet_split": splits[topo],
               "is_species_tree": topo == species_topo,
               "has_species_tree_split": splits[topo] == species_split}
        for k in ("all", "resolved", "ufboot95"):
            row[f"n_{k}"] = counts[k][topo]
            row[f"pct_{k}"] = round(100.0 * counts[k][topo] / totals[k], 2) if totals[k] else 0.0
        rows.append(row)
    return rows, totals, skipped


def write_tsv(rows, path):
    keys = list(dict.fromkeys(k for r in rows for k in r)) if rows else ["rooted_topology"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def run(tree, gene_trees, ids, trimmed, accounting, alignments, outgroup, ingroup,
        cf_stat, cf_branch, out_summary, out_topologies):
    items = []
    for path in (accounting, alignments):
        if path and os.path.exists(path):
            for r in csv.DictReader(open(path), delimiter="\t"):
                key = r["item"] + (f" {r['sample']}" if r.get("sample") else "")
                items.append((f"loci: {key}", r["value"]))
    aligns = sorted(glob.glob(os.path.join(trimmed, "*.trim"))) if trimmed else []
    items.append(("concatenated alignment loci", len(aligns)))
    items.append(("concatenated alignment columns", sum(fasta_width(p) for p in aligns)))
    tree_rows, topo, _ = species_tree_rows(tree, cf_stat, cf_branch, outgroup, ingroup)
    items.extend(tree_rows)
    topo_rows, totals, skipped = gene_tree_rows(gene_trees, ids, outgroup, ingroup, topo)
    for k, v in totals.items():
        items.append((f"gene trees {k}", v))
    items.append(("gene trees without all five taxa", skipped))
    conc = next((r for r in topo_rows if r["is_species_tree"]), None)
    items.append(("gene trees matching the rooted species tree (pct of all)",
                  conc["pct_all"] if conc else 0.0))
    items.append(("gene trees with the species tree quartet split (pct of all)",
                  round(sum(r["pct_all"] for r in topo_rows if r["has_species_tree_split"]), 2)))
    with open(out_summary, "w") as fh:
        fh.write("item\tvalue\n")
        for k, v in items:
            fh.write(f"{k}\t{'NA' if v is None else v}\n")
    write_tsv(topo_rows, out_topologies)
    print(f"rooted species tree: {topo}")
    for k, v in items:
        print(f"  {k}: {v}")
    for r in topo_rows[:8]:
        print(f"  {r['rooted_topology']:70s} {r['n_all']:6d} {r['pct_all']:6.2f}%  "
              f"resolved {r['n_resolved']}  UFBoot95 {r['n_ufboot95']}")


def main():
    if "snakemake" not in globals():
        sys.exit("Run via Snakemake (see the docstring).")
    sm = globals()["snakemake"]
    run(sm.input.tree, sm.input.gene_trees, sm.input.ids, sm.params.trimmed,
        sm.input.accounting, sm.input.alignments, sm.params.outgroup,
        list(sm.params.ingroup), sm.params.cf_stat, sm.params.cf_branch,
        sm.output.summary, sm.output.topologies)


if __name__ == "__main__":
    main()
