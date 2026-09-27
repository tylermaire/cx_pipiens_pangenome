#!/usr/bin/env python3
"""
anvio_compare.py - summarize the anvi'o gene cluster pangenome and compare it
with the OrthoFinder partition of the same proteins.

Reads the gene cluster table of anvi-summarize, the gene maps written by
anvio_inputs.py (gene caller id to transcript), OrthoFinder's Orthogroups.tsv
and Orthogroups_UnassignedGenes.tsv, the partition table
(results/pangenome/partitioned_orthogroups.tsv), and the gene of each
transcript with Liftoff's ORF flag. The last come from the annotation GFFs
(results/annotation/<form>_liftoff.gff3, read as cloud_composition.py and
transfer_quality.py read them) or, where the GFFs are not at hand, from a
table of the same content (--gene-table: sample, transcript, gene, valid_ORF,
reference_partial, reference_exception).

A gene cluster is in the class "all four" when it holds genes of all four
forms, "two or three" when of two or three, and "one form" when of one. The
OrthoFinder compartments are core, shell and cloud over the four forms, and
unassigned for genes OrthoFinder left out of orthogroups.

anvi'o builds its MCL graph from the DIAMOND hits alone, so a protein with no
hit, not even to itself (mostly fragments a few residues long), is left out
of the gene clusters without a warning. Such genes are reported apart, as
"not clustered", and are not counted in any class.

Writes, in --out:
    anvio_gene_clusters.tsv.gz  one row per gene: gene cluster, its class,
                                orthogroup and compartment, valid_ORF
    anvio_combinations.tsv      gene clusters and genes per combination of forms
    anvio_vs_orthofinder.tsv    genes by anvi'o class and OrthoFinder compartment
    anvio_orf_by_class.tsv      transferred models without a valid ORF by class
    anvio_one_form.tsv          genes of one form gene clusters, per form
    anvio_summary.tsv           section, item, value (the numbers in the text)
"""
import argparse
import collections
import csv
import glob
import gzip
import os
import re
import statistics
import sys

FORMS = ["Cx_quinquefasciatus", "Cx_pallens", "Cx_molestus", "Cx_pipiens"]
REF = "Cx_quinquefasciatus"
TRANSFERRED = FORMS[1:]
ABBR = {"Cx_quinquefasciatus": "qui", "Cx_pallens": "pal", "Cx_molestus": "mol",
        "Cx_pipiens": "pip"}
CLASSES = ["all four", "two or three", "one form"]
COMPARTMENTS = ["core", "shell", "cloud", "unassigned"]

ATTR_ID = re.compile(r"(?:^|;)ID=([^;]+)")
ATTR_PARENT = re.compile(r"(?:^|;)Parent=([^;,]+)")
TRANSCRIPT_TYPES = {"mRNA", "transcript"}


def open_any(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode.replace("t", ""))


def gene_cluster_table(summary):
    hits = sorted(glob.glob(os.path.join(summary, "*_gene_clusters_summary.txt*")))
    if not hits:
        sys.exit(f"no *_gene_clusters_summary.txt(.gz) in {summary}")
    with open_any(hits[0]) as fh:
        return list(csv.DictReader(fh, delimiter="\t")), hits[0]


def gene_maps(inputs):
    """{(form, gene caller id): transcript} and {(form, transcript): protein length}."""
    out, length = {}, {}
    for form in FORMS:
        with open(os.path.join(inputs, f"{form}.gene_map.tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                out[(form, int(r["gene_callers_id"]))] = r["transcript"]
        with open(os.path.join(inputs, f"{form}.gene_calls.tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                length[(form, out[(form, int(r["gene_callers_id"]))])] = len(r["aa_sequence"])
    return out, length


def orthogroups(of_dir):
    """{(form, transcript): orthogroup, or 'unassigned'}."""
    og_file = glob.glob(os.path.join(of_dir, "**", "Orthogroups.tsv"), recursive=True)
    un_file = glob.glob(os.path.join(of_dir, "**", "Orthogroups_UnassignedGenes.tsv"),
                        recursive=True)
    if not og_file or not un_file:
        sys.exit(f"no Orthogroups.tsv and Orthogroups_UnassignedGenes.tsv under {of_dir}")
    out = {}
    for path, unassigned in ((og_file[0], False), (un_file[0], True)):
        with open(path) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                for i, cell in enumerate(f[1:], 1):
                    if i < len(head) and head[i] in FORMS:
                        for g in cell.split(","):
                            if g.strip():
                                out[(head[i], g.strip())] = "unassigned" if unassigned else f[0]
    return out


def compartments(path):
    with open(path) as fh:
        return {r["Orthogroup"]: r["compartment"] for r in csv.DictReader(fh, delimiter="\t")}


def gene_table(annotation, table):
    """{(form, transcript): (gene, valid_ORF, reference model not clean)}."""
    out = {}
    gffs = {f: os.path.join(annotation, f"{f}_liftoff.gff3") for f in FORMS}
    if all(os.path.exists(p) for p in gffs.values()):
        for form, path in gffs.items():
            with open(path) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 9 or f[2] not in TRANSCRIPT_TYPES:
                        continue
                    m = ATTR_ID.search(f[8])
                    if not m:
                        continue
                    attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                    p = ATTR_PARENT.search(f[8])
                    tid = m.group(1).replace(".", "")
                    out[(form, tid)] = (p.group(1).replace(".", "") if p else tid,
                                        attrs.get("valid_ORF", ""),
                                        attrs.get("partial") == "true" or "exception" in attrs)
        return out, "annotation GFFs"
    if not table or not os.path.exists(table):
        sys.exit("need the annotation GFFs or --gene-table")
    with open(table) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[(r["sample"], r["transcript"])] = (
                r["gene"], r["valid_ORF"],
                r["reference_partial"] == "true" or r["reference_exception"] != "")
    return out, table


def cls(n_forms):
    return {4: "all four", 1: "one form"}.get(n_forms, "two or three")


def pct(k, n):
    return round(100 * k / n, 1) if n else "NA"


def write_tsv(path, header, rows):
    with open_any(path, "wt") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", default="results/anvio/summary")
    ap.add_argument("--inputs", default="results/anvio/inputs")
    ap.add_argument("--orthofinder", default="results/orthofinder/output")
    ap.add_argument("--partition", default="results/pangenome/partitioned_orthogroups.tsv")
    ap.add_argument("--annotation", default="results/annotation")
    ap.add_argument("--gene-table", default=None)
    ap.add_argument("--out", default="results/anvio")
    args = ap.parse_args()

    rows, src = gene_cluster_table(args.summary)
    gmap, aa_len = gene_maps(args.inputs)
    og = orthogroups(args.orthofinder)
    comp = compartments(args.partition)
    ann, ann_src = gene_table(args.annotation, args.gene_table)

    genes = []                                  # (form, transcript, cluster)
    for r in rows:
        genes.append((r["genome_name"], gmap[(r["genome_name"], int(r["gene_callers_id"]))],
                      r["gene_cluster_id"]))
    if len({(f, t) for f, t, _ in genes}) != len(genes):
        sys.exit("a gene appears twice in the gene cluster table")
    clustered = {(f, t) for f, t, _ in genes}
    unclustered = sorted({(f, t) for (f, _), t in gmap.items()} - clustered)
    missing = [g for g in list(clustered) + unclustered if g not in og or g not in ann]
    if missing:
        sys.exit(f"{len(missing)} genes missing from OrthoFinder or the annotation, e.g. {missing[:3]}")

    members = collections.defaultdict(list)
    for form, tx, c in genes:
        members[c].append((form, tx))
    forms_of = {c: frozenset(f for f, _ in m) for c, m in members.items()}
    cluster_of = {(f, tx): c for f, tx, c in genes}

    def klass(g):
        return cls(len(forms_of[cluster_of[g]]))

    def compartment(g):
        o = og[g]
        return "unassigned" if o == "unassigned" else comp.get(o, "missing")

    # kept models of each reference gene, per form
    kept = collections.defaultdict(set)         # gene -> forms with a kept model
    for (form, _), tx in gmap.items():
        kept[ann[(form, tx)][0]].add(form)

    write_tsv(os.path.join(args.out, "anvio_gene_clusters.tsv.gz"),
              ["form", "transcript", "gene", "gene_cluster", "cluster_forms", "cluster_genes",
               "anvio_class", "orthogroup", "orthofinder_compartment", "valid_ORF",
               "protein_length"],
              sorted([[f, tx, ann[(f, tx)][0], c, len(forms_of[c]), len(members[c]),
                       klass((f, tx)), og[(f, tx)], compartment((f, tx)), ann[(f, tx)][1],
                       aa_len[(f, tx)]] for f, tx, c in genes]
                     + [[f, tx, ann[(f, tx)][0], "", 0, 0, "not clustered", og[(f, tx)],
                         compartment((f, tx)), ann[(f, tx)][1], aa_len[(f, tx)]]
                        for f, tx in unclustered]))

    S = []

    def add(section, item, value):
        S.append((section, item, value))

    # ---- the gene clusters
    add("clusters", "input proteins", len(gmap))
    add("clusters", "proteins in gene clusters", len(genes))
    add("clusters", "gene clusters", len(members))
    add("not clustered", "proteins without a DIAMOND hit, left out of gene clusters",
        len(unclustered))
    if unclustered:
        lens = sorted(aa_len[g] for g in unclustered)
        add("not clustered", "protein length, median (residues)", statistics.median(lens))
        add("not clustered", "protein length, at most 10 residues", sum(1 for x in lens if x <= 10))
        add("not clustered", "protein length, largest (residues)", lens[-1])
        for form in FORMS:
            add("not clustered", f"{form} | proteins", sum(1 for f, _ in unclustered if f == form))
        for c in COMPARTMENTS:
            add("not clustered", f"OrthoFinder {c}", sum(1 for g in unclustered if compartment(g) == c))
        tr = [g for g in unclustered if g[0] != REF]
        add("not clustered", "transferred models without a valid ORF",
            sum(1 for g in tr if ann[g][1] == "False"))
        add("not clustered", "transferred models", len(tr))
    by_n = collections.Counter(len(forms_of[c]) for c in members)
    for k, lab in ((4, "all four forms"), (3, "three forms"), (2, "two forms"), (1, "one form")):
        add("clusters", f"clusters in {lab}", by_n.get(k, 0))
        add("clusters", f"clusters in {lab} (%)", pct(by_n.get(k, 0), len(members)))
    for k in CLASSES:
        n = sum(1 for f, tx, _ in genes if klass((f, tx)) == k)
        add("clusters", f"genes in {k} clusters", n)
        add("clusters", f"genes in {k} clusters (%)", pct(n, len(genes)))
    add("clusters", "clusters with one gene in each form",
        sum(1 for c, m in members.items() if len(m) == 4 and len(forms_of[c]) == 4))
    sizes = [len(m) for m in members.values()]
    add("clusters", "genes per cluster, median", statistics.median(sizes))
    add("clusters", "genes per cluster, largest", max(sizes))
    add("clusters", "clusters of one gene", sum(1 for s in sizes if s == 1))
    for form in FORMS:
        add("per form", f"{form} | clusters with its genes",
            sum(1 for c in members if form in forms_of[c]))
        add("per form", f"{form} | one form clusters",
            sum(1 for c in members if forms_of[c] == frozenset([form])))
        add("per form", f"{form} | clusters missing only this form",
            sum(1 for c in members if len(forms_of[c]) == 3 and form not in forms_of[c]))

    combos = collections.defaultdict(lambda: [0, 0])
    for c, m in members.items():
        combos[forms_of[c]][0] += 1
        combos[forms_of[c]][1] += len(m)
    write_tsv(os.path.join(args.out, "anvio_combinations.tsv"),
              ["combination"] + FORMS + ["n_forms", "n_clusters", "n_genes"],
              [[" + ".join(ABBR[f] for f in FORMS if f in k)] + [int(f in k) for f in FORMS]
               + [len(k), v[0], v[1]]
               for k, v in sorted(combos.items(), key=lambda kv: (-len(kv[0]), -kv[1][0]))])

    # ---- one form clusters: the same gene in another form, and broken models
    one_rows = []
    for form in FORMS:
        own = [(form, tx) for c in members if forms_of[c] == frozenset([form])
               for _, tx in members[c]]
        row = {"form": form,
               "clusters": sum(1 for c in members if forms_of[c] == frozenset([form])),
               "genes": len(own),
               "same_gene_in_another_form": sum(1 for g in own if kept[ann[g][0]] - {form}),
               "refseq_model_not_clean": sum(1 for g in own if ann[g][2])}
        for c in COMPARTMENTS:
            row[f"orthofinder_{c}"] = sum(1 for g in own if compartment(g) == c)
        row["without_valid_orf"] = ("NA" if form == REF
                                    else sum(1 for g in own if ann[g][1] == "False"))
        one_rows.append(row)
        add("one form clusters", f"{form} | genes", row["genes"])
        add("one form clusters", f"{form} | genes with a kept model of the same gene in another form",
            row["same_gene_in_another_form"])
        add("one form clusters", f"{form} | genes whose RefSeq model is partial or has an exception",
            row["refseq_model_not_clean"])
        add("one form clusters", f"{form} | genes whose RefSeq model is partial or has an exception (%)",
            pct(row["refseq_model_not_clean"], row["genes"]))
        add("one form clusters", f"{form} | genes OrthoFinder left unassigned",
            row["orthofinder_unassigned"])
        add("one form clusters", f"{form} | genes in OrthoFinder cloud orthogroups",
            row["orthofinder_cloud"])
    write_tsv(os.path.join(args.out, "anvio_one_form.tsv"), list(one_rows[0]),
              [list(r.values()) for r in one_rows])
    ref_all = [(f, t) for (f, _), t in gmap.items() if f == REF]
    add("one form clusters", f"{REF} | all kept models whose RefSeq model is partial or has an exception (%)",
        pct(sum(1 for g in ref_all if ann[g][2]), len(ref_all)))
    one = [(f, tx) for f, tx, _ in genes if klass((f, tx)) == "one form"]
    add("one form clusters", "all forms | genes", len(one))
    add("one form clusters", "all forms | genes with a kept model of the same gene in another form",
        sum(1 for g in one if kept[ann[g][0]] - {g[0]}))
    tr = [g for g in one if g[0] != REF]
    add("one form clusters", "transferred forms | genes", len(tr))
    add("one form clusters", "transferred forms | genes with a kept model of the same gene in another form",
        sum(1 for g in tr if kept[ann[g][0]] - {g[0]}))
    add("one form clusters", "transferred forms | genes without a valid ORF",
        sum(1 for g in tr if ann[g][1] == "False"))
    add("one form clusters", "transferred forms | genes without a valid ORF (%)",
        pct(sum(1 for g in tr if ann[g][1] == "False"), len(tr)))

    # ---- broken transferred models and the model of their reference gene:
    # how often each method puts them apart (a model left out of gene clusters,
    # or unassigned by OrthoFinder, counts as apart)
    ref_model = {ann[(f, t)][0]: (f, t) for (f, _), t in gmap.items() if f == REF}
    for label, flag in (("without a valid ORF", "False"), ("with a valid ORF", "True")):
        models = [(f, t) for (f, _), t in gmap.items()
                  if f != REF and ann[(f, t)][1] == flag and ann[(f, t)][0] in ref_model]
        apart_c = sum(1 for g in models
                      if g not in cluster_of or ref_model[ann[g][0]] not in cluster_of
                      or cluster_of[g] != cluster_of[ref_model[ann[g][0]]])
        apart_o = sum(1 for g in models
                      if og[g] == "unassigned" or og[ref_model[ann[g][0]]] != og[g])
        add("reference gene", f"transferred models {label}", len(models))
        add("reference gene", f"transferred models {label} | not in the gene cluster of the reference model",
            apart_c)
        add("reference gene", f"transferred models {label} | not in the gene cluster of the reference model (%)",
            pct(apart_c, len(models)))
        add("reference gene", f"transferred models {label} | not in the orthogroup of the reference model",
            apart_o)
        add("reference gene", f"transferred models {label} | not in the orthogroup of the reference model (%)",
            pct(apart_o, len(models)))
    # genes of one form clusters that are broken transfers, or reference
    # models with at least one broken transfer of the same gene
    broken_genes = {ann[(f, t)][0] for (f, _), t in gmap.items()
                    if f != REF and ann[(f, t)][1] == "False"}
    ref_one = [g for g in one if g[0] == REF]
    add("one form clusters", f"{REF} | genes with a transferred model without a valid ORF",
        sum(1 for g in ref_one if ann[g][0] in broken_genes))
    either = sum(1 for g in one if (g[0] != REF and ann[g][1] == "False")
                 or (g[0] == REF and ann[g][0] in broken_genes))
    add("one form clusters", "all forms | broken transfers or reference models with a broken transfer",
        either)
    add("one form clusters", "all forms | broken transfers or reference models with a broken transfer (%)",
        pct(either, len(one)))

    # ---- transferred models without a valid ORF, by class
    orf_rows = []
    for form in TRANSFERRED:
        for k in CLASSES:
            g = [(f, tx) for f, tx, _ in genes if f == form and klass((f, tx)) == k]
            bad = sum(1 for x in g if ann[x][1] == "False")
            orf_rows.append([form, k, len(g), bad, pct(bad, len(g))])
            add("models without a valid ORF", f"{form} | {k} (%)", pct(bad, len(g)))
    write_tsv(os.path.join(args.out, "anvio_orf_by_class.tsv"),
              ["form", "anvio_class", "n_genes", "n_without_valid_orf", "pct_without_valid_orf"],
              orf_rows)

    # ---- against OrthoFinder
    xt = collections.Counter((klass((f, tx)), compartment((f, tx))) for f, tx, _ in genes)
    if any(c == "missing" for _, c in xt):
        sys.exit("orthogroups missing from the partition table")
    write_tsv(os.path.join(args.out, "anvio_vs_orthofinder.tsv"),
              ["anvio_class"] + COMPARTMENTS + ["total"],
              [[k] + [xt.get((k, c), 0) for c in COMPARTMENTS]
               + [sum(xt.get((k, c), 0) for c in COMPARTMENTS)] for k in CLASSES])
    for k in CLASSES:
        n = sum(xt.get((k, c), 0) for c in COMPARTMENTS)
        for c in COMPARTMENTS:
            add("genes by class and compartment", f"{k} | {c}", xt.get((k, c), 0))
            add("genes by class and compartment", f"{k} | {c} (%)", pct(xt.get((k, c), 0), n))
    same = sum(n for (k, c), n in xt.items()
               if (k, c) in {("all four", "core"), ("two or three", "shell"),
                             ("one form", "cloud"), ("one form", "unassigned")})
    add("against OrthoFinder", "genes in the corresponding class and compartment", same)
    add("against OrthoFinder", "genes in the corresponding class and compartment (%)",
        pct(same, len(genes)))

    og_members = collections.defaultdict(set)
    for f, tx, _ in genes:
        if og[(f, tx)] != "unassigned":
            og_members[og[(f, tx)]].add((f, tx))
    add("against OrthoFinder", "orthogroups with ingroup genes in gene clusters", len(og_members))
    exact = sum(1 for m in members.values()
                if og[m[0]] != "unassigned" and og_members[og[m[0]]] == set(m))
    add("against OrthoFinder", "gene clusters identical to an orthogroup (ingroup genes)", exact)
    add("against OrthoFinder", "gene clusters identical to an orthogroup (%)", pct(exact, len(members)))
    add("against OrthoFinder", "orthogroups whose genes fall in more than one gene cluster",
        sum(1 for m in og_members.values() if len({cluster_of[g] for g in m}) > 1))
    add("against OrthoFinder", "gene clusters holding genes of more than one orthogroup",
        sum(1 for m in members.values()
            if len({og[g] for g in m if og[g] != "unassigned"}) > 1))
    core = [o for o in og_members if comp[o] == "core"]
    add("against OrthoFinder", "core orthogroups", len(core))
    add("against OrthoFinder", "core orthogroups whose genes are all in four form gene clusters",
        sum(1 for o in core if all(len(forms_of[cluster_of[g]]) == 4 for g in og_members[o])))
    four = [c for c in members if len(forms_of[c]) == 4]
    add("against OrthoFinder", "four form gene clusters", len(four))
    add("against OrthoFinder", "four form gene clusters whose genes are all in core orthogroups",
        sum(1 for c in four if all(compartment(g) == "core" for g in members[c])))

    vpath = os.path.join(args.out, "versions.tsv")      # written by anvio_pangenome.sh
    if os.path.exists(vpath):
        with open(vpath) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                add("versions", r["tool"], r["version"])
    add("sources", "gene cluster table", os.path.relpath(src))
    add("sources", "genes and ORF flags", ann_src if ann_src == "annotation GFFs"
        else os.path.basename(ann_src))
    write_tsv(os.path.join(args.out, "anvio_summary.tsv"), ["section", "item", "value"], S)
    for s, i, v in S:
        print(f"{s:30s} {i:78s} {v}")


if __name__ == "__main__":
    main()
