#!/usr/bin/env python3
"""
collect_manuscript_values.py - every number the manuscript cites, in one table.

Each value is read from the workflow output that produces it and written with
its source path, so the text can be checked against a single file after every
run instead of against a dozen. Missing inputs are written as NA with the path
that was expected, never skipped silently.

Output columns: section, item, sample, value, source

Snakemake provides:
    params.samples, params.ingroup, params.reference, params.parameters (dict
    of run parameters to record verbatim), params.tools (package names whose
    installed versions are read from .snakemake/conda)
    output[0]
"""

import csv
import glob
import json
import os
import re
import sys

ROWS = []


def add(section, item, value, source, sample=""):
    if value is None or (isinstance(value, float) and value != value):
        value = "NA"
    ROWS.append({"section": section, "item": item, "sample": sample,
                 "value": value, "source": source})


def read_tsv(path):
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        lines = [l for l in fh if l.strip() and not l.startswith("#")]
    if not lines:
        return []
    return list(csv.DictReader(lines, delimiter="\t"))


def flatten(section, path, rows):
    """Add every cell of a table, labelled by its first column, or by the
    first two when the first alone is not unique (pairwise tables)."""
    if not rows:
        return
    keys = list(rows[0])
    firsts = [str(r[keys[0]]) for r in rows]
    n_key = 1 if len(set(firsts)) == len(firsts) or len(keys) < 3 else 2
    for r in rows:
        label = " vs ".join(str(r[k]) for k in keys[:n_key])
        for k in keys[n_key:]:
            if r[k] == "":          # sparse tables leave inapplicable cells blank
                continue
            add(section, f"{label} | {k}", r[k], path)


def count_fasta(path):
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


# ---------------------------------------------------------------- assemblies
def quast(samples):
    path = "results/quast/report.tsv"
    rows = read_tsv(path)
    if rows is None:
        add("assembly", "quast report", None, path)
        return
    wanted = {"Total length": "total_length_bp", "# contigs": "contigs",
              "N50": "N50_bp", "GC (%)": "GC_pct", "# N's per 100 kbp": "Ns_per_100kbp"}
    for r in rows:
        key = r.get("Assembly")
        if key in wanted:
            for s in samples:
                add("assembly", wanted[key], r.get(s), path, s)


def busco(samples):
    for mode, root in (("genome", "results/busco"), ("protein", "results/busco_proteins")):
        for s in samples:
            hits = glob.glob(os.path.join(root, s, "short_summary.specific.*.json"))
            if not hits:
                add("busco", f"{mode} summary", None, os.path.join(root, s), s)
                continue
            d = json.load(open(hits[0]))
            res = d.get("results", {})
            line = res.get("one_line_summary", "")
            add("busco", f"{mode} one_line_summary", line, hits[0], s)
            for label, pat in (("complete_pct", r"C:([\d.]+)%"), ("single_pct", r"S:([\d.]+)%"),
                               ("duplicated_pct", r"D:([\d.]+)%"), ("fragmented_pct", r"F:([\d.]+)%"),
                               ("missing_pct", r"M:([\d.]+)%")):
                m = re.search(pat, line)
                add("busco", f"{mode} {label}", m.group(1) if m else None, hits[0], s)
            v = d.get("versions", {})
            add("busco", f"{mode} busco_version", v.get("busco"), hits[0], s)
            add("busco", f"{mode} lineage", d.get("lineage_dataset", {}).get("name"), hits[0], s)
            if mode == "genome":
                add("busco", f"{mode} gene_predictor",
                    d.get("parameters", {}).get("gene_predictor"), hits[0], s)


def proteins(samples):
    for s in samples:
        path = f"results/proteins/{s}.fa"
        add("annotation", "proteins_one_per_gene", count_fasta(path), path, s)
    path = "results/annotation/transfer_quality.tsv"
    rows = read_tsv(path)
    if rows is None:
        add("annotation", "transfer quality", None, path)
    else:
        for r in rows:
            for k in ("n_kept_models", "liftoff_flags", "n_invalid_orf", "pct_invalid_orf",
                      "pct_inframe_stop", "pct_missing_start", "pct_missing_stop",
                      "pct_mismatch_ref_protein", "n_ref_model_not_clean",
                      "pct_invalid_orf_clean_ref"):
                if r.get(k) == "":      # Liftoff flags do not apply to the reference
                    continue
                add("annotation", k, r.get(k), path, r.get("sample", ""))
    path = "results/annotation/transfer_quality_by_compartment.tsv"
    rows = read_tsv(path)
    if rows is None:
        add("annotation", "transfer quality by compartment", None, path)
    else:
        for r in rows:
            add("annotation", f"{r['compartment']} pct_invalid_orf",
                r.get("pct_invalid_orf"), path, r.get("sample", ""))


# ---------------------------------------------------------------- pangenome
def pangenome():
    for path, section in (("results/pangenome/pangenome_summary.tsv", "pangenome"),):
        rows = read_tsv(path)
        if rows is None:
            add(section, "summary", None, path)
            continue
        for r in rows:
            add(section, f"{r['compartment']} n_orthogroups", r.get("n_orthogroups"), path)
            add(section, f"{r['compartment']} genes_ingroup", r.get("total_genes_ingroup"), path)
    path = "results/pangenome/cloud_composition_summary.tsv"
    rows = read_tsv(path)
    if rows is None:
        add("cloud", "composition", None, path)
    else:
        for r in rows:
            for k, v in r.items():
                if k != "form":
                    add("cloud", k, v, path, r["form"])
    path = "results/validation/absence_summary.tsv"
    rows = read_tsv(path)
    if rows is None:
        add("absence", "summary", None, path)
    else:
        for r in rows:
            for k, v in r.items():
                if k != "compartment":
                    add("absence", f"{r['compartment']} {k}", v, path)


# ---------------------------------------------------------------- phylogeny
def phylogeny():
    add("phylogeny", "sco_fastas", len(glob.glob("results/phylo/sco_fastas/*.fa")),
        "results/phylo/sco_fastas")
    trims = sorted(glob.glob("results/phylo/trimmed/*.trim"))
    add("phylogeny", "trimmed_alignments", len(trims), "results/phylo/trimmed")
    sites = cells = gaps = 0
    for path in trims:
        seqs, name = {}, None
        for line in open(path):
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
        seqs = ["".join(v) for v in seqs.values()]
        if not seqs:
            continue
        sites += len(seqs[0])
        cells += sum(len(x) for x in seqs)
        gaps += sum(x.count("-") + x.upper().count("X") + x.count("?") for x in seqs)
    add("phylogeny", "supermatrix_sites", sites if trims else None, "results/phylo/trimmed")
    add("phylogeny", "supermatrix_missing_pct",
        round(100.0 * gaps / cells, 4) if cells else None, "results/phylo/trimmed")

    path = "results/phylo/concord.cf.stat"
    rows = read_tsv(path)
    if rows:
        for k in ("gCF", "gCF_N", "gDF1", "gDF1_N", "gDF2", "gDF2_N", "gN",
                  "sCF", "sDF1", "sDF2", "sN", "Length"):
            add("phylogeny", k, rows[0].get(k), path)
    else:
        add("phylogeny", "concordance", None, path)

    rooted()

    for path, section in (("results/phylo/quartet_topology_counts.tsv", "quartet"),
                          ("results/phylo/quartet_asymmetry_by_support.tsv", "quartet_support"),
                          ("results/phylo/quartet_site_patterns.tsv", "quartet_sites"),
                          ("results/phylo/quartet_robustness.tsv", "quartet_robustness"),
                          ("results/phylo/branch_length_summary.tsv", "branch_lengths"),
                          ("results/phylo/sco_pairwise_identity.tsv", "sco_identity"),
                          ("results/phylo/locus_accounting.tsv", "loci")):
        rows = read_tsv(path)
        if rows is None:
            add(section, "table", None, path)
            continue
        flatten(section, path, rows)


def rooted():
    """The outgroup analyses (V5): the rooted species tree, rooted gene tree
    topologies and the D statistics."""
    path = "results/phylo/rooted/rooted_summary.tsv"
    rows = read_tsv(path)
    if rows is None:
        add("rooted", "summary", None, path)
    else:
        for r in rows:
            add("rooted", r["item"], r["value"], path)
    for path, section in (("results/phylo/rooted/rooted_topology_counts.tsv", "rooted_topologies"),
                          ("results/phylo/dstat/d_statistics.tsv", "dstat")):
        rows = read_tsv(path)
        if rows is None:
            add(section, "table", None, path)
            continue
        flatten(section, path, rows)


# ---------------------------------------------------------------- CAFE
def cafe():
    path = "results/cafe/output/Gamma_results.txt"
    if os.path.exists(path):
        text = open(path).read()
        for label, pat in (("lambda", r"Lambda:\s*([\d.eE+-]+)"),
                           ("alpha", r"Alpha:\s*([\d.eE+-]+)"),
                           ("neg_log_likelihood", r"-lnL\):\s*([\d.eE+-]+)")):
            m = re.search(pat, text)
            add("cafe", label, m.group(1) if m else None, path)
    else:
        add("cafe", "model results", None, path)
    sig = read_tsv("results/cafe/significant_families.tsv")
    add("cafe", "significant_families", len(sig) if sig is not None else None,
        "results/cafe/significant_families.tsv")
    for path, section in (("results/cafe/branch_summary.tsv", "cafe_branches"),
                          ("results/cafe/transfer_bias_summary.tsv", "transfer_bias"),
                          ("results/cafe/transfer_bias_by_copy_number.tsv", "transfer_bias_bins"),
                          ("results/cafe/transfer_bias_by_lineage.tsv", "transfer_bias_lineage")):
        rows = read_tsv(path)
        if rows is None:
            add(section, "table", None, path)
            continue
        flatten(section, path, rows)


# ---------------------------------------------------------------- ANI, synteny, TE
def ani_synteny_te(ingroup):
    for path, section in (("results/synteny/ani_pairs.tsv", "ani"),
                          ("results/synteny/synteny_summary.tsv", "synteny"),
                          ("results/synteny/inversion_recurrence.tsv", "inversion_recurrence")):
        rows = read_tsv(path)
        if rows is None:
            add(section, "table", None, path)
            continue
        flatten(section, path, rows)
    inv = read_tsv("results/synteny/inversions_detail.tsv")
    if inv:
        spans = [int(r["span_bp"]) for r in inv]
        add("synteny", "inversions_total", len(spans), "results/synteny/inversions_detail.tsv")
        add("synteny", "inversions_at_least_1Mb", sum(s >= 1_000_000 for s in spans),
            "results/synteny/inversions_detail.tsv")
        big = max(inv, key=lambda r: int(r["span_bp"]))
        add("synteny", "largest_inversion",
            f"{big['sample1']} vs {big['sample2']} {big.get('chromosome', '')} "
            f"{big['seqid']}:{big['start_bp']}-{big['end_bp']} ({big['n_genes']} genes)",
            "results/synteny/inversions_detail.tsv")
    for s in ingroup:
        path = f"results/repeats/{s}/{s}.fasta.tbl"
        if not os.path.exists(path):
            add("repeats", "tbl", None, path, s)
            continue
        text = open(path).read()
        for label, pat in (("bases_masked_pct", r"bases masked:\s*[\d,]+ bp\s*\(\s*([\d.]+) %\)"),
                           ("interspersed_pct", r"Total interspersed repeats:\s*[\d,]+ bp\s*([\d.]+) %"),
                           ("simple_repeats_pct", r"Simple repeats:\s*\d+\s*[\d,]+ bp\s*([\d.]+) %"),
                           ("low_complexity_pct", r"Low complexity:\s*\d+\s*[\d,]+ bp\s*([\d.]+) %")):
            m = re.search(pat, text)
            add("repeats", label, m.group(1) if m else None, path, s)


# ---------------------------------------------------------------- tools and parameters
def tools(names):
    found = {}
    for meta in glob.glob(".snakemake/conda/*/conda-meta/*.json"):
        base = os.path.basename(meta)
        m = re.match(r"(.+?)-(\d[^-]*)-[^-]+\.json$", base)
        if m and m.group(1) in names:
            found.setdefault(m.group(1), set()).add(m.group(2))
    if not glob.glob(".snakemake/conda/*/conda-meta"):
        add("tools", "conda environments", "not present; see tools_observed",
            ".snakemake/conda/*/conda-meta")
        return
    for n in names:
        add("tools", n, ",".join(sorted(found[n])) if n in found else None,
            ".snakemake/conda/*/conda-meta")


def versions_from_outputs():
    """Versions and command lines the outputs themselves record. These are
    what ran, whatever the conda environment files ask for."""
    def first_match(paths, pattern, flags=0):
        for path in paths:
            if not os.path.exists(path):
                continue
            with open(path, errors="replace") as fh:
                head = fh.read(200_000)
            m = re.search(pattern, head, flags)
            if m:
                return m.group(1).strip(), path
        return None, (paths[0] if paths else "")

    checks = [
        ("iqtree", ["results/phylo/concat_tree.iqtree"], r"IQ-TREE (?:multicore version )?(\d[\w.\-]*)"),
        ("orthofinder", sorted(glob.glob("results/orthofinder/output/*/Log.txt")),
         r"OrthoFinder version (\d[\w.\-]*)"),
        ("liftoff", sorted(glob.glob("results/annotation/*_liftoff.gff3")), r"# Liftoff v(\d[\w.\-]*)"),
        ("liftoff command", sorted(glob.glob("results/annotation/*_liftoff.gff3")),
         r"^# \S*liftoff (.+)$"),
        ("repeatmasker", sorted(glob.glob("results/repeats/*/*.tbl")),
         r"RepeatMasker version (\S+)"),
        ("rmblastn", sorted(glob.glob("results/repeats/*/*.tbl")), r"rmblastn version (\S+)"),
        ("skani command", ["results/synteny/ani_pairs.tsv"], r"^# (skani .+)$"),
    ]
    for name, paths, pattern in checks:
        value, source = first_match(paths, pattern, re.M)
        add("tools_observed", name, value, source)


def parameters(params):
    for k, v in params.items():
        add("parameters", k, v, "config/config.yaml and workflow rules")


def main():
    if "snakemake" not in globals():
        sys.exit("Run via Snakemake.")
    sm = globals()["snakemake"]
    samples = list(sm.params.samples)
    ingroup = list(sm.params.ingroup)
    quast(samples)
    busco(samples)
    proteins(samples)
    pangenome()
    phylogeny()
    cafe()
    ani_synteny_te(ingroup)
    tools(list(sm.params.tools))
    versions_from_outputs()
    parameters(dict(sm.params.parameters))
    with open(sm.output[0], "w", newline="") as fh:
        w = csv.DictWriter(fh, lineterminator="\n", fieldnames=["section", "item", "sample", "value", "source"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(ROWS)
    missing = sum(1 for r in ROWS if r["value"] == "NA")
    print(f"{len(ROWS)} values written, {missing} NA")


if __name__ == "__main__":
    main()
