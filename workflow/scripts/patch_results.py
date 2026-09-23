#!/usr/bin/env python3
"""
patch_results.py - apply post-processing fixes to a finished run without
rerunning the pipeline.

Every step reads files the full workflow already produced and runs the same
script its Snakemake rule runs, with the same inputs and parameters. Nothing
upstream is touched: no Liftoff, OrthoFinder, alignment, IQ-TREE, CAFE5,
minimap2, miniprot or RepeatMasker. The full workflow stays the reference and
produces the same outputs from scratch; this exists so a fix to a summary
script does not cost a cloud rerun.

Keep the inputs and params below in step with workflow/rules/*.smk.

Steps, in default order:
  absence     classify_absence.py          (rule classify_absence)
  cloud       cloud_composition.py         (rule cloud_composition)
  quartet     quartet_asymmetry.py         (rule quartet_asymmetry)
  divergence  divergence_diagnostics.py    (rule divergence_diagnostics)
  cafe        parse_cafe.py, then cafe_transfer_bias.py
  synteny     synteny_summary.py           (rule synteny_summary)
  values      collect_manuscript_values.py (rule manuscript_values)

Usage, from the repository root with results/ holding the run:
    python workflow/scripts/patch_results.py
    python workflow/scripts/patch_results.py --steps absence,cloud
"""

import argparse
import csv
import os
import runpy
import sys
import time
from types import SimpleNamespace

SCRIPTS = os.path.join("workflow", "scripts")
STEPS = ["absence", "cloud", "quartet", "divergence", "cafe", "synteny", "values"]


def read_config(path="config/config.yaml"):
    try:
        import yaml
        return yaml.safe_load(open(path))
    except ImportError:
        pass
    # minimal reader for this file's two level layout
    cfg, section = {}, None
    for raw in open(path):
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        key, _, val = line.strip().partition(":")
        val = val.strip()
        if not raw.startswith(" "):
            section = key
            cfg[key] = _scalar(val) if val else {}
        elif isinstance(cfg.get(section), dict):
            cfg[section][key] = _scalar(val)
    return cfg


def _scalar(v):
    for cast in (int, float):
        try:
            return cast(v)
        except ValueError:
            pass
    return v


def samples(path="config/samples.tsv"):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    all_s = [r["sample"] for r in rows]
    ingroup = [r["sample"] for r in rows if r["is_outgroup"].strip().lower() == "false"]
    return all_s, ingroup


def run(script, **sections):
    sm = SimpleNamespace(**{k: SimpleNamespace(**v) if isinstance(v, dict) else v
                            for k, v in sections.items()})
    t = time.time()
    print(f"\n=== {script} ===", flush=True)
    runpy.run_path(os.path.join(SCRIPTS, script), init_globals={"snakemake": sm},
                   run_name="__main__")
    print(f"=== {script} done in {time.time() - t:.0f} s ===", flush=True)


def step_absence(cfg, all_s, ingroup, ref):
    run("classify_absence.py",
        input={"manifest": "results/validation/queries/manifest.tsv",
               "of": "results/orthofinder/output",
               "prot_pafs": [f"results/validation/paf/prot_vs_{t}.paf" for t in ingroup],
               "dna_pafs": [f"results/validation/paf/dna_vs_{t}.paf" for t in ingroup],
               "gffs": [f"results/annotation/{s}_liftoff.gff3" for s in ingroup],
               "proteins": [f"results/proteins/{s}.fa" for s in ingroup]},
        params={"ingroup": ingroup, "min_id": 0.80, "min_cov": 0.70,
                "dna_min_cov": 0.50, "same_gene_id": 0.90, "max_overlap": 10},
        output={"calls": "results/validation/absence_calls.tsv",
                "summary": "results/validation/absence_summary.tsv"})


def step_cloud(cfg, all_s, ingroup, ref):
    run("cloud_composition.py",
        input={"table": "results/pangenome/partitioned_orthogroups.tsv",
               "of": "results/orthofinder/output"},
        params={"ingroup": ingroup, "reference": ref},
        output={"per_og": "results/pangenome/cloud_composition.tsv",
                "summary": "results/pangenome/cloud_composition_summary.tsv"})


def step_quartet(cfg, all_s, ingroup, ref):
    run("quartet_asymmetry.py",
        input={"tree": "results/phylo/concat_tree.treefile",
               "concord": "results/phylo/concord.cf.tree"},
        params={"gene_trees": "results/phylo/all_gene_trees.nwk",
                "cf_stat": "results/phylo/concord.cf.stat",
                "trimmed_dir": "results/phylo/trimmed"},
        output={"topology": "results/phylo/quartet_topology_counts.tsv",
                "support": "results/phylo/quartet_asymmetry_by_support.tsv",
                "sites": "results/phylo/quartet_site_patterns.tsv"})


def step_divergence(cfg, all_s, ingroup, ref):
    run("divergence_diagnostics.py",
        input={"tree": "results/phylo/concat_tree.treefile",
               "concord": "results/phylo/concord.cf.tree",
               "sco": "results/phylo/sco_fastas"},
        params={"gene_trees": "results/phylo/all_gene_trees.nwk",
                "trimmed_dir": "results/phylo/trimmed",
                "sco_dir": "results/phylo/sco_fastas"},
        output={"branch_lengths": "results/phylo/branch_length_summary.tsv",
                "pairwise": "results/phylo/sco_pairwise_identity.tsv",
                "loci": "results/phylo/locus_accounting.tsv"})


def step_cafe(cfg, all_s, ingroup, ref):
    run("parse_cafe.py",
        input=["results/cafe/output"],
        params={"pvalue": cfg["cafe"]["pvalue_threshold"]},
        output={"significant": "results/cafe/significant_families.tsv",
                "summary": "results/cafe/branch_summary.tsv"})
    run("cafe_transfer_bias.py",
        input={"counts": "results/cafe/gene_counts_filtered.tsv",
               "sig": "results/cafe/significant_families.tsv",
               "branch": "results/cafe/branch_summary.tsv"},
        params={"reference": ref, "ingroup": ingroup},
        output={"summary": "results/cafe/transfer_bias_summary.tsv",
                "bins": "results/cafe/transfer_bias_by_copy_number.tsv",
                "lineage": "results/cafe/transfer_bias_by_lineage.tsv"})


def step_synteny(cfg, all_s, ingroup, ref):
    pafs = [f"results/synteny/paf/{a}_vs_{b}.paf"
            for i, a in enumerate(ingroup) for b in ingroup[i + 1:]]
    run("synteny_summary.py",
        input={"gffs": [f"results/annotation/{s}_liftoff.gff3" for s in ingroup],
               "pafs": pafs, "ani": "results/synteny/ani_matrix.tsv"},
        params={"samples": ingroup, "reference": ref, "n_chrom": 3,
                "min_span": 100_000, "min_genes": 3, "large_span": 1_000_000},
        output={"summary": "results/synteny/synteny_summary.tsv",
                "orientation": "results/synteny/chromosome_orientation.tsv",
                "inversions": "results/synteny/inversions_detail.tsv",
                "recurrence": "results/synteny/inversion_recurrence.tsv"})


def step_values(cfg, all_s, ingroup, ref):
    params = {
        "samples": all_s, "ingroup": ingroup, "reference": ref,
        "tools": ["liftoff", "gffread", "orthofinder", "diamond", "mafft", "trimal",
                  "iqtree", "cafe", "skani", "minimap2", "miniprot", "busco", "quast",
                  "repeatmodeler", "repeatmasker", "eggnog-mapper", "biopython"],
        "parameters": {
            "liftoff -s (min child feature identity)": cfg["liftoff"]["identity_threshold"],
            "liftoff -sc (copy identity; inert without -copies)": cfg["liftoff"]["copy_identity"],
            "liftoff -copies / -polish": "not used",
            "orthofinder": f"-M {cfg['orthofinder']['method']} -S {cfg['orthofinder']['search']} -a 4",
            "iqtree species tree": f"-p per locus partitions -m {cfg['iqtree']['model']} -bb {cfg['iqtree']['bootstrap']}",
            "iqtree gene trees": "-m MFP -bb 1000, one per trimmed SCO",
            "concordance": "--gcf all_gene_trees.nwk --scf 100",
            "cafe": f"-p -k {cfg['cafe']['n_gamma_categories']}; ML topology rooted between the sister pairs, unit tip lengths (assumed, not estimated)",
            "busco lineage": cfg["busco"]["lineage"],
            "minimap2 synteny": "-x asm20 --secondary=no on the three longest sequences",
            "synteny inversions": "runs of >= 3 anchors spanning >= 100 kb after orientation normalisation",
            "absence validation": "miniprot --outn 1; minimap2 -x asm20 -N 1; hit id >= 0.80, protein cov >= 0.70, DNA cov >= 0.50; same gene protein identity >= 0.90",
        }}
    run("collect_manuscript_values.py", params=params,
        output=["results/manuscript_values.tsv"])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--steps", default=",".join(STEPS),
                    help=f"comma separated, any of {','.join(STEPS)}")
    args = ap.parse_args()
    steps = [s.strip() for s in args.steps.split(",") if s.strip()]
    unknown = [s for s in steps if s not in STEPS]
    if unknown:
        sys.exit(f"unknown steps: {unknown}")
    if not os.path.isdir("results"):
        sys.exit("run from the repository root, with results/ holding the run")
    cfg = read_config()
    all_s, ingroup = samples()
    ref = cfg["reference"]["name"]
    for s in STEPS:
        if s in steps:
            globals()[f"step_{s}"](cfg, all_s, ingroup, ref)


if __name__ == "__main__":
    main()
