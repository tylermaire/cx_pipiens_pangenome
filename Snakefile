"""
Culex pipiens Complex Pangenome Analysis Pipeline
==================================================
Authors: Tyler Maire and collaborators
Date: 2026
License: MIT
"""
import pandas as pd
import os
import shutil
from itertools import combinations
 
# -- Load configuration --
configfile: "config/config.yaml"
 
# -- Load sample table --
samples = pd.read_csv(config["samples"], sep="\t", index_col="sample")
ALL_SAMPLES = list(samples.index)
INGROUP_SAMPLES = list(samples[samples["is_outgroup"] == False].index)
OUTGROUP_SAMPLES = list(samples[samples["is_outgroup"] == True].index)
 
# -- Pairwise combinations for synteny (ingroup only) --
PAIRS = list(combinations(INGROUP_SAMPLES, 2))
 
# -- Reference name shortcut --
REF = config["reference"]["name"]
 
# -- Include rule modules --
include: "workflow/rules/download.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/annotation.smk"
include: "workflow/rules/orthology.smk"
include: "workflow/rules/functional.smk"
include: "workflow/rules/phylogenomics.smk"
include: "workflow/rules/gene_families.smk"
include: "workflow/rules/synteny.smk"
include: "workflow/rules/repeats.smk"
include: "workflow/rules/figures.smk"
 
# -- Default target: run the entire pipeline --
rule all:
    input:
        # QC
        expand("results/busco/{s}", s=ALL_SAMPLES),
        "results/quast",
        # Annotation
        expand("results/proteins/{s}.fa", s=ALL_SAMPLES),
        expand("results/busco_proteins/{s}", s=ALL_SAMPLES),
        # Pangenome
        "results/pangenome/pangenome_summary.tsv",
        # Functional
        "results/functional/go_enrichment_results.tsv",
        # Phylogenomics
        "results/phylo/concord.cf.tree",
        # Gene families
        "results/cafe/significant_families.tsv",
        # Synteny
        expand("results/synteny/{a}_vs_{b}_syri.out",
            zip, a=[p[0] for p in PAIRS],
            b=[p[1] for p in PAIRS]),
        # Repeats
        expand("results/repeats/{s}/{s}.fasta.tbl",
            s=INGROUP_SAMPLES),
