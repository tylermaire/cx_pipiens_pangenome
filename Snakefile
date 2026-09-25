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

# Tools that take a thread count as an argument (IQ-TREE in particular) abort
# rather than scale down when asked for more threads than the machine has.
# Cap the configured value so a run on a smaller box degrades instead of dying.
import os as _os
_avail = _os.cpu_count() or 1
if config["threads"] > _avail:
    print(f"WARNING: config threads={config['threads']} exceeds {_avail} "
          f"available cores; capping to {_avail}")
    config["threads"] = _avail
INGROUP_SAMPLES = list(samples[samples["is_outgroup"] == False].index)
OUTGROUP_SAMPLES = list(samples[samples["is_outgroup"] == True].index)

# -- Annotation and download source of each sample (V5) --
# annotation: reference = the RefSeq annotation Liftoff transfers from;
#             liftoff   = models transferred from the reference;
#             native    = the sample's own gene set (the outgroup, from Ensembl).
# source:     ncbi (datasets CLI) or ensembl (FTP; genome and gene set together).
def _column(name, default):
    if name in samples.columns:
        return samples[name].fillna(default).astype(str).str.strip().to_dict()
    return {s: default for s in ALL_SAMPLES}

ANNOTATION = _column("annotation", "liftoff")
SOURCE = _column("source", "ncbi")
LIFTOFF_SAMPLES = [s for s in ALL_SAMPLES if ANNOTATION[s] == "liftoff"]
NATIVE_SAMPLES = [s for s in ALL_SAMPLES if ANNOTATION[s] == "native"]
NCBI_SAMPLES = [s for s in ALL_SAMPLES if SOURCE[s] == "ncbi"]
ENSEMBL_SAMPLES = [s for s in ALL_SAMPLES if SOURCE[s] == "ensembl"]
OUTGROUP = config["outgroup"]["name"]
if OUTGROUP_SAMPLES != [OUTGROUP]:
    raise ValueError(f"config outgroup {OUTGROUP} does not match the outgroup rows of "
                     f"{config['samples']}: {OUTGROUP_SAMPLES}")


def one_of(names):
    """Wildcard constraint matching exactly these sample names (or nothing)."""
    import re as _re
    return "|".join(_re.escape(n) for n in names) if names else "(?!x)x"

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
include: "workflow/rules/validation.smk"
include: "workflow/rules/report.smk"

# -- Default target: build everything except SyRI (strand-correction needed) --
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
        "results/pangenome/cloud_composition_summary.tsv",
        # Functional
        "results/functional/key_families_wide.tsv",
        # Phylogenomics (already built)
        "results/phylo/concord.cf.tree",
        "results/phylo/quartet_asymmetry_by_support.tsv",
        "results/phylo/branch_length_summary.tsv",
        "results/annotation/transfer_quality.tsv",
        "results/manuscript_values.tsv",
        # Gene families (already rebuilt)
        "results/cafe/significant_families.tsv",
        # Repeats - RepeatMasker tables + TE proximity
        expand("results/repeats/{s}/{s}.fasta.tbl", s=INGROUP_SAMPLES),
        "results/repeats/te_gene_proximity_summary.tsv",
        "results/validation/absence_summary.tsv",
        "results/cafe/transfer_bias_summary.tsv",
        "results/synteny/synteny_summary.tsv",
        # Rooted analyses with the outgroup (V5)
        "results/phylo/quartet_robustness.tsv",
        "results/phylo/rooted/rooted_summary.tsv",
        "results/phylo/dstat/d_statistics.tsv",
        # Manuscript figures, tables and supplement, built from results/
        "figures/revision/Figure_1_pangenome.pdf",
        "tables/manuscript_tables.json",
        "tables/Supplementary_Tables.xlsx",
