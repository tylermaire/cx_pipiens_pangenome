# -------------------------------------------------------------------------
# Publication figures
# -------------------------------------------------------------------------
# make_all_figures.R runs either as this rule or standalone (see its header).
# Figure 7 (anvi'o) is produced in the anvi'o web interface rather than by
# this workflow; the script skips it when anvio_data/ is absent, so the rule
# does not declare it as an output.

rule figures:
    """Regenerate every figure the workflow can produce."""
    input:
        busco=expand("results/busco_proteins/{s}", s=ALL_SAMPLES),
        quast="results/quast",
        pangenome="results/pangenome/pangenome_summary.tsv",
        partition="results/pangenome/partitioned_orthogroups.tsv",
        tree="results/phylo/concord.cf.tree",
        ani="results/synteny/ani_matrix.tsv",
        gffs=expand("results/annotation/{s}_liftoff.gff3", s=INGROUP_SAMPLES),
        repeats=expand("results/repeats/{s}/{s}.fasta.tbl", s=INGROUP_SAMPLES),
        families="results/functional/key_families_wide.tsv"
    output:
        "figures/Figure_1_assembly_quality.png",
        "figures/Figure_2_pangenome.png",
        "figures/Figure_3_species_tree.png",
        "figures/Figure_4_ani.png",
        "figures/Figure_5_synteny.png",
        "figures/Figure_6_te_composition.png",
        "figures/Figure_8_key_families.png"
    params:
        proj="."
    conda: "../envs/r.yaml"
    script: "../scripts/make_all_figures.R"
