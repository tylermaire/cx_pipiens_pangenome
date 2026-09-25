# -------------------------------------------------------------------------
# Manuscript figures, tables and supplement
# -------------------------------------------------------------------------
# Built from results/ by the same scripts that can be run by hand from the
# repository root (see each script's header). make_all_figures.R and
# make_figures*.py drew the figures of the first submission; they are kept for
# reference but are no longer part of the workflow. Figure 7 of that
# submission (anvi'o) is made in the anvi'o interface, outside the workflow.

REVISION_FIGURES = ["Figure_1_pangenome", "Figure_2_phylogeny", "Figure_3_copy_number",
                    "Figure_4_synteny"]

rule revision_figures:
    """Figures 1 to 4 of the revised manuscript (vector PDF, 1200 dpi TIFF and
    a PNG preview each). Figure 2 shows the species tree rooted with the
    outgroup and the D statistics when those outputs exist."""
    input:
        "results/pangenome/pangenome_summary.tsv",
        "results/pangenome/partitioned_orthogroups.tsv",
        "results/orthofinder/output",
        "results/annotation/transfer_quality_by_compartment.tsv",
        "results/validation/absence_summary.tsv",
        "results/phylo/quartet_topology_counts.tsv",
        "results/phylo/quartet_asymmetry_by_support.tsv",
        "results/phylo/quartet_site_patterns.tsv",
        "results/phylo/quartet_robustness.tsv",
        "results/phylo/rooted/rooted_tree.treefile",
        "results/phylo/rooted/rooted_concord.cf.tree",
        "results/phylo/dstat/d_statistics.tsv",
        "results/cafe/transfer_bias_by_copy_number.tsv",
        "results/synteny/chromosome_orientation.tsv",
        expand("results/annotation/{s}_liftoff.gff3", s=INGROUP_SAMPLES)
    output:
        expand("figures/revision/{f}.{ext}", f=REVISION_FIGURES, ext=["pdf", "tif"])
    conda: "../envs/figures.yaml"
    shell: "python workflow/scripts/make_revision_figures.py --dpi 1200"

rule manuscript_tables:
    """Tables 1 to 4 (JSON read by the manuscript build) and the
    supplementary tables workbook (S1 to S12)."""
    input:
        "results/manuscript_values.tsv",
        "results/pangenome/partitioned_orthogroups.tsv",
        "results/pangenome/cloud_composition.tsv",
        "results/validation/absence_calls.tsv",
        "results/annotation/transfer_quality.tsv",
        "results/phylo/locus_accounting.tsv",
        "results/phylo/sco_pairwise_identity.tsv",
        "results/phylo/rooted/rooted_summary.tsv",
        "results/phylo/rooted/rooted_topology_counts.tsv",
        "results/phylo/dstat/d_statistics.tsv",
        "results/cafe/significant_families.tsv",
        "results/cafe/transfer_bias_summary.tsv",
        "results/synteny/ani_pairs.tsv",
        "results/synteny/synteny_summary.tsv",
        "results/synteny/chromosome_orientation.tsv",
        "results/synteny/inversions_detail.tsv",
        "results/synteny/inversion_recurrence.tsv",
        "results/functional/key_families_wide.tsv",
        "results/functional/key_families_long.tsv",
        expand("results/synteny/genomes/{s}.chromosomes.fasta", s=INGROUP_SAMPLES)
    output:
        json="tables/manuscript_tables.json",
        xlsx="tables/Supplementary_Tables.xlsx"
    conda: "../envs/figures.yaml"
    shell:
        """
        python workflow/scripts/make_manuscript_tables.py --out {output.json}
        python workflow/scripts/make_supplement.py {output.xlsx}
        """
