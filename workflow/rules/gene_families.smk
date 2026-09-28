rule prepare_cafe_input:
    """Gene counts of the ingroup forms and an ultrametric tree for CAFE5.
    The topology is the species tree rooted with the outgroup (V5); earlier
    runs rooted the unrooted four taxon tree between its two pairs, which
    the outgroup now tests instead of assuming. Branch lengths are assumed."""
    input:
        counts="results/orthofinder/output",
        tree="results/phylo/concat_tree.treefile",
        rooted="results/phylo/rooted/rooted_tree.treefile"
    output:
        counts="results/cafe/gene_counts_filtered.tsv",
        tree="results/cafe/ultrametric_tree.nwk"
    params:
        ingroup=INGROUP_SAMPLES,
        outgroup=OUTGROUP_SAMPLES
    conda: "../envs/phylo.yaml"
    script: "../scripts/format_cafe_input.py"

rule cafe5:
    """Run CAFE5 gene family evolution analysis."""
    input:
        counts="results/cafe/gene_counts_filtered.tsv",
        tree="results/cafe/ultrametric_tree.nwk"
    output: directory("results/cafe/output")
    params:
        k=config["cafe"]["n_gamma_categories"]
    conda: "../envs/cafe.yaml"
    shell:
        """
        # No silent fallback: a gamma failure (usually memory) must stop
        # the run rather than quietly producing Base-model results that
        # are not comparable to anything else.
        cafe5 -i {input.counts} -t {input.tree} \
            -p -k {params.k} -o {output}
        test -s {output}/Gamma_family_results.txt || {{
            echo "ERROR: CAFE5 gamma model produced no results" >&2
            exit 1
        }}
        """

CAFE_PVALUES = str(config["cafe"].get("pvalues", "cafe")).lower()

rule cafe_exact_pvalues:
    """Exact family P values under the fitted CAFE5 model.

    CAFE5 estimates each family's P value from 1,000 families it simulates
    with an unseeded random number generator, so the estimates change from run
    to run and every family with the same counts shares one estimate. On the
    four taxon tree the distribution can be summed exactly; this rule does so
    with CAFE5's definition and model (workflow/scripts/cafe_exact_pvalues.py)."""
    input:
        counts="results/cafe/gene_counts_filtered.tsv",
        tree="results/cafe/ultrametric_tree.nwk",
        cafe_dir="results/cafe/output"
    output:
        table="results/cafe/family_pvalues.tsv"
    conda: "../envs/phylo.yaml"
    script: "../scripts/cafe_exact_pvalues.py"

rule parse_cafe_results:
    """Extract significantly evolving gene families, by the exact P values
    (config cafe: pvalues: exact) or by CAFE5's estimates ("cafe")."""
    input:
        cafe_dir="results/cafe/output",
        pvalues=(["results/cafe/family_pvalues.tsv"] if CAFE_PVALUES == "exact" else [])
    output:
        significant="results/cafe/significant_families.tsv",
        summary="results/cafe/branch_summary.tsv"
    params:
        pvalue=config["cafe"]["pvalue_threshold"],
        source=CAFE_PVALUES
    conda: "../envs/phylo.yaml"
    script: "../scripts/parse_cafe.py"

rule cafe_transfer_bias:
    """Control: is the CAFE per-lineage result separable from Liftoff bias?

    Gene models outside the reference are transferred copies, and transfer
    loses copies in proportion to family size. This rule measures that and
    reports whether the families CAFE calls significant are the same families
    that lose copies, so the confound is quantified rather than assumed."""
    input:
        counts="results/cafe/gene_counts_filtered.tsv",
        sig="results/cafe/significant_families.tsv",
        branch="results/cafe/branch_summary.tsv",
        cafe_dir="results/cafe/output"
    output:
        summary="results/cafe/transfer_bias_summary.tsv",
        bins="results/cafe/transfer_bias_by_copy_number.tsv",
        lineage="results/cafe/transfer_bias_by_lineage.tsv"
    params:
        reference=config["reference"]["name"],
        ingroup=INGROUP_SAMPLES
    conda: "../envs/phylo.yaml"
    script: "../scripts/cafe_transfer_bias.py"
