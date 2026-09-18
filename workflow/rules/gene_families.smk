rule prepare_cafe_input:
    """Filter gene counts and create ultrametric tree for CAFE5."""
    input:
        counts="results/orthofinder/output",
        tree="results/phylo/concat_tree.treefile"
    output:
        counts="results/cafe/gene_counts_filtered.tsv",
        tree="results/cafe/ultrametric_tree.nwk"
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

rule parse_cafe_results:
    """Extract significantly evolving gene families."""
    input: "results/cafe/output"
    output:
        significant="results/cafe/significant_families.tsv",
        summary="results/cafe/branch_summary.tsv"
    params:
        pvalue=config["cafe"]["pvalue_threshold"]
    conda: "../envs/phylo.yaml"
    script: "../scripts/parse_cafe.py"
