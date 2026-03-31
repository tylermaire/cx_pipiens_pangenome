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
        cafe5 -i {input.counts} -t {input.tree} \
            -p -k {params.k} -o {output} || \
        cafe5 -i {input.counts} -t {input.tree} \
            -o {output}
        """

rule parse_cafe_results:
    """Extract significantly evolving gene families."""
    input: "results/cafe/output"
    output:
        significant="results/cafe/significant_families.tsv",
        summary="results/cafe/branch_summary.tsv"
    params:
        pvalue=config["cafe"]["pvalue_threshold"]
    script: "../scripts/parse_cafe.py"
