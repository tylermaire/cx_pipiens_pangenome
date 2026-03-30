# functional.smk - Functional annotation and GO enrichment
# Owner: Person B
# Tools: eggNOG-mapper, R (topGO/clusterProfiler)
# Outputs: results/eggnog/, results/functional/
 
rule eggnog_mapper:
    input: "results/proteins/{sample}.fa"
    output: "results/eggnog/{sample}.emapper.annotations"
    params: data_dir=config.get("eggnog_data", "resources/eggnog_data")
    threads: config["threads"]
    conda: "../envs/eggnog.yaml"
    shell:
        """
        emapper.py -i {input} --output results/eggnog/{wildcards.sample} \
            --data_dir {params.data_dir} --cpu {threads} --override
        """
 
rule go_enrichment:
    input:
        annotations=expand("results/eggnog/{s}.emapper.annotations",
            s=INGROUP_SAMPLES),
        partitions="results/pangenome/partitioned_orthogroups.tsv"
    output:
        enrichment="results/functional/go_enrichment_results.tsv",
        heatmap="results/figures/go_enrichment_heatmap.pdf"
    script: "../scripts/go_enrichment.R"
