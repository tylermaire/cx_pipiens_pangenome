rule eggnog_download:
    """Download eggNOG database (one-time)."""
    output: directory("resources/eggnog_data/eggnog.db")
    conda: "../envs/eggnog.yaml"
    shell:
        """
        mkdir -p resources/eggnog_data
        download_eggnog_data.py --data_dir resources/eggnog_data/ -y
        """

rule eggnog_mapper:
    """Functional annotation with eggNOG-mapper."""
    input:
        proteins="results/proteins/{sample}.fa",
        db="resources/eggnog_data/eggnog.db"
    output: "results/eggnog/{sample}.emapper.annotations"
    threads: 8
    conda: "../envs/eggnog.yaml"
    shell:
        """
        mkdir -p results/eggnog
        emapper.py -i {input.proteins} \
            --output results/eggnog/{wildcards.sample} \
            --data_dir resources/eggnog_data/ \
            --cpu {threads} --override
        """

rule go_enrichment:
    """GO enrichment analysis: core vs shell vs cloud."""
    input:
        annotations=expand("results/eggnog/{s}.emapper.annotations",
            s=INGROUP_SAMPLES),
        partitions="results/pangenome/partitioned_orthogroups.tsv"
    output:
        enrichment="results/functional/go_enrichment_results.tsv"
    script: "../scripts/go_enrichment.R"
