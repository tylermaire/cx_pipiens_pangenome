rule busco:
    input: "resources/genomes/{sample}.fasta"
    output: directory("results/busco/{sample}")
    params:
        lineage=config["busco"]["lineage"],
        mode=config["busco"]["mode"]
    threads: config["threads"]
    conda: "../envs/busco.yaml"
    shell:
        """
        busco -i {input} -l {params.lineage} -o {wildcards.sample} \
            -m {params.mode} -c {threads} \
            --out_path results/busco/ --force
        """
 
rule quast:
    input: expand("resources/genomes/{s}.fasta", s=ALL_SAMPLES)
    output: directory("results/quast")
    params: labels=",".join(ALL_SAMPLES)
    conda: "../envs/quast.yaml"
    shell:
        """
        quast {input} -o {output} --labels {params.labels} \
            --min-contig 1000
        """
