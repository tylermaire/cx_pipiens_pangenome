REF = config["reference"]["name"]

rule liftoff:
    input:
        target="resources/genomes/{sample}.fasta",
        ref=f"resources/genomes/{REF}.fasta",
        gff=f"resources/annotations/{REF}.gff3"
    output: "results/annotation/{sample}_liftoff.gff3"
    params:
        sc=config["liftoff"]["coverage_threshold"],
        s=config["liftoff"]["identity_threshold"]
    threads: config["threads"]
    conda: "../envs/liftoff.yaml"
    shell:
        """
        mkdir -p results/annotation/liftoff_intermediates/{wildcards.sample}
        liftoff -g {input.gff} {input.target} {input.ref} \
            -o {output} -p {threads} -sc {params.sc} -s {params.s} \
            -dir results/annotation/liftoff_intermediates/{wildcards.sample}
        """

rule copy_ref_annotation:
    input: f"resources/annotations/{REF}.gff3"
    output: f"results/annotation/{REF}_liftoff.gff3"
    shell: "cp {input} {output}"

rule extract_proteins:
    input:
        fasta="resources/genomes/{sample}.fasta",
        gff="results/annotation/{sample}_liftoff.gff3"
    output: "results/proteins/{sample}.fa"
    conda: "../envs/gffread.yaml"
    shell:
        """
        gffread {input.gff} -g {input.fasta} -y {output}.tmp
        sed 's/\\.//g' {output}.tmp > {output}
        rm {output}.tmp
        """

rule busco_proteins:
    input: "results/proteins/{sample}.fa"
    output: directory("results/busco_proteins/{sample}")
    params: lineage=config["busco"]["lineage"]
    threads: config["threads"]
    conda: "../envs/busco.yaml"
    shell:
        """
        busco -i {input} -l {params.lineage} -o {wildcards.sample} \
            -m proteins -c {threads} \
            --out_path results/busco_proteins/ --force
        """
