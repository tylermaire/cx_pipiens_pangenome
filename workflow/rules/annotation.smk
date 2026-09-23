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
        cp {input.gff} results/annotation/liftoff_intermediates/{wildcards.sample}/ref.gff3
        liftoff -g results/annotation/liftoff_intermediates/{wildcards.sample}/ref.gff3 \
            {input.target} {input.ref} \
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
    output:
        fa="results/proteins/{sample}.fa",
        report="results/proteins/{sample}_isoform_report.tsv"
    conda: "../envs/gffread.yaml"
    shell:
        """
        gffread {input.gff} -g {input.fasta} -y {output.fa}.tmp
        python workflow/scripts/longest_isoform.py \
            --gff {input.gff} --fasta {output.fa}.tmp \
            --out {output.fa} --report {output.report} \
            --strip-periods
        rm {output.fa}.tmp
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

rule transfer_quality:
    """Broken protein models among those kept for analysis: internal stop
    codons, missing start methionine, missing terminal stop. A transfer that
    shifts the reading frame shows up here and nowhere else, because
    extract_proteins strips the stop symbols."""
    input:
        gff="results/annotation/{sample}_liftoff.gff3",
        fasta="resources/genomes/{sample}.fasta",
        kept="results/proteins/{sample}.fa"
    output: "results/annotation/{sample}_transfer_quality.tsv"
    conda: "../envs/gffread.yaml"
    shell:
        """
        gffread {input.gff} -g {input.fasta} -y {output}.raw.faa
        python workflow/scripts/transfer_quality.py --raw {output}.raw.faa \
            --kept {input.kept} --sample {wildcards.sample} --out {output}
        rm -f {output}.raw.faa
        """

rule transfer_quality_summary:
    input: expand("results/annotation/{s}_transfer_quality.tsv", s=ALL_SAMPLES)
    output: "results/annotation/transfer_quality.tsv"
    shell:
        "head -n 1 {input[0]} > {output} && "
        "for f in {input}; do tail -n +2 $f >> {output}; done"
