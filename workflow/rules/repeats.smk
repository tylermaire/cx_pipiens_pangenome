rule repeatmodeler:
    """Build de novo repeat library from reference genome."""
    input:
        ref=f"resources/genomes/{config['reference']['name']}.fasta"
    output: "results/repeats/custom_repeat_lib.fa"
    threads: 16
    conda: "../envs/repeatmasker.yaml"
    shell:
        """
        WORKDIR=$(pwd)
        mkdir -p results/repeats/repeatmodeler_workdir
        cd results/repeats/repeatmodeler_workdir
        BuildDatabase -name culex_db $WORKDIR/{input.ref}
        RepeatModeler -database culex_db -threads {threads} -LTRStruct || \
            RepeatModeler -database culex_db -threads {threads}
        cp culex_db-families.fa $WORKDIR/{output}
        """

rule repeatmasker:
    """Mask genome with custom repeat library."""
    input:
        fasta="resources/genomes/{sample}.fasta",
        lib="results/repeats/custom_repeat_lib.fa"
    output:
        out="results/repeats/{sample}/{sample}.fasta.out",
        tbl="results/repeats/{sample}/{sample}.fasta.tbl"
    threads: 8
    conda: "../envs/repeatmasker.yaml"
    shell:
        """
        mkdir -p results/repeats/{wildcards.sample}
        cp {input.fasta} results/repeats/{wildcards.sample}/{wildcards.sample}.fasta
        RepeatMasker -lib {input.lib} -pa {threads} -gff -xsmall -a \
            -dir results/repeats/{wildcards.sample}/ \
            results/repeats/{wildcards.sample}/{wildcards.sample}.fasta
        """


# -------------------------------------------------------------------------
# Gene-to-TE proximity analysis (core vs shell vs cloud)
# -------------------------------------------------------------------------

rule te_gene_proximity:
    """For each gene in each ingroup genome, distance to nearest TE;
    summarise per pangenome compartment."""
    input:
        gff_files=expand("results/annotation/{s}_liftoff.gff3", s=INGROUP_SAMPLES),
        repeatmasker=expand("results/repeats/{s}/{s}.fasta.out", s=INGROUP_SAMPLES),
        partitions="results/pangenome/partitioned_orthogroups.tsv",
        orthogroups_tsv="results/orthofinder/output",
    output:
        per_gene="results/repeats/te_gene_proximity_per_gene.tsv",
        summary="results/repeats/te_gene_proximity_summary.tsv",
    params:
        samples=INGROUP_SAMPLES,
        feature="gene",
    conda: "../envs/phylo.yaml"
    script: "../scripts/te_gene_proximity.py"
