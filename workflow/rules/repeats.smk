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
