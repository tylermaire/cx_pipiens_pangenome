rule nucmer_align:
    """Pairwise whole-genome alignment with nucmer."""
    input:
        ref="resources/genomes/{ref}.fasta",
        query="resources/genomes/{query}.fasta"
    output:
        delta="results/synteny/{ref}_vs_{query}.delta",
        filtered="results/synteny/{ref}_vs_{query}.filtered.delta"
    threads: 8
    conda: "../envs/mummer.yaml"
    shell:
        """
        nucmer -t {threads} -l 100 -c 500 \
            {input.ref} {input.query} \
            -p results/synteny/{wildcards.ref}_vs_{wildcards.query}
        delta-filter -m -i 90 -l 100 \
            {output.delta} > {output.filtered}
        """

rule syri:
    """Detect structural variants with SyRI from nucmer delta."""
    input:
        delta="results/synteny/{ref}_vs_{query}.filtered.delta",
        ref="resources/genomes/{ref}.fasta",
        query="resources/genomes/{query}.fasta"
    output: "results/synteny/{ref}_vs_{query}_syri.out"
    conda: "../envs/syri.yaml"
    shell:
        """
        syri -c {input.delta} -r {input.ref} -q {input.query} \
            -k --nc 3 --nosnp --no-chrmatch \
            --dir results/synteny/ \
            --prefix {wildcards.ref}_vs_{wildcards.query}_ || \
            touch {output}
        """
