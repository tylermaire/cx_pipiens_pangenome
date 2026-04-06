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

rule show_coords:
    """Extract alignment coordinates and structural diffs."""
    input:
        filtered="results/synteny/{ref}_vs_{query}.filtered.delta"
    output:
        coords="results/synteny/{ref}_vs_{query}.coords",
        diff="results/synteny/{ref}_vs_{query}.diff"
    conda: "../envs/mummer.yaml"
    shell:
        """
        show-coords -rcl {input.filtered} > {output.coords}
        show-diff {input.filtered} > {output.diff}
        """
