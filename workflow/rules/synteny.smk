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


# -------------------------------------------------------------------------
# Structural variation with SyRI (required by `rule all` -> *_syri.out)
# -------------------------------------------------------------------------

rule syri:
    """Detect structural rearrangements with SyRI from nucmer alignments."""
    input:
        ref="resources/genomes/{ref}.fasta",
        query="resources/genomes/{query}.fasta",
        filtered="results/synteny/{ref}_vs_{query}.filtered.delta",
        coords="results/synteny/{ref}_vs_{query}.coords",
    output:
        syri="results/synteny/{ref}_vs_{query}_syri.out",
        summary="results/synteny/{ref}_vs_{query}_syri.summary",
    threads: 4
    conda: "../envs/syri.yaml"
    shell:
        r"""
        WORKDIR=results/synteny/syri_{wildcards.ref}_vs_{wildcards.query}
        mkdir -p $WORKDIR
        syri -c {input.coords} -d {input.filtered} \
             -r {input.ref} -q {input.query} \
             --prefix {wildcards.ref}_vs_{wildcards.query}_ \
             --dir $WORKDIR --nc {threads} --all
        cp $WORKDIR/{wildcards.ref}_vs_{wildcards.query}_syri.out     {output.syri}
        cp $WORKDIR/{wildcards.ref}_vs_{wildcards.query}_syri.summary {output.summary}
        """


rule check_3Rb_inversion:
    """Test whether the 3Rb paracentric inversion (Ryazansky et al. 2024)
    is shared with the other Cx. pipiens forms."""
    input:
        syri_outs=expand(
            "results/synteny/{ref}_vs_{query}_syri.out",
            ref=[REF],
            query=[s for s in INGROUP_SAMPLES if s != REF],
        ),
    output:
        table="results/synteny/3Rb_inversion_summary.tsv",
        flag="results/synteny/3Rb_inversion_call.txt",
    params:
        # TODO: replace these placeholders with the real chr3 contig name and
        # 3Rb breakpoint coordinates from Ryazansky et al. 2024 (GCF_015732765.1).
        chr3_label="NC_051862.1",
        rb3_start=20_000_000,
        rb3_end=55_000_000,
        min_inv_len=100_000,
        shared_thresh=0.5,
    script: "../scripts/check_3Rb_inversion.py"
