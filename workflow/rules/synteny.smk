rule filter_chromosomes:
    """Reduce each genome FASTA to its 3 longest contigs and rename them
    chr1/chr2/chr3 by length rank, so SyRI sees matched chromosome names
    across pairs."""
    input:  "resources/genomes/{sample}.fasta"
    output: "results/synteny/genomes/{sample}.chromosomes.fasta"
    shell:
        "python workflow/scripts/filter_chromosomes.py {input} {output} 3"

rule nucmer_align:
    """Pairwise whole-genome alignment with nucmer (chromosomes only)."""
    input:
        ref="results/synteny/genomes/{ref}.chromosomes.fasta",
        query="results/synteny/genomes/{query}.chromosomes.fasta"
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
        delta-filter -m -i 90 -l 100 {output.delta} > {output.filtered}
        """

rule show_coords:
    """SyRI-compatible coordinate dump."""
    input:
        filtered="results/synteny/{ref}_vs_{query}.filtered.delta"
    output:
        coords="results/synteny/{ref}_vs_{query}.coords",
        diff="results/synteny/{ref}_vs_{query}.diff"
    conda: "../envs/mummer.yaml"
    shell:
        """
        show-coords -THrd {input.filtered} > {output.coords}
        show-diff {input.filtered} > {output.diff}
        """

rule syri:
    """Detect structural rearrangements with SyRI."""
    input:
        ref="results/synteny/genomes/{ref}.chromosomes.fasta",
        query="results/synteny/genomes/{query}.chromosomes.fasta",
        filtered="results/synteny/{ref}_vs_{query}.filtered.delta",
        coords="results/synteny/{ref}_vs_{query}.coords"
    output:
        syri="results/synteny/{ref}_vs_{query}_syri.out",
        summary="results/synteny/{ref}_vs_{query}_syri.summary"
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
        )
    output:
        table="results/synteny/3Rb_inversion_summary.tsv",
        flag="results/synteny/3Rb_inversion_call.txt"
    params:
        # In the chromosome-renamed files, chr1 = NC_051862.1 (longest contig,
        # which is chromosome 3 in NCBI naming for this assembly).
        chr3_label="chr1",
        rb3_start=108_955_000,
        rb3_end=117_355_000,
        min_inv_len=100_000,
        shared_thresh=0.5
    script: "../scripts/check_3Rb_inversion.py"

# -------------------------------------------------------------------------
# Gene-anchor synteny summary
# -------------------------------------------------------------------------
# The collinearity and inversion figures in the Results were previously
# produced outside the workflow and could not be regenerated from the repo.
# These rules replace that, and normalise chromosome orientation before
# counting: several chromosome pairs are deposited on opposite strands
# between assemblies, which otherwise makes collinear blocks read as
# whole-chromosome inversions.

rule synteny_minimap2:
    """Whole-chromosome alignment for a real (not partial) identity estimate."""
    input:
        ref="results/synteny/genomes/{ref}.chromosomes.fasta",
        query="results/synteny/genomes/{query}.chromosomes.fasta"
    output: "results/synteny/paf/{ref}_vs_{query}.paf"
    threads: config["threads"]
    conda: "../envs/minimap2.yaml"
    shell:
        "mkdir -p results/synteny/paf && "
        "minimap2 -t {threads} -x asm20 --secondary=no "
        "{input.ref} {input.query} > {output}"

rule synteny_summary:
    """Anchor-based collinearity and inversion counts across ingroup pairs."""
    input:
        gffs=expand("results/annotation/{s}_liftoff.gff3", s=INGROUP_SAMPLES),
        pafs=[f"results/synteny/paf/{a}_vs_{b}.paf"
              for i, a in enumerate(INGROUP_SAMPLES)
              for b in INGROUP_SAMPLES[i + 1:]],
        ani="results/synteny/ani_matrix.tsv"
    output:
        summary="results/synteny/synteny_summary.tsv",
        orientation="results/synteny/chromosome_orientation.tsv",
        inversions="results/synteny/inversions_detail.tsv"
    params:
        samples=INGROUP_SAMPLES,
        n_chrom=3,
        min_span=100_000,
        min_genes=3
    conda: "../envs/phylo.yaml"
    script: "../scripts/synteny_summary.py"
