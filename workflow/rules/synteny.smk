rule filter_chromosomes:
    """Reduce each genome FASTA to its 3 longest contigs and rename them
    chr1/chr2/chr3 by length rank, so SyRI sees matched chromosome names
    across pairs."""
    input:  "resources/genomes/{sample}.fasta"
    output: "results/synteny/genomes/{sample}.chromosomes.fasta"
    shell:
        "python workflow/scripts/filter_chromosomes.py {input} {output} 3"

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
