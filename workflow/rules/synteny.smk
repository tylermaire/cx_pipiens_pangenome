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

rule skani_ani:
    """Pairwise whole-genome ANI. Previously produced by hand; two versions
    existed in the repository disagreeing by 0.2-0.3% per cell with no record
    of the parameters used, so this rule replaces both."""
    input:
        genomes=expand("resources/genomes/{s}.fasta", s=ALL_SAMPLES)
    output:
        matrix="results/synteny/ani_matrix.tsv",
        long="results/synteny/ani_pairs.tsv"
    params:
        samples=ALL_SAMPLES,
        extra=["-s", "80"]     # report down to 80% ANI; see script header
    conda: "../envs/skani.yaml"
    script: "../scripts/skani_matrix.py"

rule synteny_minimap2:
    """Whole-chromosome alignment for a real (not partial) identity estimate."""
    input:
        ref="results/synteny/genomes/{ref}.chromosomes.fasta",
        query="results/synteny/genomes/{query}.chromosomes.fasta"
    output: "results/synteny/paf/{ref}_vs_{query}.paf"
    threads: config["threads"]
    # asm20 on 550 Mb chromosome sets peaks near 150 GB, so these must not
    # run concurrently: Snakemake schedules by cores alone unless told
    # otherwise, and three at once exceeds any single machine here.
    resources:
        mem_mb=160000
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
        inversions="results/synteny/inversions_detail.tsv",
        recurrence="results/synteny/inversion_recurrence.tsv"
    params:
        samples=INGROUP_SAMPLES,
        reference=REF,
        n_chrom=3,
        min_span=100_000,
        min_genes=3,
        large_span=1_000_000
    conda: "../envs/phylo.yaml"
    script: "../scripts/synteny_summary.py"
