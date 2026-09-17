"""
Validation of pangenome absence calls (cloud and shell).

Every gene model in this study was transferred from the Cx. quinquefasciatus
reference by Liftoff, so an orthogroup recorded as absent from a genome may
simply be a gene the transfer missed. This module tests each absence with two
independent probes and the target genome's own annotation.

  miniprot   representative protein vs each ingroup genome (divergence tolerant,
             but can land on a paralog)
  minimap2   that gene's genomic locus vs each ingroup genome (positionally
             stricter, degrades faster with distance)
  context    whether the hit lands on annotated ground, and if so whether that
             gene was assigned to a different orthogroup

Only absences that neither probe can find are treated as real.
"""

# Matches how partition_pangenome.py splits ingroup from outgroup.
INGROUP = [s for s in ALL_SAMPLES if "tarsalis" not in s.lower()]

rule extract_absence_queries:
    """One representative protein and genomic locus per cloud/shell orthogroup."""
    input:
        table="results/pangenome/partitioned_orthogroups.tsv",
        of="results/orthofinder/output",
        proteins=expand("results/proteins/{s}.fa", s=INGROUP),
        gffs=expand("results/annotation/{s}_liftoff.gff3", s=INGROUP)
    output:
        dir=directory("results/validation/queries"),
        manifest="results/validation/queries/manifest.tsv",
        faa="results/validation/queries/queries.faa",
        fna="results/validation/queries/queries.fna"
    params:
        ingroup=INGROUP,
        flank=0
    conda: "../envs/phylo.yaml"
    script: "../scripts/extract_absence_queries.py"

rule miniprot_index:
    input: "resources/genomes/{target}.fasta"
    output: "results/validation/index/{target}.mpi"
    threads: config["threads"]
    conda: "../envs/miniprot.yaml"
    shell: "miniprot -t {threads} -d {output} {input}"

rule absence_protein_align:
    """Protein probe: miniprot against one target genome."""
    input:
        faa="results/validation/queries/queries.faa",
        index="results/validation/index/{target}.mpi"
    output: "results/validation/paf/prot_vs_{target}.paf"
    threads: config["threads"]
    conda: "../envs/miniprot.yaml"
    shell:
        "mkdir -p results/validation/paf && "
        "miniprot -t {threads} --outn 1 {input.index} {input.faa} > {output}"

rule absence_dna_align:
    """DNA probe: minimap2 of the genomic locus against one target genome."""
    input:
        fna="results/validation/queries/queries.fna",
        genome="resources/genomes/{target}.fasta"
    output: "results/validation/paf/dna_vs_{target}.paf"
    threads: config["threads"]
    conda: "../envs/miniprot.yaml"
    shell:
        "mkdir -p results/validation/paf && "
        "minimap2 -t {threads} -x asm20 -N 1 --secondary=no "
        "{input.genome} {input.fna} > {output}"

rule classify_absence:
    """Call every absence event using both probes plus annotation context."""
    input:
        manifest="results/validation/queries/manifest.tsv",
        of="results/orthofinder/output",
        prot_pafs=expand("results/validation/paf/prot_vs_{t}.paf", t=INGROUP),
        dna_pafs=expand("results/validation/paf/dna_vs_{t}.paf", t=INGROUP),
        gffs=expand("results/annotation/{s}_liftoff.gff3", s=INGROUP)
    output:
        calls="results/validation/absence_calls.tsv",
        summary="results/validation/absence_summary.tsv"
    params:
        ingroup=INGROUP,
        min_id=0.80,
        min_cov=0.70,
        dna_min_cov=0.50,
        same_gene_id=0.90,
        max_overlap=10
    conda: "../envs/phylo.yaml"
    script: "../scripts/classify_absence.py"
