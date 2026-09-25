REF = config["reference"]["name"]

rule liftoff:
    """Transfer the reference annotation to an ingroup assembly."""
    wildcard_constraints:
        sample=one_of(LIFTOFF_SAMPLES)
    input:
        target="resources/genomes/{sample}.fasta",
        ref=f"resources/genomes/{REF}.fasta",
        gff=f"resources/annotations/{REF}.gff3"
    output: "results/annotation/{sample}_liftoff.gff3"
    params:
        sc=config["liftoff"]["copy_identity"],
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


rule native_annotation:
    """A sample's own gene set (the outgroup's Ensembl annotation), with the
    type prefixes Ensembl puts on identifiers removed. results/annotation/
    <sample>_liftoff.gff3 is the annotation used for every sample, whatever
    its source: the RefSeq annotation for the reference, Liftoff transfers
    for the other ingroup forms and the native gene set for the outgroup."""
    wildcard_constraints:
        sample=one_of(NATIVE_SAMPLES)
    input: "resources/annotations/{sample}.gff3"
    output: "results/annotation/{sample}_liftoff.gff3"
    conda: "../envs/phylo.yaml"
    shell: "python workflow/scripts/normalize_gff.py {input} {output}"

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

rule extract_cds:
    """Spliced coding sequence of each protein kept by extract_proteins, under
    the same name, for the codon alignments of the rooted analyses."""
    input:
        fasta="resources/genomes/{sample}.fasta",
        gff="results/annotation/{sample}_liftoff.gff3",
        proteins="results/proteins/{sample}.fa"
    output: "results/cds/{sample}.fa"
    conda: "../envs/gffread.yaml"
    shell:
        """
        mkdir -p results/cds
        gffread {input.gff} -g {input.fasta} -x {output}.tmp
        python workflow/scripts/filter_cds.py --cds {output}.tmp \
            --proteins {input.proteins} --out {output} --strip-periods
        rm {output}.tmp
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
    """Broken transferred models among those kept for analysis, read from
    Liftoff's own ORF checks in the GFF (valid_ORF, missing start or stop,
    in frame stop, matches_ref_protein), overall and by pangenome
    compartment. Liftoff runs without -polish and extract_proteins strips
    stop symbols, so these models are otherwise invisible downstream.
    RefSeq models that are partial or carry an exception (a corrected indel
    in the reference assembly) are counted apart, because most of their
    transfers fail the check whatever the target genome holds. The outgroup
    keeps its own gene set and has no Liftoff flags."""
    input:
        gffs=expand("results/annotation/{s}_liftoff.gff3", s=ALL_SAMPLES),
        proteins=expand("results/proteins/{s}.fa", s=ALL_SAMPLES),
        table="results/pangenome/partitioned_orthogroups.tsv",
        of="results/orthofinder/output"
    output:
        summary="results/annotation/transfer_quality.tsv",
        by_compartment="results/annotation/transfer_quality_by_compartment.tsv"
    params:
        samples=ALL_SAMPLES,
        reference=REF,
        native=NATIVE_SAMPLES
    conda: "../envs/phylo.yaml"
    script: "../scripts/transfer_quality.py"
