# orthology.smk - Ortholog inference and pangenome partitioning
# Owner: Person B
# Tools: OrthoFinder
# Outputs: results/orthofinder/output/, results/pangenome/
 
rule setup_orthofinder_input:
    """Copy protein FASTAs into OrthoFinder input directory."""
    input: expand("results/proteins/{s}.fa", s=ALL_SAMPLES)
    output: directory("results/orthofinder/input")
    run:
        os.makedirs(output[0], exist_ok=True)
        for f in input:
            sample = os.path.basename(f)
            shutil.copy(f, os.path.join(output[0], sample))
 
rule orthofinder:
    input: "results/orthofinder/input"
    output: directory("results/orthofinder/output")
    params:
        method=config["orthofinder"]["method"],
        search=config["orthofinder"]["search"]
    threads: config["threads"]
    conda: "../envs/orthofinder.yaml"
    shell:
        """
        orthofinder -f {input} -t {threads} -a 4 \
            -M {params.method} -S {params.search} -o {output}
        """
 
rule partition_pangenome:
    input: "results/orthofinder/output"
    output:
        table="results/pangenome/partitioned_orthogroups.tsv",
        summary="results/pangenome/pangenome_summary.tsv"
    params:
        ingroup=INGROUP_SAMPLES,
        outgroup=OUTGROUP_SAMPLES
    script: "../scripts/partition_pangenome.py"

rule cloud_composition:
    """What cloud orthogroups are made of: outgroup genes, and whether the
    same gene (Liftoff keeps reference IDs) carries a model in another
    ingroup genome under a different orthogroup."""
    input:
        table="results/pangenome/partitioned_orthogroups.tsv",
        of="results/orthofinder/output",
        gffs=expand("results/annotation/{s}_liftoff.gff3", s=INGROUP_SAMPLES),
        proteins=expand("results/proteins/{s}.fa", s=INGROUP_SAMPLES)
    output:
        per_og="results/pangenome/cloud_composition.tsv",
        summary="results/pangenome/cloud_composition_summary.tsv"
    params:
        ingroup=INGROUP_SAMPLES,
        reference=REF
    conda: "../envs/phylo.yaml"
    script: "../scripts/cloud_composition.py"
