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
    script: "../scripts/partition_pangenome.py"
