rule eggnog_download:
    """Download eggNOG database (one-time)."""
    output: directory("resources/eggnog_data/eggnog.db")
    conda: "../envs/eggnog.yaml"
    shell:
        """
        mkdir -p resources/eggnog_data
        download_eggnog_data.py --data_dir resources/eggnog_data/ -y
        """

rule eggnog_mapper:
    """Functional annotation with eggNOG-mapper."""
    input:
        proteins="results/proteins/{sample}.fa",
        db="resources/eggnog_data/eggnog.db"
    output: "results/eggnog/{sample}.emapper.annotations"
    threads: 8
    conda: "../envs/eggnog.yaml"
    shell:
        """
        mkdir -p results/eggnog
        emapper.py -i {input.proteins} \
            --output results/eggnog/{wildcards.sample} \
            --data_dir resources/eggnog_data/ \
            --cpu {threads} --override
        """

rule go_enrichment:
    """GO enrichment analysis: core vs shell vs cloud."""
    input:
        annotations=expand("results/eggnog/{s}.emapper.annotations",
            s=INGROUP_SAMPLES),
        partitions="results/pangenome/partitioned_orthogroups.tsv"
    output:
        enrichment="results/functional/go_enrichment_results.tsv"
    script: "../scripts/go_enrichment.R"


# -------------------------------------------------------------------------
# Key gene families (P450, GST, CCE, OR, GR, IR, OBP, CSP, immune)
# -------------------------------------------------------------------------

rule tabulate_key_families:
    """Tabulate copy numbers of key gene families across the four ingroup
    forms using eggNOG annotations layered on OrthoFinder orthogroups."""
    input:
        partitions="results/pangenome/partitioned_orthogroups.tsv",
        annotations=expand(
            "results/eggnog/{s}.emapper.annotations", s=INGROUP_SAMPLES
        ),
    output:
        long="results/functional/key_families_long.tsv",
        wide="results/functional/key_families_wide.tsv",
    params:
        families={
            "P450":   [r"cytochrome\s*p[-_]?450", r"\bcyp[0-9]"],
            "GST":    [r"glutathione\s*s[-_]?transferase", r"\bgst[a-z]?[0-9]?"],
            "CCE":    [r"carboxylesterase", r"\bcce[0-9]?"],
            "OR":     [r"\bodorant\s*receptor", r"\bor[0-9]+\b"],
            "GR":     [r"\bgustatory\s*receptor", r"\bgr[0-9]+\b"],
            "IR":     [r"\bionotropic\s*receptor", r"\bir[0-9]+\b"],
            "OBP":    [r"\bodorant[-_]?binding\s*protein"],
            "CSP":    [r"\bchemosensory\s*protein"],
            "immune": [r"toll[-_]?like", r"\bimd\b", r"defensin",
                       r"cecropin", r"attacin", r"thioester", r"\bpgrp\b"],
        },
    script: "../scripts/tabulate_key_families.py"
