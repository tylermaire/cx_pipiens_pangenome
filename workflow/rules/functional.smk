rule eggnog_download:
    """Download eggNOG database (one-time). Patches the dead eggnogdb.embl.de
    URL to the live eggnog5.embl.de mirror before downloading."""
    output: "resources/eggnog_data/eggnog.db"
    conda: "../envs/eggnog.yaml"
    shell:
        "bash workflow/scripts/run_eggnog_download.sh"


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
    conda: "../envs/r.yaml"
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
        # [weight, regex]: weight 2 = specific phrase or canonical symbol,
        # weight 1 = looser token that only wins when nothing better matched.
        # The old patterns were loose enough to collide (\\bor[0-9]+\\b caught
        # "Obp12") and some matched bare substrings ("gst", "cce").
        families={
            "P450":   [[2, r"cytochrome\s*p[-_ ]?450"],
                       [2, r"\bcyp\d{1,2}[a-z]{1,2}\d+\b"],
                       [1, r"\bp450\b"]],
            "GST":    [[2, r"glutathione\s+s[-_ ]?transferase"],
                       [2, r"\bgst[demostuz]\d*\b"],
                       [1, r"\bglutathione\s+transferase\b"]],
            "CCE":    [[2, r"carboxyl(?:esterase|/choline\s*esterase)"],
                       [2, r"\bcholinesterase\b"],
                       [2, r"\bcce\d{3,}\b"],
                       [1, r"\besterase\b"]],
            "OR":     [[2, r"\bodorant\s+receptor\b"],
                       [2, r"\bolfactory\s+receptor\b"],
                       [2, r"\bor\d{1,3}\b(?!.*binding)"],
                       [1, r"\b7tm[-_ ]?odorant"]],
            "GR":     [[2, r"\bgustatory\s+receptor\b"],
                       [2, r"\bgr\d{1,3}\b"],
                       [1, r"\btaste\s+receptor\b"]],
            "IR":     [[2, r"\bionotropic\s+(?:glutamate\s+)?receptor\b"],
                       [2, r"\bir\d{1,3}[a-z]?\b"],
                       [1, r"\bglutamate[-_ ]?gated\s+ion\s+channel\b"]],
            "OBP":    [[2, r"\bodorant[-_ ]?binding\s+protein\b"],
                       [2, r"\bobp\d{1,3}\b"],
                       [2, r"\bpbp[-_/]?gobp\b"],
                       [1, r"\bpheromone[-_ ]?binding\s+protein\b"]],
            "CSP":    [[2, r"\bchemosensory\s+protein\b"],
                       [2, r"\bcsp\d{1,3}\b"],
                       [1, r"\bosd\b"]],
            "immune": [[2, r"\btoll[-_ ]?like\s+receptor\b"],
                       [2, r"\bpeptidoglycan[-_ ]?recognition\b"],
                       [2, r"\bdefensin\b"], [2, r"\bcecropin\b"],
                       [2, r"\battacin\b"], [2, r"\bdiptericin\b"],
                       [2, r"\bgambicin\b"],
                       [2, r"\bthioester[-_ ]?containing\b"],
                       [2, r"\bpgrp[-_]?[ls][a-z]?\b"],
                       [1, r"\bimd\b"]],
        },
        min_genes=2,
        min_frac=0.25,
    script: "../scripts/tabulate_key_families.py"
