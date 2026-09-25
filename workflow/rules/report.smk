# report.smk - one table of every number the manuscript cites.

MANUSCRIPT_TOOLS = ["liftoff", "gffread", "orthofinder", "diamond", "mafft",
                    "trimal", "iqtree", "cafe", "skani", "minimap2", "miniprot",
                    "busco", "quast", "repeatmodeler", "repeatmasker",
                    "eggnog-mapper", "biopython"]

rule manuscript_values:
    """Collect every value the manuscript cites, with its source file, so the
    text can be checked against one table after each run."""
    input:
        "results/pangenome/pangenome_summary.tsv",
        "results/pangenome/cloud_composition_summary.tsv",
        "results/validation/absence_summary.tsv",
        "results/phylo/concord.cf.tree",
        "results/phylo/quartet_asymmetry_by_support.tsv",
        "results/phylo/quartet_robustness.tsv",
        "results/phylo/branch_length_summary.tsv",
        "results/phylo/rooted/rooted_summary.tsv",
        "results/phylo/dstat/d_statistics.tsv",
        "results/cafe/transfer_bias_summary.tsv",
        "results/synteny/synteny_summary.tsv",
        "results/synteny/ani_matrix.tsv",
        "results/functional/key_families_wide.tsv",
        "results/annotation/transfer_quality.tsv",
        "results/quast",
        expand("results/busco/{s}", s=ALL_SAMPLES),
        expand("results/busco_proteins/{s}", s=ALL_SAMPLES),
        expand("results/repeats/{s}/{s}.fasta.tbl", s=INGROUP_SAMPLES)
    output: "results/manuscript_values.tsv"
    params:
        samples=ALL_SAMPLES,
        ingroup=INGROUP_SAMPLES,
        reference=REF,
        tools=MANUSCRIPT_TOOLS,
        parameters={
            "liftoff -s (min child feature identity)": config["liftoff"]["identity_threshold"],
            "liftoff -sc (copy identity; inert without -copies)": config["liftoff"]["copy_identity"],
            "liftoff -copies / -polish": "not used",
            "orthofinder": f"-M {config['orthofinder']['method']} -S {config['orthofinder']['search']} -a 4",
            "iqtree species tree": f"-p per locus partitions -m {config['iqtree']['model']} -bb {config['iqtree']['bootstrap']}",
            "iqtree gene trees": "-m MFP -bb 1000, one per trimmed SCO",
            "concordance": "--gcf all_gene_trees.nwk --scf 100",
            "cafe": f"-p -k {config['cafe']['n_gamma_categories']}; ingroup counts; species tree rooted with {OUTGROUP}; branch lengths assumed (tips 1, root branches 0.5 when the root falls between the two pairs)",
            "outgroup": f"{OUTGROUP} {config['outgroup']['accession']} ({config['outgroup']['assembly']}), own Ensembl gene set",
            "rooted tree": f"five taxon single copy loci; MAFFT --auto, trimAl -automated1, alignments of at least 50 columns; iqtree -p per locus partitions -m {config['iqtree']['model']} -bb {config['iqtree']['bootstrap']} -o {OUTGROUP}; --gcf and --scf 100",
            "D statistics": f"codon alignments (trimAl columns of the protein alignment); sites with A, C, G or T in all five taxa; weighted block jackknife over {int((config.get('dstat') or {}).get('block_size', 5000000)) // 1000000} Mb windows of the reference assembly",
            "busco lineage": config["busco"]["lineage"],
            "minimap2 synteny": "-x asm20 --secondary=no on the three longest sequences",
            "synteny inversions": "runs of >= 3 anchors spanning >= 100 kb after orientation normalisation",
            "absence validation": "miniprot --outn 1; minimap2 -x asm20 -N 1; hit id >= 0.80, protein cov >= 0.70, DNA cov >= 0.50; same gene protein identity >= 0.90",
        }
    conda: "../envs/phylo.yaml"
    script: "../scripts/collect_manuscript_values.py"
