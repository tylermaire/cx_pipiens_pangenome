rule extract_sco_sequences:
    """Extract single-copy ortholog sequences for ingroup species."""
    input: "results/orthofinder/output"
    output: directory("results/phylo/sco_fastas")
    conda: "../envs/phylo.yaml"
    script: "../scripts/extract_sco.py"

rule concat_and_tree:
    """Build species tree from single-copy orthologs."""
    input: "results/phylo/sco_fastas"
    output:
        tree="results/phylo/concat_tree.treefile",
        concord="results/phylo/concord.cf.tree"
    params:
        model=config["iqtree"]["model"],
        bb=config["iqtree"]["bootstrap"]
    threads: config["threads"]
    conda: "../envs/phylo.yaml"
    shell:
        """
        mkdir -p results/phylo/alignments results/phylo/trimmed results/phylo/gene_trees

        for fa in results/phylo/sco_fastas/*.fa; do
            og=$(basename $fa .fa)
            if [ ! -f results/phylo/trimmed/$og.trim ]; then
                mafft --auto $fa > results/phylo/alignments/$og.aln 2>/dev/null
                trimal -in results/phylo/alignments/$og.aln \
                    -out results/phylo/trimmed/$og.trim -automated1 2>/dev/null || true
            fi
        done

        find results/phylo/trimmed -name "*.trim" -size -100c -delete
        NTRIM=$(ls results/phylo/trimmed/*.trim 2>/dev/null | wc -l)
        echo "Trimmed alignments: $NTRIM"

        for trim in results/phylo/trimmed/*.trim; do
            og=$(basename $trim .trim)
            if [ ! -f results/phylo/gene_trees/$og.treefile ]; then
                iqtree -s $trim -m MFP -bb 1000 -nt 1 \
                    --prefix results/phylo/gene_trees/$og -quiet 2>/dev/null || true
            fi
        done

        cat results/phylo/gene_trees/*.treefile > results/phylo/all_gene_trees.nwk 2>/dev/null || true

        iqtree -p results/phylo/trimmed/ \
            -m {params.model} -bb {params.bb} -nt {threads} \
            --prefix results/phylo/concat_tree -quiet

        iqtree -t results/phylo/concat_tree.treefile \
            --gcf results/phylo/all_gene_trees.nwk \
            -p results/phylo/trimmed/ \
            --scf 100 --prefix results/phylo/concord -quiet -nt {threads} || \
            cp results/phylo/concat_tree.treefile {output.concord}
        """
