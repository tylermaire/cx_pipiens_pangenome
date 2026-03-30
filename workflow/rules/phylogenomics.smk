# phylogenomics.smk - Species tree inference
# Owner: Person C
# Tools: MAFFT, trimAl, IQ-TREE, ASTRAL
# Outputs: results/phylo/
 
# TODO: Implement rules for:
# 1. extract_sco_sequences - pull single-copy orthologs from OrthoFinder
# 2. align_sco - MAFFT --auto per orthogroup
# 3. trim_alignment - trimAl -automated1
# 4. concatenate_alignments - build supermatrix
# 5. iqtree_concat - ML tree with MFP + UFBoot
# 6. gene_trees - per-locus IQ-TREE runs
# 7. astral - coalescent species tree
# 8. concordance_factors - IQ-TREE --gcf --scf
