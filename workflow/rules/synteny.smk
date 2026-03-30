# synteny.smk - Structural variation and synteny analysis
# Owner: Person D
# Tools: minimap2, SyRI, plotsr, MUMmer4
# Outputs: results/synteny/, results/figures/
 
# TODO: Implement rules for:
# 1. minimap2_align - pairwise whole-genome alignments (6 pairs)
# 2. syri - detect inversions, translocations, duplications
# 3. plotsr - publication-quality synteny figures
# 4. nucmer_dotplot - dotplot visualizations
# 5. check_3Rb_inversion - specific chr3 inversion analysis
