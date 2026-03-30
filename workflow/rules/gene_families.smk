# gene_families.smk - Gene family evolution and selection
# Owner: Person C
# Tools: CAFE5, PAML
# Outputs: results/cafe/
 
# TODO: Implement rules for:
# 1. prepare_cafe_input - filter gene counts, make ultrametric tree
# 2. cafe5 - run CAFE5 birth-death model
# 3. parse_cafe_results - extract significant families
# 4. paml_branch_site (optional) - positive selection on candidates
