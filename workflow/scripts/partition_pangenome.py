#!/usr/bin/env python3
"""
Partition OrthoFinder orthogroups into core, shell, and cloud pangenome compartments.
 
Input:  Orthogroups.GeneCount.tsv from OrthoFinder
Output: partitioned_orthogroups.tsv, pangenome_summary.tsv
"""
import pandas as pd
 
# TODO: Implement
# 1. Read snakemake.input[0] (gene count matrix)
# 2. For each orthogroup, count how many ingroup species have >= 1 gene
# 3. core = all 4, shell = 2-3, cloud = 1
# 4. Write partitioned table to snakemake.output.table
# 5. Write summary counts to snakemake.output.summary
