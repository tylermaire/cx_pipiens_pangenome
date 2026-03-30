# go_enrichment.R
# GO term enrichment analysis: core vs shell vs cloud gene sets
# Uses topGO or clusterProfiler
#
# Input:  snakemake@input[["annotations"]], snakemake@input[["partitions"]]
# Output: snakemake@output[["enrichment"]], snakemake@output[["heatmap"]]
 
library(tidyverse)
# library(topGO)       # uncomment when implementing
# library(clusterProfiler)
 
# TODO: Implement
