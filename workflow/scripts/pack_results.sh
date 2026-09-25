#!/usr/bin/env bash
# pack_results.sh - the outputs needed to check a run and to patch its summary
# tables without the cloud (patch_results.py), in one archive of about 1 GB.
#
#     bash workflow/scripts/pack_results.sh [results_v5_small.tar.gz]
#
# Run from the repository root after the workflow. Paths that a run did not
# produce are listed and skipped. The full results/ directory (tens of GB) is
# better archived separately, as for V4.
set -euo pipefail
OUT="${1:-results_v5_small.tar.gz}"
LIST=$(mktemp)
trap 'rm -f "$LIST"' EXIT

want() {
    for p in "$@"; do
        if compgen -G "$p" > /dev/null; then
            for m in $p; do printf '%s\n' "$m" >> "$LIST"; done
        else
            echo "not found, skipped: $p" >&2
        fi
    done
}

want run_v5.log config/config.yaml config/samples.tsv
want "resources/annotations/*.source.txt"
want results/manuscript_values.tsv
want "results/annotation/*_liftoff.gff3" "results/annotation/transfer_quality*.tsv"
want "results/proteins/*.fa" "results/proteins/*_isoform_report.tsv"
want "results/busco/*/short_summary*" "results/busco_proteins/*/short_summary*"
want results/quast/report.tsv results/quast/transposed_report.tsv
want "results/pangenome/*.tsv"
want "results/validation/*.tsv" results/validation/queries/manifest.tsv "results/validation/paf/*.paf"
want "results/orthofinder/output/*/Orthogroups/Orthogroups*.tsv"
want "results/orthofinder/output/*/Log.txt"
want "results/phylo/*.tsv" "results/phylo/*.txt" results/phylo/all_gene_trees.nwk
want "results/phylo/concat_tree.*" "results/phylo/concord.cf.*" results/phylo/trimmed
want "results/phylo/rooted/*.tsv" "results/phylo/rooted/*.txt" "results/phylo/rooted/*.nwk"
want "results/phylo/rooted/rooted_tree.*" "results/phylo/rooted/rooted_concord.cf.*"
want results/phylo/rooted/alignments/trimmed results/phylo/rooted/alignments/codon
want "results/phylo/dstat/*.tsv"
want "results/cafe/*.tsv" results/cafe/ultrametric_tree.nwk results/cafe/output
want "results/functional/*.tsv"
want "results/repeats/*.tsv" "results/repeats/*/*.tbl"
want "results/synteny/*.tsv" "results/synteny/paf/*.paf"
want figures/revision tables

# checkpoint and model files are large and not needed
tar -czf "$OUT" --exclude='*.ckp.gz' --exclude='*.model.gz' -T "$LIST"
echo "wrote $OUT ($(du -h "$OUT" | cut -f1)) from $(wc -l < "$LIST") paths"
