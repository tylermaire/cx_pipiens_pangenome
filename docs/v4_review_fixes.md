# V4 review fixes: what changed and how to rerun

Branch `v4-review-fixes`, based on tag `v2.1` (commit `b60e019`, the V4 run).

## What changed

| Commit | Change | Outputs affected |
| --- | --- | --- |
| Count hits on the query's own gene as clustered elsewhere | `classify_absence.py` checks whether the target genome holds a model of the query's gene (Liftoff keeps reference IDs) before trusting protein identity | `results/validation/absence_calls.tsv`, `absence_summary.tsv` |
| Add cloud composition analysis | new rule `cloud_composition` | `results/pangenome/cloud_composition*.tsv` |
| Remove the GO enrichment rule | the rule wrote a placeholder, never an enrichment | `results/functional/go_enrichment_results.tsv` removed |
| Name the excess discordant topology | `quartet_asymmetry.py` is a rule; adds support thresholds and site patterns | `results/phylo/quartet_*.tsv` |
| Add divergence and transfer quality diagnostics | new rules `divergence_diagnostics`, `transfer_quality` | `results/phylo/branch_length_summary.tsv`, `sco_pairwise_identity.tsv`, `locus_accounting.tsv`, `results/annotation/transfer_quality.tsv` |
| Label inversions by chromosome and test recurrence | `synteny_summary.py` adds homologous chromosome, relative position, anchor genes and recurrence clusters | `results/synteny/inversions_detail.tsv`, new `inversion_recurrence.tsv` |
| Collect every manuscript value into one table | new rule `manuscript_values` | `results/manuscript_values.tsv` |
| Separate zero copy families in the transfer bias control | zero reference copy bin, strict single vs multi copy Fisher test with risk ratio, two sided tests | `results/cafe/transfer_bias_*.tsv` |
| Name CAFE internal nodes | internal nodes named by the taxa below them | `results/cafe/branch_summary.tsv` |
| Name the Liftoff thresholds | config key renamed; values unchanged | none |

Nothing upstream changed: Liftoff, protein extraction, OrthoFinder, alignments,
IQ-TREE and CAFE5 itself do not need to run again.

## Rerun

On the machine that holds the V4 results (the unpacked `results_v4_full.tar.gz`
and its `.snakemake/` directory), from the repository root:

```bash
git fetch origin && git checkout v4-review-fixes

# 1. Dry run. The list should hold only the rules below plus the new ones.
snakemake --use-conda -c 32 -n --rerun-triggers mtime \
    --forcerun classify_absence parse_cafe_results cafe_transfer_bias synteny_summary

# 2. Real run.
snakemake --use-conda -c 32 --rerun-triggers mtime \
    --forcerun classify_absence parse_cafe_results cafe_transfer_bias synteny_summary
```

`--rerun-triggers mtime` stops Snakemake from rebuilding upstream steps because
of unrelated metadata; the four rules whose code changed are forced explicitly,
and the new rules run because their outputs do not exist yet. If the dry run
lists `orthofinder`, `liftoff`, `concat_and_tree` or `cafe5`, stop and check
that the V4 results were restored with their timestamps.

## Checks after the run

1. `results/validation/absence_summary.tsv`: the supported share of cloud
   absences should fall sharply (on the Sep 17 run, transcript IDs alone moved
   902 of 1,159 cloud paralog_only calls to artifacts).
2. `results/pangenome/cloud_composition_summary.tsv`: share of cloud
   orthogroups with a *Cx. tarsalis* gene, and with the same gene elsewhere.
3. `results/phylo/quartet_topology_counts.tsv`: the excess topology and its
   gDF label; `quartet_site_patterns.tsv` for the site level direction.
4. `results/manuscript_values.tsv`: every number for the text.
5. Anvi'o input: confirm the anvi'o protein FASTA holds the same IDs as
   `results/proteins/*.fa`, for example
   `grep -h '>' anvio_data/*.fa* | sort | md5sum` against
   `grep -h '>' results/proteins/Cx_{molestus,pallens,pipiens,quinquefasciatus}.fa | sort | md5sum`.

Then commit the regenerated small results, merge into `main`, and tag the
release the manuscript will cite:

```bash
git add results/ && git commit -m "Regenerate results after the V4 review fixes"
git tag -a v2.2 -m "Pipeline and results cited in the revised manuscript"
git push origin v4-review-fixes --tags
```
