# V4 review fixes: what changed and how the V4 results were patched

Branch `v4-review-fixes`, based on tag `v2.1` (commit `b60e019`, the V4 run).

The fixes only touch summary steps. Liftoff, protein extraction, OrthoFinder,
alignment, IQ-TREE, CAFE5, minimap2, miniprot, skani and RepeatMasker did not
change, so the V4 outputs were patched in place with
`workflow/scripts/patch_results.py` instead of a cloud rerun. The full
Snakemake workflow still produces the same outputs from scratch.

## What changed

| Change | Script or rule | Outputs |
| --- | --- | --- |
| Absence validation checks whether the target genome holds a model of the query's own gene (Liftoff keeps reference IDs) before it trusts protein identity; `check_paralogy.py` removed | `classify_absence.py` | `results/validation/absence_calls.tsv`, `absence_summary.tsv` |
| Cloud composition: outgroup genes, single gene orthogroups, same gene clustered elsewhere | new rule `cloud_composition` | `results/pangenome/cloud_composition*.tsv` |
| GO enrichment rule removed (it wrote a placeholder, never an enrichment) | `functional.smk` | `results/functional/go_enrichment_results.tsv` removed |
| Excess discordant topology named, tested at UFBoot thresholds, on resolved gene trees and on loci with four intact gene models, and compared with site patterns (summed, one vote per locus, without the 1% of loci richest in informative sites) | new rule `quartet_asymmetry` | `results/phylo/quartet_*.tsv`, including `quartet_robustness.tsv` |
| Branch lengths per locus, pairwise SCO identity, locus accounting | new rule `divergence_diagnostics` | `results/phylo/branch_length_summary.tsv`, `sco_pairwise_identity.tsv`, `locus_accounting.tsv` |
| Transferred model quality read from Liftoff's own ORF flags, with the RefSeq baseline and a split by compartment | rule `transfer_quality` (replaces two gffread rules) | `results/annotation/transfer_quality.tsv`, `transfer_quality_by_compartment.tsv` |
| Inversions labelled by homologous chromosome, position and anchor genes; recurrence across pairs; shared genes on a non homologous chromosome counted | `synteny_summary.py` | `results/synteny/inversions_detail.tsv`, new `inversion_recurrence.tsv`, `synteny_summary.tsv` |
| Only the families CAFE tested (it drops families absent from one side of the root); zero reference copy bin; one against two or more reference copies Fisher test with risk ratio; two sided tests | `cafe_transfer_bias.py` | `results/cafe/transfer_bias_*.tsv` |
| CAFE internal nodes named by the taxa below them | `parse_cafe.py` | `results/cafe/branch_summary.tsv` |
| Liftoff config key renamed to `copy_identity` (values unchanged; `-sc` is inert without `-copies`) | `config.yaml` | none |
| Every number the text cites, with its source file, plus tool versions read from the outputs | new rule `manuscript_values` | `results/manuscript_values.tsv` |
| Tracked results synced with the V4 run; June leftovers removed | | `results/` |
| Result tables written with LF line endings | four scripts | none (content unchanged) |
| Manuscript figures, Tables 1 to 4 and Supp. Tables S1 to S11 built from `results/` | new `make_revision_figures.py`, `make_manuscript_tables.py`, `make_supplement.py` (the last replaces `make_tables.py`) | `figures/revision/`, `tables/` (not tracked) |

Tests for every changed script are in `tests/` (`python3 tests/test_<name>.py`).

## How the V4 results were patched

From the repository root, with the V4 outputs in `results/`:

```bash
python3 workflow/scripts/patch_results.py                    # every step
python3 workflow/scripts/patch_results.py --steps transfer,values
```

Steps: `absence`, `cloud`, `transfer`, `quartet`, `divergence`, `cafe`,
`synteny`, `values`. Each step runs the same script as its Snakemake rule, with
the same inputs and parameters, through a stand in `snakemake` object. When a
rule changes, change the matching `step_` function too.

Inputs the steps read from the V4 run:

* `results/annotation/*_liftoff.gff3` and `results/proteins/*.fa` (all five genomes)
* `results/orthofinder/output/Results_*/Orthogroups/Orthogroups.tsv`
* `results/pangenome/partitioned_orthogroups.tsv`
* `results/validation/queries/manifest.tsv` and `results/validation/paf/*.paf`
* `results/phylo/concat_tree.treefile`, `concord.cf.tree`, `concord.cf.stat`,
  `all_gene_trees.nwk` and `trimmed/*.trim`; the directory listings
  `gene_tree_ids.txt` and `sco_fasta_ids.txt` stand in for `gene_trees/` and
  `sco_fastas/`
* `results/cafe/output/` and `results/cafe/gene_counts_filtered.tsv`
* `results/synteny/paf/*.paf` and `results/synteny/ani_matrix.tsv`

## Patched V4 results

* Absence validation: 93.8% of the 1,428 cloud absence events are annotation
  artifacts (1,300 of them the query's own gene clustered in another
  orthogroup), 3.4% are supported. Shell: 76.5% artifacts, 9.7% supported.
* Cloud: 464 of 476 cloud orthogroups (97.5%) hold a gene whose namesake (same
  gene ID) in another form sits in a different orthogroup; 449 contain a
  *Cx. tarsalis* gene.
* Transfer quality: Liftoff flags 30.1% (pallens), 34.3% (molestus) and 32.9%
  (pipiens) of the kept transferred models as lacking a valid ORF, 26.7 to
  31.1% when the 5.8% of RefSeq models that are partial or discrepant are set
  aside. By compartment: core 24 to 27%, shell 51 to 57%, cloud 93 to 100%,
  unassigned 98 to 99%. *Cx. tarsalis*: 75.5%.
* Quartets: the excess discordant topology is gDF2 (molestus + quinquefasciatus
  | pallens + pipiens), 1,752 vs 1,384 gene trees (55.9% of discordant,
  p = 5.3e-11), rising to 67.9% at UFBoot 95. It holds on the 6,383 gene
  trees with a resolved internal branch (58.6%, 1,346 vs 952; the 1,263
  unresolved trees split 425/406/432) and on loci with four intact models
  (57.2%, 1,066 vs 797; 81.6% at UFBoot 95, 155 vs 35). Summed parsimony
  informative sites lean the other way (21,304 vs 20,466, p = 4.2e-5), but
  the 91 loci richest in such sites hold 46.5% of them and 90 of the 91 hold
  a model that is not intact (34.2% of other loci); one vote per locus gives
  60.2% for gDF2 (1,093 vs 722), and without those 91 loci sites give 57.6%.
* Divergence: per locus median terminal branches 0.0033 to 0.0062 against 0.033
  to 0.056 in the concatenated tree; median pairwise SCO identity 98.5 to
  99.1%. 7,646 of 9,098 trimmed loci have gene trees; all 1,452 without one
  have fewer than four distinct sequences.
* CAFE transfer bias: CAFE tested 12,908 of the 13,885 families (it drops
  families absent from one side of the root). 14 of 11,012 one copy families
  (0.13%) are significant against 247 of 1,318 with two or more copies
  (18.7%; OR 181, RR 147); 72 of 578 zero reference copy families (12.5%).
  Both internal branches show only increases (141 and 45).
* Synteny: 536 inversions form 392 clusters; 95 recur in two or three pairs,
  and 93 of those share one genome across every pair they appear in (molestus
  49, pallens 23, quinquefasciatus 19, pipiens 2), so each looks like a
  rearrangement in, or an assembly error of, that one genome. Nine recurrent
  clusters span 1 Mb or more (pallens 4, quinquefasciatus 3, molestus 2).
  10.0 to 13.4% of shared genes on the chromosome scale sequences sit on a
  non homologous chromosome in one of the two genomes.

## Still open

1. Anvi'o figure: the V4 run did not produce anvi'o inputs, and the anvi'o
   figure in the draft was built from June proteins. Either rebuild it from
   `results/proteins/Cx_{molestus,pallens,pipiens,quinquefasciatus}.fa` or drop
   it; Figure 2 already shows the partition.
2. Push the branch, merge into `main`, and tag the release the manuscript
   cites:

```bash
git push origin v4-review-fixes
git checkout main && git merge --ff-only v4-review-fixes
git tag -a v2.2 -m "Pipeline and patched V4 results cited in the revised manuscript"
git push origin main --tags
```
