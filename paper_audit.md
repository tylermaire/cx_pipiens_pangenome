# Paper Audit — current draft vs new pipeline results

Audit of `Culex pipiens Complex Pangenome Analysis (4).docx` against the clean end-to-end pipeline run completed 2026-05-21.

---

## TL;DR — three things you must address

1. **The species tree topology changed.** The new tree groups *Cx. pallens* with *Cx. quinquefasciatus* as sisters, not *Cx. pallens* with *Cx. pipiens*. Multiple paragraphs in Sections 3.2, 4 (Discussion), and 5 (Conclusions) build their argument on the old topology and will need rewriting.
2. **Gene/site concordance is not 100/100/100 anymore.** The new tree shows 100/59.3/45.4 (bootstrap/gCF/sCF). The paper's claim of "complete concordance across all 8,616 gene trees at every node" is now wrong — there is substantial gene tree discordance at the deepest branch.
3. **Cx. tarsalis isn't actually in the analysis as you describe it.** The assembly is 51 Mb (versus ~533–573 Mb for the ingroup), BUSCO protein-mode is 0.0% complete, and it isn't in the species tree at all (the tree has 4 taxa, not 5). The Methods claim "BUSCO completeness ≥ 95%" and the framing of tarsalis as a tree-rooting outgroup are not consistent with what was actually done.

Everything else is small-number drift (±0–20 across most figures), narratively unchanged.

---

## 1. Species tree (Section 3.2) — MAJOR REWRITE NEEDED

### What the paper says
> "*Cx. pallens* and *Cx. pipiens* are sister taxa, while *Cx. molestus* is a sister group to this pair. *Cx. quinquefasciatus* is considered the most divergent taxon… All internal nodes showed 100% statistical support (100% ultrafast bootstrap and full concordance factor gCF = 100, sCF = 100)."

Tree topology stated: `(molestus,(pallens,pipiens),quinquefasciatus)`

### What the new pipeline produced

Newick from `results/phylo/concord.cf.tree`:
```
(Cx_molestus:0.0325,(Cx_pallens:0.0235,Cx_quinquefasciatus:0.0689)100/59.3/45.4:0.0141,Cx_pipiens:0.0339);
```

Topology: `(molestus,(pallens,quinquefasciatus),pipiens)` — **pallens and quinquefasciatus are sisters**, not pallens and pipiens. The internal branch has full bootstrap (100) but only 59.3% gCF and 45.4% sCF.

### Downstream paragraphs that need rewriting

- **Section 3.2** (paragraphs around line 519): "*Cx. pallens* and *Cx. pipiens* are sister taxa" is no longer true. Replace with "pallens and quinquefasciatus are sisters." Remove the claim that "All internal nodes showed 100% statistical support… full concordance factor gCF = 100, sCF = 100" — replace with the actual support values and note the substantial gene-tree discordance.
- **Section 4 Discussion** (around line 755): "the phylogenomic supermatrix places *Cx. pallens* and *Cx. pipiens* as sister taxa… with complete concordance across all 8,616 gene trees at every node" — both clauses are wrong. The pallens+quinquefasciatus sister relationship, with only 59.3% gCF, actually supports the hybrid-origin hypothesis for *Cx. pallens* even more directly (a hybrid signal would be expected to produce exactly this kind of topology incongruence). You should rewrite that paragraph to capitalize on this stronger signal rather than retract it.
- **Section 5 Conclusions** (around line 929): "resolves *Cx. pallens* and *Cx. pipiens* as sister taxa under a topology compatible with the hybrid-origin hypothesis" — needs to be flipped to pallens+quinquefasciatus.
- **Figure 2a caption** (line 561): "Node labels show bootstrap/gCF/sCF support" — fine as-is, but make sure the figure shows `100/59.3/45.4` not `100/100/100`.

---

## 2. Cx. tarsalis status — METHODS AND TABLE FIX NEEDED

### What the paper claims
- Section 2.1: "*Cx. tarsalis* (GCA_016859205.1; 50 Mb) was included as a phylogenetic outgroup for tree rooting"
- Section 2.1: BUSCO completeness ≥ 95% as a selection criterion
- Discussion treats this as a 5-species pangenome with proper outgroup

### What the data say
- QUAST: total length 50,989,260 bp (51 Mb), 32 contigs — this is **not a chromosome-scale assembly**; it is roughly 10% of the expected Culex genome size (~550 Mb)
- BUSCO protein-mode (diptera_odb10): 0.0% complete, 99.9% missing
- Concordance tree (`concord.cf.tree`) has only 4 taxa — **tarsalis is not in the tree**; it could not be used for rooting because too few single-copy orthologs were recovered

### What to do

**Confirmed via NCBI (May 2026):** The actual `GCA_016859205.1` assembly (CtarK1, Main et al. 2021, PacBio HiFi) is **790 Mb**, not 51 Mb. Genome BUSCO 84.8% complete per the original publication. So the file currently in `resources/genomes/Cx_tarsalis.fasta` (51 Mb, 32 contigs) is a truncated or wrong download — about 6.5% of the actual assembly. This is fixable: re-download from NCBI and rerun the affected downstream rules (liftoff for tarsalis, extract_proteins, busco_proteins, orthofinder, partition_pangenome, sco extraction, phylo, CAFE).

Effect on conclusions of fixing tarsalis:
- **Species tree rooting**: would resolve cleanly with tarsalis as outgroup. The unrooted ingroup topology `(molestus, (pallens, quinquefasciatus), pipiens)` would gain a defined root but the internal `(pallens, quinquefasciatus)` sister relationship will not change.
- **Pangenome partition**: uses 4 ingroup taxa only; not affected.
- **CAFE**: uses 4 ingroup taxa only; not affected.
- **Methods claim "BUSCO ≥ 95%" and "50 Mb assembly"**: both need correction regardless. The assembly is 790 Mb and the published BUSCO is 84.8%, which still fails the ≥95% criterion the paper specifies. The selection-criterion sentence in Section 2.1 will need either rewording (drop the ≥95% language) or honest acknowledgment that tarsalis was retained as outgroup despite not meeting the criterion.

This is the single most consequential issue. Recommend re-downloading and rerunning the affected steps before final submission.

---

## 3. Pangenome numbers (Section 3.1, Table 1, Discussion, Conclusions)

Tiny shifts (≤4 in any cell). Narrative unchanged, but numbers should be refreshed for accuracy.

| Item | Paper | New run | Δ |
|---|---|---|---|
| Total orthogroups | 16,369 | 16,369 | 0 |
| Total ingroup genes assigned | 92,540 | 92,539 | −1 |
| Core orthogroups | 11,653 (71.2%) | **11,651 (71.2%)** | −2 |
| Shell orthogroups | 4,205 (25.7%) | **4,206 (25.7%)** | +1 |
| Cloud orthogroups | 511 (3.1%) | **512 (3.1%)** | +1 |
| Cloud_molestus | 200 | 200 | 0 |
| Cloud_pallens | 114 | 114 | 0 |
| Cloud_pipiens | 142 | 142 | 0 |
| Cloud_quinquefasciatus | **55** | **56** | +1 |

The chi-square test result (X²=85.36, p<2e-16) will move by a hair but the narrative is identical — rerun once you have the figures done.

### Table 1 per-species partitions

Per-species gene counts in each compartment (from `partitioned_orthogroups.tsv`):

| Species | Paper Core | New Core | Paper Shell | New Shell | Paper Cloud | New Cloud |
|---|---|---|---|---|---|---|
| Cx. molestus | 18,109 | **18,108** | 3,625 | 3,625 | 848 | 848 |
| Cx. pallens | 18,118 | 18,118 | 3,992 | **3,993** | 474 | 474 |
| Cx. pipiens | 18,177 | 18,177 | 4,436 | **4,437** | 560 | 560 |
| Cx. quinquefasciatus | 18,998 | **18,999** | 5,036 | **5,032** | 164 | **168** |

### Table 1 assembly sizes — small but should be fixed

QUAST gives slightly different sizes than the paper rounding:

| Species | Paper (Mb) | QUAST (Mb) |
|---|---|---|
| Cx. quinquefasciatus | 578 | **573** |
| Cx. pallens | 568 | **566** |
| Cx. molestus | 560 | 560 |
| Cx. pipiens | 533 | 533 |
| Cx. tarsalis | 50 | 51 |

### Protein counts (Section 2.2)

| Species | Paper | New |
|---|---|---|
| Cx. quinquefasciatus | 24,199 | 24,199 ✓ |
| Cx. pipiens | 23,174 | 23,174 ✓ |
| Cx. pallens | 22,585 | **22,586** (+1) |
| Cx. molestus | 22,582 | 22,582 ✓ |

---

## 4. CAFE results (Section 3.3, Discussion)

All small drift. Narrative ("every branch shows expansion bias; quinquefasciatus the strongest") unchanged.

| Item | Paper | New |
|---|---|---|
| Significant families (p<0.05) | **596** | **615** (+19) |
| Lambda | 0.0945 | **0.0985** |
| Alpha | 0.168 | **0.165** |
| Final −lnL | (not in paper) | 49,122 |

Per-branch expansions/contractions:

| Lineage | Paper exp/con | New exp/con |
|---|---|---|
| Cx. pipiens | 3,185 / 1,073 | **3,192 / 913** |
| Cx. pallens | 3,104 / 1,376 | **3,107 / 1,287** |
| Cx. molestus | 3,050 / 1,593 | **3,030 / 1,426** |
| Cx. quinquefasciatus | 3,671 / 866 | **3,646 / 770** |

Binomial p<2e-16 holds for all branches.

---

## 5. SCO count (Section 2.4, 3.2)

- Paper: 8,616 SCOs
- New: **8,632 SCOs** (after filtering)

Trivial drift, but update the number wherever it appears (Section 2.4, Section 3.2, Figure 2 caption, Discussion).

---

## 6. TE content (Section 3.5, Discussion)

All slightly higher in new run (~+0.7% across the board). Narrative ("52–54% TE content, remarkably consistent") still holds, in fact even more tightly.

| Species | Paper TE % | New TE % |
|---|---|---|
| Cx. molestus | 52.88 | **53.66** |
| Cx. pallens | 53.91 | **54.65** |
| Cx. pipiens | 52.14 | **52.98** |
| Cx. quinquefasciatus | 53.36 | **54.00** |

New range is 52.98–54.65% (paper said 52.14–53.91%); update Section 3.5 and Discussion accordingly.

---

## 7. Things in the paper the pipeline didn't actually run

- **Synteny analysis (Section 3.4, Table 2, Figure 4)** — the pipeline has `workflow/rules/synteny.smk` and `workflow/scripts/check_3Rb_inversion.py`, but synteny outputs aren't in the `all:` target and weren't produced this run. The ANI numbers (97.9–98.3%) and synteny block counts (44,274 / 66,511 / etc.) in the paper come from somewhere else. Either find the prior output, or re-run with synteny wired in.
- **Anvi'o pangenome (Section 3.6, Figure 6)** — pipeline has no anvi'o rule. This must have been run manually outside Snakemake. Mention this in Methods or remove the section if results aren't reproducible.
- **GO enrichment** — pipeline produces `results/functional/go_enrichment_results.tsv` but the script is a stub that writes a placeholder. Paper doesn't currently reference GO enrichment explicitly, but if you want to add it (and you have eggNOG annotations for all four samples ready), I can write a real Fisher's-exact-style analysis.

---

## 8. Other text-level fixes

- **Section 2.1 duplicated sentence** (line 194): "Four chromosome-scale genome assemblies representing members of the *Culex pipiens* Four chromosome-scale genome assemblies representing members of the *Culex pipiens* species complex were retrieved…" — copy-paste duplication.
- **Section 2.5 vs 2.3 redundancy** — the CAFE5 description appears in both the highlighted block above Section 2.3 and again in 2.5; the actual prose body of 2.3 is still a placeholder `*[Describe OrthoFinder…]*`.
- **Figure 5 caption** says "11,653 core orthogroups; petals show form-specific cloud orthogroups. Shell: 4,205." — update to 11,651 core / 4,206 shell.
- **Figure 1a caption** says "The largest intersection (11,655) represents core orthogroups" — that 11,655 number appears nowhere else and disagrees with both old (11,653) and new (11,651). Worth verifying which one was actually in the UpSet plot.

---

## What I'll do next (with your sign-off)

Now that I know which numbers need to change, I'll:

1. Generate the figures from the new outputs (Figs 1–5 + key-families heatmap + supplementary tables) so you have current visuals before editing the text.
2. Once you're happy with the figures, I can also produce a tracked-changes version of the .docx that updates all the in-text numbers above (the safe edits) without touching the prose around the species-tree topology issue — that one needs your judgment, not just a number swap.

Tell me whether to proceed with the figures first, or if you want to talk through the tarsalis issue or the topology rewrite before I generate anything.
