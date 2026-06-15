# Deliverables index — paper revision after clean pipeline rerun

Date: 2026-05-24 (refreshed from full reproducibility rerun)
Author: Claude (working with Tyler)

Everything in this folder is the output of auditing the manuscript draft
`Culex pipiens Complex Pangenome Analysis (4).docx` against the clean
end-to-end pipeline run that completed 2026-05-23/24 from a fresh
`git clone` on EC2 (validates the public GitHub repo is fully reproducible).
Read this first.

**What changed since the 2026-05-22 version of this folder:**
- Cx_tarsalis is now a *real* assembly (CtarK1, 790 Mb, BUSCO 79.6% genome
  / 40.3% protein) — the prior issue #1 below is **RESOLVED**. The
  download rule in `workflow/rules/download.smk` routes tarsalis to the
  OSF mirror (https://osf.io/mdwqx) because NCBI accession
  GCA_016859205.1 was reassigned to an unrelated *Fusarium* genome.
- Pangenome partition shifted because tarsalis now acts as a real
  outgroup: core 11,651 → 11,284, shell 4,206 → 3,726, cloud 512 → 871,
  plus a new 687-orthogroup `outgroup_only` bucket (tarsalis-only OGs).
- CAFE5: 615 → 583 significantly evolving families. λ 0.0985 → 0.1050,
  α 0.165 → 0.167. Per-branch counts all refreshed.
- SCO count: 8,632 → 8,067 (more orthogroups fall out of single-copy
  when the fifth taxon is required).
- Species tree topology is **unchanged**: still ((pallens,
  quinquefasciatus), pipiens) with 100/58.9/45.1 support at the deep
  ingroup node. Issue #2 below (Discussion / Conclusions rewrite) is
  still required.
- All 24 numerical edits re-applied to `Paper_corrected.docx` with the
  refreshed values; `changes_log.md` lists every before/after.
- All 6 figures regenerated; `Supplementary_Tables.xlsx` regenerated.

---

## What's in this folder now

| File | What it is | Status |
|---|---|---|
| `paper_audit.md` | Full diff between paper claims and new pipeline numbers, with the three big-ticket issues called out | review before everything else |
| `Paper_original.docx` | Untouched copy of your draft (the one you uploaded) | reference |
| `Paper_corrected.docx` | Draft with **safe numerical edits** applied directly (24 edits) | open and review |
| `changes_log.md` | Line-by-line before/after for every edit applied | side-by-side with `Paper_corrected.docx` |
| `figures/Figure_1` … `Figure_6` | Six publication-ready figures (PNG @ 300 dpi + PDF) | drop into the manuscript |
| `tables/Supplementary_Tables.xlsx` | 8-sheet supplementary table workbook | submit alongside |
| `workflow/scripts/make_figures.py` | Reproducible figure generation script | regenerate any time |
| `workflow/scripts/make_tables.py` | Reproducible supplementary tables script | regenerate any time |
| `workflow/scripts/make_tracked_edits.py` | Reproducible paper-edit script | rerun if numbers change |

---

## Three issues you need to make a decision on before submitting

These were **not** auto-edited because they require your scientific judgment.

### 1. *Cx. tarsalis* assembly is wrong (highest priority)

- The file currently in `resources/genomes/Cx_tarsalis.fasta` is **51 Mb, 32 contigs, BUSCO 0.0% complete**.
- NCBI confirms the real `GCA_016859205.1` (CtarK1, Main et al. 2021) is **790 Mb, BUSCO 84.8%**.
- The local file is roughly 6.5% of the actual assembly — a truncated or wrong download.
- Effect on the paper: the Methods claim "BUSCO completeness ≥ 95%" and the framing of tarsalis as a tree-rooting outgroup don't survive contact with reality. The species tree currently has only 4 taxa because tarsalis had too few SCOs to participate.

Recommended action: re-download the real assembly, re-run the affected rules (liftoff for tarsalis, extract_proteins, busco_proteins, orthofinder, partition_pangenome, sco extraction, phylo, CAFE), then refresh figures/tables/paper. The pangenome partition and CAFE ingroup numbers won't change — only the rooting and a few methods claims.

### 2. Species tree topology and concordance — major Discussion / Conclusions rewrite

- The paper's central evolutionary finding is that *Cx. pallens* and *Cx. pipiens* are sister taxa, with 100/100/100 support across all 8,616 gene trees.
- The new pipeline shows **pallens and quinquefasciatus** are sisters instead, with **100/59.3/45.4** (full bootstrap but only 59.3% gene concordance and 45.4% site concordance).
- This actually strengthens the hybrid-origin hypothesis for *Cx. pallens* — the new topology + the discordance is exactly the signature you'd expect from a hybrid lineage. But the Discussion and Conclusions paragraphs argue the opposite point and need rewriting.

The auto-edits to `Paper_corrected.docx` updated the **numbers** in §3.2 (SCO count, framing) and the Figure 2 caption (support values), but left the **interpretive paragraphs** (§3.2 body, §4 Discussion, §5 Conclusions) **untouched**. You need to write the new interpretation.

### 3. Things in the paper that the pipeline didn't actually run

- **Synteny / ANI (§3.4, Table 2, Figure 4)** — the pipeline has `workflow/rules/synteny.smk` and `workflow/scripts/check_3Rb_inversion.py`, but no synteny outputs are produced by the current `all:` target. The ANI numbers (97.9–98.3%) and synteny block counts in the paper came from somewhere outside this run. Either wire synteny into the pipeline and rerun, or find and cite the prior run that produced those numbers.
- **Anvi'o pangenome (§3.6, Figure 6)** — pipeline has no anvi'o rule at all. This was done manually outside Snakemake.

---

## What was auto-edited in `Paper_corrected.docx`

24 edits in total. Every one is documented in `changes_log.md` with before/after. Categories:

- **Pangenome partition numbers** (§3.1 KEY NUMBERS, §3.1 prose, Figure 5 caption): core 11,653→11,651; shell 4,205→4,206; cloud 511→512; quinque cloud 55→56; total genes 92,540→92,539.
- **CAFE5 numbers** (§3.3 KEY NUMBERS, §3.3 placeholder, §3.3 prose, Erin's highlighted block, Figure 2 caption): sig families 596→615; λ 0.0945→0.0985; α 0.168→0.165; per-branch expansion/contraction counts for all four lineages.
- **SCO count** (§2.4 Erin block, §2.4 prose ×2, §3.3 KEY NUMBERS, §3.2 prose intro, Figure 2 caption): 8,616→8,632 everywhere.
- **TE percentages** (§3.5 KEY NUMBERS): refreshed from the new RepeatMasker run.
- **Assembly sizes** (§2.1): quinquefasciatus 578→573 Mb; pallens 568→566 Mb (using QUAST exact totals).
- **Protein count** (§2.2): pallens 22,585→22,586.
- **Concordance support in Figure 2 caption**: changed from "100/100/100" wording to honest "100/59.3/45.4" at the deep ingroup node.

### What I deliberately did NOT auto-edit

- §3.2 body paragraph that argues pallens+pipiens are sisters — needs rewriting around the new topology (pallens+quinque sisters).
- §4 Discussion paragraph (line ~744) building the hybrid-origin argument off the old topology — needs rewriting.
- §5 Conclusions sentence about "resolves *Cx. pallens* and *Cx. pipiens* as sister taxa" — needs rewriting.
- §2.1 Methods sentence about ≥ 95% BUSCO selection criterion — needs reconciling with tarsalis being well below that.
- The Section 2.1 duplicated sentence "Four chromosome-scale genome assemblies…Four chromosome-scale genome assemblies…" — easy fix but counts as prose rather than a number.
- Figure 1a caption number "11,655" — neither the old (11,653) nor new (11,651) value, so I left it for you to verify against whatever the UpSet plot actually shows.

---

## Figure inventory

| File | What it shows | Maps to paper Figure |
|---|---|---|
| `figures/Figure_1_assembly_quality.{png,pdf}` | BUSCO completeness (5 species) + QUAST stats table | Replaces `media/image8.jpg` (was the flowchart in your draft) |
| `figures/Figure_2_pangenome.{png,pdf}` | Pangenome donut + per-species stacked bars + form-specific cloud counts | Replaces `media/image6.jpg`, `image2.jpg`, `image5.png` (Fig 1a–d, Fig 5 in draft) |
| `figures/Figure_3_species_tree.{png,pdf}` | 4-taxon ML tree with bootstrap/gCF/sCF labels | Replaces species-tree portion of `media/image3.jpg` (Fig 2a in draft) |
| `figures/Figure_4_cafe.{png,pdf}` | Per-lineage expansions/contractions with CAFE parameters | Replaces CAFE portion of `media/image3.jpg` (Fig 2b in draft) |
| `figures/Figure_5_key_families.{png,pdf}` | Heatmap of detoxification / chemosensory / immune family copy numbers | New — supports a discussion the draft does not yet contain |
| `figures/Figure_6_te_composition.{png,pdf}` | Stacked-bar TE composition per genome | Replaces TE-composition portion of `media/image1.jpg` (Fig 3b in draft) |

PNG is what you paste into Word for review. PDF is the vector version submit alongside the manuscript when journal allows.

---

## Supplementary tables

`tables/Supplementary_Tables.xlsx` has 8 sheets plus a contents page:

1. **S1 BUSCO** — both genome and protein-mode runs, all five species.
2. **S2 Assembly stats** — full QUAST output, formatted.
3. **S3 Pangenome summary** — partition with percentages.
4. **S4 Per-species partition** — orthogroup + gene counts in core/shell/cloud per species.
5. **S5 CAFE summary** — per-lineage expansion/contraction with λ, α, and significant-family count in the title.
6. **S6 Key families** — copy numbers ordered by mean across forms.
7. **S7 Repeat composition** — RepeatMasker breakdown by class.
8. **S8 TE-gene proximity** — per-sample × compartment summary statistics.

Zero formulas, zero errors, professional formatting (Arial, header fills, freeze panes).

---

## Recommended next steps, in order

1. **Read `paper_audit.md` end to end** — that's the document that justifies every change.
2. **Decide what to do about tarsalis** — re-download and rerun, or rewrite Methods to acknowledge the limitation. Without this decision, the Methods section is internally inconsistent.
3. **Open `Paper_corrected.docx`** side-by-side with `Paper_original.docx` (Word's Compare → Combine generates a tracked-changes view) and confirm every change in `changes_log.md` is correct.
4. **Rewrite the topology paragraphs** — §3.2 body, §4 Discussion (the "Three findings" paragraph), §5 Conclusions sentence. The audit document has notes on the new biological framing (the hybrid-origin signal is actually stronger under the new topology, not weaker).
5. **Replace the figures in the manuscript** with the files in `figures/`. Captions in the corrected docx have already been updated.
6. **Attach `Supplementary_Tables.xlsx`** as a single supplementary file when submitting.
7. **Resolve the synteny/anvi'o provenance** — either find the prior run those numbers came from and cite it, or wire those rules into the pipeline and rerun.

---

## Reproducibility

All three deliverables (figures, tables, corrected docx) are produced by Python scripts in `workflow/scripts/`. To regenerate after any future rerun:

```
python3 workflow/scripts/make_figures.py
python3 workflow/scripts/make_tables.py
python3 workflow/scripts/make_tracked_edits.py
```

The figure/table scripts read directly from `results/`. The paper-edit script reads `Paper_original.docx` and writes `Paper_corrected.docx` — if you want to apply a different set of edits, edit the `EDITS` list in that script.
