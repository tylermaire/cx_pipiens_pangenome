# bioRxiv Submission Checklist
### *Pangenome Analysis Reveals Gene Content Variation and Evolutionary Dynamics Across the* Culex pipiens *Species Complex*

Prepared for submission at **https://www.biorxiv.org/submit-a-manuscript**. Work through the sections in order. Items marked **[YOU]** require your action (login, consent, final submit) — these cannot be done on your behalf.

---

## 0. Before you start — what I can and cannot do

- **I can:** prepare/verify all files (PDF, figures, tables), draft the metadata text (abstract, author list, summary), and walk you through the form.
- **I cannot:** log into your bioRxiv/account, accept the license or terms, enter author email verifications, or click the final **Submit** button. bioRxiv requires the submitting author to do these. **[YOU]** owns every consent/submit step.

---

## 1. Account & corresponding author **[YOU]**

- [ ] Log in (or register) at bioRxiv with an **institutional email** if possible.
- [ ] Corresponding author on the manuscript is **Tyler Maire — t.maire@irmcd.org** (confirm this inbox is monitored; bioRxiv sends the confirmation + DOI here).
- [ ] Have an **ORCID iD** ready for at least the submitting author (bioRxiv prompts for it; recommended).

---

## 2. Files to upload

bioRxiv accepts a single PDF of the main text **or** a manuscript file plus separate figures. The simplest, lowest-risk path is the **single combined PDF** (already prepared).

**Main manuscript (required):**
- [ ] `Culex pipiens Pangenome Analysis.pdf` — 24 pages, figures + tables embedded, US Letter. ✅ generated from the clean draft 5. **This is the file to upload as the manuscript.**

**Figures (bioRxiv prefers high-res individual files; optional if using the combined PDF, but recommended):**
- [ ] Figures 1–9 are in `figures/` as both `.pdf` and `.png` (300 dpi).
- [ ] ⚠️ **Figure 1 caveat:** the embedded Figure 1 in the PDF is the corrected version with *diptera_odb10* in italics. The standalone `figures/Figure_1_assembly_quality.png` was regenerated to match, but **verify the standalone `Figure_1_assembly_quality.pdf` is the italic version before uploading it separately** (re-run `make_all_figures.R` if unsure, or just rely on the combined PDF). Do **not** upload `Figure_1_v2_italic.png` under a figure-1 slot if you also upload the standard one — pick one.

**Supplementary material:**
- [ ] `figures/Sup Material/Supplementary_Tables.xlsx` — Supplementary Tables S1–S8.
- [ ] `figures/Sup Material/Figure_S1_synteny_dotplots.pdf` — Supplementary Figure S1.
- [ ] `figures/Sup Material/Figure_S2_pipeline_flowchart.pdf` — Supplementary Figure S2.
- [ ] Upload these in the **Supplementary Material** section (not as main figures).

---

## 3. Metadata to enter in the form

**Title** (paste exactly):
> Pangenome Analysis Reveals Gene Content Variation and Evolutionary Dynamics Across the *Culex pipiens* Species Complex

**Authors & affiliations** (enter in this order; mark Tyler as corresponding):
1. Tyler Maire — ¹ University of Maryland Global Campus; ² Indian River Mosquito Control District, Vero Beach, FL, USA *(corresponding; t.maire@irmcd.org)*
2. Erin Maley — ¹ University of Maryland Global Campus
3. Theresa Miller — ¹ University of Maryland Global Campus
4. Kyle Kosinski — ² Indian River Mosquito Control District, Vero Beach, FL, USA

- [ ] Add each author with first/last name and affiliation number(s).
- [ ] **[YOU]** Confirm every co-author consents to posting (bioRxiv emails co-authors; they must not object).

**Abstract** — paste the abstract from the manuscript (the ~290-word paragraph beginning "The *Culex pipiens* species complex contains the principal northern-hemisphere vectors…"). *Note: bioRxiv's abstract box may strip italics; that's fine for the preprint.*

**Subject area / category** (bioRxiv requires one):
- [ ] Primary: **Genomics** *(best fit)* — alternatives: **Evolutionary Biology** or **Bioinformatics**.

**Keywords:** Culex pipiens; species complex; pangenome; phylogenomics; gene family evolution; hybrid origin; West Nile virus

---

## 4. Declarations (the form asks for each)

- **Funding:** "This research received no external funding." (matches manuscript)
- **Competing interests:** "The authors declare no competing interests." (matches manuscript)
- **Data availability:** Code + intermediate results at **https://github.com/tylermaire/cx_pipiens_pangenome** ; genome assemblies via NCBI accessions (Table 1) and the *Cx. tarsalis* CtarK1 assembly via OSF **https://osf.io/mdwqx**.
  - [ ] ⚠️ **Make the GitHub repo public** before submitting (reviewers/readers must be able to reach it). Confirm it's not private.
- **Ethics:** Not applicable (no human/animal subjects; public genome data only).

---

## 5. License (choose one) **[YOU]**

bioRxiv requires you to pick a reuse license. Common choices:
- **CC-BY 4.0** — most permissive, maximizes reuse/visibility (often required if you later submit to certain journals).
- **CC-BY-NC-ND 4.0** — no commercial use, no derivatives.
- **CC0** — public domain dedication.
- **"No reuse" / All rights reserved** — bioRxiv's most restrictive option.

- [ ] **[YOU]** Pick the license. If you may submit this to a journal later, check that journal's preprint policy first; **CC-BY** is the safe default for broad compatibility. *(I can't choose this for you.)*

---

## 6. Final pre-flight checks

- [ ] Author list and order match the manuscript (4 authors, Kyle 4th).
- [ ] Corresponding email is **.org** (fixed) and reachable.
- [ ] GitHub repo is **public** and the latest commit (with the PDF + figure fixes) is pushed.
- [ ] Combined PDF opens and all 9 figures + 2 tables render (verified: 24 pp).
- [ ] Supplementary files attached in the supplementary section, not as main figures.
- [ ] Competing-interests and funding statements entered.
- [ ] License selected.

---

## 7. Submit **[YOU]**

- [ ] Review bioRxiv's PDF proof it auto-generates (it converts your upload). Check figures/tables survived conversion.
- [ ] Click **Submit**. bioRxiv screens for non-science/plagiarism (1–2 business days) before posting; you'll get the DOI by email at t.maire@irmcd.org.

---

## Outstanding items worth deciding before you post (optional, from the editorial review)

These are *not* blockers for a preprint, but bioRxiv posts are citable and semi-permanent — consider whether to address first:
- The hybrid-origin claim still rests on concordance factors (no formal D-statistic). The text frames it as "consistent with" prior evidence, which is defensible for a preprint.
- Abstract still ~290 words and notes "first pangenome" — fine for bioRxiv; a target journal may want the abstract trimmed and the novelty claim hedged ("first chromosome-scale comparative pangenome").
- A few framing items (the "pangenome" terminology caveat, "rooting is assumed") remain in the discussion but not the abstract.

If you want any of these addressed before posting, tell me and I'll do them; otherwise the current draft 5 is submission-ready as a preprint.
