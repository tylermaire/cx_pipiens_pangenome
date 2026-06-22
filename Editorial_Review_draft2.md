# Editorial Review — *Pangenome Analysis Reveals Gene Content Variation and Evolutionary Dynamics Across the* Culex pipiens *Species Complex*

**Reviewer stance:** Read as a handling editor + methods reviewer at a journal such as *Insects*, *G3*, *BMC Genomics*, or *Genome Biology and Evolution*. Overall this is a competent, well-organized, methodologically transparent manuscript with a genuinely interesting central claim (chromosome-scale corroboration of a hybrid origin for *Cx. pallens*). The reproducible Snakemake workflow is a real strength. My recommendation would be **major revision** — not because the science is weak, but because of (1) several internal-consistency errors and leftover editing artifacts that must be fixed before review, (2) a handful of framing/scope claims that need softening or support, and (3) some stylistic patterns that read as machine-generated and will draw reviewer suspicion. None of these are fatal.

Issues are grouped by severity. Verified items were checked computationally against the manuscript's own numbers.

---

## A. Must-fix before submission — errors and artifacts (a reviewer will see these immediately)

These are the ones that would make me, as an editor, question how carefully the manuscript was prepared.

**A1. Two stray "?" editing marks left in the body text.** These look like unresolved author queries and *must* go:
- §1, ¶2 (Introduction): *"these differences have profound epidemiological consequences: **form? molestus** is implicated as a critical bridge vector…"* — should read "the *molestus* form" or "*Cx. p. molestus*".
- §3.2 (Pangenome architecture): *"92,821 genes (98.2 % of ingroup gene calls) were assigned to **15,881 ? orthogroups**…"* — delete the stray "?". (The number 15,881 is correct; see A2.)

**A2. "Core" numbers are internally inconsistent / conflated.** In §3.2:
- Text: "11,284 core … orthogroups" (this is correct; 11,284 + 3,726 + 871 = 15,881 ✓, and 71.1/23.5/5.5 % ✓).
- Same paragraph: "core orthogroup membership was essentially identical across forms (**11,275** orthogroups containing core members per species)". The 11,275 vs 11,284 discrepancy is unexplained — pick one and explain the 9-orthogroup difference, or make the wording consistent.
- Table 1 "Core genes" column lists ~18,450–19,439 per form. These are **genes**, but the text discusses **core orthogroups** (~11,275). Readers will collide these two numbers. Clarify in the Table 1 caption and §3.2 that one is gene counts and the other orthogroup counts, and make sure the relationship (≥11,275 orthogroups containing 18,476 genes in *molestus*, etc.) is stated explicitly.

**A3. Table 1 vs Methods protein-count ordering.** Methods §2.2 reports protein totals as "24,531 (*Cx. quinquefasciatus*), 23,089 (*Cx. pallens*), 23,199 (*Cx. molestus*), 23,738 (*Cx. pipiens*) and 20,026 (*Cx. tarsalis*)." Table 1's "Proteins" column matches, but §3.1 says the range runs "from 23,089 (*Cx. pallens*) to 24,531." That's fine, but double-check every place the per-form protein counts appear — there are at least four — so they're identical throughout. (I verified the values are mutually consistent; the risk is transcription drift across sections.)

**A4. Author contributions omit two of the four authors.** The CRediT statement lists only T.M., E.M., and T.Mi. With **Kyle Kosinski now added as the 4th author, K.K. must appear** (e.g., "K.K.: Writing – review & editing" or whatever is accurate). Also note E.M. and T.Mi. are each credited only with "Writing – review & editing" — journals increasingly query whether co-authors meet authorship criteria with review/editing alone. Confirm their contributions are fully and accurately described.

**A5. Correspondence email domain mismatch.** Correspondence is given as `t.maire@irmcd.com`, but the institution is the Indian River Mosquito Control District and the editor comments came from `…@irmcd.org`. Verify `.com` vs `.org` — a bounced corresponding-author email is a real submission problem.

**A6. Title has a hard line break mid-sentence.** The title is split across two paragraphs ("…Gene Content Variation and" / "Evolutionary Dynamics Across the *Culex pipiens* Species Complex"). Make sure this is a single title field, not two paragraphs, or it will mis-import into submission systems and reference managers.

**A7. Reference / citation housekeeping.**
- "RepeatMasker Open-4.0 (Smit et al., **2013–2015**)" — a date range as a citation year is unusual for most journals; many will want a single year or "Smit, Hubley & Green 2013–2015. RepeatMasker Open-4.0. http://www.repeatmasker.org" with an access date. Also confirms the text cites it as "Smit et al., 2013–2015".
- *Vinogradova 2000, 2003* and *Smith & Fonseca 2004* are cited in text and present in the list — good — but they are **not in the verified reference set I was given earlier**; double-check each against the actual source (volume/pages) since they were not part of the previously verified bibliography.
- Confirm every Supplementary item (S1–S8, Figs S1–S2) actually exists and is uploaded; the text references all of them.

---

## B. Scientific / framing issues a reviewer will push on

**B1. "Pangenome" terminology applied across species.** You pre-empt this well in §1 ("it is in this clade-level sense that we use the term") and §4.1, and that's the right instinct — but a hardline reviewer will still object that a four-genome, between-species "pangenome" is really a comparative orthology/gene-content analysis. Consider stating the caveat once more in the Abstract (a single clause) so it isn't seen as overclaiming in the part most people read.

**B2. The headline hybrid-origin claim rests on concordance factors alone.** The discordance signal (gCF 58.9 %, sCF 45.1 %, near-balanced alternative quartets) is genuinely suggestive, and you are appropriately careful ("consistent with", "molecular signature expected when"). But you explicitly note you have **not** run the formal tests that would distinguish ILS from introgression (D-statistics / ABBA-BABA, phylogenetic networks, *Dsuite*, QuIBL, HyDe). A reviewer will almost certainly ask why these weren't run, since the input data (SCO alignments, gene trees) are already in hand. Two options: (a) run at least one formal introgression test and report it, which would substantially strengthen the central claim; or (b) soften the framing further and frame the hybrid-origin support as "corroborating prior population-genetic evidence" rather than independent genomic demonstration. Option (a) is the stronger paper. With n = 4 ingroup taxa a four-taxon D-statistic with *Cx. tarsalis* as outgroup is directly feasible.

**B3. Four-taxon design and rooting.** §4.1 defends the four-genome panel honestly. The weak point you flag yourself: *Cx. tarsalis* "does not enter the concatenation SCO phylogeny at sufficient density," so the ingroup tree is effectively **unrooted** and "Tree rooting is illustrative" (Fig. 3 caption). That means the directionality of the (*pallens*+*quinquefasciatus*) vs (*pipiens*+*molestus*) split — central to the hybrid-origin narrative — is not actually established by your tree. Be explicit in §3.3/§4.2 that the rooting is assumed from prior work, not demonstrated here, and temper any language implying you *resolved* the root.

**B4. ANI / *Cx. tarsalis* omission framed as a transient glitch.** §3.5 and Methods §2.6 attribute the missing *Cx. tarsalis* ANI to "a transient genome-download failure in the analysis environment." Reviewers dislike seeing pipeline/infrastructure excuses in a results section — it reads as incomplete analysis. Since you clearly have the *Cx. tarsalis* assembly (it's in BUSCO, QUAST, and the orthology partition), just **run skani on it** and complete the matrix; it's a one-command fix and removes an obvious reviewer target. Same applies to the synteny analysis excluding *Cx. tarsalis*.

**B5. Alignment-based ANI methodology is heterogeneous.** §2.6: for three pairs ANI comes from Liftoff gene-set alignments (60–75 Mb CDS), for the other three from "partial chromosome-1 alignments (first ~2 Mb)." Mixing a 60–75 Mb basis with a ~2 Mb basis across pairs of the *same* matrix is methodologically uneven and a reviewer will flag it. Either compute all six pairs on a consistent basis or drop the alignment-based matrix to supplementary with an explicit "cross-check only, non-uniform basis" caveat (you partly do this, but the unevenness should be stated plainly).

**B6. Inversion calls are acknowledged to be noisy.** §3.5 candidly says many small inversions "likely combine bona fide micro-inversions with annotation-projection noise." Good. But you still report "258–294 candidate inversions ≥ 100 kb per pair" in the Abstract and Conclusions as if it were a clean result. Consider reporting only the high-confidence large blocks (≥ 1 Mb) as headline numbers and relegating the full counts to "candidate" status in supplementary, so the abstract isn't carrying a number you partly disavow in the results.

**B7. §3.8 "Key gene families" is keyword-match shallow.** Copy numbers are derived from "gene-name keyword matches to the Liftoff-transferred annotations." For detoxification/chemosensory/immunity families — which are exactly the families prone to mis-annotation and naming inconsistency after annotation transfer — keyword matching is a weak method and a reviewer in vector biology will say so. Either add a sentence acknowledging this limitation explicitly, or buttress with an HMM/orthogroup-based count for at least the headline families (P450, OR).

**B8. CAFE interpretation: expansion bias may be partly methodological.** Every branch shows 2.5–5.1× expansion bias. Universal, strong expansion across *all* lineages can be an artifact of (i) using the *Cx. quinquefasciatus* RefSeq as the annotation source (reference-bias inflating its counts) combined with (ii) Liftoff under-transfer in the others, or (iii) the gene-count table including many small/spurious orthogroups. You note the ranking "mirrors the order of lineage protein-count totals" — which is itself a hint that the signal may track annotation completeness rather than biology. Add a sentence addressing whether the expansion bias survives correction for annotation-source asymmetry, or temper the claim.

**B9. Abstract reports "first pangenome" — confirm novelty claim.** "First pangenome characterisation of the … *Cx. pipiens* complex" is a strong priority claim. Make sure it's still true at submission (do a final literature check), and consider "first chromosome-scale comparative pangenome" to be precise and defensible.

---

## C. AI / automation tells (will draw reviewer or screening-tool suspicion)

You asked specifically for these. None individually proves machine assistance, but the **cluster** is what reads as generated. In order of how noticeable they are:

**C1. A near-verbatim sentence is repeated three times.** "*…the largest set of form-restricted [cloud orthogroups / clusters] despite having an intermediate total protein count*" appears in §3.2, §3.7, and §3.8 (the phrase "despite having an intermediate total protein count" occurs **3×**). Human authors vary their phrasing on re-statement; verbatim repetition of a distinctive clause is one of the strongest tells here. Reword two of the three.

**C2. Heavy reliance on a small set of interpretive hedges.** "consistent with" appears **10×** (4× as "is consistent with"), plus "expected when" 3×, "suggests/suggesting that" 3×. This hedge-cluster is characteristic of generated scientific prose. Keep the caution, but vary the constructions and cut redundant hedging.

**C3. Flowery interpretive flourishes that don't match the surrounding terse methods voice.** "on its face, paradoxical"; "the molecular **fingerprint** expected"; "molecular **signature**"; "molecular **corroboration**"; "strikingly uniform" (2×); "remarkable conservation." These rhetorical phrases are a stylistic register shift from the otherwise dry, quantitative writing — a common signature of AI-assembled drafts. Trim to one such phrase per section at most.

**C4. Self-referential "§" cross-references throughout (11×), including "we return to this in §4."** Internal "§3.2"/"§4" pointers and "we develop this interpretation in the Discussion" are unusually dense and read like an outline that was filled in. Most journals also prefer "Section 3.2" or "(see Discussion)" over the "§" glyph. Reduce the density and convert "§" to the journal's house style.

**C5. Methods describe the *computing environment* rather than just the method.** Phrases like "owing to a transient genome-download failure in the analysis environment" and "for the remaining three pairs, partial chromosome-1 alignments … were generated de novo" expose the ad-hoc, iterative way the analysis was run. Real methods sections describe what was done as a finished protocol, not the troubleshooting history. Rewrite B4/B5 items as clean protocol statements.

**C6. Over-explained rationale.** §4.1 opens "the rationale is worth making explicit because it shapes the interpretive weight of every result that follows," and §2 frequently explains *why* a standard step was done (e.g., why concordance factors expose disagreement, why an outgroup makes a partition rootable). This tutorial register — explaining textbook rationale to the reader — is a mild AI tell and also costs words. Assume an expert reader; cut the justifications for standard choices.

**C7. Minor mechanical tells.** Space-before-percent ("98.2 %") appears 69× — fine if consistent with house style, but several journals set "98.2%"; pick one. British/American spelling is mixed ("characterisation", "hybridisation", "behavioural" vs "visualization", "recognized") — **standardize to one** (the journal will specify; many *Culex* venues are US English). This mix is itself a subtle tell that text came from multiple sources.

---

## D. Smaller points / polish

- **Abstract is one ~290-word block.** Many journals cap abstracts at 200–250 words and/or want it structured. Tighten.
- **"CtarK1 … OSF mirror as a special case" (§2.1)** reads as pipeline-speak ("special case"). Reframe for a reader: explain the NCBI accession was reassigned and you used the authors' OSF deposit instead.
- **Figure 3 caption** says rooting "is illustrative; *Cx. tarsalis* does not enter the supermatrix at sufficient density." Good honesty, but see B3 — make the body text equally clear.
- **Highlights** bullet "~53% transposable-element content" vs body "52.69–54.43% interspersed repeats (+~3% simple/low-complexity)". The ~53% figure conflates interspersed-only with total — make the highlight and body use the same definition.
- **"medically-important"** (hyphenated, 3×) — most style guides write "medically important" (no hyphen after an -ly adverb). Minor but a copyeditor will change every instance.
- **Define abbreviations on first use** consistently: SCO, gCF/sCF, ILS, NNI, ANI, TE, CRediT initials. Most are defined, but verify NNI (§3.3) and minbit (§3.7) are spelled out.
- **"~99 % of total assembly span" / "~50% TE content"** — make sure tilde-approximations in methods are replaced with the actual computed values where you have them.

---

## E. What's genuinely good (keep)

So the review isn't all red ink:
- The **reproducible Snakemake + Conda workflow with a tagged commit and bundled intermediates** is exactly what reviewers want and is better than most comparable papers.
- **Internal quantitative consistency is high** — I recomputed the core/shell/cloud sum (15,881 ✓), the percentages (71.1/23.5/5.5 ✓), all four CAFE expansion ratios (5.1/4.2/3.0/2.5 ✓), the ANI mean (93.87 % ✓), and the anvi'o core fraction (56.5 % ✓), and they all check out.
- The **anvi'o cross-validation of the OrthoFinder partition** (§3.7, convergent 16,673 vs 15,881) is a nice methodological control.
- The **honest limitations section** (§4.5) and the explicit defense of the small panel (§4.1) show scientific maturity — they just need the AI-register trimmed.
- The **hybrid-origin synthesis** (§4.2) is well-sourced and is the paper's strongest contribution.

---

## Suggested priority order for revision

1. Fix all artifacts and inconsistencies in **A1–A7** (fast, mandatory).
2. Complete the *Cx. tarsalis* **ANI + synteny** (B4) and run **one formal introgression test** (B2) — these two additions most strengthen the paper.
3. Revise the **AI-tell language** in **C1–C6** (reword the 3× repeated sentence first).
4. Soften/clarify rooting and "pangenome" framing (B1, B3).
5. Copyedit pass for spelling consistency, hyphenation, and abstract length (C7, D).
