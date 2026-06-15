# Numerical corrections applied to Paper_corrected.docx

Each row below shows a single text replacement applied directly
to the manuscript. The original `Paper_original.docx` is preserved
for diff purposes. Open both files side-by-side or use Word's
Compare → Combine feature to generate a tracked-changes view.

**Edits applied: 24 of 24**

---

## Applied edits

### 1. §2.1 quinquefasciatus genome size 578→573 Mb (QUAST exact total)

**Before:**  (GCF_015732765.1; 578 Mb), 

**After:**   (GCF_015732765.1; 573 Mb), 

---

### 2. §2.1 pallens genome size 568→566 Mb

**Before:**  (GCF_016801865.2; 568 Mb), 

**After:**   (GCF_016801865.2; 566 Mb), 

---

### 3. §2.2 pallens protein count 22,585→23,089 (new annotation lift)

**Before:** , 22,585 for 

**After:**  , 23,089 for 

---

### 4. §2.4 SCOs 8,616→8,067 in Erin block

**Before:** We pulled 8,616 single copy orthologs

**After:**  We pulled 8,067 single copy orthologs

---

### 5. §2.4 CAFE parameters in Erin block

**Before:** Then CAFE5 v5.1.0 (gamma model, k=3) tested which gene families are expanding or contracting. Lambda=0.0945, alpha=0.168.

**After:**  Then CAFE5 v5.1.0 (gamma model, k=3) tested which gene families are expanding or contracting. Lambda=0.1050, alpha=0.167.

---

### 6. §2.4 SCOs 8,616→8,067 in prose

**Before:** After filtering, a total of 8,616 SCOs were included.

**After:**  After filtering, a total of 8,067 SCOs were included.

---

### 7. §2.4 8,616→8,067 alignments

**Before:** Trimming resulted in a total of 8,616 alignments which were then concatenated into a supermatrix.

**After:**  Trimming resulted in a total of 8,067 alignments which were then concatenated into a supermatrix.

---

### 8. §3.1 KEY NUMBERS — refreshed pangenome partition (now includes outgroup-only category since tarsalis is a real participant)

**Before:** Total orthogroups: 16,369 | Total genes assigned: 92,540 (97.9%) | Core: 11,653 (71.2%) | Shell: 4,205 (25.7%) | Cloud: 511 (3.1%)

**After:**  Total orthogroups: 16,568 | Total genes assigned: 92,821 (98.2%) | Core: 11,284 (71.1%) | Shell: 3,726 (23.5%) | Cloud: 871 (5.5%) | Outgroup-only: 687

---

### 9. §3.1 KEY NUMBERS — all four form-specific cloud counts refreshed

**Before:** Form-specific cloud: molestus=200, pallens=114, pipiens=142, quinquefasciatus=55

**After:**  Form-specific cloud: molestus=297, pallens=203, pipiens=245, quinquefasciatus=126

---

### 10. §3.1 prose — full refresh: dropped 'Some thing like' stub language, added outgroup compartment, rebalanced percentages

**Before:** Some thing like OrthoFinder identified 16,369 orthogroups comprising 92,540 genes (97.9% assigned). The pangenome consisted of 11,653 core (71.2%), 4,205 shell (25.7%), and 511 cloud (3.1%) orthogroups.

**After:**  OrthoFinder identified 16,568 orthogroups across the four ingroup forms plus Cx. tarsalis as outgroup. Within the ingroup, 92,821 genes (98.2%) were assigned to 15,881 orthogroups, partitioned as 11,284 core (71.1%), 3,726 shell (23.5%), and 871 cloud (5.5%). An additional 687 orthogroups were present only in the Cx. tarsalis outgroup.

---

### 11. Figure 5 caption — pangenome counts refreshed

**Before:** Figure 5. Flower diagram of the Cx. pipiens complex pangenome. Central circle shows 11,653 core orthogroups; petals show form-specific cloud orthogroups. Shell: 4,205.

**After:**  Figure 5. Flower diagram of the Cx. pipiens complex pangenome. Central circle shows 11,284 core orthogroups; petals show form-specific cloud orthogroups. Shell: 3,726.

---

### 12. Figure 2 caption — SCOs, sig families, and corrected concordance support

**Before:** Figure 2. (a) Maximum likelihood species tree from 8,616 concatenated single-copy orthologs. Node labels show bootstrap/gCF/sCF support. Scale bar in substitutions per site. (b) CAFE5 gene family expansions (green) and contractions (red) per lineage. Subtitle indicates 596 significantly evolving families.

**After:**  Figure 2. (a) Maximum likelihood species tree from 8,067 concatenated single-copy orthologs. Node labels show bootstrap/gCF/sCF support (100/58.9/45.1 at the deep ingroup node — substantial gene-tree discordance). Scale bar in substitutions per site. (b) CAFE5 gene family expansions (green) and contractions (red) per lineage. Subtitle indicates 583 significantly evolving families.

---

### 13. §3.3 KEY NUMBERS SCO count

**Before:** SCOs used: 8,616 (after filtering alignments &lt;50 aa)

**After:**  SCOs used: 8,067 (after filtering alignments &lt;50 aa)

---

### 14. §3.3 KEY NUMBERS — sig families and CAFE params

**Before:** CAFE5: 596 significantly evolving families (p&lt;0.05), lambda=0.0945, alpha=0.168

**After:**  CAFE5: 583 significantly evolving families (p&lt;0.05), lambda=0.1050, alpha=0.167

---

### 15. §3.3 per-branch pipiens

**Before:** Cx. pipiens: 3,185 expansions / 1,073 contractions

**After:**  Cx. pipiens: 3,363 expansions / 801 contractions

---

### 16. §3.3 per-branch pallens

**Before:** Cx. pallens: 3,104 exp / 1,376 con

**After:**  Cx. pallens: 3,281 exp / 1,085 con

---

### 17. §3.3 per-branch molestus

**Before:** Cx. molestus: 3,050 exp / 1,593 con

**After:**  Cx. molestus: 3,220 exp / 1,263 con

---

### 18. §3.3 per-branch quinque

**Before:** Cx. quinquefasciatus: 3,671 exp / 866 con

**After:**  Cx. quinquefasciatus: 3,733 exp / 729 con

---

### 19. §3.2 intro — SCO count and concordance framing (does NOT change which taxa are sisters; topology paragraphs still need your editorial pass per audit)

**Before:** Phylogenomic analysis based on 8,616 single-copy orthologues produced a fully resolved species tree with full statistical support at all nodes (Figure 2a). According to the obtained topology,

**After:**  Phylogenomic analysis based on 8,067 single-copy orthologues produced a species tree with full bootstrap support but substantial gene- and site-tree discordance at the deep ingroup node (Figure 2a, see §3.2 below). According to the obtained topology,

---

### 20. §3.3 placeholder — refreshed numbers

**Before:** [Describe CAFE5 results. 596 significantly evolving families. All branches show significantly more expansions than contractions (binomial p&lt;2e-16). Cx. quinquefasciatus shows the strongest expansion bias (3,671 vs 866). Discuss the gamma rate heterogeneity (alpha=0.168 indicates strong rate variation among families). Reference Figure 2b.]

**After:**  [Describe CAFE5 results. 583 significantly evolving families. All branches show significantly more expansions than contractions (binomial p&lt;2e-16). Cx. quinquefasciatus shows the strongest expansion bias (3,733 vs 729). Discuss the gamma rate heterogeneity (alpha=0.167 indicates strong rate variation among families). Reference Figure 2b.]

---

### 21. §3.3 prose — 596→583 sig families

**Before:** Gene family evolutionary changes were identified through CAFE5. A total of 596 gene families showed either expansions or contractions within the Cx. pipiens complex with a p value less than 0.05.

**After:**  Gene family evolutionary changes were identified through CAFE5. A total of 583 gene families showed either expansions or contractions within the Cx. pipiens complex with a p value less than 0.05.

---

### 22. §3.3 prose — CAFE λ and α

**Before:** The rates of genes were determined to be λ = 0.0945, with the low value of α = 0.168 implying the presence of considerable differences in evolutionary rates between gene families.

**After:**  The rates of genes were determined to be λ = 0.1050, with the low value of α = 0.167 implying the presence of considerable differences in evolutionary rates between gene families.

---

### 23. §3.3 prose — per-branch numbers refreshed

**Before:** , where 3,671 gene families expanded, while only 866 contracted. A similar yet slightly less extreme pattern was observed in Cx. pipiens (3,185 vs. 1,073), Cx. pallens (3,104 vs. 1,376), and Cx. molestus (3,050 vs. 1,593).

**After:**  , where 3,733 gene families expanded, while only 729 contracted. A similar yet slightly less extreme pattern was observed in Cx. pipiens (3,363 vs. 801), Cx. pallens (3,281 vs. 1,085), and Cx. molestus (3,220 vs. 1,263).

---

### 24. §3.5 KEY NUMBERS — TE percentages

**Before:** TE content (interspersed): molestus 52.88%, pallens 53.91%, pipiens 52.14%, quinquefasciatus 53.36%

**After:**  TE content (interspersed): molestus 53.46%, pallens 54.43%, pipiens 52.69%, quinquefasciatus 53.86%

---

