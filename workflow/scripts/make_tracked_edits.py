#!/usr/bin/env python3
"""
make_tracked_edits.py — Apply numerical corrections to the paper docx.

Strategy: rather than attempt fragile XML-level tracked-changes wrapping
(which broke on multi-w:t runs), we apply each edit as a direct string
substitution inside <w:t> elements only. We also write a sidecar
changes_log.md that documents every (old, new) pair so the user has
explicit before/after transparency for review.

The original Paper_original.docx is preserved untouched. The corrected
version is written to Paper_corrected.docx via the docx skill's pack tool.
"""
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ORIG = ROOT / "Paper_original.docx"
WORKDIR = ROOT / "unpacked_clean"
OUT = ROOT / "Paper_corrected.docx"
LOG = ROOT / "changes_log.md"

UNPACK = Path("/sessions/sweet-adoring-wright/mnt/.claude/skills/docx/scripts/office/unpack.py")
PACK = Path("/sessions/sweet-adoring-wright/mnt/.claude/skills/docx/scripts/office/pack.py")

# (old_text, new_text, description)
# Updated 2026-05-24 with values from the clean reproducibility rerun
# that includes the corrected Cx_tarsalis assembly as a real outgroup.
EDITS = [
    # ---- Section 2.1 — genome size in parentheses
    (" (GCF_015732765.1; 578 Mb), ",
     " (GCF_015732765.1; 573 Mb), ",
     "§2.1 quinquefasciatus genome size 578→573 Mb (QUAST exact total)"),
    (" (GCF_016801865.2; 568 Mb), ",
     " (GCF_016801865.2; 566 Mb), ",
     "§2.1 pallens genome size 568→566 Mb"),

    # ---- Section 2.2 — pallens protein count
    (", 22,585 for ", ", 23,089 for ",
     "§2.2 pallens protein count 22,585→23,089 (new annotation lift)"),

    # ---- Section 2.4 highlighted (Erin) block — SCOs and CAFE parameters
    ("We pulled 8,616 single copy orthologs",
     "We pulled 8,067 single copy orthologs",
     "§2.4 SCOs 8,616→8,067 in Erin block"),
    ("Then CAFE5 v5.1.0 (gamma model, k=3) tested which gene families are expanding or contracting. Lambda=0.0945, alpha=0.168.",
     "Then CAFE5 v5.1.0 (gamma model, k=3) tested which gene families are expanding or contracting. Lambda=0.1050, alpha=0.167.",
     "§2.4 CAFE parameters in Erin block"),

    # ---- Section 2.4 prose body — SCO counts
    ("After filtering, a total of 8,616 SCOs were included.",
     "After filtering, a total of 8,067 SCOs were included.",
     "§2.4 SCOs 8,616→8,067 in prose"),
    ("Trimming resulted in a total of 8,616 alignments which were then concatenated into a supermatrix.",
     "Trimming resulted in a total of 8,067 alignments which were then concatenated into a supermatrix.",
     "§2.4 8,616→8,067 alignments"),

    # ---- Section 3.1 KEY NUMBERS and Section 3.1 prose summary
    ("Total orthogroups: 16,369 | Total genes assigned: 92,540 (97.9%) | Core: 11,653 (71.2%) | Shell: 4,205 (25.7%) | Cloud: 511 (3.1%)",
     "Total orthogroups: 16,568 | Total genes assigned: 92,821 (98.2%) | Core: 11,284 (71.1%) | Shell: 3,726 (23.5%) | Cloud: 871 (5.5%) | Outgroup-only: 687",
     "§3.1 KEY NUMBERS — refreshed pangenome partition (now includes outgroup-only category since tarsalis is a real participant)"),
    ("Form-specific cloud: molestus=200, pallens=114, pipiens=142, quinquefasciatus=55",
     "Form-specific cloud: molestus=297, pallens=203, pipiens=245, quinquefasciatus=126",
     "§3.1 KEY NUMBERS — all four form-specific cloud counts refreshed"),
    ("Some thing like OrthoFinder identified 16,369 orthogroups comprising 92,540 genes (97.9% assigned). The pangenome consisted of 11,653 core (71.2%), 4,205 shell (25.7%), and 511 cloud (3.1%) orthogroups.",
     "OrthoFinder identified 16,568 orthogroups across the four ingroup forms plus Cx. tarsalis as outgroup. Within the ingroup, 92,821 genes (98.2%) were assigned to 15,881 orthogroups, partitioned as 11,284 core (71.1%), 3,726 shell (23.5%), and 871 cloud (5.5%). An additional 687 orthogroups were present only in the Cx. tarsalis outgroup.",
     "§3.1 prose — full refresh: dropped 'Some thing like' stub language, added outgroup compartment, rebalanced percentages"),

    # ---- Figure captions
    ("Figure 5. Flower diagram of the Cx. pipiens complex pangenome. Central circle shows 11,653 core orthogroups; petals show form-specific cloud orthogroups. Shell: 4,205.",
     "Figure 5. Flower diagram of the Cx. pipiens complex pangenome. Central circle shows 11,284 core orthogroups; petals show form-specific cloud orthogroups. Shell: 3,726.",
     "Figure 5 caption — pangenome counts refreshed"),
    ("Figure 2. (a) Maximum likelihood species tree from 8,616 concatenated single-copy orthologs. Node labels show bootstrap/gCF/sCF support. Scale bar in substitutions per site. (b) CAFE5 gene family expansions (green) and contractions (red) per lineage. Subtitle indicates 596 significantly evolving families.",
     "Figure 2. (a) Maximum likelihood species tree from 8,067 concatenated single-copy orthologs. Node labels show bootstrap/gCF/sCF support (100/58.9/45.1 at the deep ingroup node — substantial gene-tree discordance). Scale bar in substitutions per site. (b) CAFE5 gene family expansions (green) and contractions (red) per lineage. Subtitle indicates 583 significantly evolving families.",
     "Figure 2 caption — SCOs, sig families, and corrected concordance support"),

    # ---- Section 3.3 KEY NUMBERS and per-branch CAFE numbers
    ("SCOs used: 8,616 (after filtering alignments &lt;50 aa)",
     "SCOs used: 8,067 (after filtering alignments &lt;50 aa)",
     "§3.3 KEY NUMBERS SCO count"),
    ("CAFE5: 596 significantly evolving families (p&lt;0.05), lambda=0.0945, alpha=0.168",
     "CAFE5: 583 significantly evolving families (p&lt;0.05), lambda=0.1050, alpha=0.167",
     "§3.3 KEY NUMBERS — sig families and CAFE params"),
    ("Cx. pipiens: 3,185 expansions / 1,073 contractions",
     "Cx. pipiens: 3,363 expansions / 801 contractions",
     "§3.3 per-branch pipiens"),
    ("Cx. pallens: 3,104 exp / 1,376 con",
     "Cx. pallens: 3,281 exp / 1,085 con",
     "§3.3 per-branch pallens"),
    ("Cx. molestus: 3,050 exp / 1,593 con",
     "Cx. molestus: 3,220 exp / 1,263 con",
     "§3.3 per-branch molestus"),
    ("Cx. quinquefasciatus: 3,671 exp / 866 con",
     "Cx. quinquefasciatus: 3,733 exp / 729 con",
     "§3.3 per-branch quinque"),

    # ---- Section 3.2 prose intro (SCO + topology framing)
    ("Phylogenomic analysis based on 8,616 single-copy orthologues produced a fully resolved species tree with full statistical support at all nodes (Figure 2a). According to the obtained topology,",
     "Phylogenomic analysis based on 8,067 single-copy orthologues produced a species tree with full bootstrap support but substantial gene- and site-tree discordance at the deep ingroup node (Figure 2a, see §3.2 below). According to the obtained topology,",
     "§3.2 intro — SCO count and concordance framing (does NOT change which taxa are sisters; topology paragraphs still need your editorial pass per audit)"),

    # ---- Section 3.3 placeholder block and prose
    ("[Describe CAFE5 results. 596 significantly evolving families. All branches show significantly more expansions than contractions (binomial p&lt;2e-16). Cx. quinquefasciatus shows the strongest expansion bias (3,671 vs 866). Discuss the gamma rate heterogeneity (alpha=0.168 indicates strong rate variation among families). Reference Figure 2b.]",
     "[Describe CAFE5 results. 583 significantly evolving families. All branches show significantly more expansions than contractions (binomial p&lt;2e-16). Cx. quinquefasciatus shows the strongest expansion bias (3,733 vs 729). Discuss the gamma rate heterogeneity (alpha=0.167 indicates strong rate variation among families). Reference Figure 2b.]",
     "§3.3 placeholder — refreshed numbers"),
    ("Gene family evolutionary changes were identified through CAFE5. A total of 596 gene families showed either expansions or contractions within the Cx. pipiens complex with a p value less than 0.05.",
     "Gene family evolutionary changes were identified through CAFE5. A total of 583 gene families showed either expansions or contractions within the Cx. pipiens complex with a p value less than 0.05.",
     "§3.3 prose — 596→583 sig families"),
    ("The rates of genes were determined to be λ = 0.0945, with the low value of α = 0.168 implying the presence of considerable differences in evolutionary rates between gene families.",
     "The rates of genes were determined to be λ = 0.1050, with the low value of α = 0.167 implying the presence of considerable differences in evolutionary rates between gene families.",
     "§3.3 prose — CAFE λ and α"),
    (", where 3,671 gene families expanded, while only 866 contracted. A similar yet slightly less extreme pattern was observed in Cx. pipiens (3,185 vs. 1,073), Cx. pallens (3,104 vs. 1,376), and Cx. molestus (3,050 vs. 1,593).",
     ", where 3,733 gene families expanded, while only 729 contracted. A similar yet slightly less extreme pattern was observed in Cx. pipiens (3,363 vs. 801), Cx. pallens (3,281 vs. 1,085), and Cx. molestus (3,220 vs. 1,263).",
     "§3.3 prose — per-branch numbers refreshed"),

    # ---- Section 3.5 KEY NUMBERS — TE percentages
    ("TE content (interspersed): molestus 52.88%, pallens 53.91%, pipiens 52.14%, quinquefasciatus 53.36%",
     "TE content (interspersed): molestus 53.46%, pallens 54.43%, pipiens 52.69%, quinquefasciatus 53.86%",
     "§3.5 KEY NUMBERS — TE percentages"),
]


def main():
    if WORKDIR.exists():
        # Cannot rm -rf existing unpacked dir in this environment (read-only flags),
        # so unpack to a fresh sibling. Append a suffix until we find an unused one.
        suffix = 0
        while True:
            cand = ROOT / f"unpacked_clean_{suffix}"
            if not cand.exists():
                break
            suffix += 1
        workdir = cand
    else:
        workdir = WORKDIR
    print(f"Unpacking to {workdir}…")
    r = subprocess.run([sys.executable, str(UNPACK), str(ORIG), str(workdir)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout); print(r.stderr); sys.exit(1)
    print(r.stdout.strip())

    doc = workdir / "word" / "document.xml"
    xml = doc.read_text()

    applied = []
    missing = []
    for old, new, desc in EDITS:
        if old in xml:
            xml = xml.replace(old, new, 1)
            applied.append((old, new, desc))
        else:
            missing.append((old, new, desc))

    doc.write_text(xml)
    print(f"\nApplied {len(applied)} / {len(EDITS)} edits.")
    if missing:
        print(f"WARNING: {len(missing)} edits did not match — see changes_log.md.")

    # Write changes log
    with open(LOG, "w") as f:
        f.write("# Numerical corrections applied to Paper_corrected.docx\n\n")
        f.write("Each row below shows a single text replacement applied directly\n")
        f.write("to the manuscript. The original `Paper_original.docx` is preserved\n")
        f.write("for diff purposes. Open both files side-by-side or use Word's\n")
        f.write("Compare → Combine feature to generate a tracked-changes view.\n\n")
        f.write(f"**Edits applied: {len(applied)} of {len(EDITS)}**\n\n")
        f.write("---\n\n## Applied edits\n\n")
        for i, (old, new, desc) in enumerate(applied, start=1):
            f.write(f"### {i}. {desc}\n\n")
            f.write(f"**Before:** {old}\n\n")
            f.write(f"**After:**  {new}\n\n---\n\n")
        if missing:
            f.write("## Edits that could not be matched\n\n")
            f.write("These exact strings were not found in the document. They may\n")
            f.write("have been edited manually since the audit, or the surrounding\n")
            f.write("text contains characters that differ. Apply manually.\n\n")
            for i, (old, new, desc) in enumerate(missing, start=1):
                f.write(f"### {i}. {desc}\n\n")
                f.write(f"**Before:** {old}\n\n")
                f.write(f"**After:**  {new}\n\n---\n\n")
    print(f"Wrote changes log: {LOG}")

    # Repack
    print(f"\nPacking to {OUT}…")
    r = subprocess.run(
        [sys.executable, str(PACK), str(workdir), str(OUT),
         "--original", str(ORIG)],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        print(r.stdout); print(r.stderr); sys.exit(1)
    print(r.stdout.strip())


if __name__ == "__main__":
    main()
