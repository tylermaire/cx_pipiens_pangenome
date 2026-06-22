@echo off
REM ============================================================
REM  Commit + push all manuscript/figure/ANI updates to GitHub
REM  Run this from the project folder. Git will prompt for your
REM  GitHub credentials on push if needed.
REM ============================================================
cd /d "%~dp0"

echo.
echo === Current branch and status ===
git branch --show-current
git status --short

echo.
echo === Staging all changes (Word temp files are gitignored) ===
git add -A

echo.
echo === Files staged for commit ===
git status --short

echo.
echo === Creating commit ===
git commit -m "Manuscript draft 5: add Cx. tarsalis ANI finding, 4th author, figure 1 italics fix, editorial revisions; add PDF" -m "- Regenerated Figure 1 with italic diptera_odb10 axis label (make_all_figures.R fix)" -m "- Completed skani ANI run incl. Cx. tarsalis (below detection floor); matrix saved" -m "- Added Kyle Kosinski as 4th author (IRMCD) + author-contributions entry" -m "- Fixed artifacts, number inconsistencies, captions; standardized US spelling" -m "- Added manuscript PDF and tracked/clean docx drafts; added editorial review"

echo.
echo === Pushing to origin/main (you may be prompted to log in) ===
git push origin main

echo.
echo === Done. Final status: ===
git status --short
echo.
echo If the push succeeded, your commit is now on GitHub.
pause
