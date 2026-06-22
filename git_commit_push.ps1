# PowerShell script: commit + push manuscript/figure/ANI updates to GitHub.
# Run from the project folder:  .\git_commit_push.ps1
# Git will prompt for GitHub credentials on push if needed.

Set-Location -Path $PSScriptRoot

Write-Host "`n=== Current branch and status ===" -ForegroundColor Cyan
git branch --show-current
git status --short

Write-Host "`n=== Staging all changes (Word temp files are gitignored) ===" -ForegroundColor Cyan
git add -A
git status --short

Write-Host "`n=== Creating commit ===" -ForegroundColor Cyan
git commit -m "Manuscript draft 5: add Cx. tarsalis ANI finding, 4th author, Figure 1 italics fix, editorial revisions; add PDF" -m "- Regenerated Figure 1 with italic diptera_odb10 axis label (make_all_figures.R fix)" -m "- Completed skani ANI run incl. Cx. tarsalis (below detection floor); matrix saved" -m "- Added Kyle Kosinski as 4th author (IRMCD) + author-contributions entry" -m "- Fixed artifacts, number inconsistencies, captions; standardized US spelling" -m "- Added manuscript PDF and tracked/clean docx drafts; added editorial review"

Write-Host "`n=== Pushing to origin/main (you may be prompted to log in) ===" -ForegroundColor Cyan
git push origin main

Write-Host "`n=== Done. Final status: ===" -ForegroundColor Green
git status --short
