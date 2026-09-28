# Preprint files (superseded)

These files belong to the first version of this study, posted as a bioRxiv
preprint in July 2026 (https://doi.org/10.64898/2026.07.25.740720). They are
kept for the record and are not part of the current workflow.

That version used *Cx. tarsalis* (CtarK1) as the outgroup. The current version
(release v3.0, prepared for the Journal of Medical Entomology) uses
*Cx. perexiguus* instead, and several of the preprint's results changed or were
withdrawn: the cloud orthogroups that the preprint treated as form specific
genes turned out to be mostly annotation transfer artifacts, and the CAFE
results tracked copy loss during transfer. See the main README and
`docs/v5_run.md` for the current analyses.

| Path | What it is |
|------|------------|
| `Culex pipiens Pangenome Analysis.pdf` | the preprint manuscript |
| `figures/` | preprint Figures 1 to 9, Figures S1 and S2, and the preprint supplementary tables (`Sup Material/`) |
| `workflow_scripts/` | the R and Python scripts that drew those figures (formerly in `workflow/scripts/`) |
| `python_figures/` | a later Python figure suite for a revision that was not used |

The current figures are drawn by `workflow/scripts/make_revision_figures.py`
into `figures/revision/`. The scripts here were written to run from the
repository root with their original paths, so they will not run from this
folder without editing.
