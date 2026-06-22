# Culex pipiens complex pangenome

Snakemake pipeline and analysis code for a comparative pangenome study of the *Culex pipiens* species complex. We compare four chromosome-scale assemblies — *Cx. pipiens*, *Cx. p. pallens*, *Cx. p. molestus*, and *Cx. quinquefasciatus* — using *Cx. tarsalis* as an outgroup, and quantify gene-content variation, the species tree, gene-family turnover, whole-genome similarity, synteny, and repeat content.

Maire T., Maley E., Miller T., Kosinski K. *Pangenome analysis reveals gene content variation and evolutionary dynamics across the* Culex pipiens *species complex.* Manuscript (`Culex pipiens Pangenome Analysis.pdf`).

## Genomes

| Form | Accession | Size (Mb) |
|------|-----------|-----------|
| *Cx. quinquefasciatus* (annotation reference) | GCF_015732765.1 | 573 |
| *Cx. p. pallens* | GCF_016801865.2 | 566 |
| *Cx. p. molestus* | GCA_024516115.1 | 560 |
| *Cx. p. pipiens* | GCA_963924435.1 | 533 |
| *Cx. tarsalis* (outgroup) | CtarK1, OSF [osf.io/mdwqx](https://osf.io/mdwqx) | 790 |

The original NCBI record for *Cx. tarsalis* (`GCA_016859205.1`) was reassigned to an unrelated organism, so the pipeline pulls the CtarK1 assembly (Main et al. 2021) from OSF instead of NCBI. The download rule handles this case automatically.

## How it works

Gene models from the *Cx. quinquefasciatus* RefSeq reference are lifted onto every other assembly with Liftoff, which removes annotation-source differences between genomes. From the resulting protein sets the pipeline runs OrthoFinder for orthology, partitions orthogroups into core/shell/cloud, builds a single-copy-orthologue species tree in IQ-TREE 2 with gene- and site-concordance factors, models gene-family change with CAFE5, estimates whole-genome ANI with skani, derives gene-anchor synteny from the Liftoff coordinates, and annotates repeats with RepeatModeler2/RepeatMasker.

| Rule file | Main tools |
|-----------|-----------|
| `annotation.smk` | Liftoff, gffread |
| `orthology.smk` | OrthoFinder, DIAMOND |
| `phylogenomics.smk` | MAFFT, trimAl, IQ-TREE 2 |
| `gene_families.smk` | CAFE5 |
| `synteny.smk` | skani, minimap2 |
| `repeats.smk` | RepeatModeler2, RepeatMasker |

Two analyses sit outside the main workflow. The anvi'o gene-cluster pangenome (Figure 8) was run by hand; the commands are in the Methods. `synteny.smk` still contains an older nucmer/SyRI path that was not used for the reported figures — ANI and synteny in the paper come from skani and the Liftoff projections.

A note on the outgroup ANI: *Cx. tarsalis* sits in subgenus *Culex* outside the *Cx. pipiens* complex, and it is too divergent from the ingroup for skani to score — their genome-wide identity is below the ~80% floor at which skani reports ANI, so it returns no value for any *Cx. tarsalis*–ingroup pair. That is the expected result for an outgroup at this distance, not a quirk of the data. The reported matrix is therefore ingroup-only.

## Running it

```bash
git clone https://github.com/tylermaire/cx_pipiens_pangenome.git
cd cx_pipiens_pangenome
conda install -c bioconda -c conda-forge snakemake

snakemake -n --use-conda --cores 16          # dry run
snakemake --use-conda --cores 16             # full run
```

Genomes are downloaded on demand into `resources/` and are not tracked here. Expect a long run: RepeatModeler dominates the wall time. It was developed on a 16-core / 128 GB machine, where an end-to-end run takes roughly a day and a half.

## Figures

Figures are built from the files in `results/` by `workflow/scripts/make_all_figures.R` (ggplot2, patchwork, ggtree); the anvi'o figure is made separately. Outputs land in `figures/` as PNG and PDF.

| | |
|--|--|
| Fig 1 | Assembly quality (BUSCO, QUAST) |
| Fig 2 | Pangenome composition |
| Fig 3 | Species tree with concordance factors |
| Fig 4 | CAFE5 gene-family expansions and contractions |
| Fig 5 | Whole-genome ANI heatmap |
| Fig 6 | Pairwise synteny dotplots |
| Fig 7 | Transposable-element composition |
| Fig 8 | Anvi'o gene-cluster pangenome |
| Fig 9 | Copy-number heatmap for detox/chemosensory/immunity families |
| Fig S1–S2 | Per-pair synteny dotplots; pipeline flowchart |

## Layout

```
config/      config.yaml, samples.tsv
workflow/    rules/, envs/, scripts/
results/     summary tables tracked here; large intermediates are gitignored
figures/     manuscript figures (+ Sup Material/ for S1, S2, supplementary tables)
```

## Tool citations

OrthoFinder (Emms & Kelly 2019); DIAMOND (Buchfink et al. 2015); MAFFT (Katoh & Standley 2013); trimAl (Capella-Gutiérrez et al. 2009); IQ-TREE 2 (Minh et al. 2020); CAFE5 (Mendes et al. 2020); skani (Shaw & Yu 2023); minimap2 (Li 2018); Liftoff (Shumate & Salzberg 2021); RepeatModeler2/RepeatMasker (Flynn et al. 2020); anvi'o (Eren et al. 2021); Snakemake (Mölder et al. 2021).

## License

MIT — see [LICENSE](LICENSE).

## Contact

Tyler Maire, Indian River Mosquito Control District — [ORCID 0009-0007-6451-5408](https://orcid.org/0009-0007-6451-5408)
