# Culex pipiens complex pangenome

Snakemake pipeline and analysis code for a comparative genomic study of the *Culex pipiens* complex. We compare four chromosome scale assemblies (*Cx. quinquefasciatus*, *Cx. pipiens pallens*, and *Cx. pipiens* forms *molestus* and *pipiens*) with *Cx. perexiguus* as the outgroup, and quantify gene content variation, the species tree and its discordance, gene flow between forms, gene family turnover, whole genome similarity, synteny and repeat content.

Maire T., Maley E., Miller T., Kosinski K. *Phylogenomic relationships and gene content of four forms of the* Culex pipiens *L. (Diptera: Culicidae) complex.* Manuscript.

## Genomes

| Form | Accession | Size (Mb) | Gene models |
|------|-----------|-----------|-------------|
| *Cx. quinquefasciatus* (annotation reference) | GCF_015732765.1 | 573 | RefSeq |
| *Cx. pipiens pallens* | GCF_016801865.2 | 566 | Liftoff from the reference |
| *Cx. pipiens* f. *molestus* | GCA_024516115.1 | 560 | Liftoff from the reference |
| *Cx. pipiens* f. *pipiens* | GCA_963924435.1 | 533 | Liftoff from the reference |
| *Cx. perexiguus* (outgroup) | GCA_964243045.1 (idCulPerx1.1) | 436 | its own Ensembl gene set |

*Cx. perexiguus* belongs to the Univittatus Subgroup of the Pipiens Group, outside the *Cx. pipiens* complex. Its chromosome level assembly (PacBio HiFi and Hi-C) and gene set are taken from the Ensembl FTP site so that the sequence names of the genome and the annotation match (`download.smk`, rule `download_ensembl`). Runs up to V4 used *Cx. tarsalis* CtarK1, a contig level assembly whose transferred gene models were mostly broken; that version is on branch `v4-review-fixes` (tag `v2.1` for the V4 run).

## How it works

Gene models from the *Cx. quinquefasciatus* RefSeq annotation are transferred to the other three ingroup assemblies with Liftoff, so the ingroup gene sets share one annotation source; the outgroup keeps its own annotation. One protein per gene goes to OrthoFinder, and orthogroups are partitioned into core, shell and cloud over the four ingroup forms.

| Rule file | What it does | Main tools |
|-----------|--------------|------------|
| `download.smk` | assemblies and annotations (NCBI Datasets; Ensembl FTP for the outgroup) | datasets, curl |
| `qc.smk` | assembly statistics and completeness | QUAST, BUSCO |
| `annotation.smk` | annotation transfer, proteins and coding sequences, transfer quality | Liftoff, gffread |
| `orthology.smk` | orthogroups, pangenome partition, cloud composition | OrthoFinder, DIAMOND |
| `validation.smk` | tests of every absence implied by the partition | miniprot, minimap2 |
| `phylogenomics.smk` | four taxon species tree, concordance and quartet tests; five taxon tree rooted with the outgroup; D statistics | MAFFT, trimAl, IQ-TREE |
| `gene_families.smk` | gene family change on the rooted species tree, and the annotation transfer control | CAFE5 |
| `synteny.smk` | whole genome identity and gene anchor synteny | skani, minimap2 |
| `repeats.smk` | repeat library and masking | RepeatModeler2, RepeatMasker |
| `functional.smk` | key gene families | eggNOG mapper |
| `report.smk` | every value the manuscript cites, with its source file | |
| `figures.smk` | Figures 1 to 4, Tables 1 to 4 and Supplementary Tables S1 to S12 | matplotlib, openpyxl |

The outgroup analyses (V5) use orthogroups with one gene in each of the five taxa: proteins are aligned with MAFFT and trimmed with trimAl, coding sequences are placed codon by codon on the trimmed protein alignments, the concatenated protein alignment gives a species tree rooted with *Cx. perexiguus*, and ABBA BABA tests (D statistics, block jackknife over windows of the reference assembly) ask whether one form shares more derived alleles with another than incomplete lineage sorting allows. The planned tests are set in `config/config.yaml`.

The anvi'o gene cluster pangenome is run by hand, outside the workflow.

## Running it

```bash
git clone https://github.com/tylermaire/cx_pipiens_pangenome.git
cd cx_pipiens_pangenome
conda install -c bioconda -c conda-forge snakemake

snakemake -n --use-conda --cores 16          # dry run
snakemake --use-conda --cores 16             # full run
```

`docs/v5_run.md` gives the steps used for the V5 run on AWS, including which outputs to bring back. The RepeatModeler2 library depends only on the *Cx. quinquefasciatus* assembly, so V5 reuses the V4 library in `data/repeat_library/` (set `repeats: prebuilt_library: ""` in `config/config.yaml` to build it again, about 18 hours on 32 cores). Genomes are downloaded into `resources/` and are not tracked.

## Figures and tables

`workflow/scripts/make_revision_figures.py`, `make_manuscript_tables.py` and `make_supplement.py` build the figures, tables and supplementary workbook from `results/`; the workflow runs them at the end (`figures.smk`), and they can be run by hand from the repository root. Outputs go to `figures/revision/` and `tables/`.

## Layout

```
config/      config.yaml, samples.tsv
data/        the RepeatModeler2 library reused by V5
workflow/    rules/, envs/, scripts/
results/     summary tables tracked here; large intermediates are gitignored
tests/       python3 tests/test_<name>.py
docs/        notes on each version of the analysis
```

## Tool citations

OrthoFinder (Emms & Kelly 2019); DIAMOND (Buchfink et al. 2015); MAFFT (Katoh & Standley 2013); trimAl (Capella-Gutiérrez et al. 2009); IQ-TREE 2 (Minh et al. 2020); CAFE5 (Mendes et al. 2020); skani (Shaw & Yu 2023); minimap2 (Li 2018); miniprot (Li 2023); Liftoff (Shumate & Salzberg 2021); RepeatModeler2/RepeatMasker (Flynn et al. 2020); BUSCO (Manni et al. 2021); eggNOG mapper (Cantalapiedra et al. 2021); anvi'o (Eren et al. 2021); Snakemake (Mölder et al. 2021).

## License

MIT, see [LICENSE](LICENSE).

## Contact

Tyler Maire, Indian River Mosquito Control District, [ORCID 0009-0007-6451-5408](https://orcid.org/0009-0007-6451-5408)
