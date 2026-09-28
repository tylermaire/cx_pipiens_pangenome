# Culex pipiens complex pangenome

Snakemake workflow and analysis code for a comparative genomic study of the *Culex pipiens* complex. We compare four chromosome scale assemblies (*Cx. quinquefasciatus*, *Cx. pipiens pallens*, and *Cx. pipiens* forms *molestus* and *pipiens*), rooted with *Cx. perexiguus* as the outgroup, and measure gene content, the species tree and its discordance, gene flow between forms, gene family change, genome identity, synteny and repeat content. Release `v3.0` holds the code and results of the version prepared for the *Journal of Medical Entomology*.

Maire T., Maley E., Miller T., Kosinski K. *Phylogenomic relationships and gene content of four forms of the* Culex pipiens *L. (Diptera: Culicidae) complex.* Manuscript in preparation for the *Journal of Medical Entomology*.

An earlier version of the study was posted as a bioRxiv preprint (https://doi.org/10.64898/2026.07.25.740720). Its manuscript, figures and figure scripts are in [`archive/preprint/`](archive/preprint/).

## Main results

- Rooted with *Cx. perexiguus*, the four forms fall into two lineages, (*molestus*, *pipiens*) and (*pallens*, *quinquefasciatus*). The rooted tree comes from the 4,592 five taxon loci whose gene models are intact.
- Gene tree discordance is asymmetric. Gene trees pairing *molestus* with *quinquefasciatus* and *pallens* with *pipiens* make up 55.8% of discordant four taxon gene trees (binomial *P* = 1.0 × 10⁻¹⁰).
- ABBA BABA tests show *molestus* sharing excess derived alleles with *quinquefasciatus* (D = −0.140), and *pallens* with *pipiens* (D = −0.233) and *molestus* (D = −0.137). This fits gene flow between the two lineages and *Cx. pipiens* ancestry in *pallens*.
- Transferred gene models that share an error create false shared derived alleles. On all loci they move the root to the *quinquefasciatus* branch, which is why the rooted analyses use loci with intact models.
- Most apparent differences in gene content are annotation transfer artifacts. Of the absences implied by cloud orthogroups, 91.4% are the same gene clustered elsewhere or a locus without a usable model, and the families CAFE calls significant track copy loss during transfer. An anvi'o pangenome of the same proteins agrees with OrthoFinder on the core.
- The forms are close at the protein level (median identity 98.5 to 99.1%), and synteny is conserved (83.1 to 88.1% of the genes they share on homologous chromosomes are collinear).

`results/manuscript_values.tsv` lists every number in the paper with the file it comes from.

## Genomes

| Form | Accession | Size (Mb) | Gene models |
|------|-----------|-----------|-------------|
| *Cx. quinquefasciatus* (annotation reference) | GCF_015732765.1 | 573 | RefSeq |
| *Cx. pipiens pallens* | GCF_016801865.2 | 566 | Liftoff from the reference |
| *Cx. pipiens* f. *molestus* | GCA_024516115.1 | 560 | Liftoff from the reference |
| *Cx. pipiens* f. *pipiens* | GCA_963924435.1 | 533 | Liftoff from the reference |
| *Cx. perexiguus* (outgroup) | GCA_964243045.1 (idCulPerx1.1) | 436 | its own Ensembl gene set |

*Cx. perexiguus* belongs to the Univittatus Subgroup of the Pipiens Group, outside the *Cx. pipiens* complex. Its chromosome level assembly (PacBio HiFi and Hi-C) and gene set are taken from the Ensembl FTP site so that the sequence names of the genome and the annotation match (`download.smk`, rule `download_ensembl`). Earlier versions used *Cx. tarsalis* CtarK1, a contig level assembly whose transferred gene models were mostly broken.

## Versions

| Version | Where | Outgroup | Status |
|---------|-------|----------|--------|
| V5 run, release `v3.0` | tag `v3.0`; branches `main` and `v5-perexiguus` | *Cx. perexiguus* | current; the JME manuscript |
| V4 run | tag `v2.1`; branch `v4-review-fixes` | *Cx. tarsalis* | superseded |
| earlier clean reruns | tags `v2.0-rerun`, `v1.0-rerun` | *Cx. tarsalis* | superseded |
| preprint | `archive/preprint/` | *Cx. tarsalis* | bioRxiv, July 2026 |

`docs/v5_run.md` records how the V5 run was done, the problems met and the checks made outside the workflow; `docs/v4_review_fixes.md` does the same for V4.

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
| `gene_families.smk` | gene family change on the rooted species tree (family P values computed exactly under the fitted model), and the annotation transfer control | CAFE5 |
| `synteny.smk` | whole genome identity and gene anchor synteny | skani, minimap2 |
| `repeats.smk` | repeat library and masking | RepeatModeler2, RepeatMasker |
| `functional.smk` | key gene families | eggNOG mapper |
| `report.smk` | every value the manuscript cites, with its source file | |
| `figures.smk` | Figures 1 to 4, Tables 1 to 4 and Supplementary Tables S1 to S13 (Figure 5 and S14 from the anvi'o tables, when present) | matplotlib, openpyxl |

The outgroup analyses use orthogroups with one gene in each of the five taxa. Proteins are aligned with MAFFT and trimmed with trimAl, and coding sequences are placed codon by codon on the trimmed protein alignments. The concatenated protein alignment of the loci whose four ingroup gene models are intact gives the species tree rooted with *Cx. perexiguus*, and ABBA BABA tests (D statistics, block jackknife over windows of the reference assembly) ask whether one form shares more derived alleles with another than incomplete lineage sorting allows. The planned tests are set in `config/config.yaml`.

The anvi'o gene cluster pangenome runs by hand, outside the workflow, with `workflow/scripts/anvio_pangenome.sh` (anvi'o 9, on the same proteins as OrthoFinder). `workflow/scripts/anvio_compare.py`, called at its end, compares the gene clusters with the OrthoFinder partition and writes the tables in `results/anvio` that Figure 5 and Supplementary Table S14 are drawn from. See `docs/v5_run.md` for the tools and versions.

## Running it

```bash
git clone https://github.com/tylermaire/cx_pipiens_pangenome.git
cd cx_pipiens_pangenome
git checkout v3.0                              # the release cited in the paper
conda install -c bioconda -c conda-forge snakemake

snakemake -n --use-conda --cores 16          # dry run
snakemake --use-conda --cores 16             # full run
```

`docs/v5_run.md` gives the steps used for the V5 run on AWS, including which outputs to bring back. The RepeatModeler2 library depends only on the *Cx. quinquefasciatus* assembly, so V5 reuses the V4 library in `data/repeat_library/` (set `repeats: prebuilt_library: ""` in `config/config.yaml` to build it again, about 18 hours on 32 cores). Genomes are downloaded into `resources/` and are not tracked. Large outputs (transferred annotations, alignments, gene trees, CAFE output and the anvi'o databases) are not tracked either; they will be archived on Zenodo with the release.

## Figures and tables

`workflow/scripts/make_revision_figures.py`, `make_manuscript_tables.py` and `make_supplement.py` build the figures, tables and supplementary workbook from `results/`. The workflow runs them at the end (`figures.smk`), and they can also be run by hand from the repository root. Figures go to `figures/revision/` as vector PDF, 1,200 dpi TIFF and a PNG preview, and are not tracked.

| Figure | What it shows |
|--------|---------------|
| 1 | pangenome partition, transferred models without a valid ORF, and tests of the cloud |
| 2 | species tree rooted with *Cx. perexiguus*, gene tree discordance and ABBA BABA tests |
| 3 | copy retention after annotation transfer and CAFE significance by reference copy number |
| 4 | pairwise synteny of the four ingroup assemblies |
| 5 | anvi'o gene clusters compared with the OrthoFinder partition |

The supplementary workbook, `tables/Supplementary_Tables.xlsx` (Supp. Tables S1 to S14), is tracked; each sheet names the result files it was built from.

## Layout

```
config/      config.yaml, samples.tsv
data/        the RepeatModeler2 library reused by V5
workflow/    rules/, envs/, scripts/
results/     summary tables tracked here; large intermediates are gitignored
tables/      Supplementary_Tables.xlsx (Tables 1 to 4 are built here too)
tests/       python3 tests/test_<name>.py
docs/        notes on each version of the analysis
archive/     files of the preprint version, kept for the record
```

## Tool citations

Liftoff 1.6.3 (Shumate & Salzberg 2021); minimap2 (Li 2018); miniprot (Li 2023); gffread (Pertea & Pertea 2020); QUAST 5.3.0 (Gurevich et al. 2013); BUSCO 6.1.0 (Tegenfeldt et al. 2025); OrthoFinder 3.1.5 (Emms et al. 2026); DIAMOND (Buchfink et al. 2021); MAFFT (Katoh & Standley 2013); trimAl (Capella-Gutiérrez et al. 2009); IQ-TREE 3.1.3 (Wong et al. 2026) with ModelFinder (Kalyaanamoorthy et al. 2017), UFBoot (Hoang et al. 2018) and concordance factors (Minh et al. 2020); ABBA BABA tests (Green et al. 2010; Durand et al. 2011) with a weighted block jackknife (Busing et al. 1999); CAFE 5 (Mendes et al. 2020); skani (Shaw & Yu 2023); RepeatModeler 2.0.7 (Flynn et al. 2020); RepeatMasker 4.1.7 (Smit et al. 2013–2015); eggNOG mapper 2.1.12 (Cantalapiedra et al. 2021) with eggNOG 5.0 (Huerta-Cepas et al. 2019); anvi'o 9 (Eren et al. 2021; Delmont & Eren 2018) with MCL (van Dongen & Abreu-Goodger 2012), FAMSA (Deorowicz et al. 2016) and the minbit criterion (Benedict et al. 2014); Snakemake (Mölder et al. 2021).

## License

MIT, see [LICENSE](LICENSE).

## Contact

Tyler Maire, Indian River Mosquito Control District, [ORCID 0009-0007-6451-5408](https://orcid.org/0009-0007-6451-5408)
