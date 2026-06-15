# Culex pipiens Complex Pangenome Analysis

A reproducible Snakemake pipeline for the first pangenome analysis of the *Culex pipiens* species complex, comparing four chromosome-scale genome assemblies with *Cx. tarsalis* as an outgroup.

**Manuscript:** Maire T., Maley E., Miller T. Pangenome analysis reveals gene content variation and evolutionary dynamics across the *Culex pipiens* species complex. *In preparation* for *Insect Biochemistry and Molecular Biology*.

## Genome Assemblies

| Form | Accession | Size (Mb) | Level |
|------|-----------|-----------|-------|
| *Cx. quinquefasciatus* (reference) | GCF_015732765.1 | 573 | Chromosome |
| *Cx. pipiens* f. *pallens* | GCF_016801865.2 | 566 | Chromosome |
| *Cx. pipiens* f. *molestus* | GCA_024516115.1 | 560 | Chromosome |
| *Cx. pipiens* f. *pipiens* | GCA_963924435.1 | 533 | Chromosome |
| *Cx. tarsalis* (outgroup) | CtarK1 (OSF: osf.io/mdwqx) | 790 | Chromosome |

> **Note on the outgroup:** NCBI accession `GCA_016859205.1` has been reassigned to an unrelated organism, so the *Cx. tarsalis* assembly (CtarK1; Main et al.) is retrieved from the Open Science Framework mirror as a special case in the workflow.

## Pipeline Overview

![Pipeline Flowchart](figures/Sup%20Material/Figure_S2_pipeline_flowchart.png)

| Rule file | Tools | Purpose |
|-----------|-------|---------|
| `annotation.smk` | Liftoff, gffread | Annotation transfer and protein extraction |
| `orthology.smk` | OrthoFinder, DIAMOND | Orthology inference and pangenome partitioning |
| `phylogenomics.smk` | MAFFT, trimAl, IQ-TREE 2 | SCO alignment, species tree, concordance factors |
| `gene_families.smk` | CAFE5 | Gene family expansion/contraction analysis |
| `synteny.smk` | skani, minimap2 (nucmer rule retained but unused) | Whole-genome ANI and pairwise synteny |
| `repeats.smk` | RepeatModeler2, RepeatMasker | De novo TE library and genome annotation |

> **ANI / synteny methods:** The results reported in the manuscript use **skani v0.2** for whole-genome ANI (primary) with a **minimap2** alignment-based cross-check, and gene-anchor synteny derived from the Liftoff coordinate projections. A legacy `nucmer`/`delta-filter`/`show-coords` path remains in `synteny.smk` but was **not** used to produce the reported figures.

The anvi'o gene-cluster pangenome (Fig 8) was run **manually outside Snakemake**; commands are documented in the manuscript Methods.

## Quick Start

```bash
# Clone the repository
git clone https://github.com/tylermaire/cx_pipiens_pangenome.git
cd cx_pipiens_pangenome

# Install Snakemake
conda install -c bioconda -c conda-forge snakemake

# Dry run
snakemake -n --use-conda --cores 16

# Full run (requires ~128 GB RAM, ~40 hours on 16 cores)
snakemake --use-conda --cores 16 \
  results/pangenome/pangenome_summary.tsv \
  results/phylo/concord.cf.tree \
  results/cafe/significant_families.tsv \
  results/synteny/ani_matrix.tsv \
  results/synteny/synteny_summary.tsv \
  results/repeats/Cx_quinquefasciatus/Cx_quinquefasciatus.fasta.tbl \
  results/repeats/Cx_pallens/Cx_pallens.fasta.tbl \
  results/repeats/Cx_molestus/Cx_molestus.fasta.tbl \
  results/repeats/Cx_pipiens/Cx_pipiens.fasta.tbl
```

## Key Results

| Analysis | Result |
|----------|--------|
| Total orthogroups (5 taxa) | 16,568 |
| Ingroup genes partitioned | 92,821 (98.2% of ingroup gene calls) |
| Core genome | 11,284 orthogroups (71.1%) |
| Shell genome | 3,726 orthogroups (23.5%) |
| Cloud genome | 871 orthogroups (5.5%) |
| Outgroup-only (Cx. tarsalis) | 687 orthogroups |
| Form-specific cloud | molestus 297, pipiens 245, pallens 203, quinquefasciatus 126 (chi-sq = 71.9, p = 1.7e-15) |
| Species tree | (molestus,(pallens,quinquefasciatus),pipiens) |
| Deep ingroup node support | 100 / 58.9 / 45.1 (UFBoot / gCF / sCF) |
| CAFE5 significant families | 617 (p < 0.05); lambda = 0.105, alpha = 0.167 |
| ANI range (skani) | 93.07 - 94.72% |
| Synteny | 1:1 chromosome-scale; 258-294 inversions >= 100 kb per pair |
| TE content | 52.69 - 54.43% interspersed repeats |
| Anvi'o gene clusters | 16,673 across 4 genomes |

## Figures

Publication figures are generated with the R scripts in `workflow/scripts/` (ggplot2, patchwork, ggtree, UpSetR). The anvi'o figure is produced separately (see Methods).

| Figure | Description |
|--------|-------------|
| Fig 1 | Genome assembly quality (BUSCO, QUAST) |
| Fig 2 | Pangenome composition (core/shell/cloud, per-form, form-specific cloud) |
| Fig 3 | Concatenation ML species tree with gCF/sCF concordance |
| Fig 4 | CAFE5 per-lineage gene-family expansions/contractions |
| Fig 5 | Whole-genome ANI heatmap (skani) |
| Fig 6 | Pairwise gene-anchor synteny dotplots (composite) |
| Fig 7 | Transposable element composition |
| Fig 8 | Anvi'o gene-cluster pangenome (UpSet, per-genome totals, cluster sizes) |
| Fig 9 | Per-form copy-number heatmap (detoxification, chemosensory, immunity) |
| Fig S1 | Per-pair synteny dotplots (basis for composite Fig 6) |
| Fig S2 | Computational pipeline flowchart |

## Repository Structure

```
cx_pipiens_pangenome/
├── Snakefile                    # Main workflow entry point
├── config/
│   ├── config.yaml              # Pipeline parameters
│   └── samples.tsv              # Genome accessions and metadata
├── workflow/
│   ├── rules/                   # Snakemake rule files
│   │   ├── annotation.smk
│   │   ├── orthology.smk
│   │   ├── phylogenomics.smk
│   │   ├── gene_families.smk
│   │   ├── synteny.smk
│   │   └── repeats.smk
│   ├── envs/                    # Conda environment definitions
│   └── scripts/                 # Python/R helper scripts
├── results/                     # Key result summaries (tracked)
│   ├── pangenome/               # pangenome_summary.tsv, partitioned_orthogroups.tsv
│   ├── phylo/                   # species tree + concordance factor outputs
│   ├── cafe/                    # significant_families.tsv, branch summaries
│   ├── synteny/                 # ani_matrix.tsv, synteny_summary.tsv
│   └── functional/              # key gene-family tables
├── figures/                     # Main manuscript figures (PNG + PDF)
│   └── Sup Material/            # Supplementary figures (S1, S2) + Supplementary_Tables.xlsx
├── Culex_pipiens_pangenome_IBMB.docx   # Manuscript (IBMB format)
└── README.md
```

## Computational Requirements

The pipeline was developed and validated on AWS EC2 (r6i.4xlarge: 16 vCPU, 128 GB RAM, 300 GB EBS). RepeatModeler is the bottleneck at ~32 hours. Total runtime is approximately 40 hours on 16 cores.

## Citation

If you use this pipeline, please cite the tools used:

- **OrthoFinder**: Emms and Kelly (2019) Genome Biol. 20:238
- **DIAMOND**: Buchfink et al. (2015) Nat. Methods 12:59
- **MAFFT**: Katoh and Standley (2013) Mol. Biol. Evol. 30:772
- **trimAl**: Capella-Gutierrez et al. (2009) Bioinformatics 25:1972
- **IQ-TREE 2**: Minh et al. (2020) Mol. Biol. Evol. 37:1530
- **CAFE5**: Mendes et al. (2020) Bioinformatics 36:5516
- **skani**: Shaw and Yu (2023) Nat. Methods 20:1661
- **minimap2**: Li (2018) Bioinformatics 34:3094
- **Liftoff**: Shumate and Salzberg (2021) Bioinformatics 37:1639
- **RepeatModeler2/RepeatMasker**: Flynn et al. (2020) PNAS 117:9451
- **Anvi'o**: Eren et al. (2021) Nat. Microbiol. 6:3
- **Snakemake**: Molder et al. (2021) F1000Research 10:33

## License

MIT License. See [LICENSE](LICENSE) for details.

## Contact

Tyler Maire | Research Entomologist, Indian River Mosquito Control District
ORCID: [0009-0007-6451-5408](https://orcid.org/0009-0007-6451-5408)
