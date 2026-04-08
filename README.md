# Culex pipiens Complex Pangenome Analysis

A reproducible Snakemake pipeline for pangenome analysis of the *Culex pipiens* species complex, comparing four chromosome-scale genome assemblies. This is the first pangenome analysis of this globally important mosquito vector group.

**Manuscript:** Maire T. Pangenome analysis reveals gene content variation and evolutionary dynamics across the *Culex pipiens* species complex. *In preparation* for Insects (MDPI).

## Genome Assemblies

| Form | Accession | Size (Mb) | Level |
|------|-----------|-----------|-------|
| *Cx. quinquefasciatus* (reference) | GCF_015732765.1 | 578 | Chromosome |
| *Cx. pipiens* f. *pallens* | GCF_016801865.2 | 568 | Chromosome |
| *Cx. pipiens* f. *molestus* | GCA_024516115.1 | 560 | Chromosome |
| *Cx. pipiens* f. *pipiens* | GCA_963924435.1 | 533 | Chromosome |
| *Cx. tarsalis* (outgroup) | GCA_016859205.1 | 50 | Scaffold |

## Pipeline Overview

```
NCBI GenBank assemblies (.fasta + .gff3)
        |
    Liftoff v1.6.3 (annotation transfer)
        |
    gffread v0.12.7 (protein extraction)
        |
   _____|___________________
  |           |              |
OrthoFinder  nucmer       RepeatModeler
  v2.5.5    MUMmer4         v2.0.7
  |          v4.0.0           |
  |           |           RepeatMasker
  |        delta-filter      v4.1.7
  |           |
  |        show-coords
  |        show-diff
  |___________|
  |
  |--- Partitioning (core/shell/cloud)
  |--- SCO extraction (8,616 single-copy orthologs)
  |        |
  |    MAFFT v7.520 + trimAl v1.4.1
  |        |
  |    IQ-TREE v2.2.6 (ML tree, 1000 UFBoot, gCF/sCF)
  |
  |--- CAFE5 v5.1.0 (gene family evolution, gamma k=3)
```

| Rule file | Tools | Purpose |
|-----------|-------|---------|
| `annotation.smk` | Liftoff, gffread | Annotation transfer and protein extraction |
| `orthology.smk` | OrthoFinder | Orthology inference and pangenome partitioning |
| `phylogenomics.smk` | MAFFT, trimAl, IQ-TREE | SCO alignment, species tree, concordance factors |
| `gene_families.smk` | CAFE5 | Gene family expansion/contraction analysis |
| `synteny.smk` | nucmer, delta-filter, show-coords | Pairwise whole-genome synteny and ANI |
| `repeats.smk` | RepeatModeler, RepeatMasker | De novo TE library and genome annotation |

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
  results/synteny/Cx_quinquefasciatus_vs_Cx_pallens.coords \
  results/synteny/Cx_quinquefasciatus_vs_Cx_molestus.coords \
  results/synteny/Cx_quinquefasciatus_vs_Cx_pipiens.coords \
  results/synteny/Cx_pallens_vs_Cx_molestus.coords \
  results/synteny/Cx_pallens_vs_Cx_pipiens.coords \
  results/synteny/Cx_molestus_vs_Cx_pipiens.coords \
  results/repeats/Cx_quinquefasciatus/Cx_quinquefasciatus.fasta.tbl \
  results/repeats/Cx_pallens/Cx_pallens.fasta.tbl \
  results/repeats/Cx_molestus/Cx_molestus.fasta.tbl \
  results/repeats/Cx_pipiens/Cx_pipiens.fasta.tbl
```

## Key Results

| Analysis | Result |
|----------|--------|
| Total orthogroups | 16,369 (92,540 genes, 97.9% assigned) |
| Core genome | 11,653 orthogroups (71.2%) |
| Shell genome | 4,205 orthogroups (25.7%) |
| Cloud genome | 511 orthogroups (3.1%) |
| Species tree | (molestus,(pallens,pipiens),quinquefasciatus), 100% bootstrap |
| Concordance | 100/100/100 gCF/sCF at all nodes |
| CAFE5 significant families | 596 (p < 0.05) |
| ANI range | 97.9 - 98.3% |
| TE content | 52.1 - 53.9% interspersed repeats |
| Anvi'o gene clusters | 32,304 across 4 genomes |

## Figures

Publication figures are generated with the R script `scripts/Cx_pipiens_pangenome_figures.R`. Requires ggplot2, patchwork, ggtree, UpSetR, plotrix, and cowplot.

| Figure | Description |
|--------|-------------|
| Fig 1a | UpSet plot of orthogroup intersections |
| Fig 1b-d | Accumulation curve, composition bars, form-specific counts |
| Fig 2a-b | Phylogenomic tree with concordance + CAFE5 expansions/contractions |
| Fig 3a-c | ANI heatmap + TE composition + predicted protein counts |
| Fig 4a-f | Synteny dot plots (6 pairwise comparisons) |
| Fig 5 | Pangenome flower diagram |
| Fig 6 | Anvi'o interactive pangenome visualization |
| Fig S1 | Computational pipeline flowchart |

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
│   └── scripts/                 # Python helper scripts
│       ├── extract_sco.py
│       ├── partition_pangenome.py
│       ├── format_cafe_input.py
│       └── parse_cafe.py
├── scripts/
│   └── Cx_pipiens_pangenome_figures.R
├── results/                     # Key result summaries (tracked)
│   ├── pangenome/
│   ├── phylo/
│   ├── cafe/
│   ├── synteny/                 # .coords files (not tracked, regenerate with pipeline)
│   └── repeats/                 # .tbl files (not tracked, regenerate with pipeline)
├── figures/
│   └── pipeline_flowchart.pdf
└── README.md
```

## Computational Requirements

The pipeline was developed and validated on AWS EC2 (r6i.4xlarge: 16 vCPU, 128 GB RAM, 300 GB EBS). RepeatModeler is the bottleneck at ~32 hours. Total runtime is approximately 40 hours on 16 cores.

## Citation

If you use this pipeline, please cite the tools used:

- **OrthoFinder**: Emms & Kelly (2019) Genome Biology 20:238
- **IQ-TREE**: Minh et al. (2020) Molecular Biology and Evolution 37:1530
- **CAFE5**: Mendes et al. (2020) Bioinformatics 36:5516
- **MUMmer4**: Marcais et al. (2018) PLoS Computational Biology 14:e1005944
- **RepeatModeler/RepeatMasker**: Flynn et al. (2020) PNAS 117:9451
- **Anvi'o**: Eren et al. (2021) Nature Microbiology 6:3
- **Snakemake**: Molder et al. (2021) F1000Research 10:33

## License

MIT License. See [LICENSE](LICENSE) for details.

## Contact

Tyler Maire | Research Entomologist, Indian River Mosquito Control District
ORCID: [0009-0007-6451-5408](https://orcid.org/0009-0007-6451-5408)
