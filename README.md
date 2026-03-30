# Culex pipiens Complex Pangenome Analysis
 
A reproducible Snakemake pipeline for pangenome analysis of the *Culex pipiens*
species complex, comparing four chromosome-scale genome assemblies representing
*Cx. quinquefasciatus*, *Cx. pipiens* f. *pipiens*, *Cx. pipiens* f. *pallens*,
and *Cx. pipiens* f. *molestus*.
 
## Quick Start
 
```bash
# Clone the repository
git clone https://github.com/tylermaire/cx_pipiens_pangenome.git
cd cx_pipiens_pangenome
 
# Install Snakemake (if not already installed)
conda install -c bioconda -c conda-forge snakemake mamba
 
# Dry run (preview what will execute)
snakemake -n --use-conda
 
# Full run
snakemake --use-conda --cores 8
 
# Generate report
snakemake --report report.html
```
 
## Pipeline Overview
 
| Module | Purpose |
|--------|---------|
| download.smk | Fetch assemblies from NCBI |
| qc.smk | BUSCO + QUAST quality assessment |
| annotation.smk | Liftoff annotation transfer |
| orthology.smk | OrthoFinder + pangenome partitioning |
| functional.smk | eggNOG-mapper + GO enrichment |
| phylogenomics.smk | IQ-TREE + ASTRAL species tree |
| gene_families.smk | CAFE5 gene family evolution |
| synteny.smk | SyRI structural variants |
| repeats.smk | RepeatMasker TE analysis |
| figures.smk | All publication figures |
 
## Genome Assemblies
 
| Form | Accession | Source |
|------|-----------|--------|
| *Cx. quinquefasciatus* | GCF_015732765.1 | Ryazansky et al. 2024 |
| *Cx. pipiens* f. *pallens* | GCF_016801865.2 | Liu et al. 2023 |
| *Cx. pipiens* f. *molestus* | GCA_024516115.1 | Liu et al. 2023 |
| *Cx. pipiens* f. *pipiens* | GCA_963924435.1 | Hesson et al. 2025 |
| *Cx. tarsalis* (outgroup) | GCA_016859205.1 | Main et al. 2021 |
 
## Citation
 
If you use this pipeline, please cite:
 
Maire T et al. (2026) Pangenome analysis of the *Culex pipiens* complex
reveals form-specific gene content in a globally important mosquito vector
group. *Genome Biology and Evolution*. [In preparation]
 
## License
 
MIT License. See LICENSE for details.
 
## Authors
 
- Tyler Maire - Indian River Mosquito Control District / UMGC
