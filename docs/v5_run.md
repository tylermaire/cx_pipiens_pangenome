# V5: *Cx. perexiguus* as the outgroup

Branch `v5-perexiguus`, based on `v4-review-fixes` (commit `e14394f`).

## Why

Up to V4 the outgroup was *Cx. tarsalis* CtarK1, a contig level assembly.
Liftoff transfers to it lacked a valid ORF in 75.5% of kept models, skani could
not score it against any ingroup form, and the four taxon analyses could not be
rooted, so the excess discordant topology could not be tied to particular
forms. *Cx. perexiguus* sits in the Univittatus Subgroup of the Pipiens Group,
outside the *Cx. pipiens* complex and much closer to it than *Cx. tarsalis*. Its
assembly idCulPerx1.1 (GCA_964243045.1; PacBio HiFi with Hi-C, three
chromosomes, 436 Mb) comes with its own Ensembl gene set (Ruiz-López et al.
2026, Genome Biology and Evolution 18: evaf245).

## What changed

| Change | Where |
| --- | --- |
| Outgroup *Cx. perexiguus*; sample sheet gains `annotation` (reference, liftoff, native) and `source` (ncbi, ensembl) columns | `config/samples.tsv`, `config/config.yaml`, `Snakefile` |
| Genome and gene set downloaded together from the Ensembl FTP, newest gene set release found automatically, sequence names checked against each other | rule `download_ensembl`, `check_seqids.py` |
| The outgroup keeps its own annotation (Ensembl type prefixes removed from identifiers); Liftoff only for pallens, molestus and pipiens | rules `native_annotation`, `liftoff`; `normalize_gff.py` |
| Coding sequences of the kept proteins, under the same names | rule `extract_cds`; `filter_cds.py` |
| Five taxon single copy loci, codon alignments on the trimAl columns of the protein alignments | rules `extract_sco5`, `rooted_alignments`; `extract_sco5.py`, `codon_align.py` |
| Concatenated protein tree rooted with the outgroup, from loci whose four ingroup models are intact, gene and site concordance, rooted topologies of the gene trees | rules `rooted_locus_quality`, `rooted_tree`, `rooted_summary`; `locus_quality.py`, `rooted_summary.py`, `rooting.py` |
| Derived allele patterns by gene model quality, and the root and D statistics on filtered locus sets (Supp. Table S13) | rule `rooting_diagnostics`; `rooting_diagnostics.py` |
| ABBA BABA tests: four planned tests, all sites and third codon positions, all loci and loci with intact models, block jackknife over 5 Mb windows of the reference | rule `d_statistics`; `d_statistics.py`; `dstat` in `config.yaml` |
| CAFE on the ingroup counts, on the species tree rooted with the outgroup (V4 assumed the root) | `format_cafe_input.py`, rule `prepare_cafe_input` |
| The V4 RepeatModeler library is reused (about 18 hours saved) | `data/repeat_library/custom_repeat_lib_V4.fa`, rule `repeatmodeler` |
| Figures, Tables 1 to 4 and Supp. Tables S1 to S13 built by the workflow; captions read their counts and versions from the outputs; new S12 and S13 for the outgroup analyses; Figure 2 shows the rooted tree and the D statistics | rules `revision_figures`, `manuscript_tables`; `forms.py` and the three builders |
| New results in the values table (sections `rooted`, `rooted_topologies`, `dstat`) | `collect_manuscript_values.py`, `report.smk` |
| Tests, including a Snakemake run of the outgroup rules on a simulated data set with gene flow from pallens into pipiens | `tests/test_rooted_analyses.py` |

The four taxon analyses (quartet tests, divergence diagnostics, synteny,
absence validation) are unchanged. The ingroup assemblies are the same, so the
Liftoff transfers should reproduce V4 (a useful check: compare
`results/annotation/transfer_quality.tsv` with V4). Orthogroups change,
because OrthoFinder now clusters with a different fifth proteome.

## Running it on Ubuntu

The commands are the same on an AWS instance (V4 used a c7i.24xlarge, 96
vCPU) and on your own Ubuntu machine, including Ubuntu under WSL on Windows.
`$(nproc)` is the number of cores the machine has.

One time setup, if the machine has no conda or Snakemake yet:

```bash
sudo apt update && sudo apt install -y git curl tmux
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p "$HOME/miniforge3"
source "$HOME/miniforge3/bin/activate"
conda config --set channel_priority strict
conda create -y -n snakemake -c conda-forge -c bioconda snakemake
```

Each new terminal then starts with
`source "$HOME/miniforge3/bin/activate" && conda activate snakemake`.

```bash
# clone the branch directly: some TSV files on main are stored with Windows
# line endings, and a later "git checkout" stops on them as local changes
git clone -b v5-perexiguus https://github.com/tylermaire/cx_pipiens_pangenome.git
cd cx_pipiens_pangenome

# 1. the plan: 91 jobs, no errors (always with --use-conda: without it,
#    Snakemake reports every software environment as changed)
snakemake -n --use-conda --cores $(nproc)

# 2. build the conda environments first; a failure here costs minutes
snakemake --use-conda --conda-create-envs-only --cores 8

# 3. download and check the outgroup on its own
snakemake --use-conda --cores 8 resources/genomes/Cx_perexiguus.fasta resources/annotations/Cx_perexiguus.gff3
cat resources/annotations/Cx_perexiguus.source.txt
```

Step 3 prints the two addresses it used and a line such as
`... features, ... CDS, 100.00% on sequences in the FASTA`. If it stops
instead, open https://ftp.ebi.ac.uk/pub/ensemblorganisms/Culex_perexiguus/ in a
browser, find the unmasked genome FASTA and the gene set GFF3 of
GCA_964243045.1, put the two full addresses in `config/config.yaml`
(`outgroup: genome_url` and `gff_url`) and run step 3 again.

```bash
# 4. the full run, as a system service: it keeps going when the terminal
#    closes, and one job killed for memory does not stop the others
sudo systemd-run --unit=cxv5 --same-dir --uid=$(id -u) --gid=$(id -g) \
  -E PATH="$PATH" -E HOME="$HOME" -E USER="$USER" \
  -p OOMPolicy=continue -p LimitNOFILE=65536 \
  -p StandardOutput=append:"$PWD/run_v5.log" -p StandardError=append:"$PWD/run_v5.log" \
  "$(which snakemake)" --use-conda --cores $(nproc) --keep-going --rerun-incomplete \
  --rerun-triggers mtime \
  --resources mem_mb=175000 \
  --set-resources eggnog_mapper:mem_mb=20000 repeatmasker:mem_mb=16000 orthofinder:mem_mb=30000 \
  concat_and_tree:mem_mb=30000 rooted_tree:mem_mb=30000 absence_dna_align:mem_mb=20000 \
  absence_protein_align:mem_mb=10000 busco_proteins:mem_mb=5000 skani_ani:mem_mb=10000

systemctl status cxv5 --no-pager     # active (running) while it works
tail -f run_v5.log                   # Ctrl+c stops the viewing, not the run
```

`mem_mb=175000` is for a machine with 185 GB of memory; set it about 10 GB
below the total that `free -g` reports. The rule `synteny_minimap2` declares
160 GB itself (minimap2 asm20 peaks near 150 GB), so the limit lets only one
of those two jobs run at a time. Declared memory is only enforced when
`--resources mem_mb=` is given.

Expect roughly 12 to 18 hours on 96 cores, against 31 for V4 (17.7 of those
were RepeatModeler); on 16 cores expect several days. The last lines of a
finished run read `91 of 91 steps (100%) done`; with `--keep-going`, any
failures are listed by `grep -n "Error in rule" run_v5.log`. After a stop,
`sudo systemctl reset-failed cxv5` frees the unit name, and the same command
resumes where the run left off (`--rerun-triggers mtime` keeps finished jobs
from being redone because of a changed environment file or rule).

## Problems met in the V5 run

* **Out of memory (first attempt, under tmux).** Two `synteny_minimap2` jobs
  and eggNOG mapper ran together and reached 183.6 GB. The kernel killed one
  minimap2, and systemd then stopped the whole tmux session (Ubuntu 26.04
  runs each tmux session as a systemd scope that stops when one process is
  killed for memory). Fixed by the memory limits and the service form above.
* **IQ-TREE stopped in `concat_and_tree`.** On the four taxon matrix of 9,126
  loci, ModelFinder in IQ-TREE 3.1.3 failed an internal check while fitting
  FreeRate models (`ratefree.cpp: RateFree::initFromCatMinusOne`, exit status
  134). The four taxon concatenated tree now leaves FreeRate models out
  (`iqtree: concat_mrate: E,I,G,I+G` in `config.yaml`); in V4, with them
  allowed, FreeRate models were chosen for only 8 of 9,098 loci. Gene trees and
  the rooted tree keep the full model set. The IQ-TREE log of the failure is
  kept as `results/phylo/iqtree_crash_v5.log` in the run archive, and its first
  and last lines as `results/phylo/iqtree_crash_v5_excerpt.log`.
* **Shared transfer errors moved the root.** The first rooted tree, from all
  6,831 five taxon loci with a trimmed alignment (6,846 loci in all), placed
  the root on the quinquefasciatus branch
  (UFBoot 100). Split by gene model quality, the sites carrying that signal
  (pallens, molestus and pipiens derived, quinquefasciatus and the outgroup
  ancestral) came from loci whose models are not intact: 256 per locus where
  RefSeq had corrected the reference model for an error in the reference
  genome, against 4 in intact loci. In such loci the three transferred models
  share a shifted reading frame, so their codons are 91 to 99% identical to
  each other and 40 to 82% identical to quinquefasciatus and the outgroup, and
  every misaligned codon counts as a shared derived allele; errors shared by
  two forms inflate that pair in the same way. On the 4,592 intact loci the
  sites favor a root between (molestus, pipiens) and (pallens,
  quinquefasciatus), the most frequent rooted gene tree is that topology, and
  the D statistics do not change when loci are further filtered by ingroup
  identity. The rooted tree, its concordance factors and the rooted gene tree
  counts now use the intact loci (`rooted: loci: intact`, rule
  `rooted_locus_quality`); the D statistics are reported on all loci and on
  intact loci, and `rooting_diagnostics` writes the comparison (Supp. Table
  S13). The four taxon site counts (`quartet_site_patterns.tsv`) are repeated
  on intact loci for the same reason.
  Of the 4,592 intact loci, 4,017 have a rooted gene tree; the other 575
  have at most two distinct protein sequences in their trimmed alignments
  (100 with one, 475 with two; checked on the downloaded run, Sep 27), too
  few for IQ-TREE to build a tree.
* **CAFE5's P values changed the count of significant families by chance.**
  The run called 834 families significant, against 333 in V4, with nearly the
  same model (lambda 0.0449 against 0.0432). The difference is 572 families
  with one copy in one form of each pair and none in the other two (for
  example pallens and pipiens), which all share one P value. CAFE5 estimates
  each P value from 1,000 families it simulates per root size, with an
  unseeded random number generator (`compute_pvalues(..., 1000)` in
  `execute.cpp`), so the estimate carries an error of about 0.007 near 0.05.
  It put these families at 0.042 in this run and at 0.054 in V4; the exact
  value is 0.0555. In 20 further CAFE5 runs on the same counts, with lambda
  and alpha fixed at the fitted values, the number of families with estimated
  P < 0.05 ranged from 192 to 837 (median 237; above 800 in 4 of 20), and the
  estimate for these families averaged 0.056
  (`results/cafe/cafe5_reruns.tsv`, Supp. Table S8; CAFE5 built from the hahnlab/CAFE5 repository
  at commit b9e3b2e, whose P value code is that of the 5.1.0 release the run
  used). With four ingroup forms the distribution can be summed exactly, so
  the workflow now computes the P values exactly under the fitted model, as
  CAFE5 defines them (`cafe: pvalues: exact`, rule `cafe_exact_pvalues`,
  `results/cafe/family_pvalues.tsv`), and keeps CAFE5's estimate beside each.
  With exact P values 264 families are significant (334 in V4, whose
  estimates agreed with the exact values to within 0.003 on average), and
  the transfer control reads as in V4: 12 of 11,188 single copy families
  (0.1%) against 182 of 1,108 multi copy families (16.4%), odds ratio 183.
  The model, the tree and the per branch counts are CAFE5's and do not change.
  The CAFE outputs, Figure 3, the values table and Supp. Table S8 were
  rebuilt from the run's CAFE5 output with `patch_results.py --steps
  cafe,values --values-sections cafe,parameters` and the figure and table
  scripts; nothing upstream was rerun.

On your own machine rather than AWS:

* Disk: about 200 GB free (`df -h ~`); the eggNOG database alone is about 60 GB.
* Memory: 64 GB or more is safest (`free -g`). WSL gives Ubuntu half of the
  PC's memory unless `.wslconfig` in the Windows user folder sets `memory=`.
* Under WSL, clone into the Ubuntu home folder (`~`), not into `/mnt/c` or
  `/mnt/d`: the Windows drives are many times slower from Ubuntu and conda
  environments do not work well there. Keep the PC from sleeping during the run.

## After the run

```bash
# the outputs needed to check the run and update the manuscript (about 1 GB)
bash workflow/scripts/pack_results.sh results_v5_small.tar.gz

# everything, for the record, as for V4
tar -czf results_v5_full.tar.gz results figures/revision tables run_v5.log
```

Copy `results_v5_small.tar.gz` to the Culex Pangeneome folder on your computer
(from AWS for example with `aws s3 cp` or `scp`; under WSL the Windows drives
are `/mnt/c` and `/mnt/d`, so `cp` is enough), keep the full archive, and on
AWS stop the instance. With the small archive, every number, table and figure
of the manuscript can be updated, and any summary script can be fixed and
rerun locally with `patch_results.py` (steps now include `rooted` and `dstat`)
without another run.

## Reading the D statistics

`results/phylo/dstat/d_statistics.tsv`, one row per test, locus set and site
class. D(P1, P2; P3, O) > 0 means P2 and P3 share more derived alleles than P1
and P3; D < 0 the reverse; |Z| of 3 or more is the usual threshold. The planned
tests and what gene flow between each pair would do to them are described in
`config/config.yaml` (`dstat`). The simulation in
`tests/test_rooted_analyses.py` (gene flow from the pallens lineage into
pipiens) gives D < 0 in the first test and D > 0 in the third and fourth, with
the second near 0.

## What could not be checked before the run

* The Ensembl FTP layout: automated access to ftp.ebi.ac.uk is blocked from
  the build environment, so `download_ensembl` finds the files at run time
  (tested against a local copy of the expected layout) and step 3 checks it
  before anything expensive starts.
* The conda environments: the environment files are unchanged from V4 except
  `curl` added to `phylo.yaml` and the new `figures.yaml`; step 2 builds them.
* IQ-TREE: the new rules use the options of the V4 rules and were tested with
  IQ-TREE 2.0.7; V4 ran IQ-TREE 3.1.3.

## The anvi'o pangenome

The anvi'o gene cluster pangenome runs outside the workflow, on the proteins
OrthoFinder used (one per gene, `results/proteins`) and their coding sequences
(`results/cds`):

```bash
THREADS=8 bash workflow/scripts/anvio_pangenome.sh
```

`anvio_inputs.py` makes one contig per gene with an external gene call that
carries the protein itself, so anvi'o clusters those proteins rather than its
own translations (3,258 to 3,328 transferred models per form have a coding
sequence that is not a multiple of three and are flagged partial; one
*pipiens* contig of 3 nt is padded with N to anvi'o's k-mer size of 4). The
script then builds the contigs databases, the genomes storage (with
`--gene-caller workflow`, the source named in the gene calls), and the
pangenome (DIAMOND, minbit 0.5, MCL inflation 10, FAMSA), adds the default
collection, summarizes it, and calls `anvio_compare.py`, which writes the
tables in `results/anvio` (`anvio_*.tsv`, `anvio_gene_clusters.tsv.gz` and
`versions.tsv`; the databases, inputs, pangenome and summary directories are
not tracked). Figure 5 and Supplementary Table S14 are drawn from those
tables, and `patch_results.py --steps values --values-sections anvio` puts
them in `results/manuscript_values.tsv`. `anvio_compare.py` reads the
annotation GFFs for the gene of each model and its ORF flag; without them it
takes the same content from a table (`--gene-table`: sample, transcript, gene,
valid_ORF, reference_partial, reference_exception). The V5 comparison ran
where the GFFs were not at hand, from such a table made from the V5 GFFs with
the same parsing; the per gene table `anvio_gene_clusters.tsv.gz` keeps the
gene and ORF flag of every model it used.

The V5 run (September 27) used anvi'o 9 (eunice) in a Python 3.10 virtual
environment installed from PyPI, DIAMOND 2.2.8 (the version of the
OrthoFinder run), MCL 22-282 and FAMSA 2.2.2, on two cores; the pangenome
step took about 45 minutes, most of it the homogeneity indices of the gene
clusters (`--skip-homogeneity` saves that time when the interactive display
is not needed). Two problems to know about:

* FAMSA 2.4.1 to 2.5.2 cannot read standard input, which is how anvi'o
  passes sequences to it, so every gene cluster of two or more genes is left
  unaligned; anvi-pan-genome only warns ("The alignment of sequences failed
  for ...") and goes on. FAMSA 2.2.2 works.
* anvi'o builds its MCL graph from the DIAMOND hits alone, so proteins with no
  hit, not even to themselves, are left out of the gene clusters without a
  warning: 99 of the 57,815 proteins, mostly fragments (median 12 residues),
  97 of which OrthoFinder had also left unassigned. `anvio_compare.py`
  reports them as not clustered.

Result: 16,623 gene clusters, 9,465 of them with genes of all four forms and
3,016 of one form; 97.9% of the genes in four form clusters are in core
orthogroups, and 2,440 of the 2,445 transferred models in one form clusters
lack a valid ORF.

To view the pangenome in the anvi'o interface on a machine with a browser:

```bash
anvi-display-pan -p results/anvio/pan/Culex_pipiens_complex-PAN.db \
                 -g results/anvio/CULEX-GENOMES.db     # add -I localhost under WSL
```
