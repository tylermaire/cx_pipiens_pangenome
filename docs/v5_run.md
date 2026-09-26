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
| Concatenated protein tree rooted with the outgroup, gene and site concordance, rooted topologies of the gene trees | rules `rooted_tree`, `rooted_summary`; `rooted_summary.py`, `rooting.py` |
| ABBA BABA tests: four planned tests, all sites and third codon positions, all loci and loci with intact models, block jackknife over 5 Mb windows of the reference | rule `d_statistics`; `d_statistics.py`; `dstat` in `config.yaml` |
| CAFE on the ingroup counts, on the species tree rooted with the outgroup (V4 assumed the root) | `format_cafe_input.py`, rule `prepare_cafe_input` |
| The V4 RepeatModeler library is reused (about 18 hours saved) | `data/repeat_library/custom_repeat_lib_V4.fa`, rule `repeatmodeler` |
| Figures, Tables 1 to 4 and Supp. Tables S1 to S12 built by the workflow; captions read their counts and versions from the outputs; new S12 for the outgroup analyses; Figure 2 shows the rooted tree and the D statistics | rules `revision_figures`, `manuscript_tables`; `forms.py` and the three builders |
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
  kept as `results/phylo/iqtree_crash_v5.log`.

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
