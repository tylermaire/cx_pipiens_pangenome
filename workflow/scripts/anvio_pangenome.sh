#!/usr/bin/env bash
# anvio_pangenome.sh - anvi'o gene cluster pangenome of the four ingroup forms.
#
# Run by hand, outside the Snakemake workflow, from the repository root with
# results/ holding the run (results/proteins, results/cds). Needs anvi'o 9
# with DIAMOND, MCL and FAMSA on the PATH (https://anvio.org/install/).
# The V5 analysis used anvi'o 9 (eunice), DIAMOND 2.2.8, MCL 22-282 and
# FAMSA 2.2.2. anvi'o passes the sequences to FAMSA on standard input, which
# FAMSA 2.4.1 to 2.5.2 fail to read, so that every gene cluster of more than
# one gene is left unaligned (anvi-pan-genome warns and goes on); use 2.2.2
# or check the pan genome log for that warning.
#
#   THREADS=8 bash workflow/scripts/anvio_pangenome.sh
#
# Steps:
#   1. anvio_inputs.py: one contig per gene with an external gene call that
#      carries the protein OrthoFinder used (57,815 proteins in V5)
#   2. one contigs database per form, then the genomes storage
#   3. anvi-pan-genome: DIAMOND all against all, minbit 0.5, MCL inflation 10
#      (the value anvi'o suggests for closely related genomes), FAMSA
#      alignments
#   4. a default collection and anvi-summarize (the gene cluster table)
#   5. anvio_compare.py: counts by combination of forms, singletons, and the
#      comparison with the OrthoFinder partition (results/anvio/*.tsv)
#
# The interactive view, on a machine with a browser:
#   anvi-display-pan -p results/anvio/pan/Culex_pipiens_complex-PAN.db \
#                    -g results/anvio/CULEX-GENOMES.db   (add -I localhost under WSL)
set -euo pipefail

OUT=${OUT:-results/anvio}
THREADS=${THREADS:-4}
PROJECT=Culex_pipiens_complex
FORMS="Cx_quinquefasciatus Cx_pallens Cx_molestus Cx_pipiens"

mkdir -p "$OUT/inputs" "$OUT/logs"
{
    printf "tool\tversion\n"
    printf "anvio\t%s\n" "$(anvi-pan-genome --version 2>&1 | grep -m1 "Anvi'o " | sed 's/.*: *//')"
    printf "diamond\t%s\n" "$(diamond version 2>&1 | awk 'NR==1 {print $NF}')"
    printf "mcl\t%s\n" "$(mcl --version 2>&1 | awk 'NR==1 {print $NF}')"
    printf "famsa\t%s\n" "$(famsa 2>&1 | awk '/version/ {print $2; exit}')"
} > "$OUT/versions.tsv"
python workflow/scripts/anvio_inputs.py --out "$OUT/inputs" --db-dir "$OUT"

for f in $FORMS; do
    rm -f "$OUT/$f.db"
    anvi-gen-contigs-database -f "$OUT/inputs/$f.contigs.fa" -o "$OUT/$f.db" -n "$f" \
        --external-gene-calls "$OUT/inputs/$f.gene_calls.tsv" \
        -T "$THREADS" > "$OUT/logs/contigs_db_$f.log" 2>&1
done

# the gene calls come from anvio_inputs.py, whose source column is "workflow"
rm -f "$OUT/CULEX-GENOMES.db"
anvi-gen-genomes-storage -e "$OUT/inputs/external_genomes.txt" -o "$OUT/CULEX-GENOMES.db" \
    --gene-caller workflow > "$OUT/logs/genomes_storage.log" 2>&1

anvi-pan-genome -g "$OUT/CULEX-GENOMES.db" -n "$PROJECT" -o "$OUT/pan" \
    --minbit 0.5 --mcl-inflation 10 --align-with famsa \
    --num-threads "$THREADS" --overwrite-output-destinations \
    > "$OUT/logs/pan_genome.log" 2>&1

anvi-script-add-default-collection -p "$OUT/pan/$PROJECT-PAN.db" > "$OUT/logs/collection.log" 2>&1
rm -rf "$OUT/summary"
anvi-summarize -p "$OUT/pan/$PROJECT-PAN.db" -g "$OUT/CULEX-GENOMES.db" -C DEFAULT \
    -o "$OUT/summary" > "$OUT/logs/summarize.log" 2>&1

python workflow/scripts/anvio_compare.py --summary "$OUT/summary" --inputs "$OUT/inputs" \
    --out "$OUT"
