#!/usr/bin/env python3
"""Build Tables 1 to 4 of the manuscript from the workflow outputs, so no
number in a table is typed by hand.

    python workflow/scripts/make_manuscript_tables.py [--out tables/manuscript_tables.json]

Writes one JSON object per table (title, header, rows, column widths, notes)
that the manuscript build reads. Run from the repository root after the
workflow, or after patch_results.py on a stored run.
"""
import argparse
import csv
import glob
import json
import os

import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import forms  # noqa: E402

REPO = forms.REPO
ORDER = forms.INGROUP
OUTGROUP = forms.outgroup()
NAME = forms.NAME
SHORT = forms.SHORT
ASSEMBLY = {s: f"{acc} ({name})" for s, (acc, name) in forms.ACCESSION.items()}
ASSEMBLY["Cx_tarsalis"] = "CtarK1 (osf.io/mdwqx)"
# Summed length of the three chromosome scale sequences. forms reads it from
# results/synteny/genomes/<sample>.chromosomes.fasta (ingroup) or the
# downloaded genome (outgroup); these V4 measurements of the synteny FASTA
# files are used when neither file is present.
CHROM_BP_V4 = {"Cx_quinquefasciatus": 559588684, "Cx_pallens": 549446285,
               "Cx_molestus": 530551183, "Cx_pipiens": 532519368}


def chromosome_bp(sample):
    bp = forms.chromosome_scale_bp(sample)
    return bp if bp is not None else CHROM_BP_V4.get(sample)


CHROM_BP = {s: chromosome_bp(s) for s in ORDER + [OUTGROUP]}
CHROM_BP = {s: v for s, v in CHROM_BP.items() if v}


def tsv(path, skip_comments=True):
    with open(os.path.join(REPO, path)) as fh:
        lines = [l for l in fh if l.strip() and not (skip_comments and l.startswith("#"))]
    return list(csv.DictReader(lines, delimiter="\t"))


def values():
    out = {}
    for r in tsv("results/manuscript_values.tsv"):
        out[(r["section"], r["item"], r["sample"])] = r["value"]
    return out


def fmt_int(x):
    return f"{int(round(float(x))):,}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(REPO, "tables", "manuscript_tables.json"))
    args = ap.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    V = values()

    def version(tool, fallback):
        for sec in ("tools", "tools_observed"):
            v = V.get((sec, tool, ""))
            if v and v != "NA" and not v.startswith("not present"):
                return v.split(",")[0]
        return fallback

    quast = {r["Assembly"]: r for r in tsv("results/quast/report.tsv", skip_comments=False)}
    ref_busco = json.load(open(glob.glob(os.path.join(
        REPO, f"results/busco/{ORDER[0]}/short_summary*.json"))[0]))
    busco_version = ref_busco.get("versions", {}).get("busco", "6.1.0")
    busco_lineage = ref_busco.get("lineage_dataset", {}).get("name", "diptera_odb10")
    busco_n = int(ref_busco["results"]["n_markers"])
    tables = {}

    # Table 1: assemblies
    rows = []
    for s in ORDER + [OUTGROUP]:
        size = float(quast["Total length (>= 0 bp)"][s]) / 1e6
        seqs = int(quast["# contigs (>= 0 bp)"][s])
        n50 = float(quast["N50"][s]) / 1e6
        bj = json.load(open(glob.glob(os.path.join(
            REPO, f"results/busco/{s}/short_summary*.json"))[0]))["results"]
        cn50 = float(bj["Contigs N50"]) / 1e6
        chrom = f"{100 * CHROM_BP[s] / float(quast['Total length (>= 0 bp)'][s]):.1f}" \
            if s in CHROM_BP else "–"
        busco = V[("busco", "genome complete_pct", s)]
        rep = V.get(("repeats", "bases_masked_pct", s))
        mb = lambda x: f"{x:.1f}" if x >= 10 else f"{x:.2f}" if x >= 0.1 else f"{x:.3f}"
        rows.append([NAME[s], ASSEMBLY[s], f"{size:.1f}", f"{seqs:,}", mb(n50), mb(cn50),
                     chrom, busco, f"{float(rep):.1f}" if rep else "–"])
    tables["1"] = {
        "title": "Table 1. Genome assemblies",
        "header": ["Form", "Accession (assembly)", "Size (Mb)", "Sequences",
                   "Scaffold N50 (Mb)", "Contig N50 (Mb)", "In three chromosomes (%)",
                   "Complete BUSCO, genome (%)", "Repeats masked (%)"],
        "rows": rows,
        "widths": [1750, 2500, 650, 950, 850, 700, 1150, 1000, 900],
        "notes": [f"Size, sequences and scaffold N50 from QUAST {version('quast', '5.3.0')} "
                  "(all sequences counted); contig N50 as reported by BUSCO, which splits "
                  f"scaffolds at gaps. BUSCO {busco_version} with the {busco_lineage} dataset "
                  f"({busco_n:,} genes), genome mode with miniprot. Repeats were masked with "
                  "RepeatMasker using a RepeatModeler2 library built from the "
                  f"*Cx. quinquefasciatus* assembly; the outgroup, {NAME[OUTGROUP]}, was not "
                  "masked." + (" The *Cx. tarsalis* assembly is contig level."
                               if OUTGROUP == "Cx_tarsalis" else "")
                  + " Further statistics in Supp. Table S1."],
    }

    # Table 2: gene sets
    part = tsv("results/pangenome/partitioned_orthogroups.tsv")
    per = {s: {} for s in ORDER + [OUTGROUP]}
    for r in part:
        for s in per:
            per[s][r["compartment"]] = per[s].get(r["compartment"], 0) + int(r[s])
    un_file = glob.glob(os.path.join(REPO, "results/orthofinder/output/**/"
                                     "Orthogroups_UnassignedGenes.tsv"), recursive=True)[0]
    unassigned = {s: 0 for s in per}
    with open(un_file) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            for i, x in enumerate(f[1:], 1):
                if x.strip() and head[i] in unassigned:
                    unassigned[head[i]] += 1
    tq = {r["sample"]: r for r in tsv("results/annotation/transfer_quality.tsv")}
    rows = []
    og_only_note = ""
    for s in ORDER + [OUTGROUP]:
        kept = int(tq[s]["n_kept_models"])
        inval = tq[s]["pct_invalid_orf"]
        clean = tq[s]["pct_invalid_orf_clean_ref"]
        core, shell, cloud = per[s].get("core", 0), per[s].get("shell", 0), per[s].get("cloud", 0)
        og_only = per[s].get("outgroup_only", 0)
        assert core + shell + cloud + og_only + unassigned[s] == kept, s
        name = NAME[s] + (" (reference)" if s == "Cx_quinquefasciatus" else "")
        rows.append([name, f"{kept:,}", V[("busco", "protein complete_pct", s)],
                     f"{float(inval):.1f}" if inval else "–",
                     f"{float(clean):.1f}" if clean else "–",
                     f"{core:,}", f"{shell:,}", f"{cloud:,}",
                     f"{unassigned[s]:,}" + ("^a^" if og_only else "")])
        if og_only:
            og_only_note = (f"^a^A further {og_only:,} {NAME[s]} genes are in the "
                            f"{sum(1 for r in part if r['compartment'] == 'outgroup_only'):,} "
                            "outgroup only orthogroups.")
    tables["2"] = {
        "title": "Table 2. Gene sets after annotation transfer and their place in the pangenome",
        "header": ["Form", "Protein coding genes", "Complete BUSCO, protein (%)",
                   "Without a valid ORF (%)", "Without a valid ORF, clean reference models (%)",
                   "Core", "Shell", "Cloud", "Unassigned"],
        "rows": rows,
        "widths": [2000, 950, 950, 900, 1150, 850, 800, 700, 1300],
        "notes": ["One protein per gene (longest isoform). Valid open reading frames (ORF) as "
                  "recorded by Liftoff; the last ORF column sets aside models whose "
                  "*Cx. quinquefasciatus* source model is partial or carries a RefSeq sequence "
                  f"exception ({int(tq[ORDER[0]]['n_ref_model_not_clean']):,} of the kept "
                  "reference models)."
                  + (f" {NAME[OUTGROUP]} keeps {forms.OUTGROUP_ANNOTATION.get(OUTGROUP, 'its own gene set')}, "
                     "which carries no Liftoff flags."
                     if tq[OUTGROUP].get("liftoff_flags", "").startswith("no") else "")
                  + " Core, shell and cloud give the "
                  "number of genes of each form in orthogroups of that compartment, classified "
                  "by the four ingroup forms; unassigned genes were placed in no orthogroup.",
                  og_only_note],
    }

    # Table 3: pairwise identity and synteny
    syn = {(r["sample1"], r["sample2"]): r for r in tsv("results/synteny/synteny_summary.tsv")}
    n_trimmed = int(float(V.get(("phylogeny", "trimmed_alignments", ""), "0") or 0)) or \
        max(int(r["n_loci"]) for r in tsv("results/phylo/sco_pairwise_identity.tsv"))
    sco = {}
    for r in tsv("results/phylo/sco_pairwise_identity.tsv"):
        sco[frozenset([r["taxon_a"], r["taxon_b"]])] = r
    ani_rows = tsv("results/synteny/ani_pairs.tsv")
    rows = []
    for (a, b), r in syn.items():
        pair = frozenset([a, b])
        ar = next(x for x in ani_rows if frozenset([x["sample1"], x["sample2"]]) == pair)
        af = sorted(float(ar[k]) for k in ar if "align_fraction" in k)
        rows.append([f"{SHORT[a]} vs {SHORT[b]}", f"{float(ar['ani']):.2f}",
                     f"{af[0]:.1f}–{af[1]:.1f}",
                     f"{float(r['mean_alignment_identity_pct']):.2f}",
                     f"{100 * float(sco[pair]['median_identity']):.2f}",
                     f"{int(r['n_shared_genes_chr1_3']):,}",
                     f"{float(r['pct_collinear_anchors']):.1f}",
                     f"{float(r['pct_shared_genes_other_chromosome']):.1f}",
                     r["n_inversions_intra_chr"]])
    tables["3"] = {
        "title": "Table 3. Pairwise identity and synteny between the ingroup assemblies",
        "header": ["Pair", "skani ANI (%)", "skani aligned fraction (%)",
                   "minimap2 identity (%)", "Median protein identity (%)",
                   "Anchors on homologous chromosomes", "Collinear anchors (%)",
                   "Shared genes on another chromosome (%)", "Inversions ≥ 100 kb"],
        "rows": rows,
        "widths": [2100, 800, 1050, 950, 1000, 1100, 900, 1150, 850],
        "notes": ["skani ANI and aligned fraction (both directions) on the complete "
                  "assemblies; minimap2 identity is the length weighted identity of "
                  "alignments between the three chromosome scale sequences; median protein "
                  f"identity is over the {n_trimmed:,} trimmed single copy ortholog alignments. Anchors "
                  "are genes, coding or non coding, shared by both assemblies on homologous "
                  "chromosomes. Further values in Supp. Tables S7 and S9."],
    }

    # Table 4: key families
    names = {"P450": "Cytochrome P450s", "OR": "Odorant receptors",
             "OBP": "Odorant binding proteins", "CCE": "Carboxylesterases",
             "GR": "Gustatory receptors", "IR": "Ionotropic receptors",
             "immune": "Immune genes", "GST": "Glutathione S transferases",
             "CSP": "Chemosensory proteins"}
    fam = tsv("results/functional/key_families_wide.tsv")
    fam.sort(key=lambda r: -int(r["Cx_quinquefasciatus"]))
    rows = [[names.get(r["family"], r["family"])] + [r[s] for s in ORDER] for r in fam]
    tables["4"] = {
        "title": "Table 4. Genes assigned to key gene families in each form",
        "header": ["Family"] + [NAME[s] for s in ORDER],
        "rows": rows,
        "widths": [3000, 1800, 1800, 1800, 1800],
        "notes": ["Genes in orthogroups assigned to each family from eggNOG mapper annotations. "
                  "The three transferred gene sets carry fewer copies than the reference in "
                  "most families, as expected from copy loss during transfer, so the counts "
                  "are an inventory of the gene sets rather than a comparison of the forms."],
    }
    with open(args.out, "w") as fh:
        json.dump(tables, fh, indent=1, ensure_ascii=False)
    for k, t in tables.items():
        print(t["title"])
        for r in t["rows"]:
            print("   ", " | ".join(r))


if __name__ == "__main__":
    main()
