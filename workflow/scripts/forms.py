#!/usr/bin/env python3
"""forms.py - names, accessions and roles of the samples, for the scripts
that build the manuscript tables, figures and supplement.

The outgroup is read from config/samples.tsv (is_outgroup), so the builders
follow the sample sheet: Cx. perexiguus from V5 on, Cx. tarsalis before.
"""
import csv
import os

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
INGROUP = ["Cx_quinquefasciatus", "Cx_pallens", "Cx_molestus", "Cx_pipiens"]
REFERENCE = "Cx_quinquefasciatus"
TRANSFERRED = ["Cx_pallens", "Cx_molestus", "Cx_pipiens"]

# names in the manuscript tables (markdown italics)
NAME = {"Cx_quinquefasciatus": "*Cx. quinquefasciatus*", "Cx_pallens": "*Cx. pipiens pallens*",
        "Cx_molestus": "*Cx. pipiens* f. *molestus*", "Cx_pipiens": "*Cx. pipiens* f. *pipiens*",
        "Cx_perexiguus": "*Cx. perexiguus*", "Cx_tarsalis": "*Cx. tarsalis*"}
SHORT = {"Cx_quinquefasciatus": "*quinquefasciatus*", "Cx_pallens": "*pallens*",
         "Cx_molestus": "*molestus*", "Cx_pipiens": "*pipiens*",
         "Cx_perexiguus": "*Cx. perexiguus*", "Cx_tarsalis": "*Cx. tarsalis*"}
# names in the supplement (plain text)
FORM = {"Cx_quinquefasciatus": "quinquefasciatus", "Cx_pallens": "pallens",
        "Cx_molestus": "molestus", "Cx_pipiens": "pipiens",
        "Cx_perexiguus": "Cx. perexiguus", "Cx_tarsalis": "Cx. tarsalis"}
# abbreviations in figure labels
ABBR = {"Cx_quinquefasciatus": "qui", "Cx_pallens": "pal", "Cx_molestus": "mol",
        "Cx_pipiens": "pip", "Cx_perexiguus": "per", "Cx_tarsalis": "tar"}
ACCESSION = {"Cx_quinquefasciatus": ("GCF_015732765.1", "VPISU_Cqui_1.0_pri_paternal"),
             "Cx_pallens": ("GCF_016801865.2", "TS_CPP_V2"),
             "Cx_molestus": ("GCA_024516115.1", "TS_CPM_V1"),
             "Cx_pipiens": ("GCA_963924435.1", "idCulPipi1.1"),
             "Cx_perexiguus": ("GCA_964243045.1", "idCulPerx1.1"),
             "Cx_tarsalis": ("osf.io/mdwqx", "CtarK1")}
# how the outgroup's gene set was obtained, for captions
OUTGROUP_ANNOTATION = {"Cx_perexiguus": "its own Ensembl gene set",
                       "Cx_tarsalis": "models transferred with Liftoff"}


def sample_sheet(path=None):
    """Rows of config/samples.tsv, or of the sheet named by the environment
    variable CX_SAMPLE_SHEET (to rebuild tables of an earlier run)."""
    path = path or os.environ.get("CX_SAMPLE_SHEET") or os.path.join(REPO, "config", "samples.tsv")
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def outgroup(path=None):
    rows = [r["sample"] for r in sample_sheet(path)
            if r["is_outgroup"].strip().lower() == "true"]
    if len(rows) != 1:
        raise SystemExit(f"expected one outgroup in the sample sheet, found {rows}")
    return rows[0]


def all_samples(path=None):
    """Ingroup in manuscript order, then the outgroup."""
    return INGROUP + [outgroup(path)]


def figure_name(sample):
    """'Cx. quinquefasciatus' style label for figures."""
    return "Cx. " + sample.split("_", 1)[1]


def abbr_split(split):
    """'Cx_molestus+Cx_pipiens | Cx_pallens+Cx_quinquefasciatus' -> 'mol + pip | pal + qui'."""
    sides = [side.strip() for side in split.split("|")]
    return " | ".join(" + ".join(ABBR.get(t.strip(), t.strip()) for t in side.split("+"))
                      for side in sides)


def chromosome_scale_bp(sample, n=3):
    """Summed length of the n longest sequences: from the synteny FASTA of an
    ingroup assembly, or from the downloaded genome; None when neither is
    present."""
    for p in (os.path.join(REPO, "results", "synteny", "genomes", f"{sample}.chromosomes.fasta"),
              os.path.join(REPO, "resources", "genomes", f"{sample}.fasta")):
        if os.path.exists(p):
            lengths, cur = [], None
            with open(p) as fh:
                for line in fh:
                    if line.startswith(">"):
                        if cur is not None:
                            lengths.append(cur)
                        cur = 0
                    elif cur is not None:
                        cur += len(line.strip())
            if cur is not None:
                lengths.append(cur)
            return sum(sorted(lengths, reverse=True)[:n])
    return None
