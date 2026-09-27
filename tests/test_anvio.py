"""Tests for workflow/scripts/anvio_inputs.py and anvio_compare.py
(run: python3 tests/test_anvio.py).

anvio_inputs.py must hand anvi'o the proteins as they are, flag coding
sequences that are not a multiple of three as partial and pad contigs shorter
than anvi'o's k-mer size. anvio_compare.py must class gene clusters by the
forms they hold, set apart the proteins anvi'o left out of gene clusters, and
count one form genes against OrthoFinder and the ORF flags correctly.
"""
import csv
import gzip
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "..", "workflow", "scripts")
sys.path.insert(0, SCRIPTS)

import anvio_inputs  # noqa: E402

FORMS = ["Cx_quinquefasciatus", "Cx_pallens", "Cx_molestus", "Cx_pipiens"]


def read(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def test_inputs():
    with tempfile.TemporaryDirectory() as d:
        proteins = {"t1": "MKV*", "t2": "M", "t3": "MKVL"}
        cds = {"t1": "ATGAAAGTTTAA", "t2": "ATG", "t3": "ATGAAAGTTC"}
        n, n_partial, n_padded = anvio_inputs.write_form("Cx_pallens", proteins, cds, d)
        assert (n, n_partial, n_padded) == (3, 1, 1)
        calls = read(os.path.join(d, "Cx_pallens.gene_calls.tsv"))
        assert [c["aa_sequence"] for c in calls] == ["MKV", "M", "MKVL"]   # stop removed
        assert [c["partial"] for c in calls] == ["0", "0", "1"]
        assert [c["stop"] for c in calls] == ["12", "3", "10"]            # the CDS only
        contigs = anvio_inputs.read_fasta(os.path.join(d, "Cx_pallens.contigs.fa"))
        assert contigs["pal_000002"] == "ATGN"                             # padded to k = 4
        gmap = read(os.path.join(d, "Cx_pallens.gene_map.tsv"))
        assert [(g["contig"], g["transcript"]) for g in gmap] == [
            ("pal_000001", "t1"), ("pal_000002", "t2"), ("pal_000003", "t3")]


def toy(d):
    """Four forms, gene A in all four (one gene cluster), gene B in all four
    but with the pipiens model broken and alone, gene C in qui only (its
    transfers missing), and a pallens fragment of gene D with no DIAMOND hit."""
    inputs = os.path.join(d, "inputs")
    os.makedirs(inputs)
    genes = {f: ["A", "B"] for f in FORMS}
    genes["Cx_quinquefasciatus"] += ["C"]
    genes["Cx_pallens"] += ["D"]
    table = []
    for f in FORMS:
        gm, gc = [], []
        for i, g in enumerate(genes[f]):
            tx = f"rna-{g}"
            gm.append([f"x_{i}", i, tx, 0])
            gc.append([i, f"x_{i}", 0, 30, "f", 0, 1, "workflow", "v5", "M" * (3 if g == "D" else 10)])
            valid = "" if f == "Cx_quinquefasciatus" else ("False" if (g, f) in {("B", "Cx_pipiens"),
                                                                                  ("D", "Cx_pallens")}
                                                            else "True")
            table.append([f, tx, f"gene-{g}", valid, "true" if g == "C" else "", ""])
        write(os.path.join(inputs, f"{f}.gene_map.tsv"),
              ["contig", "gene_callers_id", "transcript", "partial"], gm)
        write(os.path.join(inputs, f"{f}.gene_calls.tsv"),
              ["gene_callers_id", "contig", "start", "stop", "direction", "partial", "call_type",
               "source", "version", "aa_sequence"], gc)
    write(os.path.join(d, "genes.tsv"), ["sample", "transcript", "gene", "valid_ORF",
                                         "reference_partial", "reference_exception"], table)
    summary = os.path.join(d, "summary")
    os.makedirs(summary)
    rows, uid = [], 0
    clusters = {"GC_1": [(f, 0) for f in FORMS],                      # gene A
                "GC_2": [(f, 1) for f in FORMS[:3]],                  # gene B without pipiens
                "GC_3": [("Cx_pipiens", 1)],                          # broken B of pipiens
                "GC_4": [("Cx_quinquefasciatus", 2)]}                 # gene C
    for c, members in clusters.items():
        for f, i in members:
            uid += 1
            rows.append([uid, c, "", f, i])
    with gzip.open(os.path.join(summary, "Toy_gene_clusters_summary.txt.gz"), "wt") as fh:
        fh.write("unique_id\tgene_cluster_id\tbin_name\tgenome_name\tgene_callers_id\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    of = os.path.join(d, "of", "Results", "Orthogroups")
    os.makedirs(of)
    head = ["Orthogroup", "Cx_molestus", "Cx_pallens", "Cx_perexiguus", "Cx_pipiens",
            "Cx_quinquefasciatus"]
    write(os.path.join(of, "Orthogroups.tsv"), head,
          [["OG0", "rna-A", "rna-A", "p1", "rna-A", "rna-A"],
           ["OG1", "rna-B", "rna-B", "p2", "rna-B", "rna-B"],
           ["OG2", "", "", "p3", "", "rna-C"]])
    write(os.path.join(of, "Orthogroups_UnassignedGenes.tsv"), head,
          [["OG3", "", "rna-D", "", "", ""]])
    write(os.path.join(d, "partition.tsv"), ["Orthogroup", "compartment"],
          [["OG0", "core"], ["OG1", "core"], ["OG2", "cloud"]])
    return inputs, summary


def test_compare():
    with tempfile.TemporaryDirectory() as d:
        inputs, summary = toy(d)
        out = os.path.join(d, "out")
        os.makedirs(out)
        write(os.path.join(out, "versions.tsv"), ["tool", "version"], [["anvio", "eunice (v9)"]])
        subprocess.run([sys.executable, os.path.join(SCRIPTS, "anvio_compare.py"),
                        "--summary", summary, "--inputs", inputs,
                        "--orthofinder", os.path.join(d, "of"),
                        "--partition", os.path.join(d, "partition.tsv"),
                        "--annotation", os.path.join(d, "none"),
                        "--gene-table", os.path.join(d, "genes.tsv"), "--out", out],
                       check=True, capture_output=True)
        S = {(r["section"], r["item"]): r["value"] for r in read(os.path.join(out, "anvio_summary.tsv"))}
        assert S[("clusters", "input proteins")] == "10"
        assert S[("clusters", "proteins in gene clusters")] == "9"
        assert S[("not clustered", "proteins without a DIAMOND hit, left out of gene clusters")] == "1"
        assert S[("not clustered", "OrthoFinder unassigned")] == "1"
        assert S[("clusters", "clusters in all four forms")] == "1"
        assert S[("clusters", "clusters in three forms")] == "1"
        assert S[("clusters", "clusters in one form")] == "2"
        assert S[("versions", "anvio")] == "eunice (v9)"
        combos = {r["combination"]: r for r in read(os.path.join(out, "anvio_combinations.tsv"))}
        assert combos["qui + pal + mol"]["n_genes"] == "3"
        assert combos["pip"]["n_clusters"] == "1"
        xt = {r["anvio_class"]: r for r in read(os.path.join(out, "anvio_vs_orthofinder.tsv"))}
        assert (xt["all four"]["core"], xt["two or three"]["core"]) == ("4", "3")
        assert (xt["one form"]["core"], xt["one form"]["cloud"], xt["one form"]["total"]) == ("1", "1", "2")
        one = {r["form"]: r for r in read(os.path.join(out, "anvio_one_form.tsv"))}
        assert one["Cx_pipiens"]["same_gene_in_another_form"] == "1"      # B kept elsewhere
        assert one["Cx_pipiens"]["without_valid_orf"] == "1"
        assert one["Cx_quinquefasciatus"]["same_gene_in_another_form"] == "0"   # C nowhere else
        assert one["Cx_quinquefasciatus"]["refseq_model_not_clean"] == "1"
        assert one["Cx_quinquefasciatus"]["without_valid_orf"] == "NA"
        orf = {(r["form"], r["anvio_class"]): r for r in read(os.path.join(out, "anvio_orf_by_class.tsv"))}
        assert orf[("Cx_pipiens", "one form")]["pct_without_valid_orf"] == "100.0"
        assert orf[("Cx_pallens", "two or three")]["pct_without_valid_orf"] == "0.0"
        with gzip.open(os.path.join(out, "anvio_gene_clusters.tsv.gz"), "rt") as fh:
            per_gene = list(csv.DictReader(fh, delimiter="\t"))
        assert len(per_gene) == 10
        assert [r["anvio_class"] for r in per_gene if r["gene_cluster"] == ""] == ["not clustered"]


if __name__ == "__main__":
    test_inputs()
    test_compare()
    print("ok")
