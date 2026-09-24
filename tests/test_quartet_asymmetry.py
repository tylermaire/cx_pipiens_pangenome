"""Tests for workflow/scripts/quartet_asymmetry.py (run: python3 tests/test_quartet_asymmetry.py)."""
import collections
import importlib.util
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "quartet_asymmetry", os.path.join(HERE, "..", "workflow", "scripts", "quartet_asymmetry.py"))
qa = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qa)

TAXA = frozenset(["A", "B", "C", "D"])


def main():
    # IQ-TREE gene tree and cf.tree label styles
    taxa, pair, sup = qa.parse_tree("(A:0.01,(B:0.02,C:0.03)97:0.004,D:0.05);")
    assert taxa == TAXA and pair == frozenset("BC") and sup == 97.0
    _, pair, sup = qa.parse_tree("(A:0.01,(B:0.02,C:0.03)100/59/47.6:0.004,D:0.05);")
    assert pair == frozenset("BC") and sup == 100.0
    _, pair, sup = qa.parse_tree("((A,D),B,C);")
    assert pair == frozenset("AD") and sup is None
    assert qa.canonical(frozenset("BC"), TAXA) == frozenset("AD")

    # species tree AB|CD; AC|BD in excess over AD|BC; supports vary
    trees = ([(frozenset("AB"), 100)] * 50 + [(frozenset("AC"), 99)] * 30 +
             [(frozenset("BD"), 40)] * 10 + [(frozenset("AD"), 99)] * 10)
    cf = {"gDF1_N": "40", "gDF2_N": "10"}
    topo, support, _ = qa.analyse(trees, TAXA, ["A", "B"], cf)
    roles = {r["role"]: r for r in topo}
    assert roles["species_tree"]["n_gene_trees"] == 50
    assert roles["major_discordant"]["split"] == "A+C | B+D"
    assert roles["major_discordant"]["n_gene_trees"] == 40
    assert roles["major_discordant"]["iqtree_label"] == "gDF1"
    assert roles["minor_discordant"]["iqtree_label"] == "gDF2"
    s95 = next(r for r in support if r["min_ufboot"] == 95)
    assert (s95["n_major"], s95["n_minor"], s95["n_concordant"]) == (30, 10, 50)

    # site patterns: AABB supports A+B, ABAB supports A+C, gaps and 3 state columns ignored
    with tempfile.TemporaryDirectory() as tmp:
        with open(os.path.join(tmp, "OG1.trim"), "w") as fh:
            fh.write(">A\nKKLMW-\n>B\nKLKMWK\n>C\nRKKMRK\n>D\nRLLQRK\n")
        counts, n_files, n_cols = qa.site_patterns(tmp, TAXA)
    assert n_files == 1 and n_cols == 6
    assert counts[frozenset("AB")] == 2       # cols 1 and 5: K K R R, W W R R
    assert counts[frozenset("AC")] == 1       # col 2: K L K L
    assert counts[frozenset("AD")] == 1       # col 3: L K K L
    assert sum(counts.values()) == 4          # col 4 has three of a kind, col 6 a gap

    # one vote per locus: OG1 votes A+B, OG2 votes A+C, OG3 ties, OG4 has no
    # informative sites; the top 1% (one locus) is OG1 with 4 of 9 sites
    with tempfile.TemporaryDirectory() as tmp:
        files = {"OG1": ">A\nKKLMW-\n>B\nKLKMWK\n>C\nRKKMRK\n>D\nRLLQRK\n",
                 "OG2": ">A\nKKK\n>B\nLLL\n>C\nKKK\n>D\nLLL\n",
                 "OG3": ">A\nLK\n>B\nKL\n>C\nKK\n>D\nLL\n",
                 "OG4": ">A\nKK\n>B\nKK\n>C\nKK\n>D\nKK\n"}
        for name, text in files.items():
            with open(os.path.join(tmp, f"{name}.trim"), "w") as fh:
                fh.write(text)
        per_locus = []
        counts, n_files, _ = qa.site_patterns(tmp, TAXA, per_locus)
    assert n_files == 4 and len(per_locus) == 4
    _, _, sites = qa.analyse(trees, TAXA, ["A", "B"], cf, counts, per_locus)
    by_role = {r["role"]: r for r in sites}
    assert by_role["locus_majority_species_tree"]["n_loci"] == 1
    assert by_role["locus_majority_major_discordant_gene_trees"]["n_loci"] == 1
    assert by_role["locus_majority_minor_discordant_gene_trees"]["n_loci"] == 0
    assert by_role["locus_majority_undecided"]["n_loci"] == 2
    assert by_role["locus_majority_undecided"]["n_loci_no_informative_sites"] == 1
    assert by_role["median"]["n_informative_sites"] == 2.5
    assert by_role["concentration"]["n_loci"] == 1
    assert by_role["concentration"]["n_informative_sites"] == 4
    assert by_role["concentration"]["pct"] == round(100 * 4 / 9, 2)
    assert by_role["top_loci_species_tree"]["n_informative_sites"] == 2
    # outside the top locus: OG2 gives A+C 3 sites, OG3 gives A+C 1 and A+D 1
    rest = next(r for r in sites if r["split"].startswith("discordant sites outside"))
    assert rest["n_informative_sites"] == 5 and rest["pct"] == 80.0

    # internal branch length follows the cherry
    assert qa.internal_length("(A:0.01,(B:0.02,C:0.03)97:0.004,D:0.05);") == 0.004
    assert qa.internal_length("(A:0.01,(B:0.02,C:0.03)33:0.0000010000,D:0.05);") <= qa.MIN_BRANCH

    # robustness: loci L1 to L6; L5 has an unresolved internal branch, L6 a broken model
    by_locus = [("L1", frozenset("AB"), 100, 0.01), ("L2", frozenset("AC"), 99, 0.02),
                ("L3", frozenset("AC"), 96, 0.01), ("L4", frozenset("AD"), 97, 0.01),
                ("L5", frozenset("AD"), 30, 1e-6), ("L6", frozenset("AC"), 99, 0.01)]
    intact = {"L1": True, "L2": True, "L3": True, "L4": True, "L5": True, "L6": False}
    splits = (qa.canonical(frozenset("AB"), TAXA), qa.canonical(frozenset("AC"), TAXA),
              qa.canonical(frozenset("AD"), TAXA))
    rows = {r["subset"]: r for r in qa.robustness_rows(by_locus, TAXA, splits, intact)}
    r = rows["all gene trees"]
    assert (r["n_gene_trees"], r["n_concordant"], r["n_major"], r["n_minor"]) == (6, 1, 3, 2)
    r = rows["internal branch above the minimum length"]
    assert (r["n_gene_trees"], r["n_major"], r["n_minor"]) == (5, 3, 1)
    assert rows["internal branch at the minimum length"]["n_minor"] == 1
    r = rows["loci with four intact models"]
    assert (r["n_gene_trees"], r["n_major"], r["n_minor"]) == (5, 2, 2)
    assert rows["loci with four intact models, UFBoot >= 95"]["n_gene_trees"] == 4
    assert rows["loci with a model that is not intact"]["n_major"] == 1

    # intact loci: reference R, transferred T1 to T3; OG2 has a transferred model
    # without a valid ORF, OG3 a partial reference source, OG4 is not single copy
    with tempfile.TemporaryDirectory() as tmp:
        of = os.path.join(tmp, "Results_X", "Orthogroups")
        os.makedirs(of)
        with open(os.path.join(of, "Orthogroups.tsv"), "w") as fh:
            fh.write("Orthogroup\tR\tT1\tT2\tT3\n"
                     "OG1\trna-a1\trna-a1\trna-a1\trna-a1\n"
                     "OG2\trna-b1\trna-b1\trna-b1\trna-b1\n"
                     "OG3\trna-c1\trna-c1\trna-c1\trna-c1\n"
                     "OG4\trna-d1, rna-e1\trna-d1\trna-d1\trna-d1\n")
        gffs = {}
        for form in ("R", "T1", "T2", "T3"):
            path = os.path.join(tmp, f"{form}_liftoff.gff3")
            with open(path, "w") as fh:
                for tid in ("a.1", "b.1", "c.1", "d.1"):
                    attrs = f"ID=rna-{tid};Parent=gene-{tid[0]}"
                    if form == "R" and tid == "c.1":
                        attrs += ";partial=true"
                    if form != "R":
                        valid = "False" if (form == "T2" and tid == "b.1") else "True"
                        attrs += f";valid_ORF={valid}"
                    fh.write(f"chr1\tx\tmRNA\t1\t100\t.\t+\t.\t{attrs}\n")
            gffs[form] = path
        got = qa.intact_loci(tmp, gffs, "R")
    assert got == {"OG1": True, "OG2": False, "OG3": False}, got

    # locus quality: the top locus by informative sites is the broken one
    per_locus = [collections.Counter({splits[0]: 9}), collections.Counter({splits[1]: 1}),
                 collections.Counter({splits[2]: 1})]
    rows = qa.locus_quality_rows(per_locus, ["OG2", "OG1", "OG3"], splits,
                                 {"OG1": True, "OG2": False, "OG3": False}, top_share=0.34)
    assert (rows[0]["n_loci"], rows[0]["n_loci_not_intact"]) == (1, 1)
    assert (rows[1]["n_loci"], rows[1]["n_loci_not_intact"]) == (2, 1)
    print("quartet_asymmetry: all tests passed")


if __name__ == "__main__":
    sys.exit(main())
