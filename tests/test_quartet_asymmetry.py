"""Tests for workflow/scripts/quartet_asymmetry.py (run: python3 tests/test_quartet_asymmetry.py)."""
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
    print("quartet_asymmetry: all tests passed")


if __name__ == "__main__":
    sys.exit(main())
