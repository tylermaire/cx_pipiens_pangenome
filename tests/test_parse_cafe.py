"""Tests for workflow/scripts/parse_cafe.py (run: python3 tests/test_parse_cafe.py)."""
import csv
import os
import runpy
import sys
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "workflow", "scripts", "parse_cafe.py")

FAMILY = "#FamilyID\tpvalue\tSignificant at 0.05\nOG1\t0.001\ty\nOG2\t0.5\tn\nOG3\t0.04\ty\n"
CLADE = ("#Taxon_ID\tIncrease\tDecrease\n<6>\t10\t0\nA<1>\t2\t8\nB<2>\t5\t5\n"
         "<7>\t3\t0\nC<3>\t1\t9\nD<4>\t4\t6\n")
ASR = ("#nexus\nBEGIN TREES;\n"
       "  TREE OG1 = ((A<1>_2:1,B<2>*_3:1)<6>*_2:0.5,(C<3>_1:1,D<4>_1:1)<7>_1:0.5)<5>_2;\n"
       "END;\n")


def main():
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "cafe")
        os.makedirs(out)
        for name, text in (("Gamma_family_results.txt", FAMILY),
                           ("Gamma_clade_results.txt", CLADE), ("Gamma_asr.tre", ASR)):
            open(os.path.join(out, name), "w").write(text)
        sm = types.SimpleNamespace(
            input=[out], params=types.SimpleNamespace(pvalue=0.05),
            output=types.SimpleNamespace(significant=os.path.join(tmp, "sig.tsv"),
                                         summary=os.path.join(tmp, "branch.tsv")))
        runpy.run_path(SCRIPT, init_globals={"snakemake": sm}, run_name="__main__")
        sig = list(csv.DictReader(open(sm.output.significant), delimiter="\t"))
        branch = list(csv.DictReader(open(sm.output.summary), delimiter="\t"))
    assert [r["orthogroup"] for r in sig] == ["OG1", "OG3"]
    names = {r["node_label"]: r["taxon"] for r in branch}
    assert names["<6>"] == "A+B" and names["<7>"] == "C+D", names
    assert names["A<1>"] == "A" and all(r["taxon"] for r in branch)
    test_exact()
    print("parse_cafe: all tests passed")


def test_exact():
    """With params.source exact, the exact P decides; CAFE5's estimate is kept."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "cafe")
        os.makedirs(out)
        for name, text in (("Gamma_family_results.txt", FAMILY),
                           ("Gamma_clade_results.txt", CLADE), ("Gamma_asr.tre", ASR)):
            open(os.path.join(out, name), "w").write(text)
        table = os.path.join(tmp, "family_pvalues.tsv")
        open(table, "w").write("orthogroup\tA\tB\tC\tD\tcafe_pvalue\texact_pvalue\n"
                               "OG1\t1\t1\t1\t1\t0.001\t0.0012\n"
                               "OG2\t1\t1\t1\t1\t0.5\t0.03\n"
                               "OG3\t1\t1\t1\t1\t0.04\t0.0555\n")
        sm = types.SimpleNamespace(
            input=types.SimpleNamespace(cafe_dir=out, pvalues=[table]),
            params=types.SimpleNamespace(pvalue=0.05, source="exact"),
            output=types.SimpleNamespace(significant=os.path.join(tmp, "sig.tsv"),
                                         summary=os.path.join(tmp, "branch.tsv")))
        runpy.run_path(SCRIPT, init_globals={"snakemake": sm}, run_name="__main__")
        sig = list(csv.DictReader(open(sm.output.significant), delimiter="\t"))
        assert [r["orthogroup"] for r in sig] == ["OG1", "OG2"], sig
        assert set(sig[0]) == {"orthogroup", "pvalue", "cafe_pvalue", "exact_pvalue"}, sig[0]
        assert float(sig[1]["pvalue"]) == 0.03 and float(sig[1]["cafe_pvalue"]) == 0.5
        # the two tables must list the same families
        open(table, "a").write("OG9\t1\t1\t1\t1\t0.5\t0.5\n")
        try:
            runpy.run_path(SCRIPT, init_globals={"snakemake": sm}, run_name="__main__")
        except SystemExit as err:
            assert "different families" in str(err), err
        else:
            raise AssertionError("a family missing from CAFE5's results was accepted")


if __name__ == "__main__":
    sys.exit(main())
