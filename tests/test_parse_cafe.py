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
    print("parse_cafe: all tests passed")


if __name__ == "__main__":
    sys.exit(main())
