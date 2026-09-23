"""Tests for workflow/scripts/cafe_transfer_bias.py (run: python3 tests/test_cafe_transfer_bias.py)."""
import csv
import os
import runpy
import sys
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "workflow", "scripts", "cafe_transfer_bias.py")


def main():
    # reference R, transferred T1 and T2; 20 single copy families (1 significant),
    # 10 multi copy (5 significant), 4 with no reference copy (2 significant)
    rows = [("Family ID", "Desc", "R", "T1", "T2", "O")]
    sig = set()
    for i in range(20):
        rows.append((f"S{i}", f"S{i}", 1, 1, 1 if i else 0, 1))
    sig.add("S0")
    for i in range(10):
        rows.append((f"M{i}", f"M{i}", 3, 2, 1, 1))
        if i < 5:
            sig.add(f"M{i}")
    for i in range(4):
        rows.append((f"Z{i}", f"Z{i}", 0, 1, 0, 1))
        if i < 2:
            sig.add(f"Z{i}")
    with tempfile.TemporaryDirectory() as tmp:
        counts = os.path.join(tmp, "counts.tsv")
        with open(counts, "w") as fh:
            for r in rows:
                fh.write("\t".join(map(str, r)) + "\n")
        sigf = os.path.join(tmp, "sig.tsv")
        with open(sigf, "w") as fh:
            fh.write("orthogroup\tpvalue\n" + "".join(f"{s}\t0.01\n" for s in sorted(sig)))
        branch = os.path.join(tmp, "branch.tsv")
        with open(branch, "w") as fh:
            fh.write("taxon\tnode_label\tincrease\tdecrease\n"
                     "R\tR<1>\t5\t5\nT1\tT1<2>\t1\t8\nT2\tT2<3>\t0\t12\nT1+T2\t<5>\t3\t0\n")
        sm = types.SimpleNamespace(
            input=types.SimpleNamespace(counts=counts, sig=sigf, branch=branch),
            params=types.SimpleNamespace(reference="R", ingroup=["R", "T1", "T2"]),
            output=types.SimpleNamespace(summary=os.path.join(tmp, "summary.tsv"),
                                         bins=os.path.join(tmp, "bins.tsv"),
                                         lineage=os.path.join(tmp, "lineage.tsv")))
        runpy.run_path(SCRIPT, init_globals={"snakemake": sm}, run_name="__main__")
        summary = {r["metric"]: r["value"] for r in
                   csv.DictReader(open(sm.output.summary), delimiter="\t")}
        bins = {r["ref_copies"]: r for r in csv.DictReader(open(sm.output.bins), delimiter="\t")}
        lineage = list(csv.DictReader(open(sm.output.lineage), delimiter="\t"))
    assert bins["0"]["n_families"] == "4" and bins["0"]["n_significant"] == "2"
    assert bins["0"]["mean_retention"] in ("", "nan", "NaN")
    assert summary["n_families_single_copy"] == "20" and summary["n_families_multi_copy"] == "10"
    assert summary["n_significant_zero_ref_copy"] == "2"
    # strictly single against multi: (5/5) / (1/19) = 19
    assert abs(float(summary["odds_ratio_multi_vs_single_copy"]) - 19.0) < 1e-6
    assert abs(float(summary["risk_ratio_multi_vs_single_copy"]) - 10.0) < 1e-6
    assert {r["taxon"] for r in lineage} == {"R", "T1", "T2"}
    print("cafe_transfer_bias: all tests passed")


if __name__ == "__main__":
    sys.exit(main())
