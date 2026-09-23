"""Tests for divergence_diagnostics.py and transfer_quality.py (run: python3 tests/test_diagnostics.py)."""
import importlib.util
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, "..", "workflow", "scripts", f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


dd = load("divergence_diagnostics")
tq = load("transfer_quality")


def main():
    concat = "(A:0.04,(B:0.03,C:0.05)100:0.013,D:0.038);"
    genes = ["(A:0.004,(B:0.003,C:0.006)97:0.005,D:0.0000010000);",
             "(A:0.30,(B:0.002,D:0.006)40:0.002,C:0.005);",
             "(A:0.002,(C:0.004,D:0.004)88:0.001,B:0.001);"]
    rows, n = dd.branch_table(concat, genes)
    by = {r["branch"]: r for r in rows}
    assert n == 3 and by["A"]["concat_tree"] == 0.04
    assert by["A"]["gene_tree_median"] == 0.004 and by["A"]["share_above_0.1"] == round(1 / 3, 4)
    assert by["D"]["share_minimum"] == round(1 / 3, 4)
    assert by["internal"]["concat_tree"] == 0.013 and by["internal"]["n_gene_trees"] == 3

    with tempfile.TemporaryDirectory() as tmp:
        with open(os.path.join(tmp, "OG1.trim"), "w") as fh:
            fh.write(">A\nMKLV-A\n>B\nMKLVQA\n>C\nMRLVQA\n")
        with open(os.path.join(tmp, "OG2.trim"), "w") as fh:
            fh.write(">A\nMKK\n>B\nMKK\n>C\nMKK\n")
        pair_rows, invariant, few = dd.pairwise_table(sorted(
            os.path.join(tmp, f) for f in os.listdir(tmp)))
    ab = next(r for r in pair_rows if (r["taxon_a"], r["taxon_b"]) == ("A", "B"))
    assert ab["n_loci"] == 2 and ab["median_identity"] == 1.0
    ac = next(r for r in pair_rows if (r["taxon_a"], r["taxon_b"]) == ("A", "C"))
    assert ac["share_below_0.95"] == 0.5          # OG1: 4 of 5 shared columns
    assert invariant == {"OG2"} and few == {"OG1", "OG2"}

    raw = {"t1": "MKLV.", "t2": "MK.LV.", "t3": "KLV.", "t4": "MKLV", "t5": "MAA."}
    r = tq.assess(raw, ["t1", "t2", "t3", "t4", "t9"])
    assert r["n_models"] == 4 and r["n_internal_stop"] == 1
    assert r["n_no_start_met"] == 1 and r["n_no_terminal_stop"] == 1
    assert r["n_complete"] == 1 and r["pct_complete"] == 25.0
    print("diagnostics: all tests passed")


if __name__ == "__main__":
    sys.exit(main())
