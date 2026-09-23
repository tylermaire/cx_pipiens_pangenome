"""Smoke test for collect_manuscript_values.py against the tracked results
(run from the repository root: python3 tests/test_collect_manuscript_values.py)."""
import csv
import os
import runpy
import sys
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
SCRIPT = os.path.join(ROOT, "workflow", "scripts", "collect_manuscript_values.py")


def main():
    os.chdir(ROOT)
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "values.tsv")
        sm = types.SimpleNamespace(
            params=types.SimpleNamespace(
                samples=["Cx_quinquefasciatus", "Cx_pallens", "Cx_molestus",
                         "Cx_pipiens", "Cx_tarsalis"],
                ingroup=["Cx_quinquefasciatus", "Cx_pallens", "Cx_molestus", "Cx_pipiens"],
                reference="Cx_quinquefasciatus", tools=["iqtree", "cafe"],
                parameters={"busco lineage": "diptera_odb10"}),
            output=[out])
        runpy.run_path(SCRIPT, init_globals={"snakemake": sm}, run_name="__main__")
        rows = list(csv.DictReader(open(out), delimiter="\t"))
    get = lambda sec, item, sample="": next(
        r["value"] for r in rows if (r["section"], r["item"], r["sample"]) == (sec, item, sample))
    assert get("assembly", "total_length_bp", "Cx_tarsalis") == "789668571"
    assert get("busco", "genome busco_version", "Cx_tarsalis") == "6.0.0"
    assert get("phylogeny", "gN") == "7646"
    assert get("pangenome", "cloud n_orthogroups") == "479"
    assert get("parameters", "busco lineage") == "diptera_odb10"
    assert get("tools", "iqtree") == "NA"          # no conda envs in a fresh clone
    print(f"collect_manuscript_values: smoke test passed ({len(rows)} rows)")


if __name__ == "__main__":
    sys.exit(main())
