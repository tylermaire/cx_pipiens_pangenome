"""Smoke test for collect_manuscript_values.py against the tracked results
(run from the repository root: python3 tests/test_collect_manuscript_values.py).

The tracked results are those of the last run, so the outgroup is read from
them (the QUAST column that is not an ingroup form) rather than assumed, and
values are checked for form rather than for the numbers of one run."""
import csv
import os
import runpy
import sys
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
SCRIPT = os.path.join(ROOT, "workflow", "scripts", "collect_manuscript_values.py")
INGROUP = ["Cx_quinquefasciatus", "Cx_pallens", "Cx_molestus", "Cx_pipiens"]


def outgroup_in_results():
    with open(os.path.join(ROOT, "results", "quast", "report.tsv")) as fh:
        head = fh.readline().rstrip("\n").split("\t")
    others = [c for c in head[1:] if c not in INGROUP]
    assert len(others) == 1, head
    return others[0]


def main():
    os.chdir(ROOT)
    og = outgroup_in_results()
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "values.tsv")
        sm = types.SimpleNamespace(
            params=types.SimpleNamespace(
                samples=INGROUP + [og], ingroup=INGROUP,
                reference="Cx_quinquefasciatus", tools=["iqtree", "cafe"],
                parameters={"busco lineage": "diptera_odb10"}),
            output=[out])
        runpy.run_path(SCRIPT, init_globals={"snakemake": sm}, run_name="__main__")
        rows = list(csv.DictReader(open(out), delimiter="\t"))
    get = lambda sec, item, sample="": next(
        r["value"] for r in rows if (r["section"], r["item"], r["sample"]) == (sec, item, sample))
    assert int(get("assembly", "total_length_bp", og)) > 100_000_000
    assert get("busco", "genome busco_version", og) not in ("", "NA")
    assert int(get("phylogeny", "gN")) > 1000
    assert int(get("pangenome", "cloud n_orthogroups")) > 0
    assert get("parameters", "busco lineage") == "diptera_odb10"
    assert get("tools", "conda environments").startswith("not present")
    # the outgroup analyses are collected when present and flagged NA when not
    rooted = [r for r in rows if r["section"] == "rooted"]
    assert rooted, "no rooted section"
    if not os.path.exists(os.path.join(ROOT, "results", "phylo", "rooted", "rooted_summary.tsv")):
        assert rooted[0]["value"] == "NA"
    print(f"collect_manuscript_values: smoke test passed ({len(rows)} rows, outgroup {og})")


if __name__ == "__main__":
    sys.exit(main())
