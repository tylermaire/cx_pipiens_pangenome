"""Tests for workflow/scripts/synteny_summary.py (run: python3 tests/test_synteny_summary.py)."""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "synteny_summary", os.path.join(HERE, "..", "workflow", "scripts", "synteny_summary.py"))
ss = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ss)


def build():
    """Reference R with chromosomes NC_2 (genes a0..a19) and NC_10 (b0..b19).
    A and B both carry the a5..a10 block inverted; B's second chromosome is
    deposited reverse complemented."""
    step = 200_000
    R = {f"a{i}": ("NC_2", i * step) for i in range(20)}
    R.update({f"b{i}": ("NC_10", i * step) for i in range(20)})
    inv = list(range(5, 11))
    def with_inversion(chrom_a, chrom_b, reverse_b=False):
        g = {}
        for i in range(20):
            j = (inv[0] + inv[-1] - i) if i in inv else i
            g[f"a{i}"] = (chrom_a, j * step)
        for i in range(20):
            pos = (19 - i) * step if reverse_b else i * step
            g[f"b{i}"] = (chrom_b, pos)
        return g
    return {"R": R, "A": with_inversion("A1", "A2"),
            "B": with_inversion("B1", "B2", reverse_b=True)}


def main():
    genes = build()
    homolog = ss.homolog_names(genes, ["R", "A", "B"], "R", 2)
    # natural order: NC_2 before NC_10
    assert homolog[("R", "NC_2")] == "chr1" and homolog[("R", "NC_10")] == "chr2"
    assert homolog[("A", "A1")] == "chr1" and homolog[("B", "B2")] == "chr2"

    n, pct, inv_ra, orients = ss.analyse_pair("R", "A", genes, 2, 100_000, 3, homolog)
    assert n == 40 and len(inv_ra) == 1
    d = inv_ra[0]
    assert d["chromosome"] == "chr1" and d["n_genes"] == 6
    assert (d["first_gene"], d["last_gene"]) == ("a5", "a10")
    assert d["rel_start"] == round(5 / 19, 3)

    _, _, inv_rb, orients_b = ss.analyse_pair("R", "B", genes, 2, 100_000, 3, homolog)
    assert len(inv_rb) == 1, "orientation normalisation must not report chr2 as inverted"
    assert sum(o["orientation"] == "reverse" for o in orients_b) == 1

    _, pct_ab, inv_ab, _ = ss.analyse_pair("A", "B", genes, 2, 100_000, 3, homolog)
    assert inv_ab == [] and pct_ab == 100.0

    all_inv = inv_ra + inv_rb
    clusters = ss.recurrence_clusters(all_inv)
    assert clusters[0] == clusters[1], "the shared inversion must form one cluster"
    rows = ss.recurrence_table(all_inv, clusters, 1_000_000)
    assert len(rows) == 1 and rows[0]["n_pairs"] == 2 and rows[0]["any_at_least_large_span"]

    # a gene placed on the non-homologous chromosome is counted apart
    assert ss.count_other_chromosome(genes["R"], genes["A"], 2) == (40, 0)
    moved = dict(genes["A"])
    moved["a3"] = ("A2", 50 * 200_000)
    assert ss.count_other_chromosome(genes["R"], moved, 2) == (40, 1)
    print("synteny_summary: all tests passed")


if __name__ == "__main__":
    sys.exit(main())
