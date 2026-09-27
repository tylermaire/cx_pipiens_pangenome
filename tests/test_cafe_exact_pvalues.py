"""Tests for workflow/scripts/cafe_exact_pvalues.py
(run: python3 tests/test_cafe_exact_pvalues.py).

The exact P value must equal (1) a brute force sum over every configuration
of counts and (2) the share of simulated families at most as likely as the
observed one, which is how CAFE5 estimates it.
"""
import itertools
import os
import runpy
import sys
import tempfile
import types

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "..", "workflow", "scripts")
sys.path.insert(0, SCRIPTS)

import cafe_exact_pvalues as cx  # noqa: E402

LAM = 0.045
TREE = "((Cx_molestus:1.0,Cx_pipiens:1.0):0.5,(Cx_pallens:1.0,Cx_quinquefasciatus:1.0):0.5);"


def test_transition():
    P = cx.transition(LAM, 1.0, 60)
    a = LAM / (1 + LAM)
    assert abs(P[1, 0] - a) < 1e-12
    assert abs(P[1, 1] - (1 - a) ** 2) < 1e-12          # (1 - a)^2 for one copy
    assert np.allclose(P[:20].sum(axis=1), 1.0, atol=1e-9)
    assert P[0, 0] == 1.0 and P[0, 1:].sum() == 0.0


def test_tree():
    pairs = cx.two_pairs(TREE)
    assert [(p[0][0], p[1][0], p[2]) for p in pairs] == [
        ("Cx_molestus", "Cx_pipiens", 0.5), ("Cx_pallens", "Cx_quinquefasciatus", 0.5)]
    for bad in ("(((A:1,B:1):0.5,C:1.5):0.5,D:2);", "((A:1,B:1,C:1):1,D:2);"):
        try:
            cx.two_pairs(bad)
        except ValueError:
            continue
        raise AssertionError(f"accepted {bad}")


def joint(n, s, t_int=0.5, t_tip=1.0):
    """P(counts of the four leaves | root size s), by brute force."""
    Ti, Tt = cx.transition(LAM, t_int, n), cx.transition(LAM, t_tip, n)
    pair = np.einsum("a,ax,ay->xy", Ti[s], Tt, Tt)
    return np.einsum("xy,zw->xyzw", pair, pair)


def test_against_enumeration():
    n = 24
    d = cx.PairDistribution(LAM, 0.5, 1.0, 1.0, n)
    for s in (1, 2, 3):
        J = joint(n, s)
        for counts in ((0, 1, 0, 1), (1, 1, 1, 1), (2, 1, 1, 1), (0, 0, 1, 1), (3, 1, 0, 2)):
            brute = J[J <= J[counts] * (1 + 1e-9)].sum()
            fast = cx.p_value_at_root(d, d, counts, s)
            assert abs(brute - fast) < 1e-9, (s, counts, brute, fast)


def test_against_simulation():
    """CAFE5's estimator, with 200,000 families instead of 1,000."""
    rng = np.random.default_rng(1)
    n = 30
    Ti, Tt = cx.transition(LAM, 0.5, n), cx.transition(LAM, 1.0, n)
    d = cx.PairDistribution(LAM, 0.5, 1.0, 1.0, n)
    F = d.matrix(1)[0]
    k = 200_000

    def draw(row_p, size):
        p = row_p / row_p.sum()
        return rng.choice(len(p), size=size, p=p)

    lik = np.ones(k)
    for _ in range(2):                         # the two pairs below the root
        a = draw(Ti[1], k)                     # size at the pair's ancestor
        x = np.empty(k, dtype=int)
        y = np.empty(k, dtype=int)
        for v in np.unique(a):
            idx = np.where(a == v)[0]
            x[idx] = draw(Tt[v], len(idx))
            y[idx] = draw(Tt[v], len(idx))
        lik *= F[x, y]
    for counts in ((0, 1, 0, 1), (1, 1, 1, 0)):
        exact = cx.p_value_at_root(d, d, counts, 1)
        obs = F[counts[0], counts[1]] * F[counts[2], counts[3]]
        est = np.mean(lik <= obs * (1 + 1e-9))
        se = np.sqrt(exact * (1 - exact) / k)
        assert abs(est - exact) < 4 * se, (counts, est, exact, se)


def test_run():
    with tempfile.TemporaryDirectory() as tmp:
        cafe = os.path.join(tmp, "output")
        os.makedirs(cafe)
        open(os.path.join(cafe, "Gamma_results.txt"), "w").write(
            f"Model Gamma Final Likelihood (-lnL): 100\nLambda: {LAM}\nAlpha: 0.3\n")
        species = ["Cx_quinquefasciatus", "Cx_pallens", "Cx_molestus", "Cx_pipiens"]
        fams = {"OG1": (1, 1, 1, 1), "OG2": (1, 0, 0, 1), "OG3": (0, 1, 0, 1),
                "OG4": (5, 1, 1, 1), "OG5": (0, 0, 1, 1)}   # OG5: not at the root
        with open(os.path.join(tmp, "counts.tsv"), "w") as fh:
            fh.write("Desc\tFamily ID\t" + "\t".join(species) + "\n")
            for f, c in fams.items():
                fh.write(f"{f}\t{f}\t" + "\t".join(map(str, c)) + "\n")
        with open(os.path.join(cafe, "Gamma_family_results.txt"), "w") as fh:
            fh.write("#FamilyID\tpvalue\tSignificant at 0.05\n")
            for f, p in (("OG1", 0.999), ("OG2", 0.042), ("OG3", 0.042), ("OG4", 0.0)):
                fh.write(f"{f}\t{p}\t{'y' if p < 0.05 else 'n'}\n")
        open(os.path.join(tmp, "tree.nwk"), "w").write(TREE + "\n")
        out = os.path.join(tmp, "family_pvalues.tsv")
        sm = types.SimpleNamespace(
            input=types.SimpleNamespace(counts=os.path.join(tmp, "counts.tsv"),
                                        tree=os.path.join(tmp, "tree.nwk"), cafe_dir=cafe),
            output=types.SimpleNamespace(table=out))
        runpy.run_path(os.path.join(SCRIPTS, "cafe_exact_pvalues.py"),
                       init_globals={"snakemake": sm}, run_name="__main__")
        rows = [l.rstrip("\n").split("\t") for l in open(out)]
    head, body = rows[0], {r[0]: r for r in rows[1:]}
    assert head == ["orthogroup", "Cx_molestus", "Cx_pipiens", "Cx_pallens",
                    "Cx_quinquefasciatus", "cafe_pvalue", "exact_pvalue"], head
    assert sorted(body) == ["OG1", "OG2", "OG3", "OG4"]          # tested families only
    p = {f: float(r[-1]) for f, r in body.items()}
    assert p["OG1"] > 0.99
    assert abs(p["OG2"] - p["OG3"]) < 1e-12                      # the same pattern class
    n = 5 + max(50, 5 // 5)
    d = cx.PairDistribution(LAM, 0.5, 1.0, 1.0, n)
    assert abs(p["OG3"] - cx.p_value_at_root(d, d, (0, 1, 1, 0), 1)) < 1e-6
    # largest count 5: root sizes 1 to rint(6.25) = 6, the largest P kept
    best = max(cx.p_value_at_root(d, d, (1, 1, 1, 5), s) for s in range(1, 7))
    assert abs(p["OG4"] - best) < 1e-6
    assert p["OG4"] < 0.01


def main():
    test_transition()
    test_tree()
    test_against_enumeration()
    test_against_simulation()
    test_run()
    print("cafe_exact_pvalues: all tests passed")


if __name__ == "__main__":
    sys.exit(main())
