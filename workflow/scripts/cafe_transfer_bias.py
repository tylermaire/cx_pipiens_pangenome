#!/usr/bin/env python3
"""
cafe_transfer_bias.py - control analysis for the CAFE gene family results.

Every gene model outside the reference genome was transferred by Liftoff from
the Cx. quinquefasciatus annotation. Transfer is not uniform: a family with
several near-identical paralogues gives the aligner several chances to place a
model ambiguously, so copies are lost in proportion to how many there are.
CAFE sees that loss as contraction and cannot distinguish it from gene loss.

This script quantifies the bias from files the workflow already produces, and
writes the numbers the manuscript cites rather than leaving them asserted.

Three measurements:

  retention    mean copies in the three transferred genomes divided by copies
               in the reference, binned by reference copy number. A flat curve
               near 1.0 would mean transfer is unbiased with respect to family
               size; a declining curve means it is not.

  enrichment   whether families CAFE calls significantly evolving are the same
               families that lose copies in transfer. Reported as the
               significance rate per copy-number bin, two sided Mann-Whitney
               tests on reference copy number and on transfer deficit (families
               with at least one reference copy), and a two sided Fisher test
               with odds ratio and risk ratio for multi-copy (2 or more
               reference copies) against strictly single-copy families.

  zero copy    families with no reference copy get their own bin. Every
               transferred model comes from a reference gene, so these families
               exist only because transferred copies clustered apart from their
               source; their significance rate is reported, not pooled with
               single-copy families. (Pooling them is what produced the
               earlier odds ratio of 32.4, which was labelled single-copy.)

  per-lineage  net CAFE contraction against each genome's total copy-number
               deficit relative to the reference. With four ingroup genomes
               this is four points and is reported as supporting rather than
               primary evidence. The per-lineage increase and decrease counts
               come from Gamma_clade_results and cover ALL families, not only
               the significant ones.

Snakemake provides:
    input.counts   results/cafe/gene_counts_filtered.tsv
    input.sig      results/cafe/significant_families.tsv
    input.branch   results/cafe/branch_summary.tsv
    params.reference, params.ingroup
    output.summary, output.bins, output.lineage
"""

import sys

import numpy as np
import pandas as pd

try:
    from scipy.stats import mannwhitneyu, fisher_exact, pearsonr
except ImportError:
    sys.exit("scipy required; add it to the rule's conda environment")

BINS = [(0, 0), (1, 1), (2, 2), (3, 5), (6, 10), (11, 10_000)]


def fmt_p(p):
    """Format a p value; one that underflows to 0 is reported as a bound."""
    return "< 1e-300" if p == 0 else f"{p:.3g}"


def bin_label(lo, hi):
    return f"{lo}" if lo == hi else (f"{lo}+" if hi >= 10_000 else f"{lo}-{hi}")


def main():
    if "snakemake" not in globals():
        sys.exit("Run via Snakemake.")

    ref = snakemake.params.reference
    ingroup = list(snakemake.params.ingroup)
    lifted = [s for s in ingroup if s != ref]

    d = pd.read_csv(snakemake.input.counts, sep="\t")
    if "Family ID" in d.columns:
        d = d.rename(columns={"Family ID": "orthogroup"})
    sig = pd.read_csv(snakemake.input.sig, sep="\t")
    idcol = "orthogroup" if "orthogroup" in sig.columns else sig.columns[0]
    sig_ids = set(sig[idcol].astype(str))

    anon = sum(1 for x in sig_ids if str(x).strip() in ("", "n/a", "nan"))
    if anon:
        sys.exit(f"{anon}/{len(sig_ids)} significant families have no identifier; "
                 f"format_cafe_input.py must write the orthogroup id into Desc.")

    d["ref_copies"] = d[ref]
    d["mean_lifted"] = d[lifted].mean(axis=1)
    d["deficit"] = d["ref_copies"] - d["mean_lifted"]
    d["retention"] = np.where(d.ref_copies > 0, d.mean_lifted / d.ref_copies, np.nan)
    d["is_sig"] = d["orthogroup"].astype(str).isin(sig_ids)

    rows = []
    for lo, hi in BINS:
        m = (d.ref_copies >= lo) & (d.ref_copies <= hi)
        if not m.any():
            continue
        rows.append({
            "ref_copies": bin_label(lo, hi),
            "n_families": int(m.sum()),
            "mean_lifted_copies": round(float(d.loc[m, "mean_lifted"].mean()), 4),
            "mean_retention": round(float(d.loc[m, "retention"].mean()), 4)
                              if lo > 0 else float("nan"),
            "mean_deficit": round(float(d.loc[m, "deficit"].mean()), 4),
            "n_significant": int(d.loc[m, "is_sig"].sum()),
            "pct_significant": round(100.0 * float(d.loc[m, "is_sig"].mean()), 3),
        })
    bins_df = pd.DataFrame(rows)
    bins_df.to_csv(snakemake.output.bins, sep="\t", index=False)

    # deficit is defined against a reference copy, so the rank tests use only
    # families that have one; zero-copy families are reported separately
    has_ref = d.ref_copies >= 1
    s, ns = d[has_ref & d.is_sig], d[has_ref & ~d.is_sig]
    _, p_copies = mannwhitneyu(s.ref_copies, ns.ref_copies, alternative="two-sided")
    _, p_deficit = mannwhitneyu(s.deficit, ns.deficit, alternative="two-sided")
    multi = d.ref_copies >= 2
    single = d.ref_copies == 1
    zero = d.ref_copies == 0
    table = [[int((multi & d.is_sig).sum()), int((multi & ~d.is_sig).sum())],
             [int((single & d.is_sig).sum()), int((single & ~d.is_sig).sum())]]
    odds, p_fisher = fisher_exact(table, alternative="two-sided")
    rate_multi = float(d.loc[multi, "is_sig"].mean())
    rate_single = float(d.loc[single, "is_sig"].mean())
    risk = rate_multi / rate_single if rate_single else float("inf")

    branch = pd.read_csv(snakemake.input.branch, sep="\t")
    branch = branch[branch["taxon"].isin(ingroup)]
    totals = {s_: int(d[s_].sum()) for s_ in ingroup}
    lin = []
    for _, r in branch.iterrows():
        t = r["taxon"]
        lin.append({
            "taxon": t,
            "total_copies": totals[t],
            "copy_deficit_vs_reference": totals[ref] - totals[t],
            "increase": int(r["increase"]), "decrease": int(r["decrease"]),
            "net_contraction": int(r["decrease"]) - int(r["increase"]),
        })
    lin_df = pd.DataFrame(lin).sort_values("copy_deficit_vs_reference")
    lin_df.to_csv(snakemake.output.lineage, sep="\t", index=False)

    if len(lin_df) >= 3 and lin_df["copy_deficit_vs_reference"].nunique() > 1:
        r_lin, p_lin = pearsonr(lin_df.copy_deficit_vs_reference,
                                lin_df.net_contraction)
    else:
        r_lin, p_lin = float("nan"), float("nan")

    summary = {
        "n_families": len(d),
        "n_significant": int(d.is_sig.sum()),
        "reference": ref,
        "n_families_single_copy": int(single.sum()),
        "n_families_multi_copy": int(multi.sum()),
        "n_families_zero_ref_copy": int(zero.sum()),
        "n_significant_single_copy": int((single & d.is_sig).sum()),
        "n_significant_multi_copy": int((multi & d.is_sig).sum()),
        "n_significant_zero_ref_copy": int((zero & d.is_sig).sum()),
        "retention_single_copy": round(float(d.loc[single, "retention"].mean()), 4),
        "retention_multi_copy": round(float(d.loc[multi, "retention"].mean()), 4),
        "pct_significant_single_copy": round(100.0 * rate_single, 3),
        "pct_significant_multi_copy": round(100.0 * rate_multi, 3),
        "pct_significant_zero_ref_copy": round(100.0 * float(d.loc[zero, "is_sig"].mean()), 3)
                                         if zero.any() else float("nan"),
        "odds_ratio_multi_vs_single_copy": round(float(odds), 3),
        "risk_ratio_multi_vs_single_copy": round(float(risk), 3),
        "fisher_p_two_sided": fmt_p(p_fisher),
        "mean_deficit_significant": round(float(s.deficit.mean()), 4),
        "mean_deficit_not_significant": round(float(ns.deficit.mean()), 4),
        "mannwhitney_p_ref_copies_two_sided": fmt_p(p_copies),
        "mannwhitney_p_deficit_two_sided": fmt_p(p_deficit),
        "per_lineage_pearson_r": round(float(r_lin), 4),
        "per_lineage_p": fmt_p(p_lin),
        "n_lineages": len(lin_df),
        "per_lineage_counts_scope": "all families (Gamma_clade_results)",
    }
    pd.DataFrame([summary]).T.reset_index().to_csv(
        snakemake.output.summary, sep="\t", index=False, header=["metric", "value"])

    print(f"families {len(d)}, significant {d.is_sig.sum()}, reference {ref}\n")
    print(bins_df.to_string(index=False))
    print(f"\nsignificance rate: single-copy {summary['pct_significant_single_copy']}%"
          f"  multi-copy {summary['pct_significant_multi_copy']}%"
          f"  (OR {summary['odds_ratio_multi_vs_single_copy']},"
          f" RR {summary['risk_ratio_multi_vs_single_copy']},"
          f" two sided p {summary['fisher_p_two_sided']})")
    print(f"zero reference copy families: {summary['n_families_zero_ref_copy']},"
          f" {summary['pct_significant_zero_ref_copy']}% significant")
    print(f"mean transfer deficit (families with a reference copy): significant "
          f"{summary['mean_deficit_significant']}  vs other "
          f"{summary['mean_deficit_not_significant']}"
          f"  (two sided p {summary['mannwhitney_p_deficit_two_sided']})")
    print(f"\nper-lineage net contraction vs copy deficit: "
          f"r = {summary['per_lineage_pearson_r']} over {len(lin_df)} lineages")
    print(lin_df.to_string(index=False))
    print("\nA declining retention curve with a rising significance rate means "
          "CAFE's\nper-lineage result is confounded with annotation transfer "
          "success and\ncannot be read as gene family evolution.")


if __name__ == "__main__":
    main()
