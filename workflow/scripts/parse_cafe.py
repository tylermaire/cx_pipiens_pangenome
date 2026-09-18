#!/usr/bin/env python3
"""
parse_cafe.py - extract significantly evolving gene families and per-branch
expansion/contraction counts from CAFE5 gamma-model output.

Reads two files by name rather than globbing. The previous version globbed
"*_results.txt", which matches both Gamma_results.txt (a four-line model
summary) and Gamma_family_results.txt (the per-family table), concatenated
them, and produced a file whose header was three headers joined together and
whose row count was therefore wrong by one. Both downstream consumers compute
the significant-family count as line count minus one, so that malformed file
is the sole cause of the off-by-one.

Outputs:
  significant_families.tsv  one row per family below the p threshold, with the
                            orthogroup id carried through from the Desc field
  branch_summary.tsv        per-lineage increases, decreases, ratio and an
                            exact binomial test against equal rates

The binomial test replaces a hardcoded "p<2e-16 each branch" caption. That
claim is not true of every dataset and must be read from the data.
"""

import os
import re
import sys

import pandas as pd

try:
    from scipy.stats import binomtest
except ImportError:
    binomtest = None

cafe_dir = snakemake.input[0]
pvalue_threshold = snakemake.params.pvalue

family_file = os.path.join(cafe_dir, "Gamma_family_results.txt")
clade_file = os.path.join(cafe_dir, "Gamma_clade_results.txt")

for path in (family_file, clade_file):
    if not os.path.exists(path):
        sys.exit(f"ERROR: {path} not found. CAFE5 must be run with the gamma "
                 f"model (-p -k N); Base-model output is not accepted here.")

# ---- significant families -------------------------------------------------
fam = pd.read_csv(family_file, sep="\t")
fam.columns = [c.lstrip("#").strip() for c in fam.columns]

pcol = next((c for c in fam.columns if c.lower() == "pvalue"), None)
if pcol is None:
    pcol = next((c for c in fam.columns if "p" in c.lower() and "val" in c.lower()), None)
if pcol is None:
    sys.exit(f"ERROR: no p-value column in {family_file}: {list(fam.columns)}")

idcol = fam.columns[0]
n_anon = (fam[idcol].astype(str).str.strip().isin(["", "n/a", "nan"])).sum()
if n_anon:
    print(f"WARNING: {n_anon}/{len(fam)} families have no identifier. Check that "
          f"format_cafe_input.py wrote the orthogroup id into the Desc column.")

sig = fam[fam[pcol] < pvalue_threshold].copy()
sig = sig.rename(columns={idcol: "orthogroup"})
sig.to_csv(snakemake.output.significant, sep="\t", index=False)
print(f"Significant (p < {pvalue_threshold}): {len(sig)} of {len(fam)} families")

# ---- per-branch expansion and contraction ---------------------------------
clade = pd.read_csv(clade_file, sep="\t")
clade.columns = [c.lstrip("#").strip() for c in clade.columns]

rows = []
for r in clade.itertuples(index=False):
    raw = str(getattr(r, "Taxon_ID", getattr(r, clade.columns[0], "")))
    name = re.sub(r"<\d+>$", "", raw)
    inc, dec = int(r.Increase), int(r.Decrease)
    total = inc + dec
    ratio = inc / dec if dec else float("nan")
    if binomtest is not None and total:
        p = binomtest(inc, total, 0.5).pvalue
    else:
        p = float("nan")
    rows.append({"taxon": name, "node_label": raw, "increase": inc,
                 "decrease": dec, "total": total,
                 "ratio_inc_dec": round(ratio, 4), "binomial_p": p})

summary = pd.DataFrame(rows)
summary.to_csv(snakemake.output.summary, sep="\t", index=False)

print("\nper-branch gene family change")
for r in summary.itertuples(index=False):
    print(f"  {r.taxon:24s} {r.increase:6d} inc  {r.decrease:6d} dec  "
          f"ratio {r.ratio_inc_dec:6.2f}  p={r.binomial_p:.3g}")
if binomtest is None:
    print("  (scipy unavailable; binomial p not computed)")
