#!/usr/bin/env python3
import pandas as pd
import os, glob

cafe_dir = snakemake.input[0]
pvalue_threshold = snakemake.params.pvalue
result_files = glob.glob(os.path.join(cafe_dir, "*_results.txt"))
all_results = []
for rf in result_files:
    try:
        df = pd.read_csv(rf, sep="\t")
        all_results.append(df)
    except:
        pass

if all_results:
    results = pd.concat(all_results, ignore_index=True)
    pval_col = [c for c in results.columns if "p" in c.lower() and "val" in c.lower()]
    if pval_col:
        sig = results[results[pval_col[0]] < pvalue_threshold]
        sig.to_csv(snakemake.output.significant, sep="\t", index=False)
        print(f"Significant (p < {pvalue_threshold}): {len(sig)}")
    else:
        results.to_csv(snakemake.output.significant, sep="\t", index=False)
else:
    pd.DataFrame().to_csv(snakemake.output.significant, sep="\t", index=False)

pd.DataFrame({"note": ["See significant_families.tsv"]}).to_csv(snakemake.output.summary, sep="\t", index=False)
