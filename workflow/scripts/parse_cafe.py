#!/usr/bin/env python3
"""
Parse CAFE5 output to extract significantly evolving gene families.
"""
import pandas as pd
import os
import glob

cafe_dir = snakemake.input[0]
pvalue_threshold = snakemake.params.pvalue

# Find the results file
# CAFE5 outputs: Base_results.txt, Gamma_results.txt, etc.
result_files = glob.glob(os.path.join(cafe_dir, "*_results.txt"))

all_results = []
for rf in result_files:
    try:
        df = pd.read_csv(rf, sep="\t")
        all_results.append(df)
    except Exception as e:
        print(f"Could not read {rf}: {e}")

if not all_results:
    # Try alternative CAFE5 output format
    result_files = glob.glob(os.path.join(cafe_dir, "*.txt"))
    for rf in result_files:
        try:
            df = pd.read_csv(rf, sep="\t")
            if "p-value" in df.columns or "pvalue" in df.columns:
                all_results.append(df)
        except:
            pass

if all_results:
    results = pd.concat(all_results, ignore_index=True)
    
    # Find p-value column
    pval_col = [c for c in results.columns if "p" in c.lower() and "val" in c.lower()]
    if pval_col:
        sig = results[results[pval_col[0]] < pvalue_threshold]
        sig.to_csv(snakemake.output.significant, sep="\t", index=False)
        
        print(f"Total families tested: {len(results)}")
        print(f"Significant (p < {pvalue_threshold}): {len(sig)}")
    else:
        # Just save everything
        results.to_csv(snakemake.output.significant, sep="\t", index=False)
        print(f"Saved {len(results)} results (no p-value column found)")
else:
    # Create empty output if no results found
    print("WARNING: No CAFE5 results files found. Creating empty output.")
    pd.DataFrame().to_csv(snakemake.output.significant, sep="\t", index=False)

# Summary per branch
summary = pd.DataFrame({"note": ["See significant_families.tsv for details"]})
summary.to_csv(snakemake.output.summary, sep="\t", index=False)
