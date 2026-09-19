#!/usr/bin/env python3
"""
Tabulate copy numbers of key gene families (P450, GST, CCE, OR, GR, IR,
OBP, CSP, immune) across the four ingroup forms using eggNOG-mapper
annotations layered on top of the OrthoFinder orthogroups.

Two changes from the original, both of which altered the published counts:

1. Family assignment scores every family and takes the highest, rather than
   returning the first family whose pattern happens to match. Under the old
   first-match rule the dict order decided the outcome, so OR (defined before
   OBP, with the loose pattern \\bor[0-9]+\\b) claimed genes annotated as
   "odorant binding protein 12 | Obp12" before OBP was ever tried. Patterns
   now carry a weight so a specific phrase outranks an incidental token.

2. An orthogroup is assigned by plurality with a floor, not by requiring 50%
   of its genes to classify. The old rule discarded any orthogroup where
   eggNOG annotated fewer than half the members, which penalised sparsely
   annotated families hardest; OBP fell from 2 to 0 when the longest-isoform
   filter shrank the denominators.
"""
import re
import sys
from pathlib import Path
from collections import defaultdict
import pandas as pd


INGROUP = ("Cx_molestus", "Cx_pallens", "Cx_pipiens", "Cx_quinquefasciatus")


def read_eggnog(path: Path) -> pd.DataFrame:
    """Parse an eggNOG-mapper annotations TSV (skips '##' comment lines)."""
    with open(path) as fh:
        lines = [ln for ln in fh if not ln.startswith("##")]
    from io import StringIO
    df = pd.read_csv(StringIO("".join(lines)), sep="\t", dtype=str,
                     na_values=["-"], keep_default_na=False)
    if "#query" in df.columns:
        df = df.rename(columns={"#query": "query"})
    return df


def build_patterns(cfg):
    """Accept either [weight, regex] pairs or bare regex strings (weight 1)."""
    out = {}
    for fam, pats in cfg.items():
        compiled = []
        for entry in pats:
            if isinstance(entry, (list, tuple)) and len(entry) == 2:
                w, p = entry
            else:
                w, p = 1, entry
            compiled.append((int(w), re.compile(p, re.IGNORECASE)))
        out[fam] = compiled
    return out


def classify_family(text, family_patterns):
    """Highest-scoring family for this text, or None. Ties leave it unassigned."""
    if not text:
        return None
    scores = {}
    for fam, pats in family_patterns.items():
        s = sum(w for w, p in pats if p.search(text))
        if s:
            scores[fam] = s
    if not scores:
        return None
    best = max(scores.values())
    winners = [f for f, s in scores.items() if s == best]
    return winners[0] if len(winners) == 1 else None


def main() -> None:
    if "snakemake" not in globals():
        sys.exit("Run via Snakemake.")

    part = pd.read_csv(snakemake.input.partitions, sep="\t")      # type: ignore
    annot_paths = [Path(p) for p in snakemake.input.annotations]  # type: ignore
    families_cfg: dict = snakemake.params.families                # type: ignore
    min_genes = int(getattr(snakemake.params, "min_genes", 2))    # type: ignore
    min_frac = float(getattr(snakemake.params, "min_frac", 0.25)) # type: ignore
    out_long = Path(snakemake.output.long)                        # type: ignore
    out_wide = Path(snakemake.output.wide)                        # type: ignore

    family_patterns = build_patterns(families_cfg)

    gene_to_family = {}
    n_rows = n_classified = 0
    for p in annot_paths:
        df = read_eggnog(p)
        text_cols = [c for c in ("Description", "PFAMs", "Preferred_name",
                                 "best_OG_desc") if c in df.columns]
        if not text_cols:
            sys.exit(f"No usable description columns in {p}: {list(df.columns)}")
        for _, row in df.iterrows():
            n_rows += 1
            text = " | ".join(str(row.get(c, "") or "") for c in text_cols)
            fam = classify_family(text, family_patterns)
            if fam:
                gene_to_family[str(row["query"])] = fam
                n_classified += 1

    print(f"eggNOG rows: {n_rows}; classified to a key family: {n_classified}")
    per_fam = defaultdict(int)
    for f in gene_to_family.values():
        per_fam[f] += 1
    for f in families_cfg:
        print(f"  {f:8s} {per_fam.get(f, 0):6d} genes")
        if per_fam.get(f, 0) == 0:
            print(f"    WARNING: no genes matched {f}; check its patterns "
                  f"against the eggNOG Description/PFAMs text.")

    of_dir = Path(snakemake.input.partitions).parent.parent / "orthofinder" / "output"
    of_tables = list(of_dir.rglob("Orthogroups.tsv"))
    if not of_tables:
        sys.exit(f"Could not find Orthogroups.tsv under {of_dir}")
    of = pd.read_csv(of_tables[0], sep="\t").set_index("Orthogroup")

    long_rows = []
    for og, row in of.iterrows():
        fam_counts = defaultdict(int)
        total_genes = 0
        for form in INGROUP:
            cell = row.get(form)
            if pd.isna(cell) or not str(cell).strip():
                continue
            for g in [x.strip() for x in str(cell).split(",") if x.strip()]:
                total_genes += 1
                fam = gene_to_family.get(g)
                if fam:
                    fam_counts[fam] += 1
        if not fam_counts or total_genes == 0:
            continue
        best_fam, best_n = max(fam_counts.items(), key=lambda kv: kv[1])
        if best_n < min_genes and best_n / total_genes < min_frac:
            continue
        for form in INGROUP:
            try:
                n = int(part.loc[part["Orthogroup"] == og, form].iloc[0])
            except IndexError:
                n = 0
            long_rows.append({"orthogroup": og, "form": form,
                              "copies": n, "family": best_fam})

    long_df = pd.DataFrame(long_rows)
    out_long.parent.mkdir(parents=True, exist_ok=True)
    long_df.to_csv(out_long, sep="\t", index=False)

    wide = (
        long_df.groupby(["family", "form"])["copies"]
        .sum().unstack(fill_value=0)
        .reindex(columns=INGROUP, fill_value=0)
        .reset_index()
    )
    wide.to_csv(out_wide, sep="\t", index=False)

    n_og = long_df["orthogroup"].nunique() if len(long_df) else 0
    print(f"\nassigned {n_og} orthogroups to key families")
    if len(long_df):
        print(wide.to_string(index=False))


if __name__ == "__main__":
    main()
