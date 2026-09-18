#!/usr/bin/env python3
"""
Tabulate copy numbers of key gene families (P450, GST, CCE, OR, GR, IR,
OBP, CSP, immune) across the four ingroup forms using eggNOG-mapper
annotations layered on top of the OrthoFinder orthogroups.

Inputs
------
snakemake.input.partitions    : results/pangenome/partitioned_orthogroups.tsv
snakemake.input.annotations   : list of eggNOG-mapper annotation TSVs
                                (one per ingroup proteome, comment lines start
                                with '#'; column 'Preferred_name' or '#query'
                                holds the gene identifier; 'Description' and
                                'PFAMs' hold functional descriptors)

Params
------
snakemake.params.families     : dict mapping family-label -> list of regex
                                patterns matched (case-insensitive) against
                                the eggNOG Description / PFAMs / Preferred_name
                                fields. Example:
                                  {
                                    "P450":    [r"cytochrome\\s*p[-_]?450",
                                                r"\\bcyp[0-9]"],
                                    "GST":     [r"glutathione\\s*s[-_]?transferase",
                                                r"\\bgst[a-z]?[0-9]?"],
                                    "CCE":     [r"carboxylesterase",
                                                r"\\bcce[0-9]?"],
                                    "OR":      [r"\\bodorant\\s*receptor",
                                                r"\\bor[0-9]+\\b"],
                                    "GR":      [r"\\bgustatory\\s*receptor",
                                                r"\\bgr[0-9]+\\b"],
                                    "IR":      [r"\\bionotropic\\s*receptor",
                                                r"\\bir[0-9]+\\b"],
                                    "OBP":     [r"\\bodorant[-_]?binding\\s*protein"],
                                    "CSP":     [r"\\bchemosensory\\s*protein"],
                                    "immune":  [r"toll[-_]?like", r"\\bimd\\b",
                                                r"defensin", r"cecropin",
                                                r"attacin", r"thioester",
                                                r"\\bpgrp\\b"],
                                  }

Outputs
-------
snakemake.output.long   : tidy long-format TSV
                          orthogroup  form  copies  family
snakemake.output.wide   : per-family per-form copy-number matrix:
                          family   Cx_molestus  Cx_pallens  Cx_pipiens  Cx_quinquefasciatus
"""
import re
import sys
from pathlib import Path
from collections import defaultdict
import pandas as pd


INGROUP = ("Cx_molestus", "Cx_pallens", "Cx_pipiens", "Cx_quinquefasciatus")


def read_eggnog(path: Path) -> pd.DataFrame:
    """Parse an eggNOG-mapper annotations TSV (skips '#' comment lines)."""
    with open(path) as fh:
        # eggNOG header line is preceded by '##' but the column header line
        # itself starts with '#query' - find it.
        lines = []
        for line in fh:
            if line.startswith("##"):
                continue
            lines.append(line)
    from io import StringIO
    df = pd.read_csv(StringIO("".join(lines)), sep="\t", dtype=str,
                     na_values=["-"], keep_default_na=False)
    # Rename '#query' -> 'query' if present
    if "#query" in df.columns:
        df = df.rename(columns={"#query": "query"})
    return df


def classify_family(text: str, family_patterns: dict[str, list[re.Pattern]]):
    """Return the first family whose any-pattern matches `text`, else None."""
    if not text:
        return None
    for fam, pats in family_patterns.items():
        for p in pats:
            if p.search(text):
                return fam
    return None


def main() -> None:
    if "snakemake" not in globals():
        sys.exit("Run via Snakemake.")

    part = pd.read_csv(snakemake.input.partitions, sep="\t")     # type: ignore
    annot_paths = [Path(p) for p in snakemake.input.annotations]  # type: ignore
    families_cfg: dict = snakemake.params.families                 # type: ignore
    out_long = Path(snakemake.output.long)                         # type: ignore
    out_wide = Path(snakemake.output.wide)                         # type: ignore

    # Compile regexes once
    family_patterns = {
        fam: [re.compile(p, re.IGNORECASE) for p in pats]
        for fam, pats in families_cfg.items()
    }

    # Build a gene -> family map by concatenating Description + PFAMs + Preferred_name
    # for each eggNOG row, then classifying.
    gene_to_family: dict[str, str] = {}
    for p in annot_paths:
        df = read_eggnog(p)
        text_cols = [c for c in ("Description", "PFAMs", "Preferred_name", "best_OG_desc") if c in df.columns]
        for _, row in df.iterrows():
            text = " | ".join(str(row.get(c, "") or "") for c in text_cols)
            fam = classify_family(text, family_patterns)
            if fam:
                gene_to_family[str(row["query"])] = fam

    # For each orthogroup, count copies per form, then count how many of those
    # gene IDs map to a key family. Voting strategy: an orthogroup is assigned
    # to a family if >= 50% of its gene members classify to that family.
    long_rows = []
    of_dir = Path(snakemake.input.partitions).parent.parent / "orthofinder" / "output"  # type: ignore
    # We need per-orthogroup gene IDs, which live in Orthogroups.tsv inside the
    # OrthoFinder output. Try to locate it.
    of_tables = list(of_dir.rglob("Orthogroups.tsv"))
    if not of_tables:
        sys.exit(f"Could not find Orthogroups.tsv under {of_dir}")
    of = pd.read_csv(of_tables[0], sep="\t").set_index("Orthogroup")

    for og, row in of.iterrows():
        fam_counts = defaultdict(int)
        total_genes = 0
        for form in INGROUP:
            cell = row.get(form)
            if pd.isna(cell) or not str(cell).strip():
                continue
            genes = [g.strip() for g in str(cell).split(",") if g.strip()]
            for g in genes:
                total_genes += 1
                fam = gene_to_family.get(g)
                if fam:
                    fam_counts[fam] += 1
        # Assign the orthogroup to the family with >=50% of classified genes
        if not fam_counts or total_genes == 0:
            continue
        best_fam, best_n = max(fam_counts.items(), key=lambda kv: kv[1])
        if best_n / total_genes < 0.5:
            continue
        # Emit per-form rows
        for form in INGROUP:
            try:
                n = int(part.loc[part["Orthogroup"] == og, form].iloc[0])
            except IndexError:
                n = 0
            long_rows.append({"orthogroup": og, "form": form, "copies": n, "family": best_fam})

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


if __name__ == "__main__":
    main()
