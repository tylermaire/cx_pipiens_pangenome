#!/usr/bin/env python3
"""
Check whether the chromosome 3 (3Rb) paracentric inversion identified in
Cx. quinquefasciatus (Ryazansky et al., 2024) is shared by other forms of
the Cx. pipiens complex.

Strategy
--------
For each ingroup pairwise comparison against Cx. quinquefasciatus, parse the
SyRI output (`syri.out`) — a TSV with columns:

    RefChrom  RefStart  RefEnd  RefSeq  QrySeq  QryChrom  QryStart  QryEnd
    UniqueID  ParentID  Annotation  CopyStatus

We look for entries on the chromosome that hosts the 3Rb breakpoints in the
Cx. quinquefasciatus reference (GCF_015732765.1 chromosome "NC_051845.1" or
the equivalent name from the assembly's report), filtered to annotation type
"INV" (inversion) or "INVAL" (inverted alignment), with length >= the minimum
inversion length expected for 3Rb (~5 Mb in Ryazansky et al., 2024).

For each non-reference form we report:
    - the number of inversion blocks on chr3 that overlap the 3Rb interval
    - the cumulative inverted length
    - whether the inversion appears shared (overlap >= 50% of the reference
      3Rb span) or unique to Cx. quinquefasciatus

Inputs (via snakemake.input)
    syri_outs : list[Path]   - one syri.out per pairwise comparison
    chr_map   : Path          - 2-col TSV mapping reference contig accession ->
                                chromosome label (e.g. NC_051845.1 -> chr3)

Params (via snakemake.params)
    chr3_label    : str       - the chromosome label that hosts 3Rb
    rb3_start     : int       - reference start coordinate of 3Rb (bp)
    rb3_end       : int       - reference end coordinate of 3Rb (bp)
    min_inv_len   : int       - minimum inversion block length to include (bp)
    shared_thresh : float     - fraction overlap to call "shared" (default 0.5)

Outputs (via snakemake.output)
    table : Path  - per-comparison summary TSV
    flag  : Path  - one-line interpretation:
                    "3Rb shared with: <list>; private to Cx. quinquefasciatus: <list>"
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd


SYRI_COLUMNS = [
    "RefChrom", "RefStart", "RefEnd", "RefSeq", "QrySeq",
    "QryChrom", "QryStart", "QryEnd", "UniqueID", "ParentID",
    "Annotation", "CopyStatus",
]


def parse_syri(path: Path) -> pd.DataFrame:
    """Read a SyRI output file into a DataFrame. SyRI uses '-' as missing."""
    df = pd.read_csv(
        path, sep="\t", header=None, names=SYRI_COLUMNS,
        dtype=str, na_values=["-", ""],
    )
    for col in ("RefStart", "RefEnd", "QryStart", "QryEnd"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def overlap_bp(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def summarise_one(syri_df, chr_label, rb_start, rb_end, min_inv_len):
    """Return (n_blocks, total_inv_bp, overlap_with_rb_bp)."""
    inv_types = {"INV", "INVAL", "INVDP", "INVDPAL", "INVTR", "INVTRAL"}
    mask = (
        (syri_df["RefChrom"] == chr_label)
        & (syri_df["Annotation"].isin(inv_types))
        & ((syri_df["RefEnd"] - syri_df["RefStart"]) >= min_inv_len)
    )
    sub = syri_df[mask].copy()
    n_blocks = len(sub)
    total_inv = int((sub["RefEnd"] - sub["RefStart"]).sum()) if n_blocks else 0
    if n_blocks:
        sub["overlap"] = sub.apply(
            lambda r: overlap_bp(int(r.RefStart), int(r.RefEnd), rb_start, rb_end),
            axis=1,
        )
        rb_overlap = int(sub["overlap"].sum())
    else:
        rb_overlap = 0
    return n_blocks, total_inv, rb_overlap


def main() -> None:
    if "snakemake" not in globals():
        sys.exit("This script is intended to be run via Snakemake.")

    syri_outs   = [Path(p) for p in snakemake.input.syri_outs]   # type: ignore
    chr3_label  = snakemake.params.chr3_label                      # type: ignore
    rb3_start   = int(snakemake.params.rb3_start)                  # type: ignore
    rb3_end     = int(snakemake.params.rb3_end)                    # type: ignore
    min_inv_len = int(getattr(snakemake.params, "min_inv_len", 100_000))  # type: ignore
    shared_th   = float(getattr(snakemake.params, "shared_thresh", 0.5)) # type: ignore
    out_table   = Path(snakemake.output.table)                     # type: ignore
    out_flag    = Path(snakemake.output.flag)                      # type: ignore

    rb_span = rb3_end - rb3_start
    rows, shared, private = [], [], []
    for p in syri_outs:
        df = parse_syri(p)
        # SyRI file names follow the convention ref_vs_query.syri.out
        name = p.stem.replace(".syri", "")
        if "_vs_" in name:
            ref, qry = name.split("_vs_", 1)
        else:
            ref, qry = "ref", name
        nb, total_inv, rb_overlap = summarise_one(
            df, chr3_label, rb3_start, rb3_end, min_inv_len
        )
        frac = rb_overlap / rb_span if rb_span else 0.0
        rows.append({
            "reference": ref,
            "query": qry,
            "chr3_inv_blocks": nb,
            "chr3_inv_total_bp": total_inv,
            "rb3_overlap_bp": rb_overlap,
            "rb3_overlap_fraction": round(frac, 4),
            "shared_with_reference": bool(frac >= shared_th),
        })
        (shared if frac >= shared_th else private).append(qry)

    out_table.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_table, sep="\t", index=False)

    out_flag.write_text(
        "3Rb inversion analysis (Cx. quinquefasciatus reference)\n"
        f"  shared with: {', '.join(shared) if shared else '(none)'}\n"
        f"  private to Cx. quinquefasciatus: "
        f"{', '.join(private) if private else '(none)'}\n"
    )


if __name__ == "__main__":
    main()
