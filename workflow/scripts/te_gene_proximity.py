#!/usr/bin/env python3
"""
For each gene in each ingroup genome, compute the distance to the nearest
RepeatMasker transposable-element (TE) insertion, then compare the
distance distributions between pangenome compartments (core, shell, cloud).

This addresses the question: "Are accessory (shell/cloud) genes more often
proximal to TEs than core genes?", a standard test for TE-mediated gene
gain or loss in comparative genomics.

Inputs (via snakemake.input)
    gff_files       : list of per-genome GFF3 files (gene coordinates)
    repeatmasker    : list of per-genome RepeatMasker .out files (TE coordinates)
    partitions      : results/pangenome/partitioned_orthogroups.tsv
    orthogroups_tsv : OrthoFinder Orthogroups.tsv (gene -> orthogroup mapping)

Params (via snakemake.params)
    samples : list of sample names matched to gff_files / repeatmasker in order
    feature : GFF3 feature type to include (default 'gene')

Outputs (via snakemake.output)
    per_gene  : long TSV: sample, gene, chrom, start, end, compartment,
                nearest_te_bp, nearest_te_family
    summary   : per-compartment summary statistics (median, IQR, % within 1 kb)
"""
import sys
from pathlib import Path
from collections import defaultdict
import bisect

import pandas as pd


def parse_gff(path: Path, feature: str = "gene"):
    """Yield (chrom, start, end, gene_id) for features of the requested type."""
    with open(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != feature:
                continue
            attrs = dict(
                kv.split("=", 1) for kv in f[8].split(";") if "=" in kv
            )
            gid = attrs.get("ID") or attrs.get("Name") or attrs.get("gene_id")
            if gid:
                yield f[0], int(f[3]), int(f[4]), gid


def parse_repeatmasker_out(path: Path):
    """Yield (chrom, start, end, family) from a RepeatMasker .out file.
    The file has three header lines; data lines have whitespace-separated cols
    where columns 5/6/7 are query_name, query_start, query_end and column 11
    is the repeat class/family.
    """
    with open(path) as fh:
        for _ in range(3):
            try:
                next(fh)
            except StopIteration:
                return
        for line in fh:
            cols = line.split()
            if len(cols) < 11:
                continue
            try:
                chrom = cols[4]
                start = int(cols[5])
                end = int(cols[6])
                family = cols[10]
            except (ValueError, IndexError):
                continue
            yield chrom, start, end, family


def build_te_index(rm_path: Path):
    """Return {chrom: (sorted_starts, sorted_ends, families)} for fast lookup."""
    by_chrom = defaultdict(list)
    for chrom, s, e, fam in parse_repeatmasker_out(rm_path):
        by_chrom[chrom].append((s, e, fam))
    index = {}
    for chrom, rows in by_chrom.items():
        rows.sort()
        starts = [r[0] for r in rows]
        ends = [r[1] for r in rows]
        fams = [r[2] for r in rows]
        index[chrom] = (starts, ends, fams)
    return index


def nearest_te(chrom: str, g_start: int, g_end: int, index) -> tuple[int, str]:
    """Return (distance_bp, nearest_te_family) for the nearest TE to the gene."""
    if chrom not in index:
        return (-1, "")
    starts, ends, fams = index[chrom]
    n = len(starts)
    if n == 0:
        return (-1, "")
    # Candidate TEs: the one whose start <= g_end (predecessor) and the one
    # whose start > g_end (successor). Compare distances.
    i = bisect.bisect_left(starts, g_start)
    candidates = []
    for j in (i - 1, i, i + 1):
        if 0 <= j < n:
            ts, te, tf = starts[j], ends[j], fams[j]
            if ts <= g_end and te >= g_start:
                d = 0
            elif te < g_start:
                d = g_start - te
            else:
                d = ts - g_end
            candidates.append((d, tf))
    if not candidates:
        return (-1, "")
    return min(candidates, key=lambda x: x[0])


def main() -> None:
    if "snakemake" not in globals():
        sys.exit("Run via Snakemake.")
    gff_files = [Path(p) for p in snakemake.input.gff_files]        # type: ignore
    rm_files  = [Path(p) for p in snakemake.input.repeatmasker]      # type: ignore
    partitions = pd.read_csv(snakemake.input.partitions, sep="\t")    # type: ignore
    # OrthoFinder writes to Results_<date>/Orthogroups/Orthogroups.tsv;
    # the rule passes the parent dir so we glob for the actual file.
    import glob
    og_dir = snakemake.input.orthogroups_tsv  # type: ignore
    matches = glob.glob(f"{og_dir}/Results_*/Orthogroups/Orthogroups.tsv")
    if not matches:
        sys.exit(f"Could not find Orthogroups.tsv under {og_dir}")
    orthogroups = pd.read_csv(matches[0], sep="\t")
    samples = list(snakemake.params.samples)                           # type: ignore
    feature = getattr(snakemake.params, "feature", "gene")             # type: ignore
    out_per = Path(snakemake.output.per_gene)                          # type: ignore
    out_sum = Path(snakemake.output.summary)                           # type: ignore

    # Build transcript_id -> gene_id map from all GFFs (NCBI-style GFFs put
    # gene IDs on `gene` features and OrthoFinder uses mRNA IDs from the
    # protein FASTAs, so we need to translate before joining).
    transcript_to_gene = {}
    for gff in gff_files:
        with open(gff) as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] not in ("mRNA", "transcript"):
                    continue
                attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                tid = attrs.get("ID")
                parent = attrs.get("Parent")
                if tid and parent:
                    transcript_to_gene[tid] = parent
                    transcript_to_gene[tid.replace(".", "")] = parent

    # gene_id -> (orthogroup_id, compartment, sample)
    gene_meta = {}
    og_to_compartment = dict(zip(partitions["Orthogroup"], partitions["compartment"]))
    for _, row in orthogroups.iterrows():
        og = row["Orthogroup"]
        comp = og_to_compartment.get(og)
        if comp is None:
            continue
        for sample in samples:
            cell = row.get(sample)
            if pd.isna(cell) or not str(cell).strip():
                continue
            for gid in str(cell).split(","):
                gid = gid.strip()
                real = transcript_to_gene.get(gid, gid)
                gene_meta[(sample, real)] = (og, comp, sample)

    long_rows = []
    for sample, gff, rmf in zip(samples, gff_files, rm_files):
        index = build_te_index(rmf)
        for chrom, gs, ge, gid in parse_gff(gff, feature=feature):
            meta = gene_meta.get((sample, gid))
            if meta is None:
                continue
            og, comp, _ = meta
            d, fam = nearest_te(chrom, gs, ge, index)
            long_rows.append({
                "sample": sample, "gene": gid, "chrom": chrom,
                "start": gs, "end": ge, "compartment": comp,
                "nearest_te_bp": d, "nearest_te_family": fam,
            })

    long_df = pd.DataFrame(long_rows)
    out_per.parent.mkdir(parents=True, exist_ok=True)
    long_df.to_csv(out_per, sep="\t", index=False)

    valid = long_df[long_df["nearest_te_bp"] >= 0]
    summary = (
        valid.groupby(["sample", "compartment"])["nearest_te_bp"]
        .agg(median="median",
             q25=lambda s: s.quantile(0.25),
             q75=lambda s: s.quantile(0.75),
             within_1kb=lambda s: float((s <= 1000).mean()),
             n="size")
        .reset_index()
    )
    summary.to_csv(out_sum, sep="\t", index=False)


if __name__ == "__main__":
    main()
