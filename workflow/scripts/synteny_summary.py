#!/usr/bin/env python3
"""
synteny_summary.py - pairwise gene-anchor synteny between ingroup genomes.

Produces the collinearity and inversion figures reported in the Results. The
previous version of this table was generated outside the workflow and could
not be regenerated from the repository; this rule replaces it.

Method
------
Anchors are genes sharing an identifier between two Liftoff annotations,
restricted to the three chromosome-scale sequences of each assembly (ranked by
gene count, which pairs homologous chromosomes unambiguously in this dataset).

Orientation is normalised before anything is measured. Assemblies differ in
the strand on which a chromosome was deposited, so a chromosome written
reverse-complemented in one genome makes every collinear block appear
inverted. For each chromosome pair the script compares the longest increasing
and longest decreasing anchor runs and flips the query coordinates when the
decreasing run is longer. The number of flipped pairs is reported, so the
correction is visible rather than silent.

Collinearity is the fraction of anchors lying on the longest strictly
increasing run of query positions when ordered by reference position, computed
per chromosome pair and pooled by anchor count.

An inversion is a maximal run of at least `min_genes` consecutive anchors whose
query positions decrease monotonically, spanning at least `min_span` bp in
reference coordinates, counted after orientation normalisation.

Alignment identity comes from minimap2 PAF when supplied: length-weighted
identity over all alignment blocks, not a partial-chromosome proxy.
"""

import bisect
import collections
import itertools
import os
import re
import sys


def read_genes(gff):
    """{gene_id: (seqid, midpoint)} for gene features."""
    out = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            m = re.search(r"ID=([^;]+)", f[8])
            if not m:
                continue
            out[m.group(1)] = (f[0], (int(f[3]) + int(f[4])) / 2.0)
    return out


def top_seqids(genes, n):
    """The n sequences carrying the most genes."""
    c = collections.Counter(s for s, _ in genes.values())
    return [s for s, _ in c.most_common(n)]


def lis_length(seq):
    """Longest strictly increasing subsequence length, O(n log n)."""
    tails = []
    for x in seq:
        i = bisect.bisect_left(tails, x)
        if i == len(tails):
            tails.append(x)
        else:
            tails[i] = x
    return len(tails)


def pair_chromosomes(ga, gb, tops_a, tops_b):
    """Map each of a's top seqids to b's by shared gene count (greedy, 1:1)."""
    by_a = {s: {g for g, (sq, _) in ga.items() if sq == s} for s in tops_a}
    by_b = {s: {g for g, (sq, _) in gb.items() if sq == s} for s in tops_b}
    scored = []
    for sa in tops_a:
        for sb in tops_b:
            scored.append((len(by_a[sa] & by_b[sb]), sa, sb))
    scored.sort(reverse=True)
    used_a, used_b, mapping = set(), set(), {}
    for n, sa, sb in scored:
        if sa in used_a or sb in used_b or n == 0:
            continue
        mapping[sa] = (sb, n)
        used_a.add(sa)
        used_b.add(sb)
    return mapping


def parse_paf(path):
    """(total aligned bp on the query, length-weighted percent identity)."""
    nmatch = blocklen = span = 0
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            span += int(f[3]) - int(f[2])
            nmatch += int(f[9])
            blocklen += int(f[10])
    ident = 100.0 * nmatch / blocklen if blocklen else float("nan")
    return span, ident


def read_ani(path, samples):
    """{(a, b): ani} from a skani matrix, if the file exists."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")[1:]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            for col, val in zip(cols, f[1:]):
                try:
                    out[(f[0], col)] = float(val)
                except (TypeError, ValueError):
                    pass
    return out


def analyse_pair(a, b, genes, n_chrom, min_span, min_genes):
    """Anchors, orientation, collinearity and inversions for one genome pair."""
    ga, gb = genes[a], genes[b]
    mapping = pair_chromosomes(ga, gb,
                               top_seqids(ga, n_chrom), top_seqids(gb, n_chrom))

    n_anchor_total = n_collinear_total = 0
    inversions, orientations = [], []

    for sa, (sb, _n) in mapping.items():
        anchors = []
        for g, (sq_a, pos_a) in ga.items():
            if sq_a != sa:
                continue
            loc_b = gb.get(g)
            if loc_b and loc_b[0] == sb:
                anchors.append((pos_a, loc_b[1]))
        if len(anchors) < 3:
            continue
        anchors.sort()
        q = [y for _, y in anchors]

        # Orientation: whichever direction supports the longer monotone run.
        inc = lis_length(q)
        dec = lis_length([-y for y in q])
        flipped = dec > inc
        if flipped:
            qmax = max(q)
            anchors = [(x, qmax - y) for x, y in anchors]
            q = [y for _, y in anchors]
        orientations.append({
            "sample1": a, "sample2": b, "seqid1": sa, "seqid2": sb,
            "n_anchors": len(anchors), "lis_forward": inc, "lis_reverse": dec,
            "orientation": "reverse" if flipped else "forward",
        })

        n_anchor_total += len(anchors)
        n_collinear_total += lis_length(q)

        # Inversions: maximal strictly decreasing runs, after normalisation.
        i = 0
        while i < len(anchors) - 1:
            j = i
            while j < len(anchors) - 1 and anchors[j + 1][1] < anchors[j][1]:
                j += 1
            if j - i + 1 >= min_genes:
                span = anchors[j][0] - anchors[i][0]
                if span >= min_span:
                    inversions.append({
                        "sample1": a, "sample2": b, "seqid": sa,
                        "start_bp": int(anchors[i][0]), "end_bp": int(anchors[j][0]),
                        "span_bp": int(span), "n_genes": j - i + 1,
                    })
            i = j + 1 if j > i else i + 1

    pct = 100.0 * n_collinear_total / n_anchor_total if n_anchor_total else float("nan")
    return n_anchor_total, pct, inversions, orientations


def main():
    if "snakemake" not in globals():
        sys.exit("Run via Snakemake.")

    samples = list(snakemake.params.samples)
    n_chrom = int(getattr(snakemake.params, "n_chrom", 3))
    min_span = int(getattr(snakemake.params, "min_span", 100_000))
    min_genes = int(getattr(snakemake.params, "min_genes", 3))

    genes = {}
    for gff in snakemake.input.gffs:
        s = os.path.basename(gff).replace("_liftoff.gff3", "")
        genes[s] = read_genes(gff)
        print(f"{s}: {len(genes[s])} genes")

    paf_by_pair = {os.path.basename(p).replace(".paf", ""): p
                   for p in getattr(snakemake.input, "pafs", [])}
    ani = read_ani(getattr(snakemake.input, "ani", None), samples)

    rows, all_inv, all_or = [], [], []
    for a, b in itertools.combinations(samples, 2):
        n_anchors, pct, inv, orients = analyse_pair(
            a, b, genes, n_chrom, min_span, min_genes)
        all_inv.extend(inv)
        all_or.extend(orients)

        span = ident = float("nan")
        src = "not computed"
        for key in (f"{a}_vs_{b}", f"{b}_vs_{a}"):
            if key in paf_by_pair:
                span, ident = parse_paf(paf_by_pair[key])
                src = "minimap2 asm20, chromosome-scale sequences"
                break

        top = "; ".join(
            f"{d['seqid']}:{d['start_bp']/1e6:.2f}-{d['end_bp']/1e6:.2f}Mb({d['n_genes']}genes)"
            for d in sorted(inv, key=lambda d: -d["span_bp"])[:3]) or "none"
        n_rev = sum(1 for o in orients if o["orientation"] == "reverse")

        rows.append({
            "pair": f"{a}_vs_{b}", "sample1": a, "sample2": b,
            "n_shared_genes_chr1_3": n_anchors,
            "pct_collinear_anchors": round(pct, 2),
            "n_chrom_pairs_reverse_oriented": n_rev,
            "n_inversions_intra_chr": len(inv),
            "max_inversion_size_bp": max((d["span_bp"] for d in inv), default=0),
            "top_inversions": top,
            "total_aligned_bp": span,
            "mean_alignment_identity_pct": round(ident, 4) if ident == ident else "",
            "mean_skani_ani_pct": ani.get((a, b), ani.get((b, a), "")),
            "aln_identity_source": src,
        })
        print(f"{a} vs {b}: {n_anchors} anchors, {pct:.2f}% collinear, "
              f"{len(inv)} inversions, {n_rev}/{n_chrom} chromosome pairs reversed")

    def write_tsv(path, dicts, cols):
        with open(path, "w") as out:
            out.write("\t".join(cols) + "\n")
            for d in dicts:
                out.write("\t".join(str(d.get(c, "")) for c in cols) + "\n")

    write_tsv(snakemake.output.summary, rows, list(rows[0].keys()))
    write_tsv(snakemake.output.orientation, all_or,
              ["sample1", "sample2", "seqid1", "seqid2", "n_anchors",
               "lis_forward", "lis_reverse", "orientation"])
    write_tsv(snakemake.output.inversions, all_inv,
              ["sample1", "sample2", "seqid", "start_bp", "end_bp",
               "span_bp", "n_genes"])

    n_rev_total = sum(1 for o in all_or if o["orientation"] == "reverse")
    print(f"\n{n_rev_total}/{len(all_or)} chromosome pairs required orientation "
          f"normalisation before counting.")


if __name__ == "__main__":
    main()
