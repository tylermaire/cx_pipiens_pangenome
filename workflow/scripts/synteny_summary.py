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

Every inversion row also carries the homologous chromosome (named from the
reference assembly's chromosome sequences in natural order, chr1 to chr3),
its position as a fraction of the anchored span of that chromosome, and its
first and last anchor gene. Because Liftoff gives every genome the reference
gene IDs, inversions from different genome pairs can be compared by gene
content: two inversions that share at least half their anchor genes (Jaccard
>= 0.5) are joined into one recurrence cluster, which is how the claim that
large inversions recur across pairs is tested.
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


def count_other_chromosome(ga, gb, n_chrom):
    """(shared genes on the top sequences of both genomes, how many of them
    sit on a non-homologous chromosome in b). These anchors are left out of
    the collinearity and inversion counts, so the share is reported apart."""
    tops_a, tops_b = top_seqids(ga, n_chrom), top_seqids(gb, n_chrom)
    mapping = pair_chromosomes(ga, gb, tops_a, tops_b)
    shared = other = 0
    for g, (sa, _pos) in ga.items():
        if sa not in mapping:
            continue
        loc = gb.get(g)
        if not loc or loc[0] not in tops_b:
            continue
        shared += 1
        if loc[0] != mapping[sa][0]:
            other += 1
    return shared, other


def parse_paf(path):
    """(total aligned bp on the query, length-weighted percent identity).

    Identity comes from minimap2's divergence tag, weighted by alignment
    length: de:f: (gap-compressed) when present, else dv:f: (approximate),
    which is what minimap2 writes without base-level alignment. The obvious alternative,
    nmatch/blocklen from columns 10 and 11, is not identity for long-range
    asm20 alignments: block length counts gap positions, so a 94%-identical
    pair reads as roughly 70%.

    Falls back to nmatch/blocklen only when no de: tag is present, in which
    case the value is flagged by the caller.
    """
    nmatch = blocklen = span = 0
    div_sum = div_w = 0.0
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            alen = int(f[3]) - int(f[2])
            span += alen
            nmatch += int(f[9])
            blocklen += int(f[10])
            for tag in f[12:]:
                if tag.startswith("de:f:") or tag.startswith("dv:f:"):
                    try:
                        div_sum += float(tag.split(":")[2]) * alen
                        div_w += alen
                    except ValueError:
                        pass
                    break
    if div_w:
        return span, 100.0 * (1.0 - div_sum / div_w)
    return span, 100.0 * nmatch / blocklen if blocklen else float("nan")


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


def analyse_pair(a, b, genes, n_chrom, min_span, min_genes, homolog=None):
    """Anchors, orientation, collinearity and inversions for one genome pair.

    homolog maps (sample, seqid) to a chromosome name shared by all genomes.
    """
    ga, gb = genes[a], genes[b]
    homolog = homolog or {}
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
                anchors.append((pos_a, loc_b[1], g))
        if len(anchors) < 3:
            continue
        anchors.sort()
        q = [y for _, y, _ in anchors]
        lo_a, hi_a = anchors[0][0], anchors[-1][0]
        extent = max(hi_a - lo_a, 1.0)

        # Orientation: whichever direction supports the longer monotone run.
        inc = lis_length(q)
        dec = lis_length([-y for y in q])
        flipped = dec > inc
        if flipped:
            qmax = max(q)
            anchors = [(x, qmax - y, g) for x, y, g in anchors]
            q = [y for _, y, _ in anchors]
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
                        "chromosome": homolog.get((a, sa), ""),
                        "start_bp": int(anchors[i][0]), "end_bp": int(anchors[j][0]),
                        "span_bp": int(span), "n_genes": j - i + 1,
                        "rel_start": round((anchors[i][0] - lo_a) / extent, 3),
                        "rel_end": round((anchors[j][0] - lo_a) / extent, 3),
                        "first_gene": anchors[i][2], "last_gene": anchors[j][2],
                        "genes": [g for _, _, g in anchors[i:j + 1]],
                    })
            i = j + 1 if j > i else i + 1

    pct = 100.0 * n_collinear_total / n_anchor_total if n_anchor_total else float("nan")
    return n_anchor_total, pct, inversions, orientations


def homolog_names(genes, samples, reference, n_chrom):
    """{(sample, seqid): 'chrN'}: reference chromosomes named chr1..chrN in
    natural order of their sequence IDs, other genomes by shared gene content."""
    def natural(x):
        return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", x)]
    ref_tops = sorted(top_seqids(genes[reference], n_chrom), key=natural)
    ref_name = {sq: f"chr{i + 1}" for i, sq in enumerate(ref_tops)}
    out = {(reference, sq): name for sq, name in ref_name.items()}
    for s in samples:
        if s == reference:
            continue
        mapping = pair_chromosomes(genes[s], genes[reference],
                                   top_seqids(genes[s], n_chrom), ref_tops)
        for sq, (ref_sq, _n) in mapping.items():
            out[(s, sq)] = ref_name[ref_sq]
    return out


def recurrence_clusters(inversions, min_jaccard=0.5):
    """Join inversions from different genome pairs that share at least
    min_jaccard of their anchor genes. Returns a cluster id per inversion."""
    parent = list(range(len(inversions)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    sets = [set(d["genes"]) for d in inversions]
    for i, j in itertools.combinations(range(len(inversions)), 2):
        a, b = inversions[i], inversions[j]
        if (a["sample1"], a["sample2"]) == (b["sample1"], b["sample2"]):
            continue
        inter = len(sets[i] & sets[j])
        if inter and inter / len(sets[i] | sets[j]) >= min_jaccard:
            parent[find(i)] = find(j)
    roots = {}
    return [roots.setdefault(find(i), len(roots) + 1) for i in range(len(inversions))]


def recurrence_table(inversions, clusters, large_span):
    by = collections.defaultdict(list)
    for d, c in zip(inversions, clusters):
        by[c].append(d)
    rows = []
    for c, members in sorted(by.items()):
        pairs = sorted({f"{d['sample1']}_vs_{d['sample2']}" for d in members})
        shared = set.intersection(*({d["sample1"], d["sample2"]} for d in members))
        rows.append({
            "cluster": c, "n_inversions": len(members), "n_pairs": len(pairs),
            # with two or more pairs, the genome in every one of them; a cluster
            # in every pair that includes one genome marks a rearrangement in,
            # or an assembly error of, that genome
            "genome_in_every_pair": (shared.pop() if len(shared) == 1 else "none")
                                    if len(pairs) > 1 else "",
            "chromosomes": ",".join(sorted({d["chromosome"] for d in members})),
            "max_span_bp": max(d["span_bp"] for d in members),
            "any_at_least_large_span": any(d["span_bp"] >= large_span for d in members),
            "pairs": ";".join(pairs),
            "members": ";".join(f"{d['sample1']}_vs_{d['sample2']}:{d['seqid']}:"
                                f"{d['start_bp']}-{d['end_bp']}" for d in members),
        })
    return rows


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
    reference = getattr(snakemake.params, "reference", samples[0])
    large_span = int(getattr(snakemake.params, "large_span", 1_000_000))
    homolog = homolog_names(genes, samples, reference, n_chrom)

    rows, all_inv, all_or = [], [], []
    for a, b in itertools.combinations(samples, 2):
        n_anchors, pct, inv, orients = analyse_pair(
            a, b, genes, n_chrom, min_span, min_genes, homolog)
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
        n_shared, n_other = count_other_chromosome(genes[a], genes[b], n_chrom)

        rows.append({
            "pair": f"{a}_vs_{b}", "sample1": a, "sample2": b,
            "n_shared_genes_chr1_3": n_anchors,
            "n_shared_genes_top_sequences": n_shared,
            "n_shared_genes_other_chromosome": n_other,
            "pct_shared_genes_other_chromosome":
                round(100.0 * n_other / n_shared, 2) if n_shared else "",
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
    clusters = recurrence_clusters(all_inv)
    cluster_size = collections.Counter(clusters)
    for d, c in zip(all_inv, clusters):
        d["recurrence_cluster"] = c
        d["n_pairs_in_cluster"] = len({(x["sample1"], x["sample2"])
                                       for x, k in zip(all_inv, clusters) if k == c})
    write_tsv(snakemake.output.inversions, all_inv,
              ["sample1", "sample2", "seqid", "chromosome", "start_bp", "end_bp",
               "span_bp", "n_genes", "rel_start", "rel_end", "first_gene",
               "last_gene", "recurrence_cluster", "n_pairs_in_cluster"])
    rec = recurrence_table(all_inv, clusters, large_span)
    write_tsv(snakemake.output.recurrence, rec,
              ["cluster", "n_inversions", "n_pairs", "chromosomes", "max_span_bp",
               "genome_in_every_pair", "any_at_least_large_span", "pairs", "members"])
    large = [d for d in all_inv if d["span_bp"] >= large_span]
    large_recurrent = sum(1 for d in large if d["n_pairs_in_cluster"] > 1)
    print(f"\n{len(all_inv)} inversions, {len(large)} of at least {large_span:,} bp; "
          f"{large_recurrent} of those recur in another pair "
          f"({sum(1 for c, n in cluster_size.items() if n > 1)} multi member clusters)")

    n_rev_total = sum(1 for o in all_or if o["orientation"] == "reverse")
    print(f"\n{n_rev_total}/{len(all_or)} chromosome pairs required orientation "
          f"normalisation before counting.")


if __name__ == "__main__":
    main()
