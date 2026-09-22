#!/usr/bin/env python3
"""
classify_absence.py - decide, for every absence call in the pangenome, whether
the gene is really absent, merely unannotated, or present under another
orthogroup.

Each absence event is an (orthogroup, genome) pair the orthogroup table says is
empty. Evidence is combined from five sources:

  gene ID   whether the target genome carries a protein coding model of the
            same gene. Liftoff keeps reference IDs, so every transferred model
            of one reference gene has that gene's ID in every genome
  protein   miniprot alignment of the orthogroup's representative protein
  dna       minimap2 alignment of that gene's genomic locus
  context   every protein coding model the target annotates under the hit
  identity  pairwise protein identity between the query and EACH of those models

Calls:

  clustered_elsewhere  the target carries a model of the same gene, assigned to
                       a different orthogroup (or left unassigned). Two routes:
                         same_gene_id      the target proteome holds a model
                                           with the query's gene ID. The gene
                                           is present whatever the probes say
                         overlap_identity  a probe hit overlaps a model that is
                                           near identical to the query
  unannotated_locus    a probe finds the sequence, but no protein coding model
                       overlaps the hit: Liftoff placed nothing usable there
  paralog_only         a hit overlaps protein coding models, none of them the
                       query's gene and none near identical to it
  absent               no probe finds the locus
  weak                 a hit below threshold; reported separately, never binned

Only 'absent' and 'paralog_only' support a real presence/absence difference.

Why the gene ID test runs first: a transferred model that has diverged from its
source (a frameshifted or truncated transfer) is still the same gene. Without
the ID test such a query lands on its own source gene, falls below the
identity threshold and is counted as a paralog, which inflates the supported
fraction. The first version of this script had exactly that failure; on the
2026-09-17 run 902 of 1,159 cloud 'paralog_only' calls hit a model with the
query's own transcript ID.

Overlapping transcripts without a protein (non coding, or no valid CDS) are not
counted as annotation: a hit that lands only on those is 'unannotated_locus'.

Snakemake provides:
    input.manifest, input.prot_pafs, input.dna_pafs, input.of, input.gffs,
    input.proteins
    params.min_id, params.min_cov, params.dna_min_cov, params.same_gene_id,
    params.ingroup, params.max_overlap
    output.calls, output.summary
"""

import collections
import glob
import os
import re

import pandas as pd

try:
    from Bio import Align
    from Bio.Align import substitution_matrices
except ImportError:
    raise SystemExit("Biopython required in the phylo environment")

ATTR_ID = re.compile(r"(?:^|;)ID=([^;]+)")
ATTR_PARENT = re.compile(r"(?:^|;)Parent=([^;,]+)")
TRANSCRIPT_TYPES = {"mRNA", "transcript"}

SUPPORTED = ("absent", "paralog_only")
ARTIFACT = ("clustered_elsewhere", "unannotated_locus")
CALL_ORDER = ["absent", "paralog_only", "clustered_elsewhere",
              "unannotated_locus", "weak"]


def read_fasta(path):
    seqs, rid, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if rid is not None:
                    seqs[rid] = "".join(chunks)
                rid = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if rid is not None:
        seqs[rid] = "".join(chunks)
    return seqs


def parse_paf(path):
    """{(query, target_genome): (identity, coverage, seqid, start, end)}."""
    base = os.path.basename(path).replace(".paf", "")
    target = base.split("_vs_")[-1]
    best = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            q, qlen, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
            tname, ts, te = f[5], int(f[7]), int(f[8])
            nmatch, blocklen = int(f[9]), int(f[10])
            ident = nmatch / max(blocklen, 1)
            cov = (qe - qs) / max(qlen, 1)
            key = (q, target)
            if key not in best or ident * cov > best[key][0] * best[key][1]:
                best[key] = (ident, cov, tname, ts, te)
    return best


def read_gff_transcripts(path):
    """Transcript intervals and transcript to gene map for one GFF3.

    IDs have periods removed so they match the protein FASTA headers, which
    longest_isoform.py writes with --strip-periods.
    """
    by_seq = collections.defaultdict(list)
    tx_to_gene = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in TRANSCRIPT_TYPES:
                continue
            m = ATTR_ID.search(f[8])
            if not m:
                continue
            tid = m.group(1).replace(".", "")
            p = ATTR_PARENT.search(f[8])
            tx_to_gene[tid] = p.group(1).replace(".", "") if p else tid
            by_seq[f[0]].append((int(f[3]), int(f[4]), tid))
    for k in by_seq:
        by_seq[k].sort()
    return dict(by_seq), tx_to_gene


def load_membership(of_dir, forms):
    """{(form, transcript): orthogroup} from OrthoFinder's Orthogroups.tsv."""
    og_files = glob.glob(os.path.join(of_dir, "**", "Orthogroups.tsv"),
                         recursive=True)
    if not og_files:
        raise SystemExit(f"No Orthogroups.tsv found under {of_dir}")
    members = pd.read_csv(og_files[0], sep="\t", index_col=0).fillna("")
    gene_to_og = {}
    for og, row in members.iterrows():
        for form in forms:
            if form in members.columns:
                for g in [x.strip() for x in str(row[form]).split(",") if x.strip()]:
                    gene_to_og[(form, g)] = og
    return gene_to_og


class Classifier:
    def __init__(self, ingroup, proteomes, gffs, gene_to_og, prot, dna,
                 min_id=0.80, min_cov=0.70, dna_min_cov=0.50,
                 same_gene_id=0.90, max_overlap=10):
        self.ingroup = ingroup
        self.proteomes = proteomes
        self.gene_to_og = gene_to_og
        self.prot, self.dna = prot, dna
        self.min_id, self.min_cov = min_id, min_cov
        self.dna_min_cov = dna_min_cov
        self.same_gene_id = same_gene_id
        self.max_overlap = max_overlap
        self.intervals, self.tx_to_gene, self.coding_model = {}, {}, {}
        for form in ingroup:
            iv, t2g = gffs[form]
            self.intervals[form] = iv
            self.tx_to_gene[form] = t2g
            # gene -> the transcript that carries its protein in this genome
            models = {}
            for tid in proteomes[form]:
                models.setdefault(t2g.get(tid, tid), tid)
            self.coding_model[form] = models
        aligner = Align.PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -11
        aligner.extend_gap_score = -1
        aligner.mode = "global"
        self.aligner = aligner
        self._cache = {}

    def protein_identity(self, a, b, key):
        """Identical residues over the shorter sequence, cached by gene pair."""
        if key in self._cache:
            return self._cache[key]
        if not a or not b or min(len(a), len(b)) / max(len(a), len(b)) < 0.30:
            self._cache[key] = 0.0
            return 0.0
        aln = self.aligner.align(a, b)[0]
        same = 0
        for (s1, e1), (s2, e2) in zip(*aln.aligned):
            same += sum(1 for x, y in zip(a[s1:e1], b[s2:e2]) if x == y)
        val = same / min(len(a), len(b))
        self._cache[key] = val
        return val

    def overlapping(self, form, seqid, start, end):
        """Protein coding transcripts overlapping the interval, longest overlap first."""
        hits = []
        for s, e, tid in self.intervals[form].get(seqid, []):
            if s > end:
                break
            if e >= start and tid in self.proteomes[form]:
                hits.append((min(e, end) - max(s, start), tid))
        hits.sort(reverse=True)
        return [tid for _, tid in hits[:self.max_overlap]]

    def classify(self, og, comp, rep_gene, rep_form, target):
        qseq = self.proteomes.get(rep_form, {}).get(rep_gene)
        pi, pc, pseq, ps, pe = self.prot.get((og, target), (0.0, 0.0, "", 0, 0))
        di, dc, dseq, ds, de = self.dna.get((og, target), (0.0, 0.0, "", 0, 0))
        row = {"orthogroup": og, "compartment": comp, "absent_from": target,
               "prot_identity": round(pi, 3), "prot_coverage": round(pc, 3),
               "dna_identity": round(di, 3), "dna_coverage": round(dc, 3),
               "call": "", "best_overlapping_gene": "",
               "best_protein_identity": 0.0, "match_basis": "",
               "rep_gene": rep_gene, "rep_form": rep_form,
               "other_orthogroup": ""}

        # 1. same gene present in the target under another orthogroup
        q_gene = self.tx_to_gene.get(rep_form, {}).get(rep_gene)
        t_tx = self.coding_model[target].get(q_gene) if q_gene else None
        if t_tx:
            ident = self.protein_identity(qseq, self.proteomes[target][t_tx],
                                          (rep_gene, target, t_tx))
            row.update(call="clustered_elsewhere", match_basis="same_gene_id",
                       best_overlapping_gene=t_tx,
                       best_protein_identity=round(ident, 4),
                       other_orthogroup=self.gene_to_og.get((target, t_tx),
                                                            "unassigned"))
            return row

        # 2. probes
        prot_ok = pi >= self.min_id and pc >= self.min_cov
        dna_ok = di >= self.min_id and dc >= self.dna_min_cov
        if not (prot_ok or dna_ok):
            row["call"] = "weak" if (pi > 0 or di > 0) else "absent"
            return row

        seqid, s, e = (pseq, ps, pe) if prot_ok else (dseq, ds, de)
        genes = self.overlapping(target, seqid, s, e)
        if not genes:
            row["call"] = "unannotated_locus"
            return row

        # 3. compare against EVERY overlapping model, keep the best match
        best_id, best_gene = 0.0, ""
        for g in genes:
            val = self.protein_identity(qseq, self.proteomes[target].get(g),
                                        (rep_gene, target, g))
            if val > best_id:
                best_id, best_gene = val, g
        row.update(best_overlapping_gene=best_gene,
                   best_protein_identity=round(best_id, 4),
                   other_orthogroup=(self.gene_to_og.get((target, best_gene),
                                                         "unassigned")
                                     if best_gene else ""))
        if best_id >= self.same_gene_id:
            row.update(call="clustered_elsewhere", match_basis="overlap_identity")
        else:
            row["call"] = "paralog_only"
        return row


def summarise(calls):
    """One row per compartment: counts per call plus supported and artifact shares."""
    rows = []
    for comp in ["cloud", "shell"]:
        sub = calls[calls["compartment"] == comp]
        if sub.empty:
            continue
        n = len(sub)
        vc = sub["call"].value_counts()
        row = {"compartment": comp, "n_events": n}
        for c in CALL_ORDER:
            row[c] = int(vc.get(c, 0))
        row["clustered_elsewhere_same_gene_id"] = int(
            (sub["match_basis"] == "same_gene_id").sum())
        row["clustered_elsewhere_overlap_identity"] = int(
            (sub["match_basis"] == "overlap_identity").sum())
        row["pct_supported"] = round(100.0 * sub["call"].isin(SUPPORTED).mean(), 1)
        row["pct_artifact"] = round(100.0 * sub["call"].isin(ARTIFACT).mean(), 1)
        row["pct_weak"] = round(100.0 * (sub["call"] == "weak").mean(), 1)
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    p = snakemake.params
    ingroup = list(p.ingroup)
    proteomes = {f: read_fasta(f"results/proteins/{f}.fa") for f in ingroup}
    gffs = {f: read_gff_transcripts(f"results/annotation/{f}_liftoff.gff3")
            for f in ingroup}
    gene_to_og = load_membership(snakemake.input.of, ingroup)

    prot, dna = {}, {}
    for path in snakemake.input.prot_pafs:
        prot.update(parse_paf(path))
    for path in snakemake.input.dna_pafs:
        dna.update(parse_paf(path))

    clf = Classifier(ingroup, proteomes, gffs, gene_to_og, prot, dna,
                     min_id=float(getattr(p, "min_id", 0.80)),
                     min_cov=float(getattr(p, "min_cov", 0.70)),
                     dna_min_cov=float(getattr(p, "dna_min_cov", 0.50)),
                     same_gene_id=float(getattr(p, "same_gene_id", 0.90)),
                     max_overlap=int(getattr(p, "max_overlap", 10)))

    man = pd.read_csv(snakemake.input.manifest, sep="\t")
    rows = []
    for n, r in enumerate(man.itertuples(), 1):
        for target in str(r.absent_forms).split(","):
            target = target.strip()
            if target:
                rows.append(clf.classify(r.orthogroup, r.compartment,
                                         r.rep_gene, r.rep_form, target))
        if n % 500 == 0:
            print(f"  ...{n}/{len(man)} orthogroups", flush=True)

    calls = pd.DataFrame(rows)
    calls.to_csv(snakemake.output.calls, sep="\t", index=False)
    summary = summarise(calls)
    summary.to_csv(snakemake.output.summary, sep="\t", index=False)

    print(f"\n=== absence validation (hit id>={clf.min_id} cov>={clf.min_cov}; "
          f"same gene protein identity >={clf.same_gene_id}) ===")
    print(summary.to_string(index=False))
    print("\n'absent' and 'paralog_only' support a real presence/absence difference;")
    print("'clustered_elsewhere' and 'unannotated_locus' are artifacts.")


if __name__ == "__main__":
    main()
