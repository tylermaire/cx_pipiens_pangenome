#!/usr/bin/env python3
"""
transfer_quality.py - how many transferred gene models are broken, and where
they end up in the pangenome.

Liftoff checks every transcript it places and writes the result into the GFF:
valid_ORF=True/False, and when something is wrong missing_start_codon=True,
missing_stop_codon=True or inframe_stop_codon=True, plus matches_ref_protein.
The workflow runs Liftoff without -polish, so none of these are repaired, and
extract_proteins strips stop symbols before the proteins reach OrthoFinder,
so a broken model looks like any other protein downstream. This script reads
Liftoff's own flags for the models that were kept (one per gene) and reports:

  summary         per genome: kept models, and how many carry each flag
  by_compartment  per transferred genome: share of kept models without a valid
                  ORF in core, shell, cloud, outgroup only and unassigned
                  orthogroups

The reference annotation carries no Liftoff flags. Its row counts the kept
RefSeq models that are not a clean ORF on the reference genome themselves
(partial=true, or an exception such as "unclassified transcription
discrepancy", where RefSeq corrected an indel in the assembly). A transferred
copy of such a model usually fails Liftoff's check (87 to 89% of them in the
V4 ingroup genomes; the rest pass where the target sequence supplies a clean
ORF), so each transferred row also gives the invalid ORF share among models
whose reference model is clean.

Snakemake provides:
    input.gffs, input.proteins, input.table, input.of
    params.samples, params.reference
    output.summary, output.by_compartment
"""

import collections
import csv
import glob
import os
import re
import sys

ATTR_ID = re.compile(r"(?:^|;)ID=([^;]+)")
FLAGS = ["valid_ORF", "missing_start_codon", "missing_stop_codon",
         "inframe_stop_codon", "matches_ref_protein"]
REF_KEYS = ["partial", "exception"]      # RefSeq attributes of the reference GFF
TRANSCRIPT_TYPES = {"mRNA", "transcript"}
COMPARTMENTS = ["core", "shell", "cloud", "outgroup_only", "unassigned"]


def fasta_ids(path):
    with open(path) as fh:
        return [line[1:].split()[0] for line in fh if line.startswith(">")]


def liftoff_flags(gff):
    """{transcript (periods removed): {attribute: value}} for transcripts,
    keeping Liftoff's ORF flags and the RefSeq partial and exception keys."""
    out = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in TRANSCRIPT_TYPES:
                continue
            m = ATTR_ID.search(f[8])
            if not m:
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            out[m.group(1).replace(".", "")] = {k: attrs[k] for k in FLAGS + REF_KEYS
                                                if k in attrs}
    return out


def not_clean(attrs):
    """True when a RefSeq model is itself not a clean ORF on its genome."""
    return attrs.get("partial") == "true" or "exception" in attrs


def summarise_sample(sample, kept, flags, is_reference, ref=None):
    """One summary row. ref is the reference genome's liftoff_flags() output,
    used to set aside models whose reference model is not clean."""
    found = [flags[t] for t in kept if t in flags]
    n = len(found)
    has_flags = any("valid_ORF" in f for f in found)
    row = {"sample": sample, "n_kept_models": len(kept), "n_in_gff": n,
           "liftoff_flags": "yes" if has_flags else
           ("no (reference annotation)" if is_reference else "no")}

    def count(pred):
        return sum(1 for f in found if pred(f))

    def invalid(f):
        return f.get("valid_ORF") == "False"

    for label, pred in (
            ("invalid_orf", invalid),
            ("inframe_stop", lambda f: f.get("inframe_stop_codon") == "True"),
            ("missing_start", lambda f: f.get("missing_start_codon") == "True"),
            ("missing_stop", lambda f: f.get("missing_stop_codon") == "True"),
            ("mismatch_ref_protein", lambda f: f.get("matches_ref_protein") == "False")):
        k = count(pred) if has_flags else ""
        row[f"n_{label}"] = k
        row[f"pct_{label}"] = round(100.0 * k / n, 2) if has_flags and n else ""

    # reference baseline
    if is_reference:
        row["n_ref_model_not_clean"] = count(not_clean)
        row["n_clean_ref_models"] = row["n_invalid_orf_clean_ref"] = ""
        row["pct_invalid_orf_clean_ref"] = ""
    elif has_flags and ref:
        clean = [flags[t] for t in kept if t in flags and t in ref and not not_clean(ref[t])]
        bad = sum(1 for f in clean if invalid(f))
        row["n_ref_model_not_clean"] = sum(1 for t in kept if t in flags and t in ref
                                           and not_clean(ref[t]))
        row["n_clean_ref_models"] = len(clean)
        row["n_invalid_orf_clean_ref"] = bad
        row["pct_invalid_orf_clean_ref"] = round(100.0 * bad / len(clean), 2) if clean else ""
    else:
        for k in ("n_ref_model_not_clean", "n_clean_ref_models",
                  "n_invalid_orf_clean_ref", "pct_invalid_orf_clean_ref"):
            row[k] = ""
    return row


def compartments(table_path, of_dir, samples):
    """{(sample, transcript): compartment}; genes in no orthogroup are 'unassigned'."""
    comp = {}
    with open(table_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            comp[r["Orthogroup"]] = r["compartment"]
    files = glob.glob(os.path.join(of_dir, "**", "Orthogroups.tsv"), recursive=True)
    if not files:
        raise SystemExit(f"No Orthogroups.tsv found under {of_dir}")
    out = {}
    with open(files[0]) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            og = r["Orthogroup"]
            for s in samples:
                for g in (x.strip() for x in (r.get(s) or "").split(",")):
                    if g:
                        out[(s, g)] = comp.get(og, "unassigned")
    return out


def by_compartment(sample, kept, flags, where):
    tally = collections.defaultdict(lambda: [0, 0])
    for t in kept:
        f = flags.get(t)
        if not f or "valid_ORF" not in f:
            continue
        c = where.get((sample, t), "unassigned")
        tally[c][0] += 1
        tally[c][1] += f.get("valid_ORF") == "False"
    rows = []
    order = COMPARTMENTS + sorted(set(tally) - set(COMPARTMENTS))
    for c in order:
        n, bad = tally.get(c, (0, 0))
        if not n:
            continue
        rows.append({"sample": sample, "compartment": c, "n_models": n,
                     "n_invalid_orf": bad,
                     "pct_invalid_orf": round(100.0 * bad / n, 2)})
    return rows


def write(rows, path):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    if "snakemake" not in globals():
        sys.exit("Run via Snakemake or workflow/scripts/patch_results.py.")
    sm = globals()["snakemake"]
    samples = list(sm.params.samples)
    reference = sm.params.reference
    gff = {os.path.basename(p).replace("_liftoff.gff3", ""): p for p in sm.input.gffs}
    prot = {os.path.basename(p)[:-3]: p for p in sm.input.proteins}
    where = compartments(sm.input.table, sm.input.of, samples)

    ref = liftoff_flags(gff[reference]) if reference in gff else None
    summary, comp_rows = [], []
    for s in samples:
        kept = fasta_ids(prot[s])
        flags = ref if s == reference and ref is not None else liftoff_flags(gff[s])
        row = summarise_sample(s, kept, flags, s == reference, ref)
        summary.append(row)
        if row["liftoff_flags"] == "yes":
            comp_rows.extend(by_compartment(s, kept, flags, where))
    write(summary, sm.output.summary)
    write(comp_rows, sm.output.by_compartment)
    for r in summary:
        print(f"  {r['sample']:22s} kept {r['n_kept_models']:6d}  invalid ORF "
              f"{r['pct_invalid_orf']!s:>6}%  (clean reference model "
              f"{r['pct_invalid_orf_clean_ref']!s:>6}%)  in-frame stop "
              f"{r['pct_inframe_stop']!s:>6}%  reference not clean "
              f"{r['n_ref_model_not_clean']}")
    for r in comp_rows:
        print(f"  {r['sample']:22s} {r['compartment']:10s} {r['n_models']:6d} models, "
              f"{r['pct_invalid_orf']!s:>6}% without a valid ORF")


if __name__ == "__main__":
    main()
