#!/usr/bin/env python3
"""
Resolve the chromosome-3 RefSeq contig name and the chromosome-3 coordinate
range of the 3Rb (3qB) paracentric inversion for the Cx. quinquefasciatus
reference assembly (VPISU_Cqui_1.0_pri_paternal = GCF_015732765.1).

Ryazansky et al. (2024, BMC Biology) place the 3Rb breakpoints at
108,955,000–117,355,000 bp on chromosomal arm 3q of this assembly. The
inversion spans 8.4 Mb; resolution is ~5 kb.

Run this once on the EC2 box after the genome is downloaded:

    python workflow/scripts/resolve_3Rb_coords.py \
        --fasta   resources/genomes/Cx_quinquefasciatus.fasta \
        --report  resources/genomes/Cx_quinquefasciatus.assembly_report.txt

The assembly_report.txt is auto-downloaded next to the FASTA by
download_genome.smk (it's part of the NCBI Datasets bundle); if you don't
have it, just point --report at /dev/null and we'll skip the chr-label
disambiguation and fall back to "longest contig with chromosome-3-scale
length".

Prints three lines suitable for pasting into rule check_3Rb_inversion's
`params:` block.
"""
from __future__ import annotations
import argparse
import re
import sys
from pathlib import Path


def parse_fasta_lengths(path: Path) -> dict[str, int]:
    """Return {seqid: length} from a FASTA file without loading sequences."""
    lengths = {}
    seqid = None
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if seqid is not None:
                    lengths[seqid] = n
                seqid = line[1:].split()[0]
                n = 0
            else:
                n += len(line.rstrip())
        if seqid is not None:
            lengths[seqid] = n
    return lengths


def parse_assembly_report(path: Path) -> list[dict]:
    """Parse an NCBI assembly_report.txt; return list of rows for
    Assembled-Molecule entries with chromosome names."""
    if not path or not path.exists():
        return []
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            # Columns: Sequence-Name, Sequence-Role, Assigned-Molecule,
            # Assigned-Molecule-Location/Type, GenBank-Accn, Relationship,
            # RefSeq-Accn, Assembly-Unit, Sequence-Length, UCSC-style-name
            if len(f) < 10:
                continue
            rows.append({
                "name": f[0],
                "role": f[1],
                "assigned_mol": f[2],
                "mol_type": f[3],
                "genbank": f[4],
                "refseq": f[6],
                "length": int(f[8]) if f[8].isdigit() else 0,
            })
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta",  required=True, type=Path)
    ap.add_argument("--report", required=False, type=Path,
                    help="NCBI assembly_report.txt (optional but recommended)")
    # 3Rb coordinates on chromosomal arm 3q (Ryazansky et al. 2024)
    ap.add_argument("--rb3q-start", type=int, default=108_955_000)
    ap.add_argument("--rb3q-end",   type=int, default=117_355_000)
    args = ap.parse_args()

    lengths = parse_fasta_lengths(args.fasta)
    report_rows = parse_assembly_report(args.report) if args.report else []

    # ---- Identify chromosome 3 ---------------------------------------------
    chr3 = None
    if report_rows:
        # Find the entry whose assigned_mol is "3" (or "Chr3", "chr3", etc.)
        for r in report_rows:
            am = r["assigned_mol"].lstrip("Cc").lstrip("hr")
            if am == "3" and r["role"] == "assembled-molecule":
                chr3 = r
                break
    if chr3 is None:
        # Heuristic fallback: pick the contig whose length is consistent with
        # chr3 in this assembly. The Ryazansky paper places 3q breakpoints at
        # ~117 Mb, so chr3 must be at least ~120 Mb. Among the four assembled
        # chromosomes the longest is usually chr3 in culicines.
        candidates = sorted(lengths.items(), key=lambda kv: -kv[1])[:5]
        sys.stderr.write(
            "WARNING: no assembly_report supplied; falling back to "
            "longest-contig heuristic. Verify the chosen contig by hand.\n"
            f"Top-5 longest contigs:\n"
        )
        for sid, n in candidates:
            sys.stderr.write(f"  {sid}\t{n:,} bp\n")
        chr3 = {"name": candidates[0][0], "refseq": candidates[0][0],
                "length": candidates[0][1]}

    contig_name = chr3.get("refseq") or chr3["name"]
    contig_len = lengths.get(contig_name) or chr3["length"]

    # ---- Determine 3p / 3q orientation --------------------------------------
    # The Ryazansky paper reports breakpoints on the 3q arm at 108.955 Mb and
    # 117.355 Mb. The centromere sits between 3p and 3q. We don't have the
    # centromere coord directly, but we can bound it by the chromosome length:
    # if 3q ends at 117 Mb and chr3 is ~187 Mb, then either
    #   - 3p first (length L_p), 3q starts at L_p,  centromere ~ L_p
    #     -> chr3 coord of breakpoint = L_p + 108_955_000
    #   - 3q first, centromere ~ 117 Mb (rb3q-end), 3p ends at chr3 length
    #     -> chr3 coord of breakpoint = 108_955_000 (no offset)
    #
    # Without the centromere annotation we can't pick automatically. Report
    # BOTH candidates and let the user pick by checking which one falls
    # inside the chromosome.
    cand_a_start = args.rb3q_start
    cand_a_end   = args.rb3q_end
    cand_b_start = contig_len - args.rb3q_end
    cand_b_end   = contig_len - args.rb3q_start

    print(f"# chromosome-3 contig (RefSeq): {contig_name}")
    print(f"# chromosome-3 length          : {contig_len:,} bp")
    print()
    print(f"# 3Rb breakpoints in 3q-arm space (Ryazansky et al. 2024):")
    print(f"#   start = {args.rb3q_start:,}   end = {args.rb3q_end:,}")
    print()
    print("# Two candidate translations to chromosome-3 coordinates:")
    print(f"#   A) 3q starts at position 0 of chr3 (3q-first orientation)")
    print(f"#      chr3_label    = \"{contig_name}\"")
    print(f"#      rb3_start     = {cand_a_start:_}")
    print(f"#      rb3_end       = {cand_a_end:_}")
    print()
    print(f"#   B) 3p starts at position 0 of chr3 (3p-first orientation,")
    print(f"#      breakpoints mirrored against the chromosome end)")
    print(f"#      chr3_label    = \"{contig_name}\"")
    print(f"#      rb3_start     = {cand_b_start:_}")
    print(f"#      rb3_end       = {cand_b_end:_}")
    print()
    print("# Quick disambiguation: open the reference FASTA in IGV and look")
    print("# at gene density at each end of chr3 — the centromeric end has")
    print("# very low gene density and the breakpoints should bracket it.")
    print("# Equivalently, run `samtools faidx` and check N-content per Mb.")


if __name__ == "__main__":
    main()
