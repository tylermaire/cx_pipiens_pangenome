#!/usr/bin/env python3
"""Filter a FASTA to its N longest contigs and rename them chr1/chr2/.../chrN
by length rank (descending). Gives every species the same chromosome names so
SyRI can match them across pairwise comparisons."""
import sys
from pathlib import Path

def parse_fasta(p):
    seqs, name, chunks = {}, None, []
    for line in open(p):
        line = line.rstrip()
        if line.startswith('>'):
            if name: seqs[name] = ''.join(chunks)
            name = line[1:].split()[0]; chunks = []
        else:
            chunks.append(line)
    if name: seqs[name] = ''.join(chunks)
    return seqs

infile  = sys.argv[1]
outfile = sys.argv[2]
n_keep  = int(sys.argv[3]) if len(sys.argv) > 3 else 3

seqs   = parse_fasta(infile)
ranked = sorted(seqs.items(), key=lambda kv: -len(kv[1]))[:n_keep]

Path(outfile).parent.mkdir(parents=True, exist_ok=True)
with open(outfile, 'w') as f:
    for i, (orig_name, seq) in enumerate(ranked, 1):
        f.write(f">chr{i} original={orig_name} length={len(seq)}\n")
        for k in range(0, len(seq), 80):
            f.write(seq[k:k+80] + '\n')

for i, (orig, seq) in enumerate(ranked, 1):
    print(f"chr{i}: {orig} ({len(seq):,} bp)")
