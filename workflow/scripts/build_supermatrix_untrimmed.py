#!/usr/bin/env python3
"""Build supermatrix from untrimmed MAFFT alignments (preserves variable
columns that trimAl -automated1 was stripping out)."""
from pathlib import Path

ALN_DIR  = Path('results/phylo/alignments')
OUT_FA   = Path('results/phylo/supermatrix.faa')
OUT_PART = Path('results/phylo/supermatrix.partition.txt')
INGROUP  = ['Cx_molestus', 'Cx_pallens', 'Cx_pipiens', 'Cx_quinquefasciatus']

def parse_fasta(p):
    seqs, name, chunks = {}, None, []
    for line in open(p):
        line = line.rstrip()
        if line.startswith('>'):
            if name is not None: seqs[name] = ''.join(chunks)
            name = line[1:].split()[0]; chunks = []
        else:
            chunks.append(line)
    if name is not None: seqs[name] = ''.join(chunks)
    return seqs

def informative_cols(seqs):
    L = len(next(iter(seqs.values())))
    n = 0
    for j in range(L):
        col = {s[j] for s in seqs.values() if s[j] not in '-X*?'}
        if len(col) >= 2: n += 1
    return n

concat = {sp: [] for sp in INGROUP}
parts = []
pos = 1
kept = miss = uninf = total = 0
total_info_cols = 0

for aln in sorted(ALN_DIR.glob('*.aln')):
    total += 1
    seqs = parse_fasta(aln)
    if not all(sp in seqs for sp in INGROUP):
        miss += 1; continue
    sub = {sp: seqs[sp] for sp in INGROUP}
    info = informative_cols(sub)
    if info == 0:
        uninf += 1; continue
    total_info_cols += info
    L = len(sub[INGROUP[0]])
    for sp in INGROUP: concat[sp].append(sub[sp])
    parts.append(f'AA, {aln.stem} = {pos}-{pos+L-1}')
    pos += L; kept += 1

print(f"Total: {total} | Missing species: {miss} | Uninformative: {uninf} | Kept: {kept}")
print(f"Supermatrix length: {pos-1:,} columns")
print(f"Total parsimony-informative-ish columns across kept alignments: {total_info_cols:,}")

OUT_FA.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_FA, 'w') as f:
    for sp in INGROUP:
        f.write(f'>{sp}\n')
        s = ''.join(concat[sp])
        for i in range(0, len(s), 80): f.write(s[i:i+80] + '\n')
with open(OUT_PART, 'w') as f:
    f.write('\n'.join(parts) + '\n')

print(f"Wrote {OUT_FA}  ({OUT_FA.stat().st_size:,} bytes)")
