#!/usr/bin/env python3
"""Run all 5 figure generation scripts in order."""
from pathlib import Path
import subprocess
import sys

SCRIPT_DIR = Path(__file__).parent
REPO_ROOT  = SCRIPT_DIR.parent.parent.resolve()
OUT_DIR    = REPO_ROOT / "results" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FIGURES = [
    ("figure1.py", "figure1_pangenome"),
    ("figure2.py", "figure2_tree_cafe"),
    ("figure3.py", "figure3_ani_te_proteins"),
    ("figure4.py", "figure4_synteny"),
    ("figure5.py", "figure5_flower"),
]

for script, basename in FIGURES:
    out = OUT_DIR / basename
    print(f"\n=== {script} -> {out}.pdf ===")
    cmd = [sys.executable, str(SCRIPT_DIR / script),
           "--repo-root", str(REPO_ROOT),
           "--out", str(out)]
    print(" ".join(cmd))
    subprocess.run(cmd, check=False)
print("\nDone. Figures written to:", OUT_DIR)
