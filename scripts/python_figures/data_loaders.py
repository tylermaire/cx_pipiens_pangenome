"""
Data loaders for the Cx. pipiens pangenome figure suite.

Centralises file-path handling so the figure scripts don't repeat parsing logic.
Each loader returns either a pandas DataFrame, a dict, or a Bio.Phylo tree.

Expected paths (relative to the repo root):
  results/pangenome/partitioned_orthogroups.tsv
  results/pangenome/pangenome_summary.tsv
  results/phylo/concat_tree.treefile
  results/phylo/concord.cf.tree
  results/cafe/output/Gamma_change.tab
  results/cafe/output/Gamma_clade_results.txt
  results/synteny/{ref}_vs_{query}.coords      (show-coords -THrd format)
  results/repeats/{genome}/{genome}.fasta.tbl  (RepeatMasker summary)
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np
import re

from style import INGROUP


# ----------------------------------------------------------------------------
# Pangenome partitioning
# ----------------------------------------------------------------------------
def load_partitions(repo_root: Path) -> pd.DataFrame:
    """Per-orthogroup gene counts + compartment assignment."""
    p = repo_root / "results/pangenome/partitioned_orthogroups.tsv"
    df = pd.read_csv(p, sep="\t")
    return df


def load_pangenome_summary(repo_root: Path) -> pd.DataFrame:
    """Compartment-level summary (core / shell / cloud counts)."""
    p = repo_root / "results/pangenome/pangenome_summary.tsv"
    return pd.read_csv(p, sep="\t")


def build_membership_matrix(partitions: pd.DataFrame) -> pd.DataFrame:
    """Boolean orthogroup × genome membership for UpSet."""
    mat = (partitions[INGROUP] > 0).astype(bool)
    mat.index = partitions["Orthogroup"]
    return mat


# ----------------------------------------------------------------------------
# Phylogenetic tree
# ----------------------------------------------------------------------------
def load_tree(repo_root: Path, with_cf: bool = False):
    """Load the IQ-TREE ML tree (Newick). Returns a Bio.Phylo tree object."""
    from Bio import Phylo
    from io import StringIO
    fname = "concord.cf.tree" if with_cf else "concat_tree.treefile"
    p = repo_root / "results/phylo" / fname
    return Phylo.read(p, "newick")


def load_newick_string(repo_root: Path, with_cf: bool = False) -> str:
    fname = "concord.cf.tree" if with_cf else "concat_tree.treefile"
    p = repo_root / "results/phylo" / fname
    return p.read_text().strip()


# ----------------------------------------------------------------------------
# CAFE5 results
# ----------------------------------------------------------------------------
def load_cafe_change(repo_root: Path) -> pd.DataFrame:
    """Per-family, per-branch gene count changes (Gamma_change.tab)."""
    p = repo_root / "results/cafe/output/Gamma_change.tab"
    df = pd.read_csv(p, sep="\t")
    return df


def load_cafe_clade(repo_root: Path) -> pd.DataFrame:
    """Per-branch expansion / contraction summary (Gamma_clade_results.txt)."""
    p = repo_root / "results/cafe/output/Gamma_clade_results.txt"
    return pd.read_csv(p, sep="\t")


def compute_per_branch_counts(change_df: pd.DataFrame) -> pd.DataFrame:
    """Expansions / contractions / no-change per branch from Gamma_change.tab."""
    rows = []
    for col in change_df.columns:
        if col in ("FamilyID", "#FamilyID"):
            continue
        s = pd.to_numeric(change_df[col], errors="coerce").dropna()
        rows.append({
            "branch": col,
            "expansions": int((s > 0).sum()),
            "contractions": int((s < 0).sum()),
            "no_change": int((s == 0).sum()),
            "net_change": int(s.sum()),
        })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# ANI matrix
# ----------------------------------------------------------------------------
def compute_ani_from_coords(repo_root: Path) -> pd.DataFrame:
    """Build pairwise ANI matrix from show-coords output files.

    Falls back to the manually-curated matrix from Table 2 of the paper if the
    coord files don't yet exist.
    """
    syn = repo_root / "results/synteny"
    if not syn.exists():
        return _fallback_ani()
    coord_files = list(syn.glob("*_vs_*.coords"))
    if not coord_files:
        return _fallback_ani()

    # Build ANI as alignment-length-weighted mean identity from -THrd output.
    # Columns: S1 E1 S2 E2 LEN1 LEN2 IDY REF_TAG QRY_TAG  (no header)
    import itertools
    df = pd.DataFrame(100.0, index=INGROUP, columns=INGROUP)
    for ref, qry in itertools.permutations(INGROUP, 2):
        path = syn / f"{ref}_vs_{qry}.coords"
        if not path.exists():
            continue
        cd = pd.read_csv(path, sep="\t", header=None,
                         names=["S1","E1","S2","E2","LEN1","LEN2","IDY","R","Q"])
        if len(cd):
            ani = (cd["IDY"] * cd["LEN1"]).sum() / cd["LEN1"].sum()
            df.loc[ref, qry] = ani
    # Symmetrize by averaging
    out = (df + df.T) / 2
    np.fill_diagonal(out.values, 100.0)
    return out


def _fallback_ani() -> pd.DataFrame:
    """ANI values from Table 2 of the manuscript (manual curation)."""
    data = {
        "Cx_molestus":         {"Cx_molestus": 100.0, "Cx_pallens": 98.2, "Cx_pipiens": 98.2, "Cx_quinquefasciatus": 98.2},
        "Cx_pallens":          {"Cx_molestus": 98.2,  "Cx_pallens": 100.0, "Cx_pipiens": 98.3, "Cx_quinquefasciatus": 97.9},
        "Cx_pipiens":          {"Cx_molestus": 98.2,  "Cx_pallens": 98.3, "Cx_pipiens": 100.0, "Cx_quinquefasciatus": 98.2},
        "Cx_quinquefasciatus": {"Cx_molestus": 98.2,  "Cx_pallens": 97.9, "Cx_pipiens": 98.2, "Cx_quinquefasciatus": 100.0},
    }
    return pd.DataFrame(data)


# ----------------------------------------------------------------------------
# Transposable elements
# ----------------------------------------------------------------------------
def load_te_summary(repo_root: Path) -> pd.DataFrame:
    """Per-genome TE composition from RepeatMasker .tbl files.

    Falls back to manually-curated numbers from the paper if files don't exist.
    """
    rows = []
    for g in INGROUP:
        tbl = repo_root / f"results/repeats/{g}/{g}.fasta.tbl"
        if tbl.exists():
            rows.append(_parse_repeatmasker_tbl(tbl, g))
    if not rows:
        return _fallback_te()
    return pd.DataFrame(rows).set_index("genome")


def _parse_repeatmasker_tbl(path: Path, genome: str) -> dict:
    """Parse a RepeatMasker .tbl file for per-class TE percentages."""
    text = path.read_text()
    out = {"genome": genome}
    patterns = {
        "Retroelements":    r"Retroelements\s+\d+\s+\d+\s+bp\s+([\d.]+)\s+%",
        "DNA_transposons":  r"DNA transposons\s+\d+\s+\d+\s+bp\s+([\d.]+)\s+%",
        "Rolling_circles":  r"Rolling-circles\s+\d+\s+\d+\s+bp\s+([\d.]+)\s+%",
        "Unclassified":     r"Unclassified:\s+\d+\s+\d+\s+bp\s+([\d.]+)\s+%",
        "Simple_repeats":   r"Simple repeats:\s+\d+\s+\d+\s+bp\s+([\d.]+)\s+%",
        "Low_complexity":   r"Low complexity:\s+\d+\s+\d+\s+bp\s+([\d.]+)\s+%",
        "Interspersed":     r"Total interspersed repeats:\s+\d+\s+bp\s+([\d.]+)\s+%",
        "Total_masked":     r"bases masked:\s+\d+\s+bp\s+\(\s*([\d.]+)\s*%\)",
    }
    for k, pat in patterns.items():
        m = re.search(pat, text)
        out[k] = float(m.group(1)) if m else 0.0
    return out


def _fallback_te() -> pd.DataFrame:
    """TE percentages reported in the manuscript (Section 3.5)."""
    return pd.DataFrame({
        "Interspersed": [52.88, 53.91, 52.14, 53.36],
        "Simple_repeats": [3.5, 3.5, 3.5, 3.5],   # paper doesn't break down further; placeholder
        "Low_complexity": [0.8, 0.8, 0.8, 0.8],   # placeholder
    }, index=INGROUP).rename_axis("genome")


# ----------------------------------------------------------------------------
# Protein counts (from per-genome FASTAs)
# ----------------------------------------------------------------------------
def load_protein_counts(repo_root: Path) -> pd.Series:
    """Count proteins in each genome's gffread output, or fall back to paper numbers."""
    counts = {}
    for g in INGROUP:
        fa = repo_root / f"results/proteins/{g}.fa"
        if fa.exists():
            with open(fa) as f:
                counts[g] = sum(1 for line in f if line.startswith(">"))
    if len(counts) == len(INGROUP):
        return pd.Series(counts).reindex(INGROUP)
    return pd.Series({
        "Cx_molestus":         22582,
        "Cx_pallens":          22585,
        "Cx_pipiens":          23174,
        "Cx_quinquefasciatus": 24199,
    })


# ----------------------------------------------------------------------------
# Synteny coords (for dot plots)
# ----------------------------------------------------------------------------
def load_coords(repo_root: Path, ref: str, qry: str) -> pd.DataFrame | None:
    """Load a show-coords -THrd file for one pairwise comparison.

    Returns DataFrame with S1/E1/S2/E2/LEN1/LEN2/IDY/R/Q columns, or None
    if the file doesn't exist.
    """
    path = repo_root / f"results/synteny/{ref}_vs_{qry}.coords"
    if not path.exists():
        return None
    return pd.read_csv(path, sep="\t", header=None,
                      names=["S1","E1","S2","E2","LEN1","LEN2","IDY","R","Q"])
