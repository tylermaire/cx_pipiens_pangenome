"""
Shared style for all Cx. pipiens pangenome figures.

Conventions chosen to match MDPI Insects journal requirements (Arial sans-serif,
TrueType-embedded PDF, 600 DPI raster output) and to be colorblind-safe.
Import `apply_style()` at the top of every figure script.
"""
from __future__ import annotations
import matplotlib as mpl
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# Colour palette - Okabe-Ito (colourblind-safe) with one ingroup-per-colour
# Ordered consistently with the species tree topology recovered from the
# fixed-SCO supermatrix: pallens + quinquefasciatus are sister taxa;
# molestus + pipiens form the second bipartition.
# ----------------------------------------------------------------------------
GENOME_COLORS = {
    "Cx_quinquefasciatus": "#0072B2",  # blue
    "Cx_pallens":          "#56B4E9",  # sky blue (same bipartition as quinq)
    "Cx_molestus":         "#E69F00",  # orange
    "Cx_pipiens":          "#D55E00",  # vermillion (same bipartition as molestus)
    "Cx_tarsalis":         "#999999",  # grey (outgroup)
}
GENOME_PRETTY = {
    "Cx_quinquefasciatus": "Cx. quinquefasciatus",
    "Cx_pallens":          "Cx. pipiens f. pallens",
    "Cx_molestus":         "Cx. pipiens f. molestus",
    "Cx_pipiens":          "Cx. pipiens f. pipiens",
    "Cx_tarsalis":         "Cx. tarsalis",
}
INGROUP = ["Cx_molestus", "Cx_pallens", "Cx_pipiens", "Cx_quinquefasciatus"]

# Compartment colours (core / shell / cloud) - muted blue/orange/grey
COMPARTMENT_COLORS = {
    "core":  "#1F77B4",
    "shell": "#FF7F0E",
    "cloud": "#7F7F7F",
}

# CAFE5 expansion / contraction
EXPANSION_COLOR   = "#2CA02C"   # green
CONTRACTION_COLOR = "#D62728"   # red

# Synteny forward / inverted (mummer convention)
FORWARD_COLOR  = "#1B1B1B"      # near-black for collinear
INVERTED_COLOR = "#D62728"      # red for inversions

# Figure dimensions (MDPI Insects column widths)
WIDTH_SINGLE_MM = 85.0
WIDTH_DOUBLE_MM = 180.0
MM_TO_INCH = 1 / 25.4

def fig_size(width_mm: float, height_mm: float):
    """Convert mm dimensions to (width_in, height_in) for matplotlib."""
    return (width_mm * MM_TO_INCH, height_mm * MM_TO_INCH)


def apply_style():
    """Set matplotlib rcParams for journal-quality output."""
    mpl.rcParams.update({
        # Font
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.titlesize": 11,
        # PDF embedding (TrueType, not Type 3)
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        # Lines and spines
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        # Figure
        "figure.dpi": 100,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        # Legends
        "legend.frameon": False,
        "legend.borderpad": 0.4,
    })


def save_publication(fig, basename: str, formats=("pdf", "png"), tiff=False):
    """Save a figure in multiple formats. `basename` should have no extension.

    Outputs:
      <basename>.pdf  - vector for publication
      <basename>.png  - quick preview at 300 dpi
      <basename>.tiff - 600 dpi LZW-compressed (if tiff=True)
    """
    if "pdf" in formats:
        fig.savefig(f"{basename}.pdf", format="pdf")
    if "png" in formats:
        fig.savefig(f"{basename}.png", format="png", dpi=300)
    if tiff or "tiff" in formats:
        fig.savefig(f"{basename}.tiff", format="tiff", dpi=600,
                    pil_kwargs={"compression": "tiff_lzw"})


def panel_letter(ax, letter: str, x: float = -0.12, y: float = 1.04,
                 fontsize: float = 12, weight: str = "bold"):
    """Add a panel letter (A, B, C, …) at axis-relative coordinates."""
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=fontsize, fontweight=weight,
            va="bottom", ha="left")
