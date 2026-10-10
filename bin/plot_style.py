"""
Common look of the PNG plots in report.html (gene_set_statistics.py,
fantasia_summary.py, paludamentum_report.py).

report.html shows a PNG at its physical size (pixels / dpi, 96 CSS pixels per
inch), at most as wide as the page. Text of the same point size is then
equally large in every plot, about as large as the text of the page (15 px).
FIG_WIDTH fills the page of the report (1008 px content width, 14 px of it
padding and border of the figure).
"""
from __future__ import annotations

FIG_WIDTH = 10.3       # inches
DPI = 200

# seaborn "deep"; one colour per meaning across all plots
BLUE = "#4C72B0"
GREEN = "#55A868"
ORANGE = "#DD8452"
RED = "#C44E52"
PURPLE = "#8172B3"
GREY = "#8C8C8C"
TEXT = "#1f1f1f"
MUTED = "#52514e"
MEDIAN = "#1f1f1f"     # dashed median line

# colours of the BUSCO plot (generate_plot.py of BUSCO)
BUSCO_SINGLE = "#56B4E9"
BUSCO_DUPLICATED = "#3492C7"
BUSCO_FRAGMENTED = "#F0E442"
BUSCO_MISSING = "#F04442"

RC = {
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.titleweight": "normal",
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "legend.title_fontsize": 11,
    "legend.frameon": False,
    "figure.titlesize": 13,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": "#bdbcb6",
    "axes.labelcolor": TEXT,
    "axes.titlecolor": TEXT,
    "text.color": TEXT,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "savefig.dpi": DPI,
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
}
ANNOTATION_SIZE = 9.5   # numbers written on bars and wedges


def apply(plt) -> None:
    """Set the common rcParams on matplotlib.pyplot."""
    plt.rcParams.update(RC)


def text_colour(hex_colour: str) -> str:
    """Dark or white text, whichever reads better on the colour."""
    h = hex_colour.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return TEXT if luminance > 0.5 else "white"
