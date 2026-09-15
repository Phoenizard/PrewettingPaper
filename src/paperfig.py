"""Shared publication defaults for the T-a/T-f paper figures."""

from pathlib import Path

import matplotlib.pyplot as plt


PANEL_LABEL_SIZE = 10.0


def configure():
    """Apply typography sized for a two-column journal figure."""
    plt.rcParams.update({
        "font.size": 9.0,
        "axes.labelsize": 10.0,
        "axes.titlesize": 9.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.5,
        "lines.linewidth": 1.5,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
    })


def panel_label(ax, label):
    ax.text(-0.16, 1.04, label, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=PANEL_LABEL_SIZE, fontweight="bold", clip_on=False)


def save_pdf_png(fig, output, dpi=600):
    """Save matching vector PDF and high-resolution PNG files."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    stem = output.with_suffix("")
    pdf = stem.with_suffix(".pdf")
    png = stem.with_suffix(".png")
    fig.savefig(pdf)
    fig.savefig(png, dpi=dpi)
    return pdf, png
