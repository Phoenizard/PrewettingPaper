"""Extract prewetting-line metrics from the overlay PNGs (no numeric CSV exists
for the full sweep, only 12 scattered cases). Color-agnostic: the PW line may be
drawn as several colormap-colored branches (navy/green/maroon...), not just red,
so we mask by "saturated, not-binodal, not-guide, not-legend" instead of a fixed hue.

Plot is a Binodal/PW distribution plot on the (phi1, phi2) simplex. Axes are the
black frame: bottom = phi2=0 (phi1: x 0->1), left = phi1=0 (phi2: y 0->1). All
overlay PNGs share one size (1725x1580), so one pixel->phi calibration serves all.
"""

from pathlib import Path

import numpy as np
from PIL import Image

# calibration for the shared 2070x2220 full-overlay size (pixel -> phi on [0,1] axes),
# measured from the plot frame; verified against known cases. Use the FULL overlay
# (not *_standardized_plot, whose sky-blue binodal is not cleanly separable).
CALIB = {"x0": 209, "x1": 2143, "y0": 1739, "y1": 127}  # phi1: x0->x1=0->1 ; phi2(down): y0->y1=0->1


def _load(path):
    return np.asarray(Image.open(path).convert("RGB")).astype(int)


def pw_mask(im):
    """Boolean mask of prewetting-branch pixels (any branch color).

    PW branches use a viridis-like colormap: branch 0 is dark purple ~(68,1,84)
    with low saturation, so the saturation floor must stay low (~28)."""
    H, W = im.shape[:2]
    R, G, B = im[..., 0], im[..., 1], im[..., 2]
    sat = im.max(2) - im.min(2)
    m = sat > 28                                   # colored curve, incl. dark-purple branch 0
    m &= ~((B > 140) & (G > 110) & (R < 120))      # drop teal-blue binodal
    m &= ~((abs(R - G) < 25) & (abs(G - B) < 25))  # drop gray diagonal guide / gridlines
    m[: int(0.28 * H), int(0.55 * W):] = False     # drop legend box (top-right)
    return m


def branch_colors(im, mask):
    """Distinct branch color clusters among PW pixels -> count of branches.
    n_branch = number of separate PW-line branches (NOT number of points on a line;
    point counts depend on numerical scan settings and are not reliable)."""
    pix = im[mask]
    if len(pix) == 0:
        return 0
    q = pix // 48                                  # quantize to merge anti-aliasing
    uniq, counts = np.unique(q.reshape(-1, 3), axis=0, return_counts=True)
    return int((counts >= 40).sum())               # clusters with real line length


def metrics(path):
    """Return dict of PW-line metrics for one overlay PNG, or exists=False."""
    im = _load(path)
    m = pw_mask(im)
    ys, xs = np.where(m)
    if len(xs) == 0:
        return dict(exists=False, n_branch=0, phi2_max=np.nan,
                    phi1_min=np.nan, phi1_max=np.nan, phi2_min=np.nan,
                    pw_coverage=0)
    c = CALIB
    phi1 = (xs - c["x0"]) / (c["x1"] - c["x0"])
    phi2 = (c["y0"] - ys) / (c["y0"] - c["y1"])
    # bounding box in phi-space of the PW pixels, as a coverage proxy (relative only)
    cover = (phi1.max() - phi1.min()) * (phi2.max() - phi2.min())
    return dict(
        exists=True,
        n_branch=branch_colors(im, m),
        phi2_max=float(phi2.max()),
        phi2_min=float(phi2.min()),
        phi1_min=float(phi1.min()),
        phi1_max=float(phi1.max()),
        pw_coverage=float(cover),
    )


def find_overlay(case_dir):
    """The full overlay PNG carries all branches; pick it (not *_standardized_plot)."""
    case_dir = Path(case_dir)
    cands = [p for p in case_dir.glob("overlay_*.png")
             if "standardized" not in p.name]
    return cands[0] if cands else None
