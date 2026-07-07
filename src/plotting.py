"""Phase-map rendering: bulk binodal + pre-wetting branches in the phi1-phi2 plane.

Axis convention matches the existing result/ PNGs: x = phi1, y = phi2, with the
faint anti-diagonal phi1 + phi2 = 1 (the phi_s = 0 edge of the Gibbs triangle).
Scale/proportion are not tuned to match result/ pixel-for-pixel (per user: only
the axis *content* must agree).
"""
from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def render_phase_map(
    path,
    *,
    chi,
    surf=None,
    pw=None,
    binodal=None,
    title=None,
    params_text=None,
):
    """Turn one case's inputs — experiment variables plus computed results —
    into a phase-map PNG. General (variables + results) -> figure entry point,
    independent of any experiment driver: pass the physical parameters and the
    results, get a figure. Callers (verify, future experiments) only assemble
    inputs; they do not touch matplotlib.

    Experiment variables:
      chi     : T.Chi — interaction params; fixes the bulk binodal (required).
      surf    : T.Surf — surface energetics; only recorded in the caption, does
                not change the phase-map geometry. Optional.
    Computed results:
      pw      : (N, 2) array of pre-wetting (phi1_inf, phi2_inf) points, or None.
      binodal : precomputed BinodalResult; derived from `chi` if omitted (cheap).
    """
    if binodal is None:
        import binodal as B  # lazy: only needed when not precomputed

        binodal = B.binodal_from_hull(chi)
    branches = (
        [{"points": pw, "state": "pw", "branch_id": 3}]
        if pw is not None and len(pw)
        else None
    )
    plot_phase_map(path, binodal, pw_branches=branches, title=title,
                   params_text=params_text)


def plot_phase_map(
    path,
    binodal,
    pw_branches=None,
    title=None,
    params_text=None,
    figsize=(6.2, 6.0),
):
    """Render one phase map to `path`.

    binodal   : BinodalResult
    pw_branches: optional list of dicts {'points': (N,2) array, 'state': str,
                 'branch_id': int} — pre-wetting surface-state loci.
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Gibbs triangle boundary
    ax.plot([0, 1], [1, 0], color="0.5", lw=1.0, zorder=1)
    ax.plot([0, 1], [0, 0], color="0.5", lw=0.8, zorder=1)
    ax.plot([0, 0], [0, 1], color="0.5", lw=0.8, zorder=1)

    # tie-lines (faint) + binodal endpoints
    if len(binodal.tie_lines):
        for tl in binodal.tie_lines:
            ax.plot(tl[:, 0], tl[:, 1], color="#4C72B0", lw=0.3, alpha=0.25, zorder=2)
        pts = binodal.points
        ax.scatter(
            pts[:, 0], pts[:, 1], s=3, color="#1f4e9c", zorder=3, label="Binodal"
        )

    # three-phase triangle(s)
    for tri in binodal.three_phase:
        t = np.vstack([tri, tri[0]])
        ax.plot(t[:, 0], t[:, 1], color="green", lw=1.4, ls="--", zorder=4,
                label="3-phase")

    # pre-wetting branches
    if pw_branches:
        cmap = plt.get_cmap("tab10")
        for br in pw_branches:
            p = np.asarray(br["points"])
            if p.size == 0:
                continue
            k = br.get("branch_id", 0)
            lbl = f"PW branch {k}" + (f" ({br['state']})" if "state" in br else "")
            ax.scatter(p[:, 0], p[:, 1], s=10, color=cmap(k % 10),
                       zorder=5, label=lbl)

    ax.set_xlabel(r"$\phi_1$")
    ax.set_ylabel(r"$\phi_2$")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal")
    if title:
        ax.set_title(title)
    # de-duplicate legend labels
    h, l = ax.get_legend_handles_labels()
    seen, hh, ll = set(), [], []
    for hi, li in zip(h, l):
        if li not in seen:
            seen.add(li)
            hh.append(hi)
            ll.append(li)
    if hh:
        ax.legend(hh, ll, fontsize=7, loc="upper right", framealpha=0.9)
    if params_text:
        fig.text(0.5, 0.005, params_text, ha="center", va="bottom", fontsize=7)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(path, dpi=130)
    plt.close(fig)
