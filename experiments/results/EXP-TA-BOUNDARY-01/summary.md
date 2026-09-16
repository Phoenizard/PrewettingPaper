# EXP-TA-BOUNDARY-01

This audit used only the existing `EXP-TA-REENTRY-01` CSV outputs. It did not
run the thin/thick-film boundary-value solver and did not add wall parameters.
The final audit itself took about 9 seconds; total elapsed time, including a
correction that removed appended convex-hull edge vertices from the displayed
binodal point set, was 3 minutes 33 seconds.

## Curve decisions

| common wall affinity | points | two-phase-side points | publication status |
|---:|---:|---:|---|
| -0.50 | 2 | 0 | invalid: insufficient points for a curve |
| -0.48 | 10 | 0 | valid |
| -0.46 | 16 | 0 | valid |
| -0.44 | 35 | 0 | valid |
| -0.10 | 188 | 47 | invalid |
| -0.08 | 127 | 50 | invalid |
| -0.06 | 78 | 52 | invalid |
| -0.04 | 19 | 19 | invalid |

The existing results at -0.02, -0.01 and -0.001 contain no crossings. Following
the approved publication rule, invalid curves are excluded rather than clipped,
shifted or smoothed across the binodal. The publication heatmap records their
line length as zero and leaves distance to the binodal undefined.

Two continuous-curve candidates are supplied: one retains sparse numerical
markers and one shows curves only. Both use shape-preserving interpolation,
remain inside the sampled endpoints and show the bulk two-phase side explicitly.
