# CLAUDE.md

Guidance for working in this repository. Keep this file short; link out instead of duplicating.

## Project

Pre-wetting in a ternary mixture (solvent + solute 1 + solute 2) against a chemically
preferential rigid wall: does it nucleate a thin -> thick pre-wetting film, and where in the
`(phi_1_inf, phi_2_inf)` plane. Core goal is to VERIFY someone else's results in `result/`
(their pre-wetting phase maps) with our own independent code — not to explore freely.

Full model and solving condition: [doc/note/project_plan.md](doc/note/project_plan.md).

## Rules

- Run all Python in the `numenv` conda env: `conda run -n numenv python <script>`
  (interpreter `/opt/miniconda3/envs/numenv/bin/python`). Both are pre-approved in
  `.claude/settings.local.json`. Do not use any other interpreter.
- Not a package. The code is plain scripts plus a module directory, not an installable
  library. Do not reintroduce `__init__.py` / `__all__` or a `prewet.`-style namespace.
  Scripts in `scripts/` put `src/` on `sys.path` and import the modules directly
  (`import thermo`, `import binodal`, ...). Run scripts directly.
- Writing conventions (docs and notes): no Markdown bold (`**...**`). Keep LaTeX formulas
  pure math — no Chinese (or other prose) inside a formula; label terms in prose outside
  the math instead.
- Figures: `doc/note/figures/` is for NOTE figures only (pedagogical). Experiment /
  verification results go under `out/`, never a tmp directory. Produce phase-map / binodal
  figures with the numerical code in `src/` (the paper's method), not by hand.

## SSH compute workflow

Each session the user provides an SSH link (e.g. `ssh -p 32829 root@connect.cqa1.seetacloud.com`) —
a CPU box, free to use for heavy compute. Session-start routine on the server:

1. The project lives under `autodl-tmp/`, in a directory with the same name as the GitHub
   repo (`autodl-tmp/PrewettingPaper`). `cd` there.
2. `git pull` first, always. ALWAYS `source /etc/network_turbo` before any server-side
   github op (pull/push/clone/fetch), in the same shell — proactively, not only after a
   hang; without it the server cannot reach github.com:443 (curl 28 timeout). network_turbo
   only speeds github/hf and SLOWS conda/pip, so do not leave it sourced for those.
3. Check the `numenv` conda env exists. If it is missing, STOP and ask the user to install it —
   do not build it yourself.

## Layout

- `src/` — model modules (plain, importable via `sys.path`): `thermo.py` (bulk + surface
  thermodynamics), `binodal.py` (bulk binodal via lower convex hull), `equilibrium.py`
  (§3.3 equilibrium-DE prewetting solver via `solve_bvp`), `plotting.py` (phase-map
  rendering), `cases.py` (result/ dir-name <-> params, case iteration, output-path map).
- Verification outputs mirror `result/` 1:1 under `out/verify/<chi_dir>/<om_dir>/<chibb_dir>/`
  (reuse result/'s exact dir names). Per case: `overlay.png` + `pw_line.csv`; top-level
  `out/verify/SUMMARY.csv`. Run via `scripts/verify.py` (single case, or `--all [N]`).
- `scripts/` — runnable drivers. `binodal_check.py` computes the binodal for topologies
  T-a..T-f and saves phase maps.
- `doc/note/` — derivation and intro notes plus `figures/`.
- `doc/paper/` — reference paper (Omar, Adame, Arana 2020).
- `out/`, `result/` — generated PNGs (`result/` is a pre-existing sweep to cross-check).

## References

- [PROGRESS.md](PROGRESS.md) — running progress log (current status, next steps, dated
  entries); update it as verification advances.
- [doc/note/project_plan.md](doc/note/project_plan.md) — goals, control-variable structure,
  the 6 topologies, the E1-E5 surface-experiment matrix.
- [doc/note/ternary.md](doc/note/ternary.md) — model: `f_b`, `W`, `f_surf`, the Gibbs
  surface energy `gamma`, equilibrium conditions, and the phase-space-integral algorithm.
- [doc/note/prewetting_intro.md](doc/note/prewetting_intro.md) — pedagogical intro
  (wells, phase separation, thin/thick transition).

## Known issues

- No `environment.yml` capturing the `numenv` dependencies (numpy, scipy, matplotlib);
  the env is not reproducible from the repo yet. Deferred.
