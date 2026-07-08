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
- Not a package (约定，供将来重写时遵循). Keep the code as plain scripts plus a module
  directory, not an installable library. Do not introduce `__init__.py` / `__all__` or a
  `prewet.`-style namespace. Scripts in `scripts/` should put `src/` on `sys.path` and
  import the modules directly. Run scripts directly.
- Writing conventions (docs and notes): no Markdown bold (`**...**`). Keep LaTeX formulas
  pure math — no Chinese (or other prose) inside a formula; label terms in prose outside
  the math instead.
- Figures: `doc/note/figures/` is for NOTE figures only (pedagogical). Experiment /
  verification results go under `out/`, never a tmp directory. Produce phase-map / binodal
  figures with the numerical code (the paper's method), not by hand.

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

代码已清空（`src/`、`scripts/` 无任何 .py，历次实现已删；经过见 PROGRESS.md）。
当前仓库只有文档、参考实现、数据与已跑出的结果：

- `doc/note/` — derivation and intro notes plus `figures/`; `reference_method.md` 是复现
  参考方法的自查笔记。
- `doc/paper/` — reference paper (Omar, Adame, Arana 2020).
- `reference/` — 同组成员可运行的参考实现（只读教材，不 import、不共享）。
- `unvalidate_data/` — 同组成员的 pre-wetting 相图（旧名 result/，要独立复现的对照）。
- `result_cases.txt` — 770 行三元组（chi 目录 / om 目录 / chibb 目录）case 清单。
- `out/` — 已跑出的结果（PNG + pw_line.csv），被 .gitignore 忽略、不入 git。

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
