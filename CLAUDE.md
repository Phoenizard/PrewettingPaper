# CLAUDE.md

Guidance for working in this repository. Keep this file short; link out instead of duplicating.

## Project

Pre-wetting in a ternary mixture (solvent + solute 1 + solute 2) against a chemically
preferential rigid wall: does it nucleate a thin -> thick pre-wetting film, and where in the
`(phi_1_inf, phi_2_inf)` plane.

Current phase: DATA ANALYSIS. The verification phase is over — our independent code was
judged to reproduce the reference results in `data/` correctly. The goal now is to explore
the physical meaning of how the parameters affect pre-wetting, to serve the paper. Target is
a physics journal, so the model and the numerics do not go in the main text.

The four analysis topics (one per control variable group) are in
[doc/analysis/topiclist.md](doc/analysis/topiclist.md); each gets its own note under
`doc/analysis/`. Still in the initial exploration stage — the observables are not yet
pinned down. Observable framing: pre-wetting has two independent dimensions, strength
(the jump in surface excess between thin and thick branch, not readable off a 2D phase
map, set aside) and extent (how much of the `(phi_1_inf, phi_2_inf)` plane it occupies,
measured by the length of the pre-wetting line and its distance to the binodal).

Full model and solving condition: [doc/note/project_plan.md](doc/note/project_plan.md).

## Model: Gibbs surface free energy gamma

Fields: phi_1(z), phi_2(z) on z in [0, inf); solvent phi_s(z) = 1 - phi_1(z) - phi_2(z).
Far-field (z -> inf) reservoir composition phi_{1,inf}, phi_{2,inf}. Surface values phi_i(0).

Bulk free energy density:

$$
f_b(\phi_1,\phi_2) = \frac{k_B T}{\nu}\Big[
\frac{\phi_1}{n_1}\ln\phi_1 + \frac{\phi_2}{n_2}\ln\phi_2 + \phi_s\ln\phi_s
+ \chi_{13}\,\phi_1\phi_s + \chi_{23}\,\phi_2\phi_s + \chi_{12}\,\phi_1\phi_2 \Big]
$$

Grand-potential density (zero and minimal at phi_i = phi_{i,inf}), with
mu_{i,inf} = nu_i (partial f_b / partial phi_i) evaluated at the reservoir:

$$
W(\phi_1,\phi_2) = f_b(\phi_1,\phi_2) - f_b(\phi_{1,\infty},\phi_{2,\infty})
- \sum_{i=1,2}\frac{\mu_{i,\infty}}{\nu_i}\,(\phi_i - \phi_{i,\infty})
$$

Surface interaction energy (depends only on the wall-contact values phi_i(0)):

$$
f_{\text{surf}}(\phi_1(0),\phi_2(0)) = \frac{k_B T}{\bar\nu}\Big[
\omega_1\,\phi_1(0) + \omega_2\,\phi_2(0)
+ \chi_{bb,11}\,\phi_1(0)^2 + \chi_{bb,22}\,\phi_2(0)^2 + \chi_{bb,12}\,\phi_1(0)\phi_2(0) \Big]
$$

Gibbs surface free energy (the output target):

$$
\gamma[\phi_1,\phi_2] = \int_0^\infty dz\Big[
W(\phi_1(z),\phi_2(z)) + \tfrac{1}{2}\kappa_1(\partial_z\phi_1)^2 + \tfrac{1}{2}\kappa_2(\partial_z\phi_2)^2
\Big] + f_{\text{surf}}(\phi_1(0),\phi_2(0))
$$

Parameter roles:
- chi12, chi13, chi23: bulk Flory-Huggins interactions — chi12 solute1-solute2, chi13 solute1-solvent, chi23 solute2-solvent. They set the bulk phase topology (which stage / T-a..T-f).
- omega1, omega2 (wall affinity of solute 1, 2) and chibb11, chibb22, chibb12 (surface-enhanced interactions) enter only f_surf, i.e. the wall boundary condition.
- n1, n2 solute-to-solvent size ratios; kappa1, kappa2 gradient penalties; nu, nubar, kBT scales.
Derivation: [doc/note/ternary.md](doc/note/ternary.md).

## Control variables and their value sets

A case = one (chi topology) x (omega1, omega2) x (chibb11, chibb22, chibb12) point.
Value sets below are the union measured across all 6 stages in `data/`; a given
stage samples only a subset (older stages a regular 11x11 omega grid; newer T-b/T-d/T-f
scan finer / freer points). Encoding in dir names: `_p_` = decimal point, `_m_` = minus.

- chi (chi12, chi13, chi23) — 6 topologies (stage dirs):
  - T-a: (0.0, 2.8, 0.0)
  - T-b: (0.0, 2.8, 2.6)
  - T-c: (2.3, 2.8, 2.6)
  - T-d: (2.6, 2.6, 2.6)
  - T-e: (2.8, 2.8, 2.8)
  - T-f: (-8.5, 0.0, 0.0)
- omega1: range [-1.2, 0.28] (union 64 values; older stages step 0.02 over -0.50..-0.30).
- omega2: range [-1.2, 0.28] (union 66 values; older stages step 0.02 over -0.50..-0.30).
- chibb11: {-0.10, -0.08, -0.06, -0.05, -0.04, -0.025, -0.02, 0, 0.02, 0.025, 0.04, 0.05, 0.06, 0.08, 0.10}.
- chibb22: same set as chibb11.
- chibb12: {-1.0, -0.9, ..., -0.1} plus {-0.08, -0.06, -0.05, -0.04, -0.02, 0, 0.02, 0.04, 0.05, 0.06, 0.08} plus {0.1, 0.2, ..., 0.9}.
Most cases hold chibb = (0, 0, 0) and scan only omega; the chibb sets above appear mainly in
the chibb-sweep stage (T-f, from chi_m8500_chibb_sweep).

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

- `src/` — 七个模型模块（params/model/solver/scan/bulk/pipeline/plotting，含 logutil），
  外加 `pwpix.py`（从对照 overlay PNG 颜色无关提取 pre-wetting 线像素；对照侧只有图、
  无数值 CSV，故需像素提取）。
- `scripts/` — 直接运行的脚本（run_case / plot_case / build_summary / run_verify.sh，
  外加各分析 topic 的脚本）。
- `config/` — 单一 yaml 参数源。
- `doc/analysis/` — 数据分析阶段的记录：`topiclist.md` 是 4 个 topic 的索引，
  每个 topic 一个 note。图不进 doc/analysis/，留在 `out/`。
- `doc/note/` — derivation and intro notes plus `figures/`; `reference_method.md` 是复现
  参考方法的自查笔记。
- `doc/paper/` — reference paper (Omar, Adame, Arana 2020).
- `reference/` — 同组成员可运行的参考实现（只读教材，不 import、不共享）。
- `data/` — 同组成员的 pre-wetting 相图（旧名 result/、unvalidate_data/）。已验真为真，
  现为分析阶段的数据源。每个 case 只有 PNG，无数值 CSV。
- `result_cases.txt` — 770 行三元组（chi 目录 / om 目录 / chibb 目录）case 清单。
- `out/` — 已跑出的结果（PNG + pw_line.csv），被 .gitignore 忽略、不入 git。

## References

- [PROGRESS.md](PROGRESS.md) — running progress log (current status, next steps, dated
  entries); update it as the analysis advances.
- [doc/analysis/topiclist.md](doc/analysis/topiclist.md) — the 4 analysis topics, one per
  control-variable group; index into the per-topic notes.
- [doc/note/project_plan.md](doc/note/project_plan.md) — goals, control-variable structure,
  the 6 topologies, the E1-E5 surface-experiment matrix.
- [doc/note/ternary.md](doc/note/ternary.md) — model: `f_b`, `W`, `f_surf`, the Gibbs
  surface energy `gamma`, equilibrium conditions, and the phase-space-integral algorithm.
- [doc/note/prewetting_intro.md](doc/note/prewetting_intro.md) — pedagogical intro
  (wells, phase separation, thin/thick transition).

## Known issues

- No `environment.yml` capturing the `numenv` dependencies (numpy, scipy, matplotlib);
  the env is not reproducible from the repo yet. Deferred.
