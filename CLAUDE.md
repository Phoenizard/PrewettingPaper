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

An analysis topic is bound to one chi stage: a topic = one chi topology plus the question
asked of it. There is no fixed topic list and no mapping from control-variable groups to
topics; topics are opened one at a time as the physics warrants. Which topics exist, which
are finished and which one is in progress is recorded only in
[PROGRESS.md](PROGRESS.md) under `当前状态` — read it there and do not infer, enumerate or
invent topics beyond what it states.

Observable framing (established on the omega topic, not assumed to carry over): pre-wetting
has two independent dimensions — strength (the jump in surface excess between thin and thick
branch, not readable off a 2D phase map, set aside) and extent (how much of the
`(phi_1_inf, phi_2_inf)` plane it occupies, measured by the length of the pre-wetting line
and its distance to the binodal).

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

Source of truth for the model, the solving condition and the equilibrium equations is
[doc/note/project_plan.md](doc/note/project_plan.md). (`doc/note/ternary.md` had errors in
the details and was deleted — do not restore or cite it.) The plan writes some symbols
differently; the correspondence is chi_1s = chi13, chi_2s = chi23, omega_{b,i} = omega_i,
chibb_1 = chibb11, chibb_2 = chibb22. This file keeps the chi13 / omega1 / chibb11 spelling
because it matches the directory names and the code variables.

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
- LaTeX rendering (chat replies AND notes). The renderer is strict; these break it, so they
  are banned:
  - `\text{...}` / `\mbox{...}` containing CJK. Never put Chinese inside math. Name the
    quantity with a plain symbol and explain it in prose outside the formula.
  - `\dfrac`, `\displaystyle`, and any fraction nested inside another fraction's numerator
    or denominator. Use only flat `\frac{a}{b}`; if a term is itself a ratio, write it with
    a slash, e.g. `\partial W/\partial \phi_1`.
  - Manual spacing macros inside formulas: `\ `, `\;`, `\,`, `\!`, `\quad`, `\qquad`.
  - `\Big` / `\big` / `\left` / `\right` sizing. Use plain `(`, `[`, `|`.
  - Multiple relations or a list of definitions crammed into one display block.
  Allowed: `\frac`, `\int`, `\sum`, `\partial`, `\ln`, `\phi`, `\mu`, `\chi`, `\omega`,
  `\Delta`, `^`, `_`, `=`, `+`, `-`, plain brackets. One equation per display block.
- Figures: `doc/note/figures/` is for NOTE figures only (pedagogical). Experiment /
  verification results go under `out/`, never a tmp directory. Produce phase-map / binodal
  figures with the numerical code (the paper's method), not by hand. `out/analysis/` is
  tracked by git (analysis notes plus the figures and CSVs they reference), so keep what
  lands there presentable; the rest of `out/` stays ignored.

## SSH data-access machine (fixed)

`ssh -p 23472 root@connect.westb.seetacloud.com` — passwordless, address long-lived
(not expected to change). Purpose: data access and reading ONLY, no compute. It mounts the
shared storage container at `/root/autodl-fs/pw-space`, where teammates upload experiment
data per [doc/data_format.md](doc/data_format.md). Distinct from the per-session compute
box below.

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
- `doc/analysis/` — 早期分析笔记（`angle1_omega.md`、`measure_extent_debate.md`）。
  `topiclist.md` 已作废（那份「4 个 topic 按控制变量分完」的划分不再成立），不要引用。
  新的 topic 记录一律写在 `out/analysis/<topic 目录>/` 里，图与笔记放在一起。
- `doc/note/` — `project_plan.md`（模型与求解条件的唯一出处）、intro note 与 `figures/`；
  `reference_method.md` 是复现参考方法的自查笔记。
- `doc/paper/` — reference paper (Omar, Adame, Arana 2020).
- `manuscript/` — 论文 LaTeX 项目（独立嵌套 git 仓库，见下节）。
- `reference/` — 同组成员可运行的参考实现（只读教材，不 import、不共享）。
- `data/` — 同组成员的 pre-wetting 相图（旧名 result/、unvalidate_data/）。已验真为真，
  现为分析阶段的数据源。每个 case 只有 PNG，无数值 CSV。
- `result_cases.txt` — 770 行三元组（chi 目录 / om 目录 / chibb 目录）case 清单。
- `out/` — 已跑出的结果（PNG + pw_line.csv）。`.gitignore` 是 `out/*` 加例外
  `!out/analysis/`：只有 `out/analysis/` 入 git（分析笔记及其引用的图与 CSV），
  其余结果留本地。

## Paper manuscript (`manuscript/`)

- 独立嵌套 git 仓库，remote 是 Overleaf（git.overleaf.com），与本仓库互不可见：
  外层 `.gitignore` 忽略整个 `manuscript/`，GitHub 不追踪它；它的提交与推送
  一律由用户手动执行（与本仓库的 git 一样，Claude 不代做）。
- `main.tex` 是唯一写作中心：全部正文都写在这个文件里，不拆 method/result 等多个 tex。
- `NPJ/` 是期刊模版原件（注意：文件夹名叫 NPJ，内容实为 IOP Publishing 模版），
  只读参考，`iopjournal-template.tex` 等一律不改；根目录 `iopjournal.cls` 是供
  main.tex 使用的副本。
- 编译验证在 Overleaf 端（push 后网页编译），本地不装 LaTeX 工具链。

## References

- [PROGRESS.md](PROGRESS.md) — running progress log (current status, next steps, dated
  entries); update it as the analysis advances.
  It is also the only place that says which analysis topics exist and where each one stands.
- [doc/note/project_plan.md](doc/note/project_plan.md) — the model and solving condition
  (`f_b`, `W`, `f_surf`, `gamma`, the equilibrium equations and their boundary conditions,
  the first integral), goals, the three-layer control-variable structure, the 6 topologies,
  and the stated limitations. Answer model questions from this file.
- [doc/note/prewetting_intro.md](doc/note/prewetting_intro.md) — pedagogical intro
  (wells, phase separation, thin/thick transition).

## Known issues

- No `environment.yml` capturing the `numenv` dependencies (numpy, scipy, matplotlib);
  the env is not reproducible from the repo yet. Deferred.
