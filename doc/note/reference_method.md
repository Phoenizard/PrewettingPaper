# 参考代码方法自查清单 (reference/ternary_ref_code)

复现自查用,非教学文档。精确到公式、函数签名、参数值,能照着写等价代码。

复现目标:给定一组模型参数,在 $(\phi_{1,\infty}, \phi_{2,\infty})$ 平面找预润湿相变线。
方法 = 一维线扫描 + 迟滞(薄/厚膜双分支) + 巨势相等($\Omega$ 交叉)。
三文件分工:solver.py 单点求解;thermodynamics.py 物理核心;pipeline.py 编排。

死代码警告:pipeline.py 的 process_single_task 在 line 786 无条件 return,786 行后是旧单方向实现,忽略。复现以 189-786 的双方向实现为准。

---

## 1. 物理模型

$\phi_3 = 1 - \phi_1 - \phi_2$ (溶剂)。

体相自由能密度 (get_free_energy):

$$f_b = \frac{\phi_1\ln\phi_1}{n_1} + \frac{\phi_2\ln\phi_2}{n_2} + \frac{\phi_3\ln\phi_3}{n_3} + \chi_{12}\phi_1\phi_2 + \chi_{13}\phi_1\phi_3 + \chi_{23}\phi_2\phi_3$$

交换化学势 (get_chemical_potential / mu1_2d / mu2_2d):

$$\mu_1 = \frac{\ln\phi_1 + 1}{n_1} - \frac{\ln\phi_3 + 1}{n_3} + \chi_{13}(\phi_3 - \phi_1) + (\chi_{12}-\chi_{23})\phi_2$$

$$\mu_2 = \frac{\ln\phi_2 + 1}{n_2} - \frac{\ln\phi_3 + 1}{n_3} + \chi_{23}(\phi_3 - \phi_2) + (\chi_{12}-\chi_{13})\phi_1$$

体相二阶导 (get_bulk_derivatives,雅可比用):

$$f_{11} = \frac{1}{n_1\phi_1} + \frac{1}{n_3\phi_3} - 2\chi_{13}$$

$$f_{12} = \frac{1}{n_3\phi_3} + \chi_{12} - \chi_{13} - \chi_{23}$$

$$f_{22} = \frac{1}{n_2\phi_2} + \frac{1}{n_3\phi_3} - 2\chi_{23}$$

EL 方程:$\kappa_{ii}\,\phi_i'' = \mu_i(\phi(z)) - \mu_{i,\infty}$。
BC:$z=0$ 墙面 Neumann + 表面场;$z=L$ 体相 Dirichlet $\phi_i(L)=\phi_{i,\infty}$。

---

## 2. 单点求解器 TernaryCHSolver (solver.py)

- 离散:$z=\mathrm{linspace}(0,L,N)$,$h=L/(N-1)$。未知 $U=[\phi_1(z_0..z_{N-1}), \phi_2(...)]$,长度 $2N$。
- 矩阵 A (_build_matrix_A):块 $[[A_{11},A_{12}],[A_{21},A_{22}]]$,每块 $\kappa/h^2$ 三对角。第 0 行 Neumann:M[0,0]=2k, M[0,1]=-2k。
- 牛顿迭代 (config optimization_mode=attempt_c → _solve_opt):
  - 残差 $\mathrm{Res} = A\cdot U - \mathrm{rhs}(U)$。rhs 装 $(\mu_\infty - \mu(U))$;第 0 行减墙面场 $\mathrm{fac}\cdot hs_i$ (fac=2/h);第 N-1 行加体相 Dirichlet $k_f\cdot\kappa\cdot\phi_\infty$ (k_f=1/h^2)。
  - 墙面场 $hs_1 = \omega_1 + 2\chi_{bb,11}\phi_1(0) + \chi_{bb,12}\phi_2(0)$;$hs_2$ 对称。
  - 雅可比 = A + J_bulk(对角 f11/f22,离对角 f12) + J_boundary(墙面 chibb 项)。
  - dU = spsolve;线搜索 $\alpha=1$ 起,破坏可行性 ($\phi_1>0$ & $\phi_2>0$ & $\phi_3>0$) 就 $\alpha\mathrel{*}=0.5$,到 1e-3 止。
  - 收敛:$\max|\mathrm{Res}| < \mathrm{tol}=10^{-8}$;max_iter=500。
- calculate_surface_metrics(U) → ($\Omega$, cs1, cs2):
  - $W = f - f_\infty - \mu_{1,\infty}(\phi_1-\phi_{1,\infty}) - \mu_{2,\infty}(\phi_2-\phi_{2,\infty})$
  - $f_{grad} = \tfrac{1}{2}\kappa_{11}\phi_1'^2 + \tfrac{1}{2}\kappa_{22}\phi_2'^2 + \kappa_{12}\phi_1'\phi_2'$ (梯度用 np.gradient(., h))
  - $J_0 = \omega_1\phi_1(0) + \omega_2\phi_2(0) + \chi_{bb,11}\phi_1(0)^2 + \chi_{bb,22}\phi_2(0)^2 + \chi_{bb,12}\phi_1(0)\phi_2(0)$
  - $\Omega = \mathrm{trapz}(W + f_{grad}, z) + J_0$
  - $cs_i = \mathrm{trapz}(\phi_i - \phi_{i,\infty}, z)$
- 关键:薄/厚膜完全由 U_init 决定,求解器本身不区分。

---

## 3. 预润湿判定:迟滞 + Omega 交叉 (thermodynamics.py 核心)

_run_hysteresis_scan (固定 $\phi_2$ 扫 $\phi_1$):

1. 正扫(薄膜):$\phi_1$ 小→大。初猜 _build_branch_initial_guess(mode=thin, 振幅 dphi1=0.015, dphi2=0.004)。每步用上一收敛解续种子。记 $\Omega_{thin}$, $cs_{thin}$;不收敛记 NaN 并沿用旧种子。
2. 反扫(厚膜):$\phi_1$ 大→小。初猜 mode=thick (dphi1=0.09, dphi2=0.02)。同样续种子。末尾 reverse 回正序。
3. $\mathrm{omega\_diff} = \Omega_{thin} - \Omega_{thick}$;$\mathrm{cs\_gap} = |cs_{thick} - cs_{thin}|$。

初猜构造 _build_branch_initial_guess:layer $= 1 - z/L$;dphi_surface 符号 $= \mathrm{sign}(-\omega_i)$ ($\omega<0 \to +1$);高浓时按 free_volume 缩放避免出界;$\phi_i = \phi_{i,\infty} + d\phi_i\cdot\mathrm{scale}\cdot\mathrm{layer}$;投影回可行域。

_extract_prewetting_crossings:
- valid_mask = (cs_gap > cs_threshold) & 两分支收敛 & 有限。
- 需 sum(valid_mask) $\geq$ min_hysteresis_points (=2)。
- 在连续有效区段内找 omega_diff 变号,线性插值零点 $\phi_{1,eq}$ (_interpolate_zero_crossing)。这就是相变 $\phi_1$。
- allow_terminal_boundary_crossing:相变点落区间边界时,relax 到 terminal_relax_ratio=0.67 倍阈值处理。

find_equilibrium_prewetting 额外有 auto_expand_phi1_bounds:无交叉但有效点近边界且 |omega_diff| 小,则外扩 $\phi_1$ 范围重扫一次。最后过滤物理点 ($0<\phi_1<1$ 且 $\phi_1+\phi_2<1$)。返回 PrewettingResult。

find_equilibrium_prewetting_fix_phi1:对称(固定 $\phi_1$ 扫 $\phi_2$),无 auto_expand。

---

## 4. 体相 binodal / demixed 区 (剔除假点)

compute_bulk_phase:
- regular_grid2 非均匀网格:低浓 logspace(→1e-2) + 中段 linspace(1e-2→0.9) + 高浓 linspace(0.9→)。
- hull_array 加高能辅助顶点 (energy_offset=1000) → ConvexHull 取下包络。
- 顶点索引跳跃 > vertex_jump_threshold(25) 处识别 tie-line 端点 → binodal_points。
- tie-line 另用 5 元方程 + fsolve (solve_exact_binodal_point) 严格解,仅初始化用。

compute_bulk_demixed_map:均匀网格算 f,ConvexHull 下包络平面,$f - \mathrm{env} > \mathrm{demix\_tol}$ (1e-5) → demixed。config 分辨率 451。

is_point_in_demixed_region:k 近邻投票 (k_nearest=8, vote_ratio=0.5),多数邻居 demixed 则剔除。

---

## 5. 编排 (pipeline.py)

- 外循环 scan_method=cases:按 sweep_plan.yaml 显式 case。相同 chi 组 binodal 缓存复用 (reuse_binodal_by_chi);binodal 空的 chi 跳过 (skip_chi_without_binodal)。
- 内循环每 case 两方向:
  - 固定 $\phi_2$:phi2_list 每点起 slice → find_equilibrium_prewetting → 接受非 demixed 的 $\phi_{1,eq}$ → {phi2_inf, phi1_eq}。
  - 固定 $\phi_1$:对称 → {phi1_inf, phi2_eq}。
  - slice 间 multiprocessing.Pool (num_workers=30) 并行。
- 端点二分细化 refine_prewetting_endpoints_between_phi2/phi1:"有→无"相邻扫描线间二分 (max_iter=10, tol=2.5e-4) 定位相变线截止端点。
- 输出:每 case 一 JSON (prewetting_by_scan.fix_phi2 / fix_phi1) + index.json。

---

## 6. 参数表 (正式 config sweep_plan.yaml)

| 类别 | 参数 | 值 |
|---|---|---|
| 空间 | L, N | 10.0, 1000 |
| 梯度能 | kappa_11, kappa_22, kappa_12 | 1.0, 1.0, 0.0 |
| 相互作用 | chi_12, chi_13, chi_23 | 0.0, 2.8, 2.6 |
| 聚合度 | n1, n2, n3 | 1.0, 1.0, 1.0 |
| 墙面场 | omega_1, omega_2 | omega_1 主扫 -0.52/-0.50/-0.45,omega_2=-0.30 |
| 墙面二次 | chibb_11 / 22 / 12 | 单因子扫 ±0.025~0.20 (基线 0) |
| 扫描线 | phi_1/phi_2 range | 0.0001→0.2,400 点 linspace |
| 判据 | cs_threshold | 0.1 |
| 判据 | min_hysteresis_points | 2 |
| binodal | binodal_n_grid / bulk_demixed_grid | 801 / 451 |
| binodal | binodal_min_demixed_fraction | 1e-4 |
| 端点细化 | max_iter / tol | 10 / 2.5e-4 |
| 求解器 | optimization_mode / tol / max_iter | attempt_c / 1e-8 / 500 |

unvalidate_date/ 目录命名 (chi12_..__chi13_..__chi23_.. / om1_..__om2_.. / chibb11_..) 与 case 一一对应。

---

## 复现三大坑

1. 迟滞本质是初猜差异:薄/厚膜 = 同一方程两个初猜 (thin 小振幅正扫、thick 大振幅反扫) + 沿扫描线续种子。这是唯一能同时抓两分支的机制。
2. Omega 相等 = 相变:预润湿点是 omega_diff 变号插值零点,且 cs_gap > 0.1 且两分支都收敛。
3. demixed 剔除不能省:k 近邻投票剔除体相两相区内假点。
