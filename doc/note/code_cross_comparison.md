# 两份参考代码交叉审计（静态）

审计对象：reference/ternary_ref_code（下称 A 码）与 reference/prewetting_project_clean
（下称 B 码）。unvalidate_data/ 的相图即这两份代码的产物。本文回答三个问题：
两码解的模型是否同为 doc/note/project_plan.md 的模型；两码数值方案是否一致
（即 doc/note/reference_method.md 所记方法）；两码中间参数差在哪里。
纯静态代码审计，未做任何数值运行。所有 file:line 相对各自代码根目录。

## 0. 结论摘要

- 模型：一致。两码逐公式相同（体相自由能、化学势、Hessian、EL 方程、两端边界条件、
  巨势、表面能、表面过剩），且与 project_plan.md 的模型一致；均取无量纲形式
  $k_BT/\nu = k_BT/\bar\nu = 1$。
- 方法：主干一致。固定网格 + 解析雅可比 Newton + 可行性线搜索 + 薄膜正扫/厚膜反扫续种子
  + $\Omega_{thin}-\Omega_{thick}$ 变号线性插值 + $|cs_{thick}-cs_{thin}|>0.1$ 掩码，
  与 reference_method.md 记载一致。差异集中在辅助机制：A 码多 demixed 剔除、端点二分细化、
  terminal 边界放宽、auto expand；B 码取第一个交叉为主线、另存全部交叉、事后手工剔除。
- 参数：核心数值参数全同（L=10、N=1000、kappa=1/1/0、n=1/1/1、tol=1e-8、max_iter=500、
  cs 阈值 0.1、扫描范围 [0.0001, 0.2]、thin 初猜 0.015/0.004）。真差异两处：
  每条扫描线点数 A 码 400 vs B 码 150；thick 初猜振幅 A 码 0.09/0.02 vs B 码 0.3/0.3。
  两码生产运行的 (chi, omega) 参数点不重叠（A 码跑 T-b + chibb 扫描，B 码跑 T-a/T-c/T-e
  的 omega 网格），故互印证属于同模型同方法层面，不是同点复算。

## 1. 模型对表

约定 $\phi_3 = 1 - \phi_1 - \phi_2$（溶剂），$n_1=n_2=n_3=1$（生产值）。
两码表达式逐项相同，只列一份公式，两列给出处。

| 量 | 表达式 | A 码 | B 码 | 一致 |
|---|---|---|---|---|
| 体相自由能密度 | $f_b = \sum_i \frac{\phi_i\ln\phi_i}{n_i} + \chi_{12}\phi_1\phi_2 + \chi_{13}\phi_1\phi_3 + \chi_{23}\phi_2\phi_3$ | solver.py:267-268 | criterion.py:15-22 | 是 |
| 交换化学势 | $\mu_1 = \frac{\ln\phi_1+1}{n_1} - \frac{\ln\phi_3+1}{n_3} + \chi_{13}(\phi_3-\phi_1) + (\chi_{12}-\chi_{23})\phi_2$，$\mu_2$ 对称 | solver.py:80-81 | criterion.py:36-39 | 是 |
| Hessian | $f_{11} = \frac{1}{n_1\phi_1} + \frac{1}{n_3\phi_3} - 2\chi_{13}$，$f_{12} = \frac{1}{n_3\phi_3} + \chi_{12} - \chi_{13} - \chi_{23}$，$f_{22}$ 对称 | solver.py:87-89 | solver.py:52-54 | 是 |
| EL 方程 | $\kappa_{ii}\,\phi_i''(z) = \mu_i(\phi(z)) - \mu_{i,\infty}$，经 $F(U) = AU - b(U) = 0$ 离散 | solver.py:92-108 | solver.py:57-81 | 是 |
| 墙面边界 $z=0$ | Neumann 差分行 $M[0,0]=2k$、$M[0,1]=-2k$；$hs_i = \omega_i + 2\chi_{bb,ii}\phi_i(0) + \chi_{bb,12}\phi_j(0)$，$b_i[0] \mathrel{-}= (2/h)\,hs_i$ | solver.py:98-101 | solver.py:35-36, 68-71 | 是 |
| 远场边界 $z=L$ | Dirichlet $\phi_i(L)=\phi_{i,\infty}$，经 RHS 第 $N-1$ 行 | solver.py:103-104 | solver.py:73-75 | 是（见注 1） |
| 参考势 | $W = f_b - f_{b,\infty} - \mu_{1,\infty}(\phi_1-\phi_{1,\infty}) - \mu_{2,\infty}(\phi_2-\phi_{2,\infty})$ | solver.py:286 | solver.py:103 | 是 |
| 梯度能 | $f_{grad} = \tfrac{1}{2}\kappa_{11}\phi_1'^2 + \tfrac{1}{2}\kappa_{22}\phi_2'^2 + \kappa_{12}\phi_1'\phi_2'$（$\kappa_{12}$ 项仅 A 码有），导数 np.gradient(., h) | solver.py:287-289 | solver.py:104-106 | 见注 1 |
| 表面能 | $J_0 = \omega_1\phi_1(0) + \omega_2\phi_2(0) + \chi_{bb,11}\phi_1(0)^2 + \chi_{bb,22}\phi_2(0)^2 + \chi_{bb,12}\phi_1(0)\phi_2(0)$ | solver.py:293-299 | solver.py:110-117 | 是 |
| 巨势 | $\Omega = \mathrm{trapz}(W + f_{grad},\, z) + J_0$ | solver.py:290, 304 | solver.py:108, 118 | 是 |
| 表面过剩 | $cs_i = \mathrm{trapz}(\phi_i - \phi_{i,\infty},\, z)$，判据用总量 $cs_1 + cs_2$ | solver.py:302-303; thermodynamics.py:651, 676 | criterion.py:50-52 | 是 |

注 1：$\kappa_{12}$ 相关项（梯度能交叉项、$z=L$ RHS 中的 $\kappa_{12}\phi_{j,\infty}$ 项）
仅 A 码实现，B 码没有；两码生产配置均为 $\kappa_{12}=0$，故无数值差别。

与 project_plan.md 的对应：plan 记号 $\chi_{1s},\chi_{2s}$ 即代码的 chi_13、chi_23；
plan 的 $\gamma$ 即代码的 $\Omega$，plan 的 $f_{\text{surf}}$ 即代码的 $J_0$；
plan 公式带 $k_BT/\nu$、$k_BT/\bar\nu$ 前因子，两码取 1（无量纲化），公式结构逐项对应。
plan 的墙面条件 $\kappa_i\phi_i'(0) = \partial f_{\text{surf}}/\partial\phi_i(0)$ 正是
两码 Neumann 行 + $hs_i$ 的连续极限（$hs_i = \partial J_0/\partial\phi_i(0)$）。
结论：两码解的就是 project_plan.md 的模型。

## 2. 数值方法对表

对照 reference_method.md 的清单，逐机制标注。

| 机制 | A 码 | B 码 | 判定 |
|---|---|---|---|
| 离散 | $z=\mathrm{linspace}(0,L,N)$，未知量 $U=[\phi_1;\phi_2]$ 长 $2N$，三对角块矩阵 | 同（solver.py:14-43） | 一致 |
| Newton 残差 | $F(U) = AU - b(U)$，$b$ 装 $\mu_\infty-\mu(U)$ + 两端边界项 | 同（solver.py:57-81） | 一致 |
| 雅可比 | $A$ + 体相 Hessian 块 + 边界 $\chi_{bb}$ 项 | $A$ + 体相 Hessian 块，无边界项（solver.py:137-144） | 差异（注 2） |
| 线性解 | spsolve | 同（solver.py:147） | 一致 |
| 线搜索 | $\alpha$ 从 1 减半至 $10^{-3}$，可行性判 $\phi_1>0 \wedge \phi_2>0 \wedge \phi_3>0$ | 同（solver.py:152-165） | 一致 |
| 收敛判据 | $\max|F| < 10^{-8}$，max_iter=500 | 同（solver.py:128-135） | 一致 |
| 初猜形状 | 线性 layer $= 1 - z/L$，符号 $\omega_i<0$ 取 $+1$（thermodynamics.py:340-357） | 同形状同符号（initial_guess.py:48-70，生产 method=linear_decay，main.py:248-255 自动转 reused_linear） | 一致 |
| 初猜振幅 thin | 0.015 / 0.004（硬编码，thermodynamics.py:353-354） | 0.015 / 0.004（default.yaml:30-31） | 一致 |
| 初猜振幅 thick | 0.09 / 0.02（硬编码，thermodynamics.py:350-351） | 0.3 / 0.3（default.yaml:78-79 覆盖 initial_guess.py:62-63 的默认 0.09/0.02） | 差异 |
| 初猜高浓保护 | 富集量按 $0.8\times$ 自由体积缩放（thermodynamics.py:361-365）再投影 | 只做事后可行域投影（initial_guess.py:4-16） | 差异（小） |
| 双向迟滞 | 薄膜正扫 + 厚膜反扫，各自续种子；首点冷启 | 同（main.py:107-153） | 一致 |
| 不收敛处理 | 记 NaN、保留旧种子继续 | 记 inf/NaN、不更新种子（等价于保留旧种子）（main.py:119-142） | 一致 |
| 判据量 | $\omega_{diff} = \Omega_{thin} - \Omega_{thick}$；$cs_{gap} = |cs_{thick} - cs_{thin}|$，$cs = cs_1 + cs_2$ | 同（criterion.py:65-67） | 一致 |
| 有效掩码 | $cs_{gap} > 0.1$ 且两分支收敛且有限；有效点 $\geq 2$ | 同（cs 为 NaN 的失败点自动落出掩码；criterion.py:67, 110-114） | 一致 |
| 交叉定位 | 连续有效区段内找 $\omega_{diff}$ 变号，线性插值零点 | 掩码点子序列上找变号，线性插值（criterion.py:106-130） | 差异（注 3） |
| 多交叉取舍 | 收集全部交叉，再过滤物理点 + demixed | 主线取第一个交叉（criterion.py:121），全部交叉另存 prewetting_points_all.csv | 差异 |
| terminal 边界放宽 | 有（阈值 $\times 0.67$，thermodynamics.py:759-777） | 无 | 仅 A 码 |
| auto expand 扫描界 | 有（thermodynamics.py:872-906） | 无 | 仅 A 码 |
| demixed 剔除 | 有（k 近邻 k=8、投票 0.5，thermodynamics.py:542-568） | 无（事后手工剔除脚本 scripts/result_processing/apply_manual_removed_loops_case.py） | 仅 A 码 |
| 端点二分细化 | 有（max_iter=10、tol=2.5e-4，thermodynamics.py:1190-1294） | 无 | 仅 A 码 |
| binodal 计算 | 非均匀网格 + ConvexHull 下包络 + demixed 图 | 均匀跨段网格 + ConvexHull 下包络 + 顶点跳跃阈值 25 + 组分下限（phase_diagram.py:343-357） | 差异（注 4） |
| 双方向扫描 | fix_phi2 扫 $\phi_1$ + fix_phi1 扫 $\phi_2$ | 同（main.py:107-202） | 一致 |
| 并行 | multiprocessing.Pool，slice 间并行 | ProcessPoolExecutor，slice 间并行 | 一致（不影响结果） |
| 备用分支 | 死代码（pipeline.py:786 后） | Picard 求解器、exponential_decay 初猜、wall_jump 判据（均非生产路径） | 不入审计 |

注 2：B 码雅可比缺 $\partial hs_i/\partial\phi(0)$ 的 $\chi_{bb}$ 项。残差本身含完整
$\chi_{bb}$ 边界项，故收敛解不变，只在 $\chi_{bb}\neq 0$ 时可能影响收敛速度/稳健性。
两码生产的 omega 网格均取 $\chi_{bb}=0$，此差异无效；chibb 扫描仅 A 码跑过。

注 3：B 码在掩码点组成的子序列上找相邻变号，若掩码中间有空洞（个别点 $cs_{gap}$
掉到阈下），变号可跨空洞在不相邻两个扫描点间插值；A 码只在连续有效区段内找变号。
掩码无空洞时两者相同；有空洞时 B 码可能多出跨洞插值点、A 码丢弃该处。

注 4：binodal 只用于画对照线与（A 码）剔除假点，不进入单点求解与交叉判定；
两码 binodal 算法差异不影响 pre-wetting 点本身的位置。

## 3. 参数对表

生产配置：A 码 config/sweep_plan.yaml；B 码 configs/default.yaml 之
prewetting_region_settings + configs/prewetting_three_chi_omega_grid.yaml
（经 scripts/prewetting/run_prewetting_omega_grid.py 驱动，任务 prewetting_region）。

共享值（逐项相同）：

| 参数 | 值 | A 码出处 | B 码出处 |
|---|---|---|---|
| L, N | 10.0, 1000 | sweep_plan.yaml:6-7 | default.yaml:2-3 |
| kappa_11/22/12 | 1.0, 1.0, 0.0 | sweep_plan.yaml:12-14 | default.yaml:4-6 |
| n1/n2/n3 | 1.0, 1.0, 1.0 | sweep_plan.yaml:9-11 | default.yaml:10-12 |
| 求解器 tol / max_iter | 1e-8 / 500 | solver.py:229 默认 | default.yaml:37-38 |
| cs 阈值 | 0.1 | sweep_plan.yaml:48 | default.yaml:66 |
| 扫描范围（每线） | [0.0001, 0.2] | sweep_plan.yaml:26-46 | default.yaml:68-75 |
| thin 初猜振幅 | 0.015 / 0.004 | thermodynamics.py:353-354 | default.yaml:30-31 |

差异值：

| 参数 | A 码 | B 码 | 影响见第 4 节 |
|---|---|---|---|
| 每线扫描点数 | 400（sweep_plan.yaml:30, 41） | 150（default.yaml:71, 75） | 条目 2 |
| thick 初猜振幅 | 0.09 / 0.02 | 0.3 / 0.3（default.yaml:78-79） | 条目 1 |
| 生产 chi 拓扑 | (0.0, 2.8, 2.6) 即 T-b | (0,2.8,0)/(2.3,2.8,2.6)/(2.8,2.8,2.8) 即 T-a/T-c/T-e（prewetting_three_chi_omega_grid.yaml:1-10） | 条目 6 |
| 生产 omega 点 | om1 ∈ {-0.52,-0.50,-0.45}, om2=-0.30 + chibb 单因子扫 | om1, om2 各 -0.3..-1.2 步长 -0.1 网格 | 条目 6 |
| 并行 workers | 30 | 12 | 无（不改结果） |

对 project_plan.md 的偏差（两码共同）：plan 第 2.3 节要求 $\phi_{1,\infty},\phi_{2,\infty}$
扫到 0.4、步长 0.01；两码实际都只扫到 0.2（plan 自身已标注后期扩全区间的 TODO）。
即 unvalidate_data 的相图只覆盖 $\phi_{i,\infty} \le 0.2$ 的窗口。

## 4. 差异影响评估

1. thick 初猜振幅（0.09/0.02 vs 0.3/0.3）：只作用于厚膜反扫首点（此后续种子）。
   更大的初始富集更容易落入厚膜盆地；若厚膜支在扫描线高浓端存在且吸引域不太窄，
   两种子收敛到同一支，交叉点相同；若厚支吸引域窄，B 码的大振幅更稳、A 码可能
   首点落回薄支导致整条反扫丢厚支。这是最可能造成两码结果不同的差异。
2. 每线点数（400 vs 150）：步长约 5.0e-4 vs 1.34e-3。影响交叉点线性插值精度、
   窄迟滞带（有效点 $\geq 2$ 才认）能否被采样到、以及 pre-wetting 线端点的截断位置。
   不改变物理分支，属分辨率差异。
3. 雅可比缺边界 $\chi_{bb}$ 项（B 码）：不改变收敛解；生产 $\chi_{bb}=0$ 时两码雅可比
   完全相同。
4. $\kappa_{12}$ 项缺失（B 码）：生产 $\kappa_{12}=0$，无影响。
5. 交叉定位与多交叉取舍：掩码带空洞或一线多交叉时两码取点可不同（B 码取首个、
   可跨洞插值；A 码分区段、全收集再过滤）。影响限于迟滞带边缘的个别点。
6. demixed 剔除 / 端点细化 / terminal 放宽 / auto expand（仅 A 码）：剔除机制去假点、
   细化机制提端点精度，均不移动已找到的交叉点本身；A 码产出的相图端点更精、
   假点更少，B 码相图可能含体相两相区内假点（其工作流用事后手工剔除补偿）。
7. 两码生产参数点不重叠：A 码产出 T-b（及 chibb 扫描），B 码产出 T-a/T-c/T-e 的
   omega 网格。因此本审计的互印证含义是：两份独立实现、同一模型、同一数值主干，
   各自产出的相图在方法学上可信；但没有任何一个 (chi, omega, chibb) 点被两码
   同时算过。

## 5. 局限

- 静态审计不构成数值等价证明：浮点求值顺序、clip 时机（A 码在 surface metrics 内
  对剖面先 sanitize，B 码对 $W$ 用原始剖面、仅 $f/\mu$ 内部 clip）等微差未逐行穷尽，
  只能断言公式与算法一致，不能断言输出逐位一致。
- 结论只覆盖两码的生产路径；备用分支（B 码 Picard/exponential_decay/wall_jump、
  A 码 pipeline.py:786 后死代码）未审计。
- 要把互印证升级为同点复算，需在同一组 (chi, omega, chibb) 与同一扫描配置下
  各跑一条线逐点对比 $\Omega_{thin}/\Omega_{thick}$ 与交叉点（本轮明确不做）。
