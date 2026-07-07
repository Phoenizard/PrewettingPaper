# PROGRESS

本文件记录项目重要进度，最新在上。维护约定：条目简洁，不用粗体；「当前状态」就地更新，「进度日志」只追加。

## 当前状态

- 阶段：推倒重来。此前 Claude Code 写的全部代码已删除，从头重新审视问题。
- 决定（2026-07-08）：之前所有代码已不可再维护，在其上的后续改进也无意义，故删除代码库
  所有代码脚本（src/ 全部 .py + scripts/ 全部 .py/.sh，含 diagnostics）。commit f27dcb5。
  reference/（同组成员的可运行参考实现）保留、未动，将作为后续从头实现的对照基准。
- 代码：src/、scripts/ 已清空，无任何可运行代码。
- GitHub：https://github.com/Phoenizard/PrewettingPaper （public）
- 服务器：AutoDL；当前已由用户关闭。
- 为什么推倒：核心求解一直建在 scipy solve_bvp 上，自适应网格撑不住低 phi2 的宽厚膜（爆
  max_nodes、丢厚支），围绕它叠了多层补丁仍解决不了 om2=-0.40 线延到 phi2→0 的问题；
  后期虽认识到应改固定网格并采纳参考方法，但实现方式是"凭理解重写 + 试错"、还自造了参考里
  没有的判据（accept 远场拒绝）与自己搏斗，越走越乱。结论：不值得继续维护，清零重来。
- 教训（已记 memory mirror-reference-code-dont-reinvent）：有可运行参考实现时，应逐项对齐其
  确切数值/逻辑（L、N、初值、扫描范围、判据…）再验证一致，而不是重写后试错。

## 下一步

从头重新审视问题。代码已清零，先想清楚再写。

- 先定路线：以 reference/（可运行的参考实现）为基准，逐项对齐其确切数值与逻辑（固定网格
  L/N、初值构造、双向扫描范围、crossing/迟滞判据、Newton 细节），而不是重写后试错。
  独立性 = 我们自己的实现，但数值配置照参考对齐、先复现其结果，再谈适配。
- 由用户主导决定从哪里起步、用什么骨架，Claude 不擅自开写。
- 之前的诊断结论仍可参考（见进度日志）：低 phi2 有薄/中/厚三个表面态，真 pre-wetting 是
  薄↔厚等深转变；参考图正确，线应延到 phi2→0。这些是"要复现的目标"，不是"要保的代码"。

（以下为历史遗留 TODO，代码已删，仅存档：全量 770 case 扫描、结果对比 UI、
environment.yml 固化依赖、phi_inf 扫描扩到全区间。重来后按新代码再规划。）

## 结果对比 UI 可行性

目标：快速浏览每个 stage（chi 拓扑目录，3 个）下每种配置（om × chibb）本地 vs result 的差异。

已具备（低成本）：
- 配对已解决，直接复用 cases.py：iter_cases 枚举全部 leaf；result_overlay(rel) 给 result 图，
  verify_dir(rel)/overlay.png 给本地图。UI 不需新写配对逻辑。
- 数量：3 chi ×（11×11 omega 网格 + chibb 扫描）≈ 最多 770 个 leaf；result 主图 overlay_omega1_omega2.png 存在于 363 个。

关键约束：
- result 侧只有 PNG（无数值 pw 数据），本地侧有 overlay.png + pw_line.csv。
- 两侧图来自不同渲染器/尺寸（result ~197 KB，本地 ~32 KB，坐标轴/排版不同）→ 像素级 diff 无意义，
  只能做「并排视觉对比」；现成图像 diff 工具（ImageMagick compare 等）不适用。
- 依赖顺序：本地图目前只有 1 个 case，UI 要等服务器批量跑完 out/verify/ 才有用。

方案与工作量：
- A 静态 HTML 画廊（推荐）：一个 Python 脚本遍历 iter_cases，生成 index.html——按 stage 分组、
  om/chibb 可筛选，每格 result 图与本地图并排，<img> 懒加载指向已有 PNG。离线、无服务器，
  约半天；缺图的 leaf 跳过。
- B 本地 web app（Streamlit/Flask）：下拉选 chi/om/chibb，并排 + 可选叠加滑块，交互更强但要起服务、
  加依赖，约 1-2 天。
- 建议：先做 A；若需交互式叠加再上 B。

已确认：「stage」= chi 拓扑目录（3 个），「配置」=(om, chibb) 组合；画廊按 stage 分组、om/chibb 可筛选。

## 进度日志

### 2026-07-08（推倒重来：删除全部代码）

- 决定删除代码库所有代码脚本（src/ 全部 .py + scripts/ 全部 .py/.sh，含 diagnostics），
  从头重新审视问题。commit f27dcb5（19 文件、2743 行删除）。reference/ 保留未动。
- 原因：核心一直建在 scipy solve_bvp，撑不住低 phi2 宽厚膜，多轮补丁未解决 om2=-0.40 线延到
  phi2→0；后期改固定网格 + 采纳参考方法时，实现走成"凭理解重写 + 试错 + 自造判据"，越改越乱。
  判定：已不可维护、后续改进无意义，故清零。
- 本 session 有价值的诊断结论（供重来时参考，非要保的代码）：
  - 低 phi2 出现薄/中/厚三个表面态（墙面 phi1≈0.34/0.68/0.92）；真 pre-wetting 是薄↔厚等深
    转变，中间态是势垒；参考图正确，线应延到 phi2→0。
  - solve_bvp 丢厚支的机制：厚膜变宽爆 max_nodes；真厚态存在的 phi1 区间太窄，等 gamma 交叉前
    厚态已消失。
  - 参考方法（固定网格 Newton + 双向扫描 + crossing）机制正确：厚支在高 phi1（phi1 富集角，
    phi1_scan 到 0.2）出生、warm 连续进来；phi1_scan=linspace(0.0001,0.2,400)、逐点更新
    phi1_inf、首点冷启后 warm（thermodynamics.py:667-680）。
- 教训记入 memory：mirror-reference-code-dont-reinvent（有参考实现就逐项对齐其数值，别重写试错）。

### 2026-07-07（夜：结果管理 + 画廊,发现 PW 线两端漏点）

- 结果管理落地:scripts/pull_results.sh(rsync 服务器 out/verify → 本地 $RESULTS_DB/verify,并集 merge;
  RESULTS_DB 可配置,默认 ./database,以后移动硬盘改路径即可)、build_gallery.py + serve_gallery.sh
  (离线静态画廊,本地 http 服务;下拉选 chi/om/chibb → ours 与 result 并排 + pw_line.csv 路径)。
  不用网盘(坚果云有存储/流量限制)。file:// 直开会因浏览器安全策略不显示 result 图,故走本地 http。
- 用画廊肉眼比对,发现遗留问题(见「下一步」):我们的 PW 线两端漏点,只覆盖 phi2=0.03–0.06,
  而 result 覆盖 phi2≈0–0.09。已登记,暂不处理。

### 2026-07-07（傍晚：验证优化不改结果 + 工作流成型）

- 纠正验证方法：删掉 bench/regression 这类"重算式测试"脚本，改为"产出 → 事后对比"：优化代码跑 10 基线
  → out/verify_opt（VERIFY_OUT 覆盖输出根，基线 out/verify 未动），再单独用 compare.py 逐点比。
- 结果：compare 10/10 PASS，max_dev 1e-17~5e-13，远低于 1e-3 → 优化确认不改结果，可上全量 770。
- 产出批次 opt10：16 核并行 10 case，墙钟 1022s（约 17 分钟）；10 个 case pre-wetting 点数全与基线一致。
- 踩坑与规范（详见 doc/note/Workflow.md，不入 git）：dry-run 要"短"（PW_MAX_PHI2 截断扫描）；
  一实验一 screen、按名字杀（不 pkill -f）；github 操作前先 source network_turbo；conda run 缓冲输出、
  改用直接解释器 /root/miniconda3/envs/numenv/bin/python；产出与对比分离；实验前给矩阵并等确认；
  每条 bash 先说明；git 指令单独执行；正文写全称 pre-wetting。均已入 memory。
- 新增 check-progress skill（只读、渲染实验矩阵）。

### 2026-07-07（下午：计算加速 + 工作流规范）

- 计算/绘图解耦：verify.py 只算写 CSV、不碰 matplotlib（根除并发画图竞态）；新增 plot.py（CSV->PNG）、
  plotting.render_phase_map（(实验变量+结果)->图 的通用包装）。10 个基线 CSV 已补出 overlay.png。
- 计算加速（藏开关后）：equilibrium.py 加解析 fun_jac/bc_jac（来自 hessian_fb，与 res 无关）、
  mu(res)/f_b(res) 缓存、Profile 带 sol_x/sol_y、warm-start（守卫 + 冷启动回退）；thermo.py 加
  reservoir_potentials 与可选 res_mu/res_fb；verify.py prewetting_line 串 warm + 进度日志。
- bench（T-a 模板）：Layer1=852s、both=834s、baseline>1046s(被手动停)。→ 雅可比约 1.3-1.4×；
  warm-start 空转（solve_bvp 次数 both 2290 vs l1 2322，几乎无减），守卫因常见 3 分支不通过。
- 新增 bench.py（测速 + solve_bvp 计数，按开关分档）、regression.py（内存重算基线、断言 ≤1e-3、不写 out/）、
  run_regression.sh（xargs 并行、日志镜像结果树 out/logs/regression/<...>/run.log）。
- 踩坑：nproc 读 OMP_NUM_THREADS=1 → run_regression 误以 cores=1 串行启动；已改 nproc --all（未带修复重跑）。
- 工作流规则大量确立（见「当前状态」末条，均入 memory）：dry-run 先行、实验前出矩阵并等确认、不擅自停实验、
  screen 一实验一会话只给 `screen -r`、git 单独执行、日志分类镜像结果树；新增 check-progress skill。
- 未完：regression 门未跑通确认；服务器上一个串行 reg 运行待用户决定停/留。

### 2026-07-07

- 服务器并行部署。git init + 首次提交；用 gh 建 public 仓库 Phoenizard/PrewettingPaper 并推送。
  result/（242MB PNG）不上传，改提交 result_cases.txt（770 行三元组清单）当服务器 case 来源。
- 新增 scripts/run_parallel.sh（xargs -P 按 case 并行）；verify.py 加 --no-summary / --rebuild-summary。
- 服务器：network_turbo 加速 clone，conda 建 numenv。跑 16-case 首测，调出两个 bug：
  1. BLAS 超额订阅：每 worker 默认开 16 线程，16×16 抢 16 核，load~48、各 54% CPU。
     修复：run_parallel.sh export *_NUM_THREADS=1（单线程），修后各 100% CPU、16 核干净铺满。
  2. matplotlib 并发首次导入竞态：16 进程同时首次 import 崩 "No module named backend_agg"
     （计算已成、pw_line.csv 已写，仅画图崩）。修复：并行前单进程预热缓存；verify.py 的
     --skip-existing 改以 pw_line.csv 为准、缺图只重画，抢救已算的 10 个 case。
- 经验记录：network_turbo 只对 github/hf 加速、会拖慢 conda/pip；tmux 保活；避免 pkill -f 自匹配。
- 待办：带修复重跑 16-case 确认加速比（今日到此为止，未重跑）。

### 2026-07-06

- 结果对比 UI 可行性分析（详见上节）。结论：配对可复用 cases.py，两侧图不可像素 diff、只能并排，
  推荐先做静态 HTML 画廊（约半天），依赖服务器批量结果先到位。已记入下一步 TODO。
- 全量成本估算。空间：33 KB/case（overlay.png + pw_line.csv）× 770 ≈ 25 MB，忽略不计。
  时间瓶颈在 find_states（多起点 solve_bvp）≈ 3.9 s/次调用；每 case 扫描 ~300-380 网格点
  → 约 20-25 min/case；770 case 串行单核约 320 h ≈ 10-13 天。binodal 计算 0.18 s/次可忽略。
- 结论：串行不现实，需部署服务器 + 并行（见「下一步」TODO）。并行放启动脚本、不进代码。
- 建立 PROGRESS.md。
- 现状：src/、scripts/ 代码完整无 stub；binodal 6 拓扑已出图；result/ 含 3 个 chi 目录约 770 个 case。
- 验证仅跑通 1 个 case（out/verify/SUMMARY.csv 仅 1 行），待批量执行。
