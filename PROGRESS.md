# PROGRESS

本文件记录项目重要进度，最新在上。维护约定：条目简洁，不用粗体；「当前状态」就地更新，「进度日志」只追加。

## 当前状态

- 阶段：验证 result/ 结果（同组成员的求解，做独立复核）；服务器并行部署链路已打通，调试收尾中
- 代码：src/ 5 模块 + scripts/ 3 驱动（含 run_parallel.sh），均可运行
- GitHub：https://github.com/Phoenizard/PrewettingPaper （public）
- 服务器：AutoDL，16 核 / 1TB 内存 / Ubuntu 20.04；代码在 /root/autodl-tmp/PrewettingPaper，numenv 已建
- binodal 诊断：6 种拓扑 T-a..T-f 已出图（out/_binodal_check/）
- 验证进度：本地 1 / ~770；服务器有 10 个 pw_line.csv 基线（已解耦出图）
- 计算/绘图已解耦：verify.py 只算写 pw_line.csv（不 import matplotlib），plot.py 单独 CSV->PNG
- 计算加速已落地并验证不改结果（均藏在开关后，USE_ANALYTIC_JAC / USE_WARM_START）：
  - Layer 1 解析雅可比 + mu(res)/f_b(res) 缓存 —— 约 1.3-1.4×（雅可比为主）
  - Layer 2 warm-start —— 目前空转（守卫要求恰好 2 分支，转变线附近常 3 分支，总退回冷启动）；改进搁置
- 验证已通过（方法：产出 → 事后对比，不写重算式测试）：优化代码跑 10 基线 → out/verify_opt，
  compare.py 逐点比基线，10/10 PASS，max_dev ≤ 5e-13（远低于 1e-3）。结论：优化不改结果，可用于全量 770。
  - 工具：verify.py 加 VERIFY_OUT/VERIFY_PROGRESS/PW_MAX_PHI2；plot.py 认 VERIFY_OUT；run_parallel.sh 加
    MANIFEST 覆盖 + nproc --all + 直接解释器默认；新增 compare.py；删 bench/regression（重算式测试逻辑错误）
- 工作流规则（新，均记入 memory）：ssh 服务器工作流；实验只在服务器跑、禁本地 smoke；每个实验一个 screen + 只给 `screen -r`；实验前给矩阵并等确认 + 先 dry-run；不擅自停/重启在跑实验；git 指令单独执行；日志按类型进 out/logs/，正式逐-case 实验日志镜像结果树；新增 check-progress skill

## 下一步

- [遗留问题, 暂不处理] pre-wetting 转变线两端漏点:画廊对比 case chi13=2.8 / om1=om2=-0.30 发现,
  result 参考的 PW 线覆盖 phi2≈0–0.09(贴 binodal 左翼一整条),我们只得 phi2=0.03–0.06 四个点,两端都缺。
  不是密度、也不是扫描范围问题(phi2 从 0.01 起步长 0.01,本就扫过 0.01/0.02/0.07/0.08/0.09,却没记到点)。
  根因方向:两端每个 phi2 上 find_states/变号夹逼没抓到转变(可能没同时得到薄/厚两支,或 phi1 内层网格
  没框住转变点)。属复核代码在 PW 线两端的检测缺陷,待后续查(先登记不改)。
- [done] 10 基线验证通过（compare 10/10 ≤5e-13）→ 优化代码可用于全量。
- 用优化代码跑全量 770（分批：产出 → 对比）：run_parallel.sh 走 result_cases.txt，产出到 out/verify；
  逐 chi 目录扩规模；每批先 dry-run、给矩阵等确认再上（见 doc/note/Workflow.md）。
- [搁置] warm-start 守卫改进（放宽到"≥2 分支、按 thin/thick 就近匹配"），若要再提速再做，之后仍走产出→对比验证。
- TODO: 部署服务器无人值守跑全部 770 case（串行单核约 10-13 天，见下方成本估算）。
  - 并行由启动脚本 scripts/run_parallel.sh 处理（xargs -P，不写进 Python 代码）；case 独立。
  - [done] verify.py 已加 --skip-existing / --no-summary / --rebuild-summary（并行无竞态 + 断点续跑）。
  - [done] 单 case 现实耗时：服务器单线程满速约 17-25 min/case（与本地估算同量级）。
  - 算法加速（可与并行叠乘）：find_states 沿 phi1 扫描改暖启动/延拓，减少 solve_bvp 重解。
- 比对 out/verify/ 与 result/ 的 pre-wetting 转变线，记录一致/差异
- TODO: 结果对比 UI（快速看每个 stage / 每种配置下 本地 vs result 的差异）。可行性见下。
- phi_inf 扫描区间后期扩到全区间（现为 [0,0.4]^2，步长 0.01）
- 补 environment.yml，固化 numenv 依赖

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
