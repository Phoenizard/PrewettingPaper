# PROGRESS

本文件记录项目重要进度，最新在上。维护约定：条目简洁，不用粗体；「当前状态」就地更新，「进度日志」只追加。

## 当前状态

- 阶段：验证 result/ 结果（同组成员的求解，做独立复核）
- 代码：src/ 5 模块 + scripts/ 2 驱动，均可运行
- binodal 诊断：6 种拓扑 T-a..T-f 已出图（out/_binodal_check/）
- 验证进度：1 / ~770 case

## 下一步

- TODO: 部署服务器无人值守跑全部 770 case（串行单核约 10-13 天，见下方成本估算）。
  - 并行由启动脚本（.sh，如 GNU parallel / 作业分发）处理，不写进 Python 代码；case 之间相互独立。
  - [done] verify.py 已加 --skip-existing：overlay.png + pw_line.csv 都在则跳过重算、从 csv 补 SUMMARY。
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
