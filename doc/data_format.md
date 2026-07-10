# 数据上传格式规范

组员实验数据统一按本规范上传到线上存储容器，目录结构与本仓库 `data/` 一一对应。

## 目录结构

四层：stage（chi 拓扑）→ omega → chibb → case 文件。路径模板：

```
chi12_<v>__chi13_<v>__chi23_<v>/om1_<v>__om2_<v>/chibb11_<v>__chibb22_<v>__chibb12_<v>/<文件>
```

命名编码：`p` = 小数点，`m` = 负号，同层参数间用双下划线 `__` 分隔，小数位数按需
（-0.30 → `m0p30`，-0.345 → `m0p345`，2.8 → `2p8`）。

真实例子：

```
chi12_2p3__chi13_2p8__chi23_2p6/om1_m0p30__om2_m0p44/chibb11_0p00__chibb22_0p00__chibb12_0p00/
```

## case 内必须的文件（数值数据）

每个 case 目录必须包含以下两个 CSV，能据此直接重绘 2D 相图：

| 文件 | 内容 | 列 |
|---|---|---|
| `pw_line.csv` | pre-wetting 线 | `source, phi1_inf, phi2_inf` |
| `binodal.csv` | binodal 线 | `phi1, phi2` |

`source` 取 `fix_phi1` 或 `fix_phi2`，标记该点来自哪个扫描方向。
无 pre-wetting 点的 case 也要交 `pw_line.csv`（只含表头）。

## case 内可选的文件

- `overlay.png` — 由上述数值生成的 2D 图。
- 注意：交图的组员之间必须统一绘图格式（坐标轴、范围、图例、配色）；
  未统一前只交数值不交图。

## 注意

- 不要上传 `.DS_Store` 等系统杂物文件。
