# omega topic 度量说明

数据源：pw-space/data 全库 1112 case。总表 `measures.csv`。本页对应 omega topic：T-a（chi12=0, chi13=2.8, chi23=0），chibb=0，11×11 的 (omega1, omega2) 网格，各 ∈[-0.50,-0.30] 步长 0.02，共 121 case。

## 问题

墙对两溶质的亲附力 omega1、omega2 变化时，pre-wetting 的 extent 如何变、何处熄灭。extent 分两维：pre-wetting line 的长度、line 到 binodal 的距离。

## 指标

每 case 合并全部 pre-wetting 点 $\{p_k=(\phi_{1,\infty},\phi_{2,\infty})_k\}$（不分branch），binodal 点集 $B$。均在原始平面，不归一化。长度 $L$（MST 弧长，剪空隙）：

$$
E = \mathrm{MST}\big(\{p_k\}\big),\qquad
L = \sum_{(i,j)\in E,\ \lVert p_i-p_j\rVert \le \tau} \lVert p_i - p_j \rVert,
\qquad \tau = 0.015
$$

MST（minimum spanning tree，最小生成树）：把所有点两两可连，选出一组边，把全部点连成一棵无环树，且总边长最小。直观上它沿着点的走向逐点相连、把点串成线，故其边长之和即这条线的弧长（沿弯折走，不受点疏密影响）。这里再把长度 $>\tau$ 的边剪掉——它们是跨越空隙连接分离簇的桥边，剪后剩余边长之和即"段内弧长之和、空隙不计"。

例（`mst_example.png`）：5 个点排成一条折线 + 1 个远点。实线（蓝）是 4 条段内边0.030+0.030+0.032+0.032，其和 $L=0.123$；虚线（红）0.130 是连向远点的桥边，长$>\tau=0.05$ 被剪、不计入 $L$；剪后两个连通分量，$n_{seg}=2$。

![MST 示例](../mst_example.png)

到 binodal 距离 $\bar d$（逐点最近距离取平均）：
$$
\bar d = \frac{1}{N}\sum_{k=1}^{N} \min_{b\in B}\lVert p_k - b\rVert
$$

诊断量：$n_{seg}$（剪 $\tau$ 后 MST 连通分量数）

$$
\mathrm{gap} = \sum_{\lVert p_i-p_j\rVert>\tau}\lVert p_i-p_j\rVert
$$

$$
d_{\min}=\min_k\min_{b}\lVert p_k-b\rVert
$$

$$
d_{\max}=\max_k\min_b\lVert p_k-b\rVert
$$

## 热力图

文件 `omega_length_dist.png`，两联（同一 (omega1, omega2) 网格）：左为长度 $L$，
右为到 binodal 距离 $\bar d$。无 pre-wetting 的 case 画为白格。

![omega 热力图](omega_length_dist.png)

## 结论

体相拓扑 $\chi_{12}=0,\ \chi_{13}=2.8,\ \chi_{23}=0$。

1. 两幅热力图的颜色都沿 $\omega_1$ 方向（横向）变化、沿 $\omega_2$ 方向（纵向）几乎不变：$L$ 从 $\omega_1$ 每列的 $0.010$ 到 $0.081$、$\bar d$ 从 $0.021$ 到 $0.004$，而固定 $\omega_1$ 改变 $\omega_2$ 时 $L$ 仅由 $0.027$ 缓变到 $0.053$、$\bar d$ 恒在 $0.012$ 附近。由此可见 pre-wetting 由溶质 1 的墙面亲附力 $\omega_1$ 单独控制，溶质 2 的 $\omega_2$ 几乎不起作用。

2. 沿减弱吸引方向（$\omega_1$ 由 $-0.50$ 增至 $-0.30$），左图 $L$ 增大而右图 $\bar d$ 减小（$L$：$0.010\to0.081$，$\bar d$：$0.021\to0.004$）。二者反号说明：墙对溶质 1 的偏好越弱，可触发薄膜到厚膜一级转变的远场组成范围越宽（$L$ 大），但转变越贴近体相共存（$\bar d$ 小）——弱偏好下墙只在体相已接近共存时才诱导厚膜，却能在紧贴 binodal 的一条较宽的带上诱导。反之，偏好越强，可触发范围越窄，但转变被推入离共存更远的深单相区，墙在体相远未共存时便已诱导厚膜。

3. 两幅图各有两个白格（无 pre-wetting），均位于 $\omega_1=-0.50$ 的最强吸引角，且逼近该角时 $L$ 平滑趋于 0（$\omega_1=-0.50$ 一列由 $0.015$ 降到熄灭）。由此可见 pre-wetting 在最强吸引处熄灭：过强的墙面吸引把薄膜与厚膜之间的一级跳变推过其临界点、转为连续吸附，pre-wetting line 随之连续收缩至消失。
