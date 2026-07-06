# 三元混合物 Pre-wetting 项目计划

## 0. 目的

图 2.4 的 6 种体相拓扑上加一堵有化学偏好的刚性壁面($f_{\text{surf}}$),求:$\phi_\infty=(\phi_{1,\infty},\phi_{2,\infty})$ 在 binodal 外侧单相区、接近共存时,墙催生的膜何时发生薄膜 → 厚膜一级转变。

## 1. 符号

- $\phi_1(z),\phi_2(z)$;$\phi_s=1-\phi_1-\phi_2$。
- 墙处 $\phi_1(0),\phi_2(0)$;远端 $\phi_\infty=(\phi_{1,\infty},\phi_{2,\infty})$。
- $n_1=n_2=1$;$\chi_{1s},\chi_{2s},\chi_{12}$($\chi_c=2$ 正、$-8$ 负);$\kappa_1,\kappa_2$。

$$
f_{\text{surf}}=\frac{k_BT}{\bar\nu}\Big[\omega_{b,1}\phi_1(0)+\omega_{b,2}\phi_2(0)+\chi_{bb,1}\phi_1(0)^2+\chi_{bb,2}\phi_2(0)^2+\chi_{bb,12}\,\phi_1(0)\phi_2(0)\Big]
$$

## 2. 三层控制变量

| 层级 | 对象 | 决定 |
|---|---|---|
| 0 | $(\chi_{1s},\chi_{2s},\chi_{12})$ | 体相拓扑、$W$ 地形 |
| 1 | $(\phi_{1,\infty},\phi_{2,\infty})$ | 欠饱和度、$W$ 零点 |
| 2 | $(\omega_{b,1},\omega_{b,2},\chi_{bb,\ast})$ | 薄膜/厚膜、是否 pre-wetting |

### 2.1 层级 0:6 种拓扑

| 设置 | $\chi_{1s}$ | $\chi_{2s}$ | $\chi_{12}$ | 拓扑 |
|---|---|---|---|---|
| T-a | 2.8 | 0 | 0 | 单个两相区从 $\phi_1$ 轴长出,1 临界点 |
| T-b | 2.8 | 2.6 | 0 | 两相区横跨两轴,无临界点 |
| T-c | 2.8 | 2.6 | 2.3 | 两个不相连两相区,各带临界点 |
| T-d | 2.6 | 2.6 | 2.6 | 三个不相连区 + 3 临界点 |
| T-e | 2.8 | 2.8 | 2.8 | 三相共存区 |
| T-f | 0 | 0 | −8.5 | 悬浮闭环区,2 临界点 |

### 2.2 层级 2

扫若干组 $(\omega_{b,1},\omega_{b,2},\chi_{bb,1},\chi_{bb,2},\chi_{bb,12})$,各求发生 pre-wetting 的 $(\phi_{1,\infty},\phi_{2,\infty})$ 点对。

### 2.3 层级 1

- 扫 $\phi_\infty$:$\phi_{1,\infty}\in[0,0.4]$、$\phi_{2,\infty}\in[0,0.4]$,$\phi_{1,\infty}+\phi_{2,\infty}\le 1$,步长 $0.01$。
- 每点判定薄膜侧 / 厚膜侧 / 转变($\gamma_{\text{thin}}=\gamma_{\text{thick}}$),连成转变线。
- ⚠ TODO:后期扩到全 $\phi_\infty$ 区间。

## 3. 求解

### 3.1 自由能

$$
\gamma[\phi_1(z),\phi_2(z)]=\int_0^\infty\Big[\,W+\frac{\kappa_1}{2}\Big(\frac{d\phi_1}{dz}\Big)^2+\frac{\kappa_2}{2}\Big(\frac{d\phi_2}{dz}\Big)^2\,\Big]dz+f_{\text{surf}}(\phi_1(0),\phi_2(0))
$$

$$
f_b=\frac{k_BT}{\nu}\Big[\frac{\phi_1}{n_1}\ln\phi_1+\frac{\phi_2}{n_2}\ln\phi_2+\phi_s\ln\phi_s+\chi_{1s}\phi_1\phi_s+\chi_{2s}\phi_2\phi_s+\chi_{12}\phi_1\phi_2\Big]
$$

$$
W(\phi_1,\phi_2)=f_b(\phi_1,\phi_2)-f_b(\phi_{1,\infty},\phi_{2,\infty})-\sum_{i=1,2}\frac{\partial f_b}{\partial\phi_i}\Big|_{\phi_\infty}\big(\phi_i-\phi_{i,\infty}\big)
$$

### 3.2 条件

固定全部参数与 $\phi_\infty$,最小化 $\gamma$ 得薄膜、厚膜两条剖面。

$$\gamma_{\text{thin}}=\gamma_{\text{thick}}$$

即 pre-wetting 转变。目标:扫参数与 $\phi_\infty$,连出 $(\phi_{1,\infty},\phi_{2,\infty})$ 平面上的转变线及临界终点。

### 3.3 平衡态(两条 DE)

$z>0$:

$$
\kappa_i\,\phi_i''(z)=\frac{\partial W}{\partial\phi_i},\qquad i=1,2
$$

$z=0$:

$$
\kappa_i\,\phi_i'(0)=\frac{\partial f_{\text{surf}}}{\partial\phi_i(0)},\qquad i=1,2
$$

$z\to\infty$:$\phi_i\to\phi_{i,\infty}$,$\phi_i'\to 0$。首次积分:

$$
\sum_{i=1,2}\frac{\kappa_i}{2}\phi_i'^2=W
$$

## 4. 局限

1. 线性路径 $\vec\phi(t)=\vec\phi_\infty+t(\vec\phi(0)-\vec\phi_\infty),\ t\in[0,1]$ 仅参数对称时严格;不对称(T-a/b/c)需真解 ODE 校验。
2. $\phi_\infty$ 须在单相区且贴近 binodal。
