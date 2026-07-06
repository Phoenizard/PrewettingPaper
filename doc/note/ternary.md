# 1. 系统变量定义 (System Variables)

沿用原文的体积分数符号 $\phi$，并使用下标区分组分：

- **体积分数 (Volume Fractions)**:
  - 溶质 1 (Solute 1): $\phi_1(z)$
  - 溶质 2 (Solute 2): $\phi_2(z)$
  - 溶剂 (Solvent): $\phi_s(z) = 1 - \phi_1(z) - \phi_2(z)$
- **储库状态 (Reservoir)** (对应原文 $\phi_\infty$):
  - $\phi_{1,\infty}, \phi_{2,\infty}$
- **化学势 (Chemical Potentials)** (对应原文 $\mu_\infty$):
  - $\mu_{1,\infty}, \mu_{2,\infty}$

---

# 2. 体相自由能密度 $f_b$ (Bulk Free Energy Density)

我们将原文 Eq. (15a) 扩展为三元混合物形式。

**原文符号**:

- $\nu$: 溶剂分子体积
- $n_b$: 溶质分子相对于溶剂的体积比 ($n_b = \nu_b/\nu$)
- $\chi_b$: 相互作用参数

**新模型 ($f_b$)**:

$$f_b(\phi_1,\phi_2) = \frac{k_B T}{\nu}\left[\frac{\phi_1}{n_1}\ln\phi_1 + \frac{\phi_2}{n_2}\ln\phi_2 + (1-\phi_1-\phi_2)\ln(1-\phi_1-\phi_2) + f_{\text{int}}\right]$$

其中相互作用项 $f_{\text{int}}$ (Flory-Huggins Interaction):

$$f_{\text{int}} = \chi_{1s}\phi_1(1-\phi_1-\phi_2) + \chi_{2s}\phi_2(1-\phi_1-\phi_2) + \chi_{12}\phi_1\phi_2$$

- **符号对应**:
  - $n_1, n_2$: 对应原文的 $n_b$。
  - $\chi_{1s}, \chi_{2s}$: 对应原文的 $\chi_b$。
  - $\chi_{12}$: 溶质间的新相互作用参数。

---

# 3. 广义势密度 $W$ (Grand Potential Density)

这是描述体相偏离平衡态的能量代价。我们依据原文 Eq. (10d) 和 Eq. (A3) 进行扩展。

**原文定义**: $W(\phi) = f_b(\phi) - \frac{\mu_\infty}{\nu_b}\phi + \Pi_\infty$。

**新模型 ($W$)**:

$$W(\phi_1,\phi_2) = f_b(\phi_1,\phi_2) - f_b(\phi_{1,\infty},\phi_{2,\infty}) - \sum_{i=1,2}\frac{\mu_{i,\infty}}{\nu_i}(\phi_i - \phi_{i,\infty})$$

- **注意**: 这里的 $\mu_{i,\infty}$ 定义为 $\nu_i \frac{\partial f_b}{\partial \phi_i}$，严格遵守原文 Eq. (102) 的量纲定义（即 $\mu/\nu$ 为能量密度）。
- **物理意义**: 当 $\phi_i = \phi_{i,\infty}$ 时，$W = 0$ 且为极小值。

---

# 4. 表面相互作用能 $f_{\text{surf}}$ (Surface Interaction Energy)

这是最关键的一步。我们需要替换掉原文中的 $J(\phi|_0, \phi_m)$ 和 $f_m(\phi_m)$，因为我们去掉了膜结合 ($\phi_m$)。

我们将 $J$ 中的 $\phi_m$ 相关项移除，保留并扩展原文用于描述"体相分子在表面性质"的参数 $\omega_b$ 和 $\chi_{bb}$。

**原文符号**:

- $\omega_b$: 体相分子在表面的内部自由能 (吸引/排斥)。
- $\chi_{bb}$: 表面增强相互作用参数。
- $\bar\nu$: 表面层的特征体积/面积尺度。

**新模型 ($f_{\text{surf}}$)** (仅依赖表面处的体相浓度 $\phi_i|_0$):

$$f_{\text{surf}}(\phi_1|_0,\phi_2|_0) = \frac{k_BT}{\bar\nu}\Big[\underbrace{\omega_{b,1}\phi_1|_0+\omega_{b,2}\phi_2|_0}_{\text{表面亲和能 (Surface Affinity)}} + \underbrace{\chi_{bb,1}(\phi_1|_0)^2+\chi_{bb,2}(\phi_2|_0)^2+\chi_{bb,12}\phi_1|_0\phi_2|_0}_{\text{表面增强效应 (Surface Enhancement)}}\Big]$$

- **符号对应**:
  - $\omega_{b,1},\omega_{b,2}$: 对应原文的 $\omega_b$。
  - $\chi_{bb,1},\chi_{bb,2}$: 对应原文的 $\chi_{bb}$。

---

# 5. 吉布斯表面能泛函 $\gamma$ (Gibbs Surface Free Energy)

这是你导师要求的最终输出目标。我们基于原文的吉布斯表面能定义 Eq. (6) 和泛函形式 Eq. (3)，结合上面的新定义。

**原文逻辑**: $\gamma_s = \int dz[\text{Gradient} + W] + \text{Surface Terms}$。

**新模型 ($\gamma$)**:

$$\gamma[\phi_1,\phi_2] = \int_0^\infty dz\left[W(\phi_1(z),\phi_2(z)) + \frac{1}{2}\kappa_1(\partial_z\phi_1)^2 + \frac{1}{2}\kappa_2(\partial_z\phi_2)^2\right] + f_{\text{surf}}(\phi_1|_0,\phi_2|_0)$$

- **符号对应**:
  - $\kappa_1,\kappa_2$: 对应原文的梯度惩罚系数 $\kappa$。
  - 这里为了简化，忽略了交叉梯度项 $(\nabla\phi_1\cdot\nabla\phi_2)$，这在经典预润湿理论中是常见的处理。

---

# 6. 平衡条件 (Equilibrium Conditions)

要找到稳定的预润湿状态，你需要求解使 $\delta\gamma = 0$ 的分布。对应原文 **Eq. (10a)** 和 Eq. (10b)。

1. **体相方程** ($z>0$):

$$\kappa_i \frac{d^2\phi_i}{dz^2} = \frac{\partial W}{\partial \phi_i} \quad (i=1,2)$$

2. **边界条件** ($z=0$):

$$-\kappa_i \frac{d\phi_i}{dz}\Big|_{z=0} + \frac{\partial f_{\text{surf}}}{\partial \phi_i|_0} = 0 \quad (i=1,2)$$

---

# 三元混合物预润湿理论：物理模型、数学推导与数值实现

本文档旨在建立一个无膜结合（Rigid Wall）、三元混合物（溶剂 + 溶质1 + 溶质2）的预润湿（Pre-wetting）理论模型。文档涵盖了从热力学基础概念到数值求解算法的完整逻辑链条。

## 1. 系统定义 (System Definition)

### 1.1 物理模型

我们研究一个半无限大的流体系统 ($z \in [0,\infty)$)，该流体与位于 $z=0$ 处的**刚性壁面**接触。

与原文 [cite: 92] 中动态结合的膜不同，这里的壁面是静态的，仅通过短程势场对流体施加吸附作用。

### 1.2 系统变量

为了描述三元混合物的组分分布，我们定义以下场变量：

- **体积分数 (Volume Fractions)**:
  - $\phi_1(z)$: 溶质 1 的局部体积分数。
  - $\phi_2(z)$: 溶质 2 的局部体积分数。
  - $\phi_s(z) = 1-\phi_1(z)-\phi_2(z)$: 溶剂的体积分数（不可压缩条件）。
- **储库状态 (Reservoir State)**: $\phi_{1,\infty},\phi_{2,\infty}$。这是无穷远处 ($z\to\infty$) 的平衡浓度，充当系统的"粒子源"和"能量基准"。

---

## 2. 体相热力学：自由能与化学势

### 2.1 概念：体相自由能密度 ($f_b$)

**定义**: $f_b$ 描述了在没有梯度和外场的情况下，均匀混合流体内部的自由能密度。

**动机**: 基于 Flory-Huggins 平均场理论，自由能由两部分竞争组成：

1. **混合熵 ($f_{\text{entropy}}$)**：总是倾向于让物质均匀混合，由 $-T\Delta S$ 贡献。
2. **相互作用焓 ($f_{\text{int}}$)**：描述分子间的吸引或排斥（$\chi$ 参数），驱动相分离。

**数学形式**:

$$f_b(\phi_1,\phi_2) = \frac{k_BT}{\nu}\Big[\underbrace{\frac{\phi_1}{n_1}\ln\phi_1+\frac{\phi_2}{n_2}\ln\phi_2+\phi_s\ln\phi_s}_{\text{混合熵}} + \underbrace{\chi_{1s}\phi_1\phi_s+\chi_{2s}\phi_2\phi_s+\chi_{12}\phi_1\phi_2}_{\text{相互作用焓}}\Big]$$

### 2.2 推导：化学势 ($\mu$)

**动机**: 化学势是粒子数变化的能量代价。为了计算系统相对于储库的势能差（广义势），我们需要知道每个组分的化学势。

**定义**: $\mu_i = \nu_i\frac{\partial f_b}{\partial \phi_i}$（相对于溶剂的交换化学势）。

**推导细节**:

我们需要计算 $\frac{\partial f_b}{\partial \phi_1}$，利用 $\phi_s = 1-\phi_1-\phi_2$ 及 $\frac{\partial \phi_s}{\partial \phi_1}=-1$。

1. **熵项求导**:

   利用 $\frac{d(x\ln x)}{dx} = \ln x+1$:

   $$\frac{\partial f_{\text{entropy}}}{\partial \phi_1} = \frac{1}{n_1}(\ln\phi_1+1)+(-1)(\ln\phi_s+1)$$

   $$= \frac{1}{n_1}\ln\phi_1 - \ln(1-\phi_1-\phi_2)+\left(\frac{1}{n_1}-1\right)$$

2. **相互作用项求导**:

   $$f_{\text{int}} = \chi_{1s}\phi_1\phi_s+\chi_{2s}\phi_2\phi_s+\chi_{12}\phi_1\phi_2$$

   $$\frac{\partial f_{\text{int}}}{\partial \phi_1} = \chi_{1s}[1\cdot\phi_s+\phi_1(-1)]+\chi_{2s}\phi_2(-1)+\chi_{12}\phi_2$$

   $$= \chi_{1s}(1-2\phi_1-\phi_2)-\chi_{2s}\phi_2+\chi_{12}\phi_2$$

3. **结果**:

   $$\frac{\mu_1}{\nu_1} = \frac{k_BT}{\nu}\left[\frac{1}{n_1}\ln\phi_1-\ln\phi_s+\left(\frac{1}{n_1}-1\right)+\chi_{1s}(1-2\phi_1-\phi_2)+(\chi_{12}-\chi_{2s})\phi_2\right]$$

---

## 3. 广义势密度 ($W$)：驱动力

### 3.1 概念：为什么引入 $W$？

**动机**: 我们的系统是开放的，与一个无穷大的储库相连。在平衡态时，任何位置的流体都试图达到与储库相同的化学势。

**定义**: $W(\phi)$ 是**巨正则势密度 (Grand Potential Density)** 的差值。它衡量了在局域浓度 $\phi(z)$ 偏离储库浓度 $\phi_\infty$ 时，单位体积需要付出的额外能量代价。

**性质**: 当 $\phi=\phi_\infty$ 时，$W=0$ 且为极小值。这保证了远离壁面处流体回归储库状态。

### 3.2 数学形式

依据原文 Eq. (10d) [cite: 183] 扩展：

$$W(\phi_1,\phi_2) = f_b(\phi_1,\phi_2)-f_b(\phi_{1,\infty},\phi_{2,\infty})-\sum_{i=1,2}\frac{\mu_{i,\infty}}{\nu_i}(\phi_i-\phi_{i,\infty})$$

---

## 4. 表面相互作用能 ($f_{\text{surf}}$)

### 4.1 概念：边界的作用

动机：由于去掉了原文中的膜结合模型 [cite: 98]，我们需要用一个有效势场来描述壁面。

**定义**: $f_{\text{surf}}$ 是仅作用于 $z=0$ 处流体层 ($\phi|_0$) 的自由能。它包含两个物理效应：

1. **表面亲和 (Surface Affinity, $\omega$)**：壁面"喜欢"某种组分的程度（线性项）。
2. **表面增强 (Surface Enhancement, $\chi_{bb}$)**：壁面对局部有序度或相互作用的修正（二次项）。

### 4.2 数学形式

$$f_{\text{surf}}(\phi_1|_0,\phi_2|_0) = \frac{k_BT}{\bar\nu}\Big[\underbrace{\omega_{b,1}\phi_1|_0+\omega_{b,2}\phi_2|_0}_{\text{吸附项}}+\underbrace{\chi_{bb,1}(\phi_1|_0)^2+\chi_{bb,2}(\phi_2|_0)^2+\chi_{bb,12}\phi_1|_0\phi_2|_0}_{\text{相互作用增强项}}\Big]$$

> **注**: 此处保留了溶质间的交叉表面增强项 $\chi_{bb,12}\phi_1|_0\phi_2|_0$，与第 4 节顶部定义一致（早期版本此式漏写了该项，现已更正）。

---

## 5. 吉布斯表面能泛函 ($\gamma$)：目标函数

### 5.1 概念：过剩自由能

**动机**: 这是我们整个求解的核心目标。$\gamma$ 代表了系统的**总过剩能量**。系统会自发调整其浓度分布 $\phi(z)$ 以最小化这个能量。

**组成**:

1. **梯度能**: $\frac{\kappa}{2}(\nabla\phi)^2$。形成界面需要付出代价（类似于表面张力），这惩罚了剧烈的浓度变化。
2. **体相势能**: $W(\phi)$。偏离储库状态的代价。
3. **表面能**: $f_{\text{surf}}$。与壁面接触的收益或代价。

### 5.2 数学形式

$$\gamma[\phi] = \int_0^\infty\left[W(\phi_1,\phi_2)+\sum_{i=1,2}\frac{\kappa_i}{2}\left(\frac{d\phi_i}{dz}\right)^2\right]dz + f_{\text{surf}}(\phi_i|_0)$$

---

## 6. 平衡态求解：变分法与欧拉-拉格朗日方程

**动机**: 根据热力学第二定律，平衡态对应于自由能极小值。我们通过变分法 ($\delta\gamma=0$) 寻找最优的浓度分布 $\phi(z)$。

### 6.1 变分推导

引入微扰 $\phi_i \to \phi_i+\delta\phi_i$:

$$\delta\gamma = \int_0^\infty\sum_{i=1,2}\left[\frac{\partial W}{\partial \phi_i}\delta\phi_i+\kappa_i\phi_i'\delta\phi_i'\right]dz+\sum_{i=1,2}\frac{\partial f_{\text{surf}}}{\partial \phi_i|_0}\delta\phi_i|_0$$

### 6.2 关键步骤：分部积分

处理梯度项 $\phi_i'\delta\phi_i'$ 是得到微分方程的关键：

$$\int_0^\infty \kappa_i\phi_i'\frac{d}{dz}(\delta\phi_i)\,dz = \underbrace{[\kappa_i\phi_i'\delta\phi_i]_0^\infty}_{\text{边界项}} - \int_0^\infty \frac{d}{dz}(\kappa_i\phi_i')\delta\phi_i\,dz$$

- $z\to\infty$: $\delta\phi_i=0$（储库固定），项消失。
- $z=0$: $\delta\phi_i \ne 0$，项保留为 $-\kappa_i\phi_i'(0)\delta\phi_i(0)$。

### 6.3 平衡方程

将结果整理，为了使 $\delta\gamma=0$ 恒成立，积分内部项和边界项必须分别通过零：

1. **体相方程 (ODE)**：控制内部浓度分布。

   $$\kappa_i\frac{d^2\phi_i}{dz^2} = \frac{\partial W}{\partial \phi_i} \quad (z>0)$$

2. **自然边界条件 (BC)**：控制表面处的吸附平衡。

   $$\kappa_i\frac{d\phi_i}{dz}\Big|_{z=0} = \frac{\partial f_{\text{surf}}}{\partial \phi_i|_0}$$

---

## 7. 能量守恒 (First Integral)：数学降阶技巧

### 7.1 概念：平移对称性

**动机**: 求解二阶耦合 ODE 非常困难。由于系统在 $z$ 方向上性质均匀（方程中不显含 $z$），根据 Noether 定理，必然存在一个守恒量（类似于力学中的"能量"）。这个守恒量可以将微分方程降阶。

### 7.2 推导证明

构造类哈密顿量 $H(z) = \sum_i \frac{\kappa_i}{2}(\phi_i')^2 - W$。对 $z$ 求导:

$$\frac{dH}{dz} = \sum_i \phi_i'\left[\kappa_i\phi_i''-\frac{\partial W}{\partial \phi_i}\right]$$

括号内即为 ODE，恒为 0。因此 $H(z)=\text{const}$。

在无穷远处 $\phi'\to 0$ 且 $W\to 0$，故常数为 0。

**结论公式**:

$$\sum_{i=1,2}\frac{\kappa_i}{2}\left(\frac{d\phi_i}{dz}\right)^2 = W(\phi_1,\phi_2)$$

这个公式告诉我们：**某一点的梯度能密度等于该点的广义势能密度**。

---

## 8. 算法实现原理：相空间积分法

**动机**: 为了计算 Pre-wetting 相变，我们需要比较两个不同状态（薄膜 vs 厚膜）的能量 $\gamma$。直接解 ODE 容易遇到刚性问题。利用第 7 节的能量守恒，我们可以将**空间积分转化为浓度路径积分**，从而极其高效地算出能量。

### 8.1 积分变换

原积分：$\gamma = \int_0^\infty(\text{Gradient}+W)\,dz + f_{\text{surf}}$。

代入守恒律：Gradient $=W$，被积函数变为 $2W$。

变量代换 $dz \to d\phi$（利用 $\frac{d\phi}{dz}=-\sqrt{2W/\kappa}$）：

$$\gamma = f_{\text{surf}} + \int_{\text{bulk}}^{\text{surf}} 2W \frac{1}{-\sqrt{2W/\kappa_{\text{eff}}}}\,d\phi = f_{\text{surf}} + \int_{\phi_\infty}^{\phi|_0}\sqrt{2\kappa_{\text{eff}}W}\,d\phi$$

### 8.2 线性路径假设 (Linear Path Ansatz)

在三元体系中，路径 $\vec\phi(z)$ 未知。对于参数对称的情况，我们假设 $\phi_1$ 和 $\phi_2$ 同步变化：

$$\vec\phi(t) = \vec\phi_\infty + t\cdot(\vec\phi_0-\vec\phi_\infty)$$

这也是代码中 `calculate_exact_gamma` 函数的核心数学依据。

### 8.3 预润湿判据

通过联立能量守恒和边界条件，我们得到寻找表面态 $\phi|_0$ 的代数方程：

$$\underbrace{W(\phi|_0)}_{\text{来自体相}} - \sum_i \frac{1}{2\kappa_i}\underbrace{\left(\frac{\partial f_{\text{surf}}}{\partial \phi|_0}\right)^2}_{\text{来自表面}} = 0$$

该方程的多个根对应不同的表面态（Thin/Thick）。当 $\gamma_{\text{thin}}=\gamma_{\text{thick}}$ 时，即为 Pre-wetting 相变点。
