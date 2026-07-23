# NJP 论文分节结构调研 + 本文大纲方案

目的：确定投 New Journal of Physics（NJP）时的分节风格，并据此为本文设计标题/分节。
锚点：用户引用的 Zhao 等 2021 论文。样本：NJP 15 篇 + 对照期刊（PRE/PRL/Soft Matter）5 篇 + NJP/IOP 作者指南。
本文内容线（依 CLAUDE.md）：模型 -> 现象 -> 物理解释。

## 一、核心结论

- NJP 不强制 IMRaD。IOP 作者指南原文是 “should follow the Introduction, Methods,
  Results and Discussion system, and usually consist of ...”，措辞为 “usually/推荐”，
  且明确允许作者自定分节标题，只要求可读、可复现。
  指南：https://publishingsupport.iopscience.iop.org/journals/new-journal-of-physics/
- 实测样本里几乎没有纯 IMRaD。NJP 20 篇（含对照）中，主流是内容自定式
  （content-derived）：Introduction 与 Conclusion/Discussion 作首尾书挡，中间各节
  用描述该节具体物理内容的标题命名（模型节、相图节、机制节）。少数是 hybrid（在
  Model/Results 轻骨架上挂具体内容），纯实验向论文才用 generic-IMRaD。
- 锚点论文正是内容自定式。二三级标题甚至把结论写进标题
  （如 “Membrane binding favors complete wetting and reduces the contact angle”）。
- 模型与方法的位置有一致惯例：模型定义放正文的一个描述性命名节（不叫 “Methods”，
  而叫 “Thermodynamics of ...” / “A model for ...”）；冗长的数值方法/推导下沉附录
  （NJP 2021+ 还多一个 Data availability statement 尾节）。对照期刊（PRE/PRL）更精简、
  多不设附录、方法直接留正文。
- 对本文的直接含义：可以采用锚点论文式的内容自定分节，把模型压到必要深度、数值方法
  下沉附录，正文只讲物理——这既符合 NJP 惯例，也符合 CLAUDE.md “模型与数值不入正文”。

## 二、锚点论文结构（Zhao et al., NJP 23, 123003, 2021）

DOI 10.1088/1367-2630/ac320b｜https://iopscience.iop.org/article/10.1088/1367-2630/ac320b

- 1. Introduction
- 2. Thermodynamics of phase separation with membrane binding
  - 2.1 Thermodynamics of a semi-infinite system
  - 2.2 Free energy minimization in a semi-infinite system
  - 2.3 Surface phase diagrams obtained via graphical construction
  - 2.4 Bulk and membrane free energies
- 3. Effects of membrane binding on wetting, prewetting and surface phase transitions
  - 3.1 Wetting and prewetting without phase separation in the membrane
    - 3.1.1 Membrane binding favors complete wetting and reduces the contact angle
    - 3.1.2 Prewetting transition shifts to lower concentrations
    - 3.1.3 Prewetting transition persists below the bulk critical point
  - 3.2 Wetting and prewetting with phase separation in the membrane
    - 3.2.1 Phase transitions between four distinct surface states
    - 3.2.2 Correlated and anti-correlated jumps of surface order parameters
  - 3.3 Wetting, prewetting and surface phase transitions with related bulk and surface interactions
- 4. Discussion
- Acknowledgments｜Data availability statement
- Appendix A Graphical construction ... ｜ Appendix B Condition of critical prewetting ｜ Appendix C 与既往工作对比

骨架读法：Introduction -> 模型节（描述性命名，含图解构造方法） -> 效应节（按现象/机制
逐一命名，二级标题即结论） -> Discussion；模型细节与数值方法在正文命名节 + 附录。

## 三、样本逐篇结构表

NJP（相分离/软物质/统计物理/生物物理理论）：

- Droplet ripening in concentration gradients（2017, aa6b84）：Introduction -> Local
  regulation -> Spatial organization -> Dynamics of a single drop -> Ripening of
  multiple drops -> Conclusion and outlook。content-derived。方法在正文。
- Discontinuous switching of position of two coexisting phases（2018, aad173）：
  Introduction -> Equilibrium model -> Concentration profiles -> Discontinuous
  transition ... -> ... periodic domain ... -> Experimental verification and outlook。
  content-derived。附录含 ternary Flory-Huggins 梯度项推导。
- Stress granule formation via ATP depletion（2018, aab549）：Intro -> Minimal models
  -> Steady-state ... -> Model selection -> Model comparison ... -> Summary and
  discussion。content-derived。方法下沉附录 A/B/C。
- Theory of wetting dynamics with surface binding（2024, ad80bb, 同作者谱系）：
  Introduction -> Non-equilibrium thermodynamics of wetting with surface binding ->
  Wetting dynamics with surface binding -> Conclusions and outlook。content-derived。
  附录 A-E。
- Wetting and pattern formation in non-reciprocal ternary phase separation（2025,
  ae2883, 题材最近）：Introduction -> A model for non-reciprocal ternary phase
  separation -> Quasi-static behaviour ... -> Pattern formation ... -> Conclusion
  and outlook。content-derived。附录含 bulk 自由能推导 + 数值格式。
- Thermodynamically consistent flocking（2024, ad4dd6）：Introduction -> Thermodynamic
  consistency ... -> Phase diagram for monotonic ... -> ... non-monotonic ...。
  content-derived。数值细节在附录。
- Species interconversion ... transient phase separation（2025, adccf1）：Introduction
  -> Particle deformation and species populations -> Transient phase separation ->
  Coarse-graining from particles to fields -> Discussion。content-derived。
- Bose-Einstein-like condensation ... （2020, ab90d8）：Introduction -> Scalar active
  matter ... -> Characterization of the condensation transition -> Generalized
  thermodynamics -> Degeneracy and stability ... -> Concluding remarks。content-derived。
- Shape-driven emergent behavior in active particle mixtures（2022, ac7161）：
  Introduction -> Methods -> Emergent behavior diagram -> Conclusion。hybrid。
  模拟协议在附录。
- Systematically extending classical nucleation theory（2018, aae174）：Introduction
  -> Theory -> CNT and generalizations -> Conclusions。hybrid（二级标题内容化）。
- Probing phase transitions ... single spin quantum sensor（2019, ab482d, 实验向）：
  Introduction -> Materials and methods -> Results -> Discussion -> Conclusions。
  generic-IMRaD。
- Periodic patterns displace active phase separation（2021, abe814）：Introduction ->
  Generic model for a conserved order parameter field -> Phase separation and its
  transition ... -> Multistability ... -> Stable secondary hybrid states -> Summary
  and discussions。content-derived。模拟方法在附录。
- Deforming polar active matter in a scalar field gradient（2023, acb2e5）：
  Introduction -> Model -> Scalar field only -> Polar field only -> Scalar and polar
  field -> Discussion。hybrid（第 2 节通用 Model，结果节按情形命名）。附录放推导。
- Uncovering novel phase transitions in dense dry polar active fluids（2021, abd8c0）：
  Introduction -> Dry polar active fluids -> MIPS -> Contact inhibition ... -> Phase
  behaviour -> Hydrodynamic description ... -> ... -> Conclusion and outlook。
  content-derived。数值方法整体下沉附录 A/B/C。

对照期刊（PRE/PRL/Soft Matter，wetting/prewetting/表面相变）：

- Binding potential and wetting behavior of binary liquid mixtures on surfaces
  （PRE 109, 024801, 2024；题材最贴近）：Introduction -> Surface thermodynamics and
  binding potentials -> Lattice DFT for binary mixtures -> Bulk mixture phase
  behaviour -> Density profiles at the free interface -> DFT results for miscible ...
  -> Results for immiscible ... -> Binding potentials and profiles ... -> Droplets at
  walls -> Concluding remarks。content-derived，10 节，无附录，方法全在正文。
- Scaling of wetting and prewetting transitions on nanopatterned walls（PRE 100,
  032801, 2019）：Introduction -> Density functional model -> Results: scaling of the
  pre-wetting lines -> Concluding remarks。hybrid，方法在正文。
- Complete Wetting and Drying at Sinusoidal Walls（PRE 112, 015502, 2025）：
  Introduction -> Microscopic model -> Mesoscopic approach -> Numerical results ->
  Conclusion。hybrid。
- Crossover scaling of apparent first-order wetting ...（PRE 93, 062802, 2016）：
  Introduction -> Scaling and fluctuation regimes ... -> Apparent first-order
  behaviour ... -> Rounded pre-wetting transitions ... -> Conclusions。content-derived。
- Apparent First-Order Wetting ... 2D Ising（PRL 116, 046101, 2016）：无分节，连续正文。

风格分布（不含无分节/纯实验）：content-derived 约 12 篇，hybrid 约 4 篇，generic-IMRaD
1 篇。结论稳健：内容自定式是本领域理论论文的主流，纯定式是少数。

## 四、本文内容 -> 分节映射

本文只做 omega topic 单 case 分类下的分析（范围已收窄）。三段内容：

- 模型：三元混合物（溶剂 + 溶质 1 + 溶质 2）对化学偏好硬墙的 Gibbs 表面自由能框架
  gamma = ∫(W + 梯度项) + f_surf；墙面亲附力 omega1、omega2 只进 f_surf。正文给必要
  深度，完整推导与数值放附录。
- 现象：pre-wetting line 在 (phi1_inf, phi2_inf) 平面的 extent，用长度 L（MST 弧长）
  与到 binodal 距离 dist_mean 刻画，随 omega1、omega2 变化——omega1 主控、omega2 近乎
  惰性；弱吸引则线长且贴 binodal，强吸引则线短伸深、在强吸引角熄灭。
- 物理解释：W 地形机制——∂²W/∂φ1² = 1/φ1 + 1/φs - 2χ13 可变负 -> 沿 φ1 出第二阱；
  ∂²W/∂φ2² = 1/φ2 + 1/φs 恒正 -> 沿 φ2 无阱；溶剂熵经非对角曲率 1/φs - χ13 弱耦合两
  溶质，故 omega2 有二阶效应（近乎惰性而非完全解耦）。

## 五、大纲方案

定位（用户 2026-07-21 定）：体相三元 Flory-Huggins 自由能与 binodal 均引自 Omar
Adame-Arana 2019 博士论文《Chemical control of liquid phase separation in the cell》
(TU Dresden)，是背景、不是本文卖点。卖点是三元体系下的 pre-wetting 条件——墙面亲附力
如何决定 pre-wetting 是否发生、在 (phi1_inf, phi2_inf) 平面占多大地。binodal 只作参考
线（度量 pre-wetting line 到它的距离）。故分节标题重心放 pre-wetting 与墙面亲附，模型
与 binodal 压成一小段引用背景、连同数值一起下沉附录，不设招牌“模型节”。

层级订正（用户 2026-07-21）：第三章是容器，对不同 case（chi 拓扑等）各挑一个值得
注意的分析点，暂定三个并列。omega 主控（“Wall affinity controls the extent of
prewetting”）不是章级标题，而是 (0, 2.8, 0) case（T-a 拓扑）这一个 case 的分析点，
收进 3.1。每个 case 子节自成“现象 + 机制”小闭环（含各自的 W 地形解释），不再单列
跨 case 机制章。另两个 case 及其分析点待定，留占位。

方案 A（推荐，锚点论文式内容自定 + 模型/数值/binodal 下沉附录）

- 1. Introduction
  （末尾一小段交代背景：体相三元 FH 与 binodal 取自 Adame-Arana 2019，作为给定参考）
- 2. Prewetting in a ternary mixture at a preferential wall
  （把 pre-wetting 条件讲清：墙面表面能 f_surf 与 omega 如何引出 thin->thick 跳变；
  定义 pre-wetting line 及其 extent——长度 L、到 binodal 参考线的距离；一张代表相图。
  体相模型只用一两句引用带过。此节为全章共用的框架。）
- 3. Prewetting across representative bulk topologies（容器，暂定三个 case 并列）
  - 3.1 Case (0, 2.8, 0)（T-a）：Wall affinity controls the extent of prewetting
    - 现象：omega1 主控、omega2 近乎惰性；弱吸引线长且贴 binodal；强吸引线短伸深、
      在强吸引角熄灭。
    - 机制：W 地形——沿 φ1 出第二阱、沿 φ2 无阱；溶剂熵非对角弱耦合给 omega2 二阶效应。
  - 3.2 Case [待定]：[分析点待定]
  - 3.3 Case [待定]：[分析点待定]
- 4. Conclusion and outlook
- Appendix A 体相三元 FH 自由能与 binodal（引 Adame-Arana）+ 表面能完整式与无量纲化 ｜
  Appendix B 数值方法（固定网格 Newton + 双向扫描 + crossing；extent 度量 MST 定义）｜
  Data availability statement

理由：最贴合锚点论文与 NJP 主流；标题即结论，正文纯讲 pre-wetting 物理；借来的模型/
binodal 与数值都在附录，不喧宾夺主，也满足 CLAUDE.md “模型与数值不入正文”。第 2 节立
框架，第 3 节按 case 逐一给分析点、每个 case 自带现象与机制。

第 3 章容器标题：Prewetting across representative bulk topologies（已定）。

方案 B（更紧凑，PRE 对照式，pre-wetting 条件与结果各一节，方法仍进附录）

- 1. Introduction
- 2. Prewetting condition at a preferential wall
- 3. Effect of wall affinity on the extent of prewetting
- 4. Free-energy-landscape mechanism
- 5. Conclusions

理由：若希望更短或改投 PRE。仍保持模型/binodal 为背景、不设招牌模型节。

方案 C（合并现象与机制，四节极简）

- 1. Introduction
- 2. Prewetting in a ternary mixture at a preferential wall
- 3. Wall affinity and the extent of prewetting
  （现象 + 机制交织，一节讲完 omega1 主控 + W 地形解释）
- 4. Discussion

理由：内容量偏小时的保底结构。风险：现象与机制挤在一节，说理密度过高、不易读。

推荐 A：与目标期刊惯例最合、与 CLAUDE.md 约束最契合，且把借来的模型/binodal 摆正为
背景、正文重心落在 pre-wetting 条件与墙面亲附。标题措辞可再打磨。
