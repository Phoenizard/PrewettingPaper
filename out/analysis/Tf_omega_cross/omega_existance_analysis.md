# 仅溶质吸附作用下的Prewetting存在性分析

![omega_cross_2x5](./omega_cross_2x5.png)

### Prewetting Line vs Omega

在本情况下，我们分别固定吸附力项进行分析，注意在本case中，除了墙面吸附之外，Bulk Energy中只有两个溶质的吸引，除此之外溶液没有对溶质的任何倾向。所以能产生prewetting line的现象的吸附力设置是对称的，当前我们取(0.25, -0.375)为中心去分析prewetting line的变化，事实上(-0.375， 0.25) 也会存在对应的prewetting。

当固定墙面对溶质1的吸附力，随着对溶质2的吸附力逐渐变弱（绝对值变小）可以发生prewetting的远程项的浓度要求更大，但是相应可以发生prewetting的环境更多了（体现在line更长），更加靠近Binodal的曲线。当我们固定墙面对溶质2的吸附力，随着对溶质1的排斥力越强，和上述变化一致，即可以发生prewetting的远程项的浓度要求更大，但是相应可以发生prewetting的环境更多了（体现在line更长），更加靠近Binodal的曲线。减弱对溶质2的吸附和加强对溶质1的排斥在结果上呈现出一样的变化。

![profile](./profiles/tf_profiles.png)

### Profile形状分析

在当前stage中，profile呈现出如图形状。薄膜态与厚膜态都是墙面上一高一低地吸附两种溶质。墙对两种溶质一个吸引、一个排斥，但溶质之间的吸引力很强，被排斥的溶质浓度也会被另一溶质的容度抬升起来。所以两个态里被排斥的溶质都呈现先上升后渐小的走向：紧贴墙的地方受墙的排斥，体积分数被压低；稍稍离开墙面，排斥迅速衰减，而此处被吸引的那一种溶质已经富集得很高，靠溶质间的吸引把它拉起来，形成一个不在墙上、而在墙外一小段距离处的峰值；再往外，被吸引的溶质自身降回远场组成，它失去牵引，也随之一起落下。

---

## 论文正文草稿（逐段对照，待审）

面向 main.tex 第 3 章的一个 case 子节。内容一律取自上面两块中文记录，只做语言调整与精确化，不添内容；例外是经裁定加入的三句：能出线的吸附力设置很少（用户补写的中文，已译）、线只出现在闭环 binodal 左下的一小段上、以及段 2 末尾对「线更长、更靠近 binodal」的物理释义。写作口径按物理表达：物种按角色称呼，墙用行为描述，参数值只出现在定义句的括号与图注里，不出现长度、距离的符号，也不出现扫描方向。

拟用小节标题：Prewetting in a mixture of two mutually attracting solutes

### 段 1（对应中文第 1 段：设定与对称性）

出处：「除了墙面吸附之外，Bulk Energy中只有两个溶质的吸引，除此之外溶液没有对溶质的任何倾向」「所以能产生prewetting line的现象的吸附力设置是对称的」「取(0.25, -0.375)为中心」「(-0.375, 0.25) 也会存在对应的prewetting」。

We next consider a mixture in which the only interaction in the bulk is an attraction between the two solutes, the solvent having no preference for either of them ($(\chi_{12}, \chi_{13}, \chi_{23}) = (-8.5, 0, 0)$). The wall carries no surface-enhanced interactions ($\chi_{bb,11} = \chi_{bb,22} = \chi_{bb,12} = 0$), so it acts on the mixture through its affinity for the two solutes alone. The mixture itself therefore does not distinguish the two solutes, and the wall is the only thing that does, so the wall affinities that produce a prewetting line come in symmetric pairs: a wall that repels the first solute and attracts the second gives the same prewetting as a wall that attracts the first and repels the second, with the two solutes interchanged. The conditions under which prewetting occurs are nevertheless restrictive: of the wall affinities surveyed, we find that only a small number produce a prewetting line at all. We take one wall of such a pair, which repels the first solute and attracts the second ($\omega_1 = 0.25$, $\omega_2 = -0.375$), and vary the two affinities about it one at a time.

### 段 2（对应中文第 2 段：线随两个亲附的变化）

出处：「固定墙面对溶质1的吸附力，随着对溶质2的吸附力逐渐变弱…浓度要求更大…line更长…更加靠近Binodal的曲线」「固定墙面对溶质2的吸附力，随着对溶质1的排斥力越强，和上述变化一致」「减弱对溶质2的吸附和加强对溶质1的排斥在结果上呈现出一样的变化」。加入一句经裁定的观察：线只出现在闭环 binodal 左下的一小段上。

Figure~\ref{fig:tf-cross} shows the prewetting line for these walls. The binodal of this mixture is a closed loop, and the line is found only along a short stretch of its lower left. Holding the affinity of the wall for the first solute fixed and weakening its attraction to the second, prewetting sets in only at reservoirs of higher solute content, but at more of them: the line becomes longer, and it lies closer to the binodal. Holding instead the affinity for the second solute fixed and strengthening the repulsion of the first changes the line in the same way. Weakening the attraction to the one solute and strengthening the repulsion of the other act alike. Prewetting therefore occurs over a wider span of reservoir compositions, so it is met under less restrictive conditions, but it is confined ever closer to the compositions at which the mixture demixes on its own, so that the surface transition is less clearly separated from bulk coexistence.

### 段 3（对应中文 Profile 形状分析）

出处：整段。

Figure~\ref{fig:tf-profiles} shows the composition profiles of the thin and the thick state. In both states the wall holds a large content of one solute and a small content of the other. The wall attracts one solute and repels the other, but the attraction between the solutes is strong, so the content of the repelled solute is raised by the content of the solute that gathers at the wall. In both states the repelled solute therefore rises before it decays: right against the wall its content is held down by the repulsion of the wall; a short distance away the repulsion has decayed, while the attracted solute is there already strongly enriched and draws the repelled solute up through the attraction between them, so that it reaches its largest content not at the wall but a short distance away from it; further out the attracted solute returns to its reservoir content, the pull is lost, and the repelled solute falls with it.

### 与上一版的出入

- 删除：墙与体相对立的机制、厚膜态比薄膜态携带更多溶质且层更厚、最强排斥处线贴合 binodal 全长、「越过线吸附从薄膜跳到厚膜」（第 2 节已交代）。都是中文记录里没有的。
- 改写：不再用「上排 / 下排、角色对调」的对比句式（那是用来讲两个亲附作用不同的，这里两者作用相同）；不再说线位于 binodal 的某一侧（闭环没有侧）。

### 尚未落笔的两件事

- 图注：两张图的 caption 未写。3.1 节的做法是把全部参数值放进图注，正文不带。
- 交叉引用与图文件名：草稿里用了 fig:tf-cross 与 fig:tf-profiles 两个标签，图还没进 manuscript/figures/，按你的裁定这轮不做图。