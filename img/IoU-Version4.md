# 第 1 章 采样域、TV 误差与扫描—混合模板

本章建立全文通用的采样组织方式与误差度量方式。

- **采样组织方式**：把复杂目标分布写成“可扫描索引上的混合”，先用一次扫描统计每个索引的质量（Pass 1），再为每个输出位置独立指派索引并分组（Planning），最后用第二次扫描执行各索引下的条件采样并回填（Pass 2），从而得到 $t$ 个 i.i.d.（有放回）样本。
- **误差度量方式**：统一用总变差距离（Total Variation, TV）度量单次边缘分布的偏差，并给出在“混合—条件—回填”结构下的可组合规则。
- **几何分解接口**：沿 $x$ 轴把平面划分为不交竖直 slab，使二维交面积可写为 slab 分量之和；在可分离权重结构下，每个 slab 的质量可写为一维积分 $\int F_i(y)G_i(y)\,dy$，从而把质量统计与条件采样都落到扫描线可维护的接口上。

---

## 1.1 采样对象与目标分布

### 1.1.1 输入：两类半开矩形

给定两类二维轴对齐半开矩形集合 $R_c$ 与 $R_{\bar c}$，总规模 $N=|R_c|+|R_{\bar c}|$。每个矩形 $r$ 表示为
$$
r=[L_x(r),R_x(r))\times [L_y(r),R_y(r)),
\qquad L_x(r)<R_x(r),\ \ L_y(r)<R_y(r).
$$

半开语义的含义是：边界接触不产生正面积交叠。例如若 $R_x(r)=L_x(s)$，则 $r$ 与 $s$ 在 $x$ 方向交叠长度为 $0$。这一约定能让扫描线分解与事件顺序具有严格一致的“端点归属”，避免端点处的双计数。

矩形面积定义为
$$
A(r)=\mathrm{Area}(r)=(R_x(r)-L_x(r))(R_y(r)-L_y(r))>0.
$$

---

### 1.1.2 交面积、并面积与 IoU 权重

对一维半开区间 $[a,b)$，定义其长度
$$
\mathrm{len}([a,b))=\max\{0,b-a\}.
$$

对任意两矩形 $r,s$，定义 $x$、$y$ 方向交叠长度
$$
\Delta_x(r,s)=\mathrm{len}\!\Big([\max\{L_x(r),L_x(s)\},\ \min\{R_x(r),R_x(s)\})\Big),
$$
$$
\Delta_y(r,s)=\mathrm{len}\!\Big([\max\{L_y(r),L_y(s)\},\ \min\{R_y(r),R_y(s)\})\Big).
$$

交面积与并面积分别为
$$
I(r,s)=\Delta_x(r,s)\Delta_y(r,s),\qquad 
U(r,s)=A(r)+A(s)-I(r,s).
$$

IoU 权重定义为
$$
W_{\mathrm{IoU}}(r,s)=\frac{I(r,s)}{U(r,s)}
=\frac{I(r,s)}{A(r)+A(s)-I(r,s)}.
$$

当 $I(r,s)>0$ 时有 $W_{\mathrm{IoU}}(r,s)\in(0,1]$。这是因为 $I(r,s)\le \min\{A(r),A(s)\}\le \frac{A(r)+A(s)}{2}$，从而
$A(r)+A(s)-I(r,s)\ge I(r,s)$，因此 $0<W_{\mathrm{IoU}}(r,s)\le 1$。

---

### 1.1.3 采样域、支撑集与目标分布

本书始终将样本空间写为跨类有序对集合
$$
\Omega = R_c\times R_{\bar c}.
$$
权重函数 $W_{\mathrm{IoU}}$ 视为定义在全域 $\Omega$ 上；当 $I(r,s)=0$ 时，$W_{\mathrm{IoU}}(r,s)=0$，因此不会对任何归一化和或概率造成影响。

定义正权重支撑集（也即“有效采样域”）
$$
J=\{(r,s)\in \Omega:\ I(r,s)>0\}.
$$

定义总权重（归一化常数）
$$
Z_{\mathrm{IoU}}=\sum_{(r,s)\in\Omega} W_{\mathrm{IoU}}(r,s)
=\sum_{(r,s)\in J} W_{\mathrm{IoU}}(r,s).
$$

当 $Z_{\mathrm{IoU}}>0$ 时，目标分布为
$$
\mathcal D_{\mathrm{IoU}}(r,s)=\frac{W_{\mathrm{IoU}}(r,s)}{Z_{\mathrm{IoU}}},\qquad (r,s)\in\Omega.
$$
由于 $\mathcal D_{\mathrm{IoU}}$ 在 $\Omega\setminus J$ 上恒为 $0$，因此也可等价地将其看作定义在 $J$ 上的分布。

**约定（零质量实例）**：若 $Z_{\mathrm{IoU}}=0$（等价于 $J=\varnothing$），则不存在正权重对象，采样任务输出空序列。

---

## 1.2 总变差距离与 i.i.d. 输出约束

### 1.2.1 TV 距离与采样目标

在有限样本空间 $\Omega$ 上，两分布 $\mathcal P,\mathcal Q$ 的总变差距离定义为
$$
d_{\mathrm{TV}}(\mathcal P,\mathcal Q)=\frac12\sum_{x\in\Omega}\big|\mathcal P(x)-\mathcal Q(x)\big|.
$$
并且等价地
$$
d_{\mathrm{TV}}(\mathcal P,\mathcal Q)=\max_{A\subseteq\Omega}\big|\mathcal P(A)-\mathcal Q(A)\big|.
$$

给定误差参数 $\varepsilon\in(0,1)$ 与样本数 $t$，目标是输出 $t$ 个 **相互独立且同分布（i.i.d.）** 的样本 $Z_1,\dots,Z_t\in\Omega$（等价地 $Z_j\in J$），使其单次边缘分布 $\tilde{\mathcal D}$ 满足
$$
d_{\mathrm{TV}}(\tilde{\mathcal D},\mathcal D_{\mathrm{IoU}})\le \varepsilon.
$$

---

### 1.2.2 TV 的三个基础工具

**引理 1.1（三角不等式）**  
对任意分布 $\mathcal P,\mathcal Q,\mathcal R$，有
$$
d_{\mathrm{TV}}(\mathcal P,\mathcal R)\le d_{\mathrm{TV}}(\mathcal P,\mathcal Q)+d_{\mathrm{TV}}(\mathcal Q,\mathcal R).
$$

*证明*：由 $\ell_1$ 范数三角不等式，且 $d_{\mathrm{TV}}=\frac12\|\cdot\|_1$。$\square$

---

**引理 1.2（同权混合的凸性：误差不放大）**  
设 $\{\lambda_i\}_{i=1}^m$ 为概率向量（$\lambda_i\ge0,\ \sum_i\lambda_i=1$），且 $\mathcal P_i,\mathcal Q_i$ 为同一有限域上的分布，则
$$
d_{\mathrm{TV}}\Big(\sum_{i=1}^m \lambda_i\mathcal P_i,\ \sum_{i=1}^m \lambda_i\mathcal Q_i\Big)
\le \sum_{i=1}^m \lambda_i\, d_{\mathrm{TV}}(\mathcal P_i,\mathcal Q_i).
$$
特别地，若对所有 $i$ 有 $d_{\mathrm{TV}}(\mathcal P_i,\mathcal Q_i)\le \varepsilon$，则混合后仍满足
$$
d_{\mathrm{TV}}\Big(\sum_{i=1}^m \lambda_i\mathcal P_i,\ \sum_{i=1}^m \lambda_i\mathcal Q_i\Big)\le \varepsilon.
$$

*证明*：
$$
\left\|\sum_i \lambda_i(\mathcal P_i-\mathcal Q_i)\right\|_1
\le \sum_i \lambda_i\|\mathcal P_i-\mathcal Q_i\|_1
=2\sum_i\lambda_i\,d_{\mathrm{TV}}(\mathcal P_i,\mathcal Q_i).
$$
两边除以 $2$ 即得。$\square$

---

**引理 1.3（乘积分布的 TV 线性放大界）**  
设 $\mathcal P,\mathcal Q$ 是同一有限域上的分布，则对任意整数 $t\ge 1$，
$$
d_{\mathrm{TV}}(\mathcal P^{\otimes t},\mathcal Q^{\otimes t})\le t\cdot d_{\mathrm{TV}}(\mathcal P,\mathcal Q).
$$

*证明*：用归纳。对 $t+1$，
$$
\mathcal P^{\otimes (t+1)}-\mathcal Q^{\otimes (t+1)}
=(\mathcal P^{\otimes t}-\mathcal Q^{\otimes t})\otimes \mathcal P
+\mathcal Q^{\otimes t}\otimes(\mathcal P-\mathcal Q).
$$
取 $\ell_1$ 范数并用 $\|\mu\otimes\nu\|_1=\|\mu\|_1\|\nu\|_1$ 以及 $\|\mathcal P\|_1=\|\mathcal Q\|_1=1$，得
$$
\|\mathcal P^{\otimes (t+1)}-\mathcal Q^{\otimes (t+1)}\|_1
\le \|\mathcal P^{\otimes t}-\mathcal Q^{\otimes t}\|_1+\|\mathcal P-\mathcal Q\|_1.
$$
两边除以 $2$ 并归纳即得结论。$\square$

---

## 1.3 混合恒等式与扫描—混合—回填模板

本节把“按权重采样”统一写为混合形式，并给出一个与扫描线实现高度兼容的 i.i.d. 输出模板。核心思想是：**只要目标分布能够写成索引上的混合，且每个索引下可独立做条件采样，就能用两遍扫描一次性生成 $t$ 个 i.i.d. 样本**。

---

### 1.3.1 块混合：由空间划分诱导的混合

设有限样本空间 $\Omega$ 被划分为不交并
$$
\Omega=\biguplus_{b\in\mathcal B}\Omega_b.
$$
给定非负权重函数 $W:\Omega\to \mathbb R_{\ge0}$，定义
$$
M_b=\sum_{x\in\Omega_b}W(x),\qquad
M=\sum_{b\in\mathcal B} M_b.
$$
当 $M>0$ 时定义目标分布
$$
\mathcal D(x)=\frac{W(x)}{M},\qquad x\in\Omega,
$$
并对 $M_b>0$ 的块定义块内条件分布
$$
\mathcal D_b(x)=\frac{W(x)}{M_b},\qquad x\in\Omega_b.
$$

**定理 1.4（块混合两阶段采样）**  
若按如下过程生成随机变量 $X$：

1. 先抽块 $B$，满足 $\Pr[B=b]=M_b/M$；
2. 条件于 $B=b$，从 $\Omega_b$ 按 $\mathcal D_b$ 抽取 $X$；

则 $X$ 的边缘分布为 $\mathcal D$。

*证明*：对任意 $x\in\Omega$，其属于唯一块 $b(x)$，于是
$$
\Pr[X=x]=\Pr[B=b(x)]\Pr[X=x\mid B=b(x)]
=\frac{M_{b(x)}}{M}\cdot \frac{W(x)}{M_{b(x)}}=\frac{W(x)}{M}.
$$
$\square$

---

### 1.3.2 分量混合：由权重可加分解诱导的混合

很多情形下样本空间不易划分，但权重易写成非负分量之和：
$$
W(x)=\sum_{b\in\mathcal B} W_b(x),\qquad W_b(x)\ge 0.
$$
定义分量质量与总质量
$$
M_b=\sum_{x\in\Omega} W_b(x),\qquad
M=\sum_{b\in\mathcal B} M_b=\sum_{x\in\Omega}W(x).
$$
对 $M_b>0$ 的分量定义条件分布
$$
\mathcal D^{(b)}(x)=\frac{W_b(x)}{M_b},\qquad x\in\Omega.
$$

**定理 1.5（分量混合两阶段采样）**  
若按如下过程生成随机变量 $X$：

1. 先抽分量 $B$，满足 $\Pr[B=b]=M_b/M$；
2. 条件于 $B=b$，从 $\Omega$ 按 $\mathcal D^{(b)}$ 抽取 $X$；

则 $X$ 的边缘分布为 $\mathcal D(x)=W(x)/M$。

*证明*：对任意 $x\in\Omega$，
$$
\Pr[X=x]=\sum_{b\in\mathcal B}\Pr[B=b]\Pr[X=x\mid B=b]
=\sum_{b\in\mathcal B}\frac{M_b}{M}\cdot \frac{W_b(x)}{M_b}
=\frac{1}{M}\sum_{b\in\mathcal B}W_b(x)=\frac{W(x)}{M}.
$$
$\square$

---

### 1.3.3 扫描—混合—回填：Two-Pass Scan–Mix–Backfill

许多几何结构天然提供一个可顺序遍历的索引集合 $\mathcal B$（例如按坐标排序后的 slab）。本节给出一个模板，只依赖两个接口：

- `Mass(b)`：返回索引 $b$ 的质量 $M_b\ge 0$；
- `CondSample(b,k)`：当 $M_b>0$ 时，返回 $k$ 个相互独立的样本，其单次边缘分布为索引 $b$ 下的条件分布 $\mathcal D_b$（或 $\mathcal D^{(b)}$）。

其中 $\mathcal D_b$ 可以来自块混合或分量混合；只要满足
$$
\tilde{\mathcal D}(x)=\sum_{b\in\mathcal B}\frac{M_b}{M}\,\mathcal D_b(x)
$$
就是模板最终生成的单次边缘分布。

---

**算法 1.1（Two-Pass Scan–Mix–Backfill）**  
输入：索引集合 $\mathcal B$（给定固定遍历顺序）、接口 `Mass/CondSample`、样本数 $t$。  
输出：$t$ 个 i.i.d. 样本（有放回）或空序列。

1. **Pass 1（统计质量）**：顺序遍历 $b\in\mathcal B$，计算 $M_b=\mathrm{Mass}(b)$；令 $M=\sum_{b}M_b$。若 $M=0$，返回空序列。
2. **Planning（独立指派并分组）**：对每个位置 $j=1,\dots,t$，独立抽取索引 $B_j$，满足 $\Pr[B_j=b]=M_b/M$；把位置 $j$ 追加到列表 $L[b]$ 中。
   - 实现层面可用前缀和数组做二分抽样（$O(\log|\mathcal B|)$/次），或用 alias 表实现 $O(1)$/次抽样；模板本身只要求分布正确且各 $B_j$ 相互独立。
3. **Pass 2（条件采样并回填）**：再次顺序遍历 $b\in\mathcal B$，令 $k=|L[b]|$。若 $k>0$，调用 `CondSample(b,k)` 得到 $k$ 个相互独立样本，并按 $L[b]$ 指定位置回填到输出数组中。

---

**定理 1.6（Two-Pass Scan–Mix–Backfill 的 i.i.d. 正确性）**  
假设满足以下独立性条件：

- (i) Planning 阶段的 $B_1,\dots,B_t$ 相互独立，且 $\Pr[B_j=b]=M_b/M$；
- (ii) 对每个 $b$，`CondSample(b,k)` 返回 $k$ 个相互独立、边缘分布为 $\mathcal D_b$ 的样本；
- (iii) 不同 $b$ 的 `CondSample` 调用所使用的随机性彼此独立，且与 Planning 阶段的随机性独立；

则算法 1.1 输出 $Z_1,\dots,Z_t$ 相互独立同分布，且单次边缘分布为
$$
\tilde{\mathcal D}(x)=\sum_{b\in\mathcal B}\frac{M_b}{M}\,\mathcal D_b(x).
$$

*证明*：对任意位置 $j$ 与任意 $x$，由全概率公式
$$
\Pr[Z_j=x]=\sum_{b\in\mathcal B}\Pr[B_j=b]\Pr[Z_j=x\mid B_j=b]
=\sum_{b\in\mathcal B}\frac{M_b}{M}\,\mathcal D_b(x).
$$
独立性方面：每个位置 $j$ 仅依赖其独立的指派索引 $B_j$ 以及该索引下条件采样器产生的一个样本；不同索引 $b$ 的条件采样随机性相互独立；回填是对独立随机变量族的确定性重标号，不引入新的相关性。因此 $Z_1,\dots,Z_t$ 相互独立同分布。$\square$

---

### 1.3.4 混合结构下的 TV 误差组合

在后续构造中，常见情形是：混合权重保持精确，但索引内的条件分布采用近似实现；或混合权重本身也存在近似。下面给出两条可直接调用的 TV 上界。

设精确混合为 $\mathcal M=\sum_b \lambda_b \mathcal D_b$，其中 $\lambda_b=M_b/M$。

**引理 1.7（条件分布近似：同权混合不放大）**  
若对每个 $b$ 有
$$
d_{\mathrm{TV}}(\tilde{\mathcal D}_b,\mathcal D_b)\le \varepsilon_b,
$$
则在相同混合系数 $\lambda_b$ 下，
$$
d_{\mathrm{TV}}\Big(\sum_b \lambda_b\tilde{\mathcal D}_b,\ \sum_b\lambda_b\mathcal D_b\Big)
\le \sum_b \lambda_b \varepsilon_b
\le \max_b \varepsilon_b.
$$

*证明*：由引理 1.2。$\square$

---

**引理 1.8（混合系数近似：系数偏差的 $\ell_1/2$ 上界）**  
设 $\lambda,\hat\lambda$ 为 $\mathcal B$ 上两概率向量，且对每个 $b$ 给定分布 $\mathcal D_b$。定义
$$
\mathcal M=\sum_b\lambda_b\mathcal D_b,\qquad 
\hat{\mathcal M}=\sum_b\hat\lambda_b\mathcal D_b.
$$
则
$$
d_{\mathrm{TV}}(\mathcal M,\hat{\mathcal M})
\le \frac12\sum_{b\in\mathcal B}|\lambda_b-\hat\lambda_b|.
$$

*证明*：
$$
\|\mathcal M-\hat{\mathcal M}\|_1
=\sum_{x}\left|\sum_b(\lambda_b-\hat\lambda_b)\mathcal D_b(x)\right|
\le \sum_b|\lambda_b-\hat\lambda_b|\sum_x \mathcal D_b(x)
=\sum_b|\lambda_b-\hat\lambda_b|.
$$
两边除以 $2$ 得结论。$\square$

---

## 1.4 采样域的 slab 分解与质量积分化

本节沿 $x$ 轴把平面分解为不交 slab，并在可分离权重结构下把每个 slab 的质量写成一维积分。该分解将二维几何权重的“总量统计”与“条件采样”都转化为扫描线可维护的对象。

---

### 1.4.1 slab 的定义、覆盖指示与事件顺序

收集所有输入矩形的 $x$ 端点
$$
X=\{L_x(r),R_x(r)\mid r\in R_c\cup R_{\bar c}\},
$$
去重排序得到
$$
x_0<x_1<\cdots<x_m,\qquad m\le 2N-1.
$$
定义第 $i$ 个竖直 slab
$$
S_i=[x_i,x_{i+1}),\qquad \Delta x_i=x_{i+1}-x_i>0,\qquad i=0,\dots,m-1.
$$

对任意矩形 $r$，定义其对 slab 的覆盖指示
$$
\chi_r(i)=\mathbf 1\big[S_i\subseteq [L_x(r),R_x(r))\big].
$$

**引理 1.9（slab 内覆盖一致性）**  
对任意矩形 $r$ 与 slab $S_i$，要么 $\chi_r(i)=1$，要么 $S_i\cap [L_x(r),R_x(r))=\varnothing$；不存在“只覆盖 slab 的一部分”的情形。

*证明*：若 $S_i$ 与 $[L_x(r),R_x(r))$ 相交，则存在 $x\in[x_i,x_{i+1})$ 使 $L_x(r)\le x<R_x(r)$。由于 $(x_i,x_{i+1})$ 内没有任何端点，而 $L_x(r),R_x(r)$ 都属于端点集合 $X$，因此必有 $L_x(r)\le x_i$ 且 $R_x(r)\ge x_{i+1}$，即 $S_i\subseteq[L_x(r),R_x(r))$。$\square$

**事件顺序（与半开语义对齐）**：按 $x$ 从小到大扫描；在同一 $x$ 处，先处理 $\mathrm{END}$（删除 $R_x(\cdot)=x$ 的矩形），再处理 $\mathrm{START}$（插入 $L_x(\cdot)=x$ 的矩形）。这样在进入 slab $S_i=[x_i,x_{i+1})$ 时，活跃集合恰好对应覆盖该 slab 的矩形集合。

---

### 1.4.2 交面积的 slab 求和形式

定义两矩形在 $y$ 方向的交叠长度
$$
h(r,s)=\mathrm{len}\!\Big([\max\{L_y(r),L_y(s)\},\ \min\{R_y(r),R_y(s)\})\Big).
$$
注意 $h(r,s)$ 与 $x$ 无关。

定义两矩形在 slab $S_i$ 上的交面积贡献
$$
I_i(r,s)=\Delta x_i\cdot \chi_r(i)\chi_s(i)\cdot h(r,s).
$$

**定理 1.10（交面积的 slab 分解）**  
对任意两矩形 $r,s$，有
$$
I(r,s)=\sum_{i=0}^{m-1} I_i(r,s).
$$

*证明*：由引理 1.9，在每个 slab 上，$r$ 与 $s$ 在 $x$ 方向要么同时覆盖整个 slab（交叠长度为 $\Delta x_i$），要么至少一个不覆盖（交叠长度为 $0$）。因此
$$
\Delta_x(r,s)=\sum_{i=0}^{m-1}\Delta x_i\,\chi_r(i)\chi_s(i).
$$
再乘以 $h(r,s)=\Delta_y(r,s)$ 即得
$I(r,s)=\Delta_x(r,s)h(r,s)=\sum_i I_i(r,s)$。$\square$

---

### 1.4.3 可分离权重与 slab 质量的积分表达

后续将反复使用如下可分离权重结构：对跨类有序对 $(r,s)\in\Omega$ 的权重写为
$$
W(r,s)=I(r,s)\cdot u(r)\cdot v(s),
$$
其中 $u:R_c\to\mathbb R_{\ge0}$ 与 $v:R_{\bar c}\to\mathbb R_{\ge0}$ 为非负的单矩形因子。

对 slab $S_i$ 定义分量权重
$$
W_i(r,s)=I_i(r,s)\,u(r)\,v(s),
$$
并定义 slab 质量
$$
M_i=\sum_{(r,s)\in\Omega} W_i(r,s).
$$
由定理 1.10 可得总权重分解
$$
\sum_{(r,s)\in\Omega}W(r,s)=\sum_{i=0}^{m-1} M_i.
$$

为了把 $M_i$ 写成可维护的一维积分，引入 slab 内的活跃集合
$$
A_c(i)=\{r\in R_c:\chi_r(i)=1\},\qquad 
A_{\bar c}(i)=\{s\in R_{\bar c}:\chi_s(i)=1\}.
$$
并在 $y$ 轴上定义加权覆盖函数
$$
F_i(y)=\sum_{r\in A_c(i)} u(r)\cdot \mathbf 1\big[y\in[L_y(r),R_y(r))\big],
$$
$$
G_i(y)=\sum_{s\in A_{\bar c}(i)} v(s)\cdot \mathbf 1\big[y\in[L_y(s),R_y(s))\big].
$$

为明确积分的有限范围，定义所有 $y$ 端点集合
$$
Y=\{L_y(r),R_y(r)\mid r\in R_c\cup R_{\bar c}\},
$$
去重排序为
$$
y_0<y_1<\cdots<y_q.
$$
由于 $F_i,G_i$ 仅可能在这些端点处发生跳变，并且在 $[y_0,y_q)$ 之外恒为 $0$，因此对任意可积表达式可将积分区间取为 $[y_0,y_q)$（与在 $\mathbb R$ 上积分完全等价）。

**引理 1.11（slab 质量的乘积积分形式）**  
对任意 slab $S_i$，
$$
M_i=\Delta x_i\cdot \int_{y_0}^{y_q} F_i(y)G_i(y)\,dy.
$$

*证明*：展开乘积并交换有限求和与积分：
$$
\int_{y_0}^{y_q}F_i(y)G_i(y)\,dy
=\sum_{r\in A_c(i)}\sum_{s\in A_{\bar c}(i)} u(r)v(s)
\int_{y_0}^{y_q}\mathbf 1[y\in r^y]\mathbf 1[y\in s^y]\,dy,
$$
其中 $r^y=[L_y(r),R_y(r))$。积分项等于交叠长度
$$
\int_{y_0}^{y_q}\mathbf 1[y\in r^y]\mathbf 1[y\in s^y]\,dy
=\mathrm{len}(r^y\cap s^y)=h(r,s).
$$
因此
$$
\int_{y_0}^{y_q}F_iG_i
=\sum_{r\in A_c(i)}\sum_{s\in A_{\bar c}(i)} u(r)v(s)\,h(r,s).
$$
两边乘以 $\Delta x_i$，并注意 $W_i(r,s)=\Delta x_i\,\chi_r(i)\chi_s(i)\,h(r,s)\,u(r)v(s)$，即得
$M_i=\sum_{(r,s)\in\Omega}W_i(r,s)$。$\square$

---

### 1.4.4 以 slab 为索引的几何子问题接口

由定理 1.10 与引理 1.11，权重 $W(r,s)=I(r,s)u(r)v(s)$ 可写成 slab 分量之和：
$$
W(r,s)=\sum_{i=0}^{m-1} W_i(r,s),
\qquad 
M_i=\sum_{(r,s)\in\Omega}W_i(r,s).
$$
当 $M_i>0$ 时定义 slab 条件分布
$$
\mathcal D^{(i)}(r,s)=\frac{W_i(r,s)}{M_i},\qquad (r,s)\in\Omega.
$$

于是采样任务可规约为一个以 slab 为索引的分量混合问题：

- 索引集合：$\mathcal B=\{0,1,\dots,m-1\}$（按 $x$ 顺序扫描）；
- 质量：$M_i=\Delta x_i\int_{y_0}^{y_q}F_i(y)G_i(y)\,dy$；
- 条件分布：$\mathcal D^{(i)}$。

只要能够在扫描过程中实现如下两个接口：
- `Mass(i)`：返回 $M_i$；
- `CondSample(i,k)`：返回 $k$ 个相互独立样本，其单次边缘分布为 $\mathcal D^{(i)}$；

就可以将 $\mathcal B$ 代入算法 1.1（Two-Pass Scan–Mix–Backfill），从而生成 $t$ 个 i.i.d. 样本，其单次边缘分布为
$$
\tilde{\mathcal D}(r,s)=\sum_{i=0}^{m-1}\frac{M_i}{\sum_j M_j}\,\mathcal D^{(i)}(r,s)
=\frac{W(r,s)}{\sum_{(u,v)\in\Omega}W(u,v)}.
$$

这一接口的关键点在于：Pass 1 只需要得到每个 slab 的质量 $M_i$；Pass 2 只需要在扫描状态下完成 slab 内条件采样。后续章节将围绕这两个接口分别构造可维护的数据结构与确定的最坏时间复杂度界。

---

# 第 2 章 IoU 权重的可采样分解：常数因子提案与正系数核展开

本章讨论如下加权采样问题：在跨集合有序对空间 $\Omega=R_c\times R_{\bar c}$ 上，以 IoU 权重
$$
W_{\mathrm{IoU}}(r,s)=\frac{I(r,s)}{A(r)+A(s)-I(r,s)}
$$
为未归一化密度，生成 $t$ 个 i.i.d.（有放回）样本。直接按 $W_{\mathrm{IoU}}$ 采样的困难在于：分母包含交面积 $I(r,s)$，使得权重难以写成适配扫描线的可分离形式。

本章给出一个可组合的规约，分为两步：

1. **常数因子提案与接受校正**：构造提案权重
   $$
   \rho(r,s)=\frac{I(r,s)}{A(r)+A(s)},
   $$
   并证明 $W_{\mathrm{IoU}}$ 与 $\rho$ 在支撑集上相差至多常数 2；因此可以从 $\rho$ 诱导的提案分布采样，再用一个接受概率把提案校正为 IoU 分布。为获得确定最坏时间，引入“有界拒绝”并在总变差距离（TV）下量化其偏差。

2. **正系数核展开**：在有限区间 $[S_{\min},S_{\max}]$ 上，用正系数指数和相对逼近 $1/x$，从而把 $\rho(r,s)$ 的分母 $1/(A(r)+A(s))$ 写成有限项 $\sum_{\ell=1}^L \omega_\ell e^{-\alpha_\ell(A(r)+A(s))}$。由此得到
   $$
   \rho(r,s)\ \approx\ \sum_{\ell=1}^L I(r,s)\,w_\ell(r)\,v_\ell(s),
   $$
   其中 $w_\ell$ 与 $v_\ell$ 为非负的单矩形因子。这样，提案采样被规约为 $L$ 个可分离几何子问题的混合。

本章只建立理论构造、可计算的系数生成方式，以及端到端的 TV 误差闭环；每个分量的几何采样器将在后续章节给出。

---

## 2.1 采样对象与 IoU 目标分布

### 2.1.1 记号：面积、交面积与支撑集

给定两类矩形集合 $R_c$ 与 $R_{\bar c}$。对任意矩形 $r$，记其面积为
$$
A(r)=\mathrm{Area}(r)>0.
$$
对任意跨类有序对 $P=(r,s)\in \Omega$，记
- $A(P)=A(r)$，
- $B(P)=A(s)$，
- $I(P)=I(r,s)=\mathrm{Area}(r\cap s)$。

定义正权重支撑集（有效采样域）
$$
J=\{(r,s)\in \Omega:\ I(r,s)>0\}.
$$
在 $J$ 上，IoU 权重定义为
$$
W_{\mathrm{IoU}}(P)=\frac{I(P)}{A(P)+B(P)-I(P)}\in(0,1].
$$

当 $J=\varnothing$ 时，所有权重为 $0$，采样任务输出空序列。以下默认 $J\neq \varnothing$。

### 2.1.2 IoU 目标分布与 TV 距离

定义总权重
$$
Z_{\mathrm{IoU}}=\sum_{Q\in J} W_{\mathrm{IoU}}(Q)>0,
$$
目标分布为
$$
\mathcal D_{\mathrm{IoU}}(P)=\frac{W_{\mathrm{IoU}}(P)}{Z_{\mathrm{IoU}}},\qquad P\in J.
$$

在有限样本空间 $J$ 上，两分布 $\mathcal P,\mathcal Q$ 的总变差距离为
$$
d_{\mathrm{TV}}(\mathcal P,\mathcal Q)=\frac12\sum_{P\in J}\big|\mathcal P(P)-\mathcal Q(P)\big|.
$$

---

## 2.2 常数因子提案 $\rho=\frac{I}{A+B}$

### 2.2.1 提案权重与范围

定义提案权重
$$
\rho(P)=\frac{I(P)}{A(P)+B(P)},\qquad P\in J.
$$

**引理 2.1（$\rho$ 的范围）**  
对任意 $P\in J$，
$$
0<\rho(P)\le \frac12.
$$

**证明**  
因 $P\in J$，有 $I(P)>0$ 且 $A(P)+B(P)>0$，故 $\rho(P)>0$。又有
$$
I(P)\le \min\{A(P),B(P)\}\le \frac{A(P)+B(P)}{2},
$$
因此 $\rho(P)\le 1/2$。证毕。

### 2.2.2 IoU 与 $\rho$ 的常数夹逼

由定义可得恒等式
$$
W_{\mathrm{IoU}}(P)
=\frac{I(P)}{A(P)+B(P)-I(P)}
=\frac{\rho(P)}{1-\rho(P)}.
$$

**推论 2.2（IoU 与 $\rho$ 的常数因子关系）**  
对任意 $P\in J$，
$$
\rho(P)\le W_{\mathrm{IoU}}(P)\le 2\rho(P).
$$

**证明**  
由引理 2.1，$1-\rho(P)\in[1/2,1)$，于是
$$
\frac{1}{1-\rho(P)}\in(1,2].
$$
代入 $W_{\mathrm{IoU}}=\rho/(1-\rho)$ 即得。证毕。

推论 2.2 表明：$2\rho$ 是 $W_{\mathrm{IoU}}$ 的输入无关常数因子包络，后续可用拒绝校正把 $\rho$ 的采样映射为 IoU 采样。

---

## 2.3 接受校正算子与 TV 稳定性

### 2.3.1 提案分布与接受概率

令 $\mathcal D_\rho$ 表示由 $\rho$ 归一化得到的提案分布：
$$
\mathcal D_\rho(P)=\frac{\rho(P)}{\sum_{Q\in J}\rho(Q)},\qquad P\in J.
$$

用 $2\rho$ 作包络，定义接受概率
$$
a(P)=\frac{W_{\mathrm{IoU}}(P)}{2\rho(P)}
=\frac{1}{2(1-\rho(P))}.
$$
由引理 2.1 得到
$$
a(P)\in\left[\frac12,1\right]\quad\text{对所有 }P\in J.
$$

### 2.3.2 重加权算子

对任意分布 $\mathcal P$（定义在 $J$ 上）与函数 $a:J\to(0,\infty)$，定义重加权算子
$$
\mathcal T_a(\mathcal P)(P)=\frac{\mathcal P(P)\,a(P)}{\sum_{Q\in J}\mathcal P(Q)\,a(Q)}.
$$

**定理 2.3（从 $\rho$ 校正到 IoU）**  
有
$$
\mathcal T_a(\mathcal D_\rho)=\mathcal D_{\mathrm{IoU}}.
$$

**证明**  
对任意 $P\in J$，
$$
\mathcal T_a(\mathcal D_\rho)(P)\propto \mathcal D_\rho(P)a(P)
\propto \rho(P)\cdot \frac{W_{\mathrm{IoU}}(P)}{2\rho(P)}
\propto W_{\mathrm{IoU}}(P),
$$
归一化后即为 $\mathcal D_{\mathrm{IoU}}$。证毕。

### 2.3.3 $\mathcal T_a$ 在 TV 下的 Lipschitz 性

后续会用近似提案分布替代 $\mathcal D_\rho$。为控制误差传播，需要 $\mathcal T_a$ 的稳定性。

**引理 2.4（重加权算子的 TV-Lipschitz 性）**  
设 $a(P)\in[a_{\min},a_{\max}]$ 对所有 $P\in J$ 成立，且 $a_{\min}>0$。则对任意两分布 $\mathcal P,\mathcal Q$（定义在 $J$ 上）有
$$
d_{\mathrm{TV}}\big(\mathcal T_a(\mathcal P),\mathcal T_a(\mathcal Q)\big)
\le \frac{2a_{\max}}{a_{\min}}\cdot d_{\mathrm{TV}}(\mathcal P,\mathcal Q).
$$
特别地，当 $a(P)\in[1/2,1]$ 时，
$$
d_{\mathrm{TV}}\big(\mathcal T_a(\mathcal P),\mathcal T_a(\mathcal Q)\big)
\le 4\, d_{\mathrm{TV}}(\mathcal P,\mathcal Q).
$$

**证明**  
令未归一化向量
$$
\mu_{\mathcal P}(P)=\mathcal P(P)a(P),\qquad \mu_{\mathcal Q}(P)=\mathcal Q(P)a(P),
$$
归一化常数
$$
S_{\mathcal P}=\sum_{P\in J}\mu_{\mathcal P}(P)=\mathbb E_{\mathcal P}[a],\qquad
S_{\mathcal Q}=\mathbb E_{\mathcal Q}[a].
$$
则 $\mathcal T_a(\mathcal P)=\mu_{\mathcal P}/S_{\mathcal P}$，$\mathcal T_a(\mathcal Q)=\mu_{\mathcal Q}/S_{\mathcal Q}$。

对任意非负向量 $u,v$ 且 $\|u\|_1,\|v\|_1>0$，有标准不等式
$$
\left\|\frac{u}{\|u\|_1}-\frac{v}{\|v\|_1}\right\|_1
\le \frac{2\|u-v\|_1}{\min\{\|u\|_1,\|v\|_1\}}.
$$
应用于 $u=\mu_{\mathcal P},v=\mu_{\mathcal Q}$，并用 $S_{\mathcal P},S_{\mathcal Q}\ge a_{\min}$ 得
$$
\|\mathcal T_a(\mathcal P)-\mathcal T_a(\mathcal Q)\|_1
\le \frac{2\|\mu_{\mathcal P}-\mu_{\mathcal Q}\|_1}{a_{\min}}.
$$
又因为 $a(P)\le a_{\max}$，
$$
\|\mu_{\mathcal P}-\mu_{\mathcal Q}\|_1
=\sum_{P\in J} a(P)\,|\mathcal P(P)-\mathcal Q(P)|
\le a_{\max}\sum_{P\in J}|\mathcal P(P)-\mathcal Q(P)|
=2a_{\max}\,d_{\mathrm{TV}}(\mathcal P,\mathcal Q).
$$
合并并除以 2 即得结论。证毕。

---

## 2.4 有界拒绝：确定最坏时间与 TV 误差

经典拒绝采样“直到接受”为止能得到精确的 $\mathcal T_a(\mathcal P)$，但其运行时间是随机变量。本节引入一个确定最坏时间的版本：每个输出最多尝试 $M$ 次，若全部拒绝则走固定失败分支。

### 2.4.1 有界拒绝算子

给定提案分布 $\mathcal P$ 与接受函数 $a(\cdot)\in(0,1]$，定义 $\mathcal T_{a,M}(\mathcal P)$ 为如下随机过程的输出分布：

- 独立采样 $P_1,\dots,P_M\sim\mathcal P$；
- 独立采样 $U_1,\dots,U_M\sim\mathrm{Unif}(0,1)$；
- 输出 $P_\tau$，其中
  - 若存在最小 $\tau$ 使 $U_\tau\le a(P_\tau)$，则输出 $P_\tau$；
  - 否则输出 $P_M$（固定失败分支）。

该定义保证：每个输出所消耗的提案数严格不超过 $M$。

### 2.4.2 有界拒绝的 TV 误差界

**引理 2.5（有界拒绝误差）**  
若存在 $a_{\min}>0$ 使 $a(P)\ge a_{\min}$ 对所有 $P$ 成立，则对任意提案分布 $\mathcal P$，
$$
d_{\mathrm{TV}}\big(\mathcal T_{a,M}(\mathcal P),\,\mathcal T_a(\mathcal P)\big)\le (1-a_{\min})^M.
$$
在本章 IoU 校正中 $a_{\min}=1/2$，因此
$$
d_{\mathrm{TV}}\big(\mathcal T_{a,M}(\mathcal P),\,\mathcal T_a(\mathcal P)\big)\le 2^{-M}.
$$

**证明**  
令事件 $\mathsf{Fail}$ 表示前 $M$ 次尝试均拒绝。无界拒绝（直到接受）诱导分布为 $\mathcal T_a(\mathcal P)$。有界过程输出可写成凸组合
$$
\mathcal T_{a,M}(\mathcal P)
=\Pr[\neg \mathsf{Fail}]\cdot \mathcal T_a(\mathcal P)
+\Pr[\mathsf{Fail}]\cdot \mathcal Q,
$$
其中 $\mathcal Q$ 为失败分支（输出 $P_M$）诱导的分布。于是
$$
d_{\mathrm{TV}}(\mathcal T_{a,M}(\mathcal P),\mathcal T_a(\mathcal P))
\le \Pr[\mathsf{Fail}]\cdot d_{\mathrm{TV}}(\mathcal Q,\mathcal T_a(\mathcal P))
\le \Pr[\mathsf{Fail}].
$$
每次尝试的接受概率至少为 $a_{\min}$，且各次尝试独立，因此
$$
\Pr[\mathsf{Fail}]\le (1-a_{\min})^M.
$$
证毕。

### 2.4.3 外层闭环形式

对任意近似提案分布 $\tilde{\mathcal D}_\rho$，理想输出为 $\mathcal D_{\mathrm{IoU}}=\mathcal T_a(\mathcal D_\rho)$，实际输出为 $\mathcal T_{a,M}(\tilde{\mathcal D}_\rho)$。由三角不等式、引理 2.5 与引理 2.4（在 $a\in[1/2,1]$ 下）得到
$$
d_{\mathrm{TV}}\big(\mathcal T_{a,M}(\tilde{\mathcal D}_\rho),\,\mathcal D_{\mathrm{IoU}}\big)
\le 2^{-M}+4\,d_{\mathrm{TV}}(\tilde{\mathcal D}_\rho,\mathcal D_\rho).
$$

因此端到端误差被规约为两项：有界拒绝项 $2^{-M}$ 与提案误差 $d_{\mathrm{TV}}(\tilde{\mathcal D}_\rho,\mathcal D_\rho)$。

---

## 2.5 正系数指数和逼近 $1/x$

提案权重为
$$
\rho(P)=\frac{I(P)}{A(P)+B(P)}.
$$
若能在 $x=A(P)+B(P)$ 的取值范围内相对逼近 $1/x$，则 $\rho(P)$ 可写为有限项非负可分离结构之和。本节给出一个构造性结论：在区间 $[S_{\min},S_{\max}]$ 上，$1/x$ 可以用正系数指数和相对逼近，且项数为多对数级。

### 2.5.1 面积和的区间

定义两类矩形面积的范围：
$$
A_{\min}^c=\min_{r\in R_c}A(r),\quad A_{\max}^c=\max_{r\in R_c}A(r),
$$
$$
A_{\min}^{\bar c}=\min_{s\in R_{\bar c}}A(s),\quad A_{\max}^{\bar c}=\max_{s\in R_{\bar c}}A(s).
$$
矩形非退化给出 $A_{\min}^c>0$ 与 $A_{\min}^{\bar c}>0$。令
$$
S_{\min}=A_{\min}^c+A_{\min}^{\bar c},\qquad
S_{\max}=A_{\max}^c+A_{\max}^{\bar c},\qquad
\kappa=\frac{S_{\max}}{S_{\min}}\ge 1.
$$
则对任意 $P=(r,s)\in J$，有 $A(P)+B(P)\in[S_{\min},S_{\max}]$。

### 2.5.2 正系数指数和相对逼近

**定理 2.6（$1/x$ 的正指数和相对逼近）**  
给定 $\delta\in(0,1/2)$。存在 $L$ 组正参数 $\{(\omega_\ell,\alpha_\ell)\}_{\ell=1}^L$，满足 $\omega_\ell>0,\alpha_\ell>0$，使得对所有 $x\in[S_{\min},S_{\max}]$，
$$
\left|\frac{1}{x}-\sum_{\ell=1}^L \omega_\ell e^{-\alpha_\ell x}\right|\le \frac{\delta}{x}.
$$
并且可以取
$$
L=O\!\left((\log\kappa+\log(1/\delta))\cdot \log(1/\delta)\right).
$$

下面给出一个自含构造与证明。

---

#### 证明：拉普拉斯表示 + 对数变量梯形求积 + 双侧截断

**缩放**  
令 $y=x/S_{\min}$，则 $y\in[1,\kappa]$ 且
$$
\frac{1}{x}=\frac{1}{S_{\min}}\cdot \frac{1}{y}.
$$
因此只需在 $y\in[1,\kappa]$ 上构造对 $1/y$ 的相对逼近。

**(1) 拉普拉斯表示与变量替换**  
对任意 $y>0$，
$$
\frac{1}{y}=\int_0^\infty e^{-yt}\,dt.
$$
令 $t=e^u$，则 $dt=e^u\,du$，得到
$$
\frac{1}{y}=\int_{-\infty}^{+\infty} f_y(u)\,du,\qquad
f_y(u)=e^u e^{-y e^u}.
$$
注意 $f_y(u)\ge 0$，并且该积分收敛。

**(2) 无穷梯形求积与离散化误差**  
取步长 $h>0$，定义无穷梯形近似
$$
I(y)=\int_{-\infty}^{+\infty} f_y(u)\,du,\qquad
I_h(y)=h\sum_{k\in\mathbb Z} f_y(kh).
$$

**引理 2.6.1（梯形求积误差界）**  
取任意 $d\in(0,\pi/2)$，令 $q=e^{-2\pi d/h}$。则对任意 $y\ge 1$，
$$
|I(y)-I_h(y)|\le \frac{2}{y\cos d}\cdot \frac{q}{1-q}.
$$

**证明**  
记傅里叶变换 $\widehat f_y(\xi)=\int_{-\infty}^{+\infty} f_y(u)e^{-i\xi u}\,du$。由泊松求和公式，
$$
h\sum_{k\in\mathbb Z} f_y(kh)=\sum_{n\in\mathbb Z}\widehat f_y\!\left(\frac{2\pi n}{h}\right).
$$
其中 $n=0$ 项为 $\widehat f_y(0)=I(y)$，因此
$$
I_h(y)-I(y)=\sum_{n\neq 0}\widehat f_y\!\left(\frac{2\pi n}{h}\right),
\qquad
|I_h(y)-I(y)|
\le \sum_{n\neq 0}\left|\widehat f_y\!\left(\frac{2\pi n}{h}\right)\right|.
$$

由于 $f_y$ 是整函数，可将积分路径平移到 $u\mapsto u\pm id$。取 $s=-\mathrm{sign}(\xi)\,d$，则
$$
\widehat f_y(\xi)=\int_{-\infty}^{+\infty} f_y(u+is)e^{-i\xi(u+is)}\,du
=e^{-|\xi|d}\int_{-\infty}^{+\infty} f_y(u+is)e^{-i\xi u}\,du,
$$
从而
$$
|\widehat f_y(\xi)|\le e^{-|\xi|d}\int_{-\infty}^{+\infty}|f_y(u+is)|\,du.
$$
当 $|s|=d<\pi/2$ 时 $\cos d>0$，且
$$
|f_y(u+is)|=\left|e^{u+is}\right|\cdot \left|e^{-y e^{u+is}}\right|
= e^u\cdot e^{-y e^u\cos d}.
$$
于是
$$
\int_{-\infty}^{+\infty}|f_y(u+is)|\,du
=\int_{-\infty}^{+\infty} e^u e^{-y e^u\cos d}\,du.
$$
令 $t=y\cos d\cdot e^u$，则 $dt=y\cos d\cdot e^u\,du$，上式变为
$$
\int_0^\infty \frac{1}{y\cos d}e^{-t}\,dt=\frac{1}{y\cos d}.
$$
因此
$$
|\widehat f_y(\xi)|\le \frac{1}{y\cos d}\,e^{-|\xi|d}.
$$
代入 $\xi=2\pi n/h$ 并求和得
$$
|I_h(y)-I(y)|
\le \frac{1}{y\cos d}\sum_{n\neq 0}e^{-2\pi|n|d/h}
= \frac{2}{y\cos d}\sum_{n\ge 1}q^n
= \frac{2}{y\cos d}\cdot \frac{q}{1-q}.
$$
证毕。

固定 $d=\pi/3$（此时 $\cos d=1/2$），并选取
$$
h=\min\left\{1,\ \frac{2\pi d}{\ln(24/\delta)}\right\}
=\min\left\{1,\ \frac{2\pi^2}{3\ln(24/\delta)}\right\}.
$$
则 $q=e^{-2\pi d/h}\le \delta/24$，且在 $\delta<1/2$ 时有 $q\le 1/48$，因此 $q/(1-q)\le 2q\le \delta/12$。代入引理 2.6.1 得对所有 $y\in[1,\kappa]$，
$$
|I(y)-I_h(y)|
\le \frac{2}{y(1/2)}\cdot \frac{\delta}{12}
=\frac{\delta}{3y}.
$$

**(3) 双侧截断与尾项误差**  
定义有限截断版本
$$
I_{h,K_-,K_+}(y)=h\sum_{k=-K_-}^{K_+} f_y(kh)
=h\sum_{k=-K_-}^{K_+} e^{kh}e^{-y e^{kh}}.
$$

- 左尾：对 $k\le -K_--1$，有 $e^{-y e^{kh}}\le 1$，故 $f_y(kh)\le e^{kh}$，于是
  $$
  h\sum_{k=-\infty}^{-K_--1}f_y(kh)
  \le h\sum_{k=-\infty}^{-K_--1}e^{kh}
  =\frac{h e^{-(K_-+1)h}}{1-e^{-h}}
  \le 2e^{-K_-h},
  $$
  其中用到 $h\le 1$ 时 $h/(1-e^{-h})<2$。

- 右尾：当 $y\ge 1$ 且 $u\ge 0$ 时，$f_y'(u)=f_y(u)(1-y e^u)\le 0$，因此 $f_y$ 在 $[0,\infty)$ 上单调不增。于是
  $$
  h\sum_{k=K_++1}^{\infty} f_y(kh)\le \int_{K_+h}^{\infty} f_y(u)\,du.
  $$
  对积分令 $t=y e^u$，得到
  $$
  \int_{K_+h}^{\infty} f_y(u)\,du
  =\frac{1}{y}\int_{y e^{K_+h}}^\infty e^{-t}\,dt
  \le \frac{1}{y}e^{-y e^{K_+h}}
  \le \frac{1}{y}e^{-e^{K_+h}}.
  $$

选择
$$
K_-=\left\lceil \frac{1}{h}\ln\frac{6\kappa}{\delta}\right\rceil,\qquad
K_+=\left\lceil \frac{1}{h}\ln\ln\frac{3}{\delta}\right\rceil.
$$
则对所有 $y\in[1,\kappa]$，
- 左尾误差 $\le 2e^{-K_-h}\le \delta/(3\kappa)\le \delta/(3y)$；
- 右尾误差 $\le \frac{1}{y}e^{-e^{K_+h}}\le \delta/(3y)$。

**(4) 合并误差并写成正指数和**  
由三角不等式，对任意 $y\in[1,\kappa]$，
$$
\left|\frac{1}{y}-I_{h,K_-,K_+}(y)\right|
\le |I(y)-I_h(y)| + |I_h(y)-I_{h,K_-,K_+}(y)|
\le \frac{\delta}{3y}+\frac{\delta}{3y}+\frac{\delta}{3y}
=\frac{\delta}{y}.
$$
而
$$
I_{h,K_-,K_+}(y)
=h\sum_{k=-K_-}^{K_+} e^{kh}e^{-y e^{kh}}
=\sum_{k=-K_-}^{K_+} \underbrace{(h e^{kh})}_{>0}\cdot e^{-\underbrace{(e^{kh})}_{>0}\,y}.
$$
令每个 $\ell$ 对应一个 $k$，定义
$$
\omega_\ell=h e^{kh}>0,\qquad \alpha_\ell=e^{kh}>0,
$$
即可得到对 $y$ 的正系数指数和逼近。

最后将 $y=x/S_{\min}$ 代回，并令
$$
\omega_\ell \leftarrow \frac{\omega_\ell}{S_{\min}},\qquad
\alpha_\ell \leftarrow \frac{\alpha_\ell}{S_{\min}},
$$
得到对 $x\in[S_{\min},S_{\max}]$ 的结论。

项数方面，$L=K_-+K_++1$。由于 $h=\Theta(\min\{1,1/\ln(1/\delta)\})$，有
$$
K_-=O\!\left(\frac{\ln(\kappa/\delta)}{h}\right),\qquad
K_+=O\!\left(\frac{\ln\ln(1/\delta)}{h}\right),
$$
从而
$$
L=O\!\left((\log\kappa+\log(1/\delta))\cdot \log(1/\delta)\right).
$$
定理得证。

---

### 2.5.3 系数生成算法

下面把上述构造写成显式流程。所有系数严格为正，便于后续作为“混合分量”使用。

```text
Algorithm 2.1  Positive Exponential-Sum Approximation for 1/x on [S_min, S_max]

Input:
  S_min > 0,  S_max ≥ S_min
  δ ∈ (0, 1/2)

Output:
  {(ω_ℓ, α_ℓ)}_{ℓ=1}^L  with ω_ℓ>0, α_ℓ>0 such that
  | 1/x - Σ_{ℓ=1}^L ω_ℓ e^{-α_ℓ x} | ≤ δ/x  for all x ∈ [S_min, S_max]

κ ← S_max / S_min
d ← π/3
h ← min{ 1, 2πd / ln(24/δ) }            // = min{1, 2π^2/(3 ln(24/δ))}

K_- ← ceil( (1/h) * ln(6κ/δ) )
K_+ ← ceil( (1/h) * ln ln(3/δ) )

for k = -K_- .. K_+:
    α ← exp(kh) / S_min
    ω ← (h * exp(kh)) / S_min
    append (ω, α) to the list

return list
```

**实现备注（数值稳定性）**  
当 $A(r)$ 或 $\alpha_\ell$ 较大时，$e^{-\alpha_\ell A(r)}$ 可能下溢；实际实现通常在对数域存储 $\log w_\ell(r)=-\alpha_\ell A(r)$ 与 $\log v_\ell(s)=\log\omega_\ell-\alpha_\ell A(s)$，并在需要比较/采样时再做受控的指数或对数差运算。

---

## 2.6 从核展开到可分离子权重混合

### 2.6.1 近似提案的未归一化权重

对任意 $P=(r,s)\in J$，令 $x=A(P)+B(P)$。由定理 2.6，
$$
\frac{1-\delta}{x}\le \sum_{\ell=1}^L \omega_\ell e^{-\alpha_\ell x}\le \frac{1+\delta}{x}.
$$
两边乘以 $I(P)$ 得到点态夹逼
$$
(1-\delta)\rho(P)\le \tilde w_\rho(P)\le (1+\delta)\rho(P),
$$
其中定义近似提案的未归一化权重
$$
\tilde w_\rho(P)=I(P)\sum_{\ell=1}^L \omega_\ell e^{-\alpha_\ell(A(P)+B(P))}.
$$

### 2.6.2 可分离子权重

将指数项分离到单矩形因子上。对每个 $\ell$ 定义
$$
w_\ell(r)=e^{-\alpha_\ell A(r)},\qquad
v_\ell(s)=\omega_\ell e^{-\alpha_\ell A(s)}.
$$
并定义子权重
$$
W_\ell(r,s)=I(r,s)\cdot w_\ell(r)\cdot v_\ell(s)\ge 0.
$$
则
$$
\tilde w_\rho(r,s)=\sum_{\ell=1}^L W_\ell(r,s).
$$

令
$$
Z_\ell=\sum_{(r,s)\in\Omega} W_\ell(r,s),\qquad
Z=\sum_{\ell=1}^L Z_\ell.
$$
当 $Z_\ell>0$ 时定义分量条件分布
$$
\mathcal D_\ell(r,s)=\frac{W_\ell(r,s)}{Z_\ell}.
$$
当 $Z>0$ 时定义混合系数
$$
\pi_\ell=\frac{Z_\ell}{Z},\qquad \sum_{\ell=1}^L \pi_\ell=1,
$$
从而近似提案分布写成分量混合：
$$
\tilde{\mathcal D}_\rho(P)
=\frac{\tilde w_\rho(P)}{\sum_{Q\in J}\tilde w_\rho(Q)}
=\sum_{\ell=1}^L \pi_\ell\,\mathcal D_\ell(P).
$$

---

## 2.7 权重乘性误差、混合误差与 TV 上界

### 2.7.1 乘性权重误差 $\Rightarrow$ TV

**引理 2.7（乘性权重误差到 TV）**  
设 $w,\tilde w$ 为定义在有限集合 $\Omega$ 上的非负权重，诱导分布
$$
\mathcal D(x)=\frac{w(x)}{\sum_{y\in\Omega} w(y)},\qquad
\tilde{\mathcal D}(x)=\frac{\tilde w(x)}{\sum_{y\in\Omega} \tilde w(y)}.
$$
若存在 $\delta\in(0,1)$ 使对所有 $x\in\Omega$，
$$
(1-\delta)w(x)\le \tilde w(x)\le (1+\delta)w(x),
$$
则
$$
d_{\mathrm{TV}}(\mathcal D,\tilde{\mathcal D})\le \frac{\delta}{1-\delta}.
$$

**证明**  
记 $W=\sum_x w(x)$，$\tilde W=\sum_x\tilde w(x)$。由夹逼得
$$
(1-\delta)W\le \tilde W\le (1+\delta)W.
$$
对任意 $x$，
$$
\left|\tilde{\mathcal D}(x)-\mathcal D(x)\right|
=\left|\frac{\tilde w(x)}{\tilde W}-\frac{w(x)}{W}\right|
\le \frac{|\tilde w(x)-w(x)|}{\tilde W}
+w(x)\left|\frac{1}{\tilde W}-\frac{1}{W}\right|.
$$
对 $x$ 求和并分别估计两项。第一项：
$$
\sum_x \frac{|\tilde w(x)-w(x)|}{\tilde W}
\le \frac{1}{\tilde W}\sum_x \delta w(x)=\frac{\delta W}{\tilde W}
\le \frac{\delta}{1-\delta}.
$$
第二项：
$$
\sum_x w(x)\left|\frac{1}{\tilde W}-\frac{1}{W}\right|
= W\cdot \frac{|W-\tilde W|}{W\tilde W}
=\frac{|W-\tilde W|}{\tilde W}.
$$
又
$$
|W-\tilde W|\le \sum_x |\tilde w(x)-w(x)|\le \delta W,
$$
因此第二项也不超过 $\delta/(1-\delta)$。于是
$$
\|\tilde{\mathcal D}-\mathcal D\|_1\le \frac{2\delta}{1-\delta}
\quad\Rightarrow\quad
d_{\mathrm{TV}}(\tilde{\mathcal D},\mathcal D)\le \frac{\delta}{1-\delta}.
$$
证毕。

将引理 2.7 应用于 $\Omega=J$、$w=\rho$、$\tilde w=\tilde w_\rho$，并使用 2.6.1 的点态夹逼，得到
$$
d_{\mathrm{TV}}(\tilde{\mathcal D}_\rho,\mathcal D_\rho)\le \frac{\delta}{1-\delta}.
$$

### 2.7.2 混合系数扰动 $\Rightarrow$ TV

实际装配中，$\pi_\ell$ 可能用近似值 $\hat\pi_\ell$ 替代（例如来自近似质量估计）。定义
$$
\varepsilon_{\mathrm{mix}}=\frac12\sum_{\ell=1}^L |\pi_\ell-\hat\pi_\ell|.
$$

**引理 2.8（混合系数扰动的 TV 上界）**  
给定同一组分量分布 $\{\mathcal D_\ell\}_{\ell=1}^L$，令
$$
\mathcal M=\sum_{\ell=1}^L \pi_\ell \mathcal D_\ell,\qquad
\hat{\mathcal M}=\sum_{\ell=1}^L \hat\pi_\ell \mathcal D_\ell.
$$
则
$$
d_{\mathrm{TV}}(\mathcal M,\hat{\mathcal M})\le \varepsilon_{\mathrm{mix}}.
$$

**证明**  
对任意样本点 $x$，
$$
\mathcal M(x)-\hat{\mathcal M}(x)=\sum_{\ell=1}^L (\pi_\ell-\hat\pi_\ell)\mathcal D_\ell(x).
$$
因此
$$
\|\mathcal M-\hat{\mathcal M}\|_1
\le \sum_{\ell=1}^L |\pi_\ell-\hat\pi_\ell|\sum_x \mathcal D_\ell(x)
=\sum_{\ell=1}^L |\pi_\ell-\hat\pi_\ell|.
$$
两边除以 2 得结论。证毕。

---

## 2.8 本章收束：IoU 采样规约定理

本节把 2.2–2.7 的结论汇总为一个可直接调用的规约：只要能在 $L$ 个分量 $\mathcal D_\ell$ 下生成提案（精确或近似），并配合外层有界拒绝校正，就可以得到 TV 可控、最坏时间可控的 IoU 采样器。

### 2.8.1 规约输入：分量混合提案

设实际提案分布为
$$
\widehat{\mathcal D}_\rho=\sum_{\ell=1}^L \hat\pi_\ell\,\mathcal D_\ell,
$$
其中 $\hat\pi$ 为某组混合系数（可为精确或近似）。对提案样本 $P$ 计算
$$
\rho(P)=\frac{I(P)}{A(P)+B(P)},\qquad
a(P)=\frac{1}{2(1-\rho(P))}.
$$

### 2.8.2 规约输出：有界拒绝后的 IoU 样本

生成一个 IoU 输出样本的过程如下（抽象接口形式）：

1. 独立生成提案 $P_1,\dots,P_M\sim \widehat{\mathcal D}_\rho$；
2. 独立生成 $U_1,\dots,U_M\sim\mathrm{Unif}(0,1)$；
3. 输出第一个满足 $U_m\le a(P_m)$ 的提案；若不存在则输出 $P_M$。

### 2.8.3 端到端 TV 误差闭环

**定理 2.9（规约：$\varepsilon$-TV 近似 IoU 采样）**  
固定 $\delta\in(0,1/2)$，按算法 2.1 生成正参数 $\{(\omega_\ell,\alpha_\ell)\}_{\ell=1}^L$ 并构造分量 $\{\mathcal D_\ell\}$ 与精确混合系数 $\{\pi_\ell\}$。令
$$
\varepsilon_{\mathrm{mix}}=\frac12\sum_{\ell=1}^L |\pi_\ell-\hat\pi_\ell|.
$$
使用提案分布 $\widehat{\mathcal D}_\rho=\sum_{\ell=1}^L \hat\pi_\ell \mathcal D_\ell$，并对每个输出位置执行有界拒绝（最多 $M$ 次，失败分支输出第 $M$ 次提案），接受概率为
$$
a(P)=\frac{1}{2(1-\rho(P))},\qquad \rho(P)=\frac{I(P)}{A(P)+B(P)}.
$$
则单次输出边缘分布 $\tilde{\mathcal D}$ 满足
$$
d_{\mathrm{TV}}(\tilde{\mathcal D},\mathcal D_{\mathrm{IoU}})
\le 2^{-M}+4\left(\frac{\delta}{1-\delta}+\varepsilon_{\mathrm{mix}}\right).
$$
此外，若不同输出位置的所有随机性相互独立，则输出序列为 i.i.d.（有放回）。

**证明**  

- 由 2.6.1 的点态夹逼与引理 2.7，
  $$
  d_{\mathrm{TV}}(\tilde{\mathcal D}_\rho,\mathcal D_\rho)\le \frac{\delta}{1-\delta},
  $$
  其中 $\tilde{\mathcal D}_\rho=\sum_\ell \pi_\ell\mathcal D_\ell$ 为使用精确混合系数的近似提案分布。

- 由引理 2.8（混合系数扰动），
  $$
  d_{\mathrm{TV}}(\widehat{\mathcal D}_\rho,\tilde{\mathcal D}_\rho)\le \varepsilon_{\mathrm{mix}}.
  $$
  因此三角不等式给出
  $$
  d_{\mathrm{TV}}(\widehat{\mathcal D}_\rho,\mathcal D_\rho)
  \le \frac{\delta}{1-\delta}+\varepsilon_{\mathrm{mix}}.
  $$

- 记无界校正分布为 $\mathcal T_a(\cdot)$，有界校正为 $\mathcal T_{a,M}(\cdot)$。由引理 2.5（此处 $a_{\min}=1/2$）与引理 2.4（Lipschitz 常数 4）并用三角不等式，
  $$
  \begin{aligned}
  d_{\mathrm{TV}}(\tilde{\mathcal D},\mathcal D_{\mathrm{IoU}})
  &=d_{\mathrm{TV}}\big(\mathcal T_{a,M}(\widehat{\mathcal D}_\rho),\,\mathcal T_a(\mathcal D_\rho)\big)\\
  &\le d_{\mathrm{TV}}\big(\mathcal T_{a,M}(\widehat{\mathcal D}_\rho),\,\mathcal T_a(\widehat{\mathcal D}_\rho)\big)
     +d_{\mathrm{TV}}\big(\mathcal T_a(\widehat{\mathcal D}_\rho),\,\mathcal T_a(\mathcal D_\rho)\big)\\
  &\le 2^{-M}+4\,d_{\mathrm{TV}}(\widehat{\mathcal D}_\rho,\mathcal D_\rho)\\
  &\le 2^{-M}+4\left(\frac{\delta}{1-\delta}+\varepsilon_{\mathrm{mix}}\right).
  \end{aligned}
  $$

i.i.d. 性质来自：每个输出位置使用独立的提案样本与独立的均匀随机数，且接受/失败分支是这些随机变量的确定函数。证毕。

---

## 2.9 参数选择建议

给定目标误差 $\varepsilon\in(0,1)$，一组直接可用的分配为
$$
\delta\le \frac{\varepsilon}{16},\qquad
\varepsilon_{\mathrm{mix}}\le \frac{\varepsilon}{16},\qquad
M=\left\lceil \log_2\frac{4}{\varepsilon}\right\rceil.
$$
则由定理 2.9 得
$$
d_{\mathrm{TV}}(\tilde{\mathcal D},\mathcal D_{\mathrm{IoU}})\le \varepsilon.
$$

---

## 2.10 小结与后续接口

本章完成两项规约：

- 通过 $\rho=I/(A+B)$ 与接受概率 $a(P)=1/(2(1-\rho(P)))$，把 IoU 采样化为“提案采样 + 有界拒绝校正”，并给出 $2^{-M}$ 级别的确定最坏时间误差项以及 TV 传播的常数 4。

- 通过正系数指数和相对逼近 $1/x$，把提案权重写成 $L$ 个非负可分离子权重 $W_\ell(r,s)=I(r,s)w_\ell(r)v_\ell(s)$ 的和，从而把提案采样化为 $L$ 个几何子问题的混合。

后续章节将分别实现：
- 固定 $\ell$ 的几何分量采样器（输入为 $W_\ell$ 的可分离结构）；
- 多分量混合与外层有界拒绝的端到端装配，并在同一 TV 账本下汇总所有误差项。

---

# 第 3 章 动态加权 stabbing 采样器：确定性桶更新与 TV-近似桶内抽样

本章给出一个一维动态采样原语：维护带权半开区间集合，在查询点 $y$ 处从所有覆盖 $y$ 的活跃区间中按权重比例抽样。该原语将在后续二维几何条件采样中被直接调用：在固定 slab 内先按密度采样 $y$，再对红/蓝两侧各做一次 stabbing 抽样，从而把二维权重采样降解到两个一维模块。

本章关注的设计目标同时覆盖三条“硬约束”：

- **确定性最坏更新**：每次区间插删最坏 $O(\log N)$；其中每个 canonical 归属（membership）在桶内的插删必须是最坏 $O(1)$；
- **确定性最坏查询**：每次 `Sample(y)` 的步数必须存在确定上界（不允许“直到接受”为止的随机循环）；
- **TV 误差可组合**：桶内允许近似，但误差以总变差距离（TV）显式量化，并在 “Path 混合骨架” 下不被放大，便于在整书的误差账本中直接加总。

本章的核心结构分为两层：

1. **段树分桶骨架**：把任意查询点 $y$ 的覆盖集合 $\mathcal I(y)$ 表示为根到叶路径上的桶集合的不交并，并据此把 stabbing 采样写成“先选桶、再桶内采样”的混合；
2. **桶内实现**：用 *dyadic 包络 + 平均阈值整层裁剪 + 单桶窗口扫描 + 小 alias + 有界拒绝* 在一个桶内完成近似加权抽样，同时保证更新最坏 $O(1)$、查询最坏有界、并给出桶内 TV 上界 $\varepsilon_{\mathrm{clip}}+2^{-M_{\mathrm{bucket}}}$。

---

## 3.1 抽象问题：动态加权 stabbing

设 $Y=\{y_0<y_1<\cdots<y_q\}$ 为离散端点集合（由端点压缩得到），基本段为 $I_j=[y_j,y_{j+1})$。本章所有区间均采用半开语义 $[a,b)$（$a,b\in Y$ 且 $a<b$），以避免端点处的双计数。

### 定义 3.1（动态加权 stabbing）

维护动态活跃区间集合 $\mathcal I$。每个区间 $I\in\mathcal I$ 具有静态权重 $w(I)>0$（在本章的抽象原语中视为固定）。支持操作：

- `Insert(I)`：将区间 $I$ 置为活跃；
- `Delete(I)`：将区间 $I$ 从活跃集中删除；
- `Sample(y)`：给定查询点 $y\in[y_0,y_q)$，从覆盖 $y$ 的活跃区间集合
  $$
  \mathcal I(y)=\{I\in\mathcal I:\ y\in I\}
  $$
  中按权重比例抽样并返回一个区间；若 $\mathcal I(y)=\varnothing$，则返回 $\bot$。

为统一表述，将输出域扩展为
$$
\Omega_y=\mathcal I(y)\cup\{\bot\}.
$$
定义目标分布 $\mathcal P_y$ 为：
- 若 $\mathcal I(y)\neq\varnothing$，则对 $I\in\mathcal I(y)$
  $$
  \mathcal P_y(I)=\frac{w(I)}{\sum_{J\in\mathcal I(y)}w(J)},\qquad \mathcal P_y(\bot)=0;
  $$
- 若 $\mathcal I(y)=\varnothing$，则 $\mathcal P_y(\bot)=1$。

### 定义 3.2（TV-近似 stabbing）

给定误差参数 $\varepsilon_{\mathrm{stab}}\in(0,1)$。称一个 stabbing 采样器为 $\varepsilon_{\mathrm{stab}}$-近似，若对任意 $y$，其输出分布 $\widetilde{\mathcal P}_y$（定义在同一输出域 $\Omega_y$ 上）满足
$$
d_{\mathrm{TV}}(\widetilde{\mathcal P}_y,\mathcal P_y)\le \varepsilon_{\mathrm{stab}}.
$$

---

## 3.2 段树分桶骨架：canonical buckets 与 Path 不交并

本节给出 stabbing 采样的统一混合骨架：把覆盖集合写成路径桶的不交并，使得整体采样可分解为“路径选桶 + 桶内条件采样”。

### 3.2.1 段树与规范覆盖

在线段 $[y_0,y_q)$ 上建立一棵平衡二叉段树 $\mathcal T_y$。每个节点 $v$ 对应区间 $\mathrm{seg}(v)=[y_L,y_R)$，其左右子节点对应将该区间二分的子区间；叶子节点对应基本段 $I_j=[y_j,y_{j+1})$。

对任意区间 $I=[a,b)$，定义其 **规范覆盖（canonical cover）** $\mathcal C(I)$ 为段树中的一组节点集合，使得
$$
I=\biguplus_{v\in\mathcal C(I)}\mathrm{seg}(v),
\qquad
|\mathcal C(I)|\le 2\lceil \log_2 q\rceil=O(\log N).
$$
（$\mathcal C(I)$ 可由标准递归分解得到。）

对每个节点 $v$ 维护一个桶 $B(v)$，存放所有活跃区间 $I$ 中满足 $v\in\mathcal C(I)$ 的那些区间。

### 引理 3.3（Path 唯一归属与不交并）

对任意查询点 $y$，令 $\mathrm{Path}(y)$ 为从根到包含 $y$ 的叶子节点的路径节点集合。则对任意活跃区间 $I$：

- 若 $y\notin I$，则 $I\notin B(v)$ 对所有 $v\in\mathrm{Path}(y)$；
- 若 $y\in I$，则存在唯一节点 $v^\star\in \mathrm{Path}(y)$ 使得 $v^\star\in\mathcal C(I)$，等价地 $I\in B(v^\star)$。

因此覆盖集合满足不交并分解
$$
\mathcal I(y)=\biguplus_{v\in\mathrm{Path}(y)} B(v).
$$

**证明**  
$\mathcal C(I)$ 中节点区间两两不交，且并为 $I$。路径 $\mathrm{Path}(y)$ 的节点区间两两嵌套，并包含叶子区间 $\mathrm{seg}(\mathrm{leaf}(y))\ni y$。若 $y\in I$，则 $\mathrm{seg}(\mathrm{leaf}(y))\subseteq I$ 必被 $\mathcal C(I)$ 中唯一某个节点区间覆盖，而该覆盖节点必须位于根到该叶子的路径上。唯一性来自 $\mathcal C(I)$ 的不交性。证毕。

### 3.2.2 Path 混合采样骨架

对每个桶 $v$ 维护桶权重和
$$
T(v)=\sum_{I\in B(v)} w(I),
\qquad
m(v)=|B(v)|.
$$
对查询点 $y$，路径权重和
$$
T_{\mathrm{path}}(y)=\sum_{v\in\mathrm{Path}(y)}T(v)=\sum_{I\in\mathcal I(y)}w(I).
$$
若 $T_{\mathrm{path}}(y)=0$，则 $\mathcal I(y)=\varnothing$，返回 $\bot$。

当 $T_{\mathrm{path}}(y)>0$ 时，定义混合过程：

1. 以概率 $T(v)/T_{\mathrm{path}}(y)$ 从 $\mathrm{Path}(y)$ 选取桶 $v$；
2. 条件于选中的桶 $v$，从 $B(v)$ 内按权重比例抽取一个区间并返回。

由引理 3.3，上述混合若桶内条件抽样是精确的，则整体分布精确为 $\mathcal P_y$。本章的任务是在保持混合系数完全精确（即 $T(v)$ 精确维护）的同时，实现一个桶内近似抽样器，并给出 TV 误差的可组合上界。

---

## 3.3 dyadic 包络：二幂层、统一提案与接受率下界

固定一个桶 $v$，其活跃元素（区间）集合记为 $B(v)$。对任意元素 $e\in B(v)$，权重记为 $w(e)>0$。

### 定义 3.4（dyadic level 与包络）

定义 dyadic level
$$
k(e)=\lceil \log_2 w(e)\rceil \in\mathbb Z,
$$
并定义二幂包络与接受率
$$
q(e)=2^{k(e)},\qquad a(e)=\frac{w(e)}{q(e)}.
$$

### 引理 3.5（包络性质与接受率下界）

对任意 $e$，有
$$
w(e)\le q(e) < 2w(e),
\qquad
a(e)\in\left(\tfrac12,1\right]\ \text{且}\ a(e)\ge \tfrac12.
$$

**证明**  
由 $k(e)=\lceil\log_2 w(e)\rceil$ 可得 $2^{k(e)-1}<w(e)\le 2^{k(e)}=q(e)$，因此 $w(e)\le q(e)$ 且 $q(e)<2w(e)$。于是
$$
a(e)=\frac{w(e)}{q(e)}>\frac12,\quad\text{并因而 }a(e)\ge\frac12.
$$
证毕。

该引理意味着：若在桶内先按 $q(e)$ 作为提案权重抽样，再以概率 $a(e)=w(e)/q(e)$ 接受，则被接受样本的权重正比于 $q(e)\cdot a(e)=w(e)$。

---

## 3.4 平均阈值整层裁剪：TV 可控与有效层窗口

桶内按 level 做精确动态加权选择通常需要复杂结构（前驱/平衡树/动态 alias 等）。本章采用一条更直接的路线：允许在桶内丢弃一小部分总权重（从而产生可控 TV 误差），以换取可扫描的 level 窗口。

### 3.4.1 裁剪定义

对桶 $v$，记
$$
m=m(v)=|B(v)|,\qquad T=T(v)=\sum_{e\in B(v)}w(e).
$$
当 $m=0$ 时桶为空；当 $m>0$ 时 $T>0$。

给定参数 $\varepsilon_{\mathrm{clip}}\in(0,1)$，定义阈值
$$
\tau(v)=\varepsilon_{\mathrm{clip}}\cdot \frac{T(v)}{m(v)}\quad(m(v)>0).
$$

按 *$q$ 整层* 进行裁剪：
$$
C(v)=\{e\in B(v):\ q(e)<\tau(v)\},
\qquad
K(v)=B(v)\setminus C(v).
$$
记裁剪质量与保留质量
$$
T_{\mathrm{clip}}(v)=\sum_{e\in C(v)} w(e),
\qquad
T_{\mathrm{keep}}(v)=\sum_{e\in K(v)} w(e)=T(v)-T_{\mathrm{clip}}(v).
$$

### 引理 3.6（裁剪质量上界）

若 $m(v)>0$，则
$$
T_{\mathrm{clip}}(v)\le \varepsilon_{\mathrm{clip}}\cdot T(v).
$$

**证明**  
对任意 $e\in C(v)$，由 $w(e)\le q(e)<\tau(v)$ 得 $w(e)<\tau(v)$。因此
$$
T_{\mathrm{clip}}(v)=\sum_{e\in C(v)}w(e)
< |C(v)|\cdot \tau(v)
\le m(v)\cdot \tau(v)
=\varepsilon_{\mathrm{clip}}\cdot T(v).
$$
证毕。

**推论 3.7（保留集非空）**  
若 $m(v)>0$ 且 $\varepsilon_{\mathrm{clip}}<1$，则 $K(v)\neq\varnothing$，且
$$
T_{\mathrm{keep}}(v)\ge (1-\varepsilon_{\mathrm{clip}})T(v)>0.
$$

### 3.4.2 裁剪对应的 TV 误差

定义桶内精确目标分布
$$
\mathcal P_v(e)=\frac{w(e)}{T(v)}\quad(e\in B(v)).
$$
定义裁剪后的条件分布（重归一到 $K(v)$）
$$
\mathcal P_{v,\mathrm{keep}}(e)=
\begin{cases}
\dfrac{w(e)}{T_{\mathrm{keep}}(v)}, & e\in K(v),\\[6pt]
0, & e\in C(v).
\end{cases}
$$

### 引理 3.8（删去质量等于 TV）

设 $\Omega$ 为有限集合，$\mathcal P$ 为其上分布，$S\subseteq \Omega$。令 $\mathcal P(\cdot\mid S)$ 为把 $\mathcal P$ 限制在 $S$ 上并重归一的条件分布，则
$$
d_{\mathrm{TV}}\big(\mathcal P,\mathcal P(\cdot\mid S)\big)=\mathcal P(\Omega\setminus S).
$$

**证明**  
记 $p=\mathcal P(S)\in(0,1]$，令 $\mathcal Q=\mathcal P(\cdot\mid S)$。则对 $x\in S$，$\mathcal Q(x)=\mathcal P(x)/p$；对 $x\notin S$，$\mathcal Q(x)=0$。于是
$$
\|\mathcal P-\mathcal Q\|_1
=\sum_{x\in S}\left|\mathcal P(x)-\frac{\mathcal P(x)}{p}\right|+\sum_{x\notin S}\mathcal P(x)
= (1-p) + (1-p)
=2(1-p),
$$
从而 $d_{\mathrm{TV}}(\mathcal P,\mathcal Q)=\tfrac12\|\mathcal P-\mathcal Q\|_1=1-p=\mathcal P(\Omega\setminus S)$。证毕。

由引理 3.6 与引理 3.8 得：

**推论 3.9（桶内裁剪 TV 上界）**
当 $m(v)>0$ 时，
$$
d_{\mathrm{TV}}(\mathcal P_v,\mathcal P_{v,\mathrm{keep}})
=\mathcal P_v(C(v))
=\frac{T_{\mathrm{clip}}(v)}{T(v)}
\le \varepsilon_{\mathrm{clip}}.
$$

### 3.4.3 保留层窗口长度

为了让查询时扫描代价可控，需要证明保留元素的 level 落在一个长度为 $O(\log(m/\varepsilon_{\mathrm{clip}}))$ 的窗口里。

定义窗口端点：
$$
k_{\min}(v)=\left\lceil \log_2\tau(v)\right\rceil,
\qquad
k_{\text{high}}(v)=\left\lceil \log_2(2T(v))\right\rceil.
$$

### 引理 3.10（有效层窗口）

设 $m(v)>0$。则任意保留元素 $e\in K(v)$ 的 level 满足
$$
k_{\min}(v)\le k(e)\le k_{\text{high}}(v).
$$
并且窗口长度满足
$$
k_{\text{high}}(v)-k_{\min}(v)+1
\le 3+\left\lceil \log_2\frac{2m(v)}{\varepsilon_{\mathrm{clip}}}\right\rceil
=O\!\left(\log\frac{m(v)}{\varepsilon_{\mathrm{clip}}}\right).
$$

**证明**  

- 上界：对任意 $e\in B(v)$ 有 $w(e)\le T(v)$，由引理 3.5 得 $q(e)<2w(e)\le 2T(v)$，即 $k(e)=\log_2 q(e) < \log_2(2T(v))$，从而 $k(e)\le \lceil \log_2(2T(v))\rceil=k_{\text{high}}(v)$。
- 下界：对任意 $e\in K(v)$ 有 $q(e)\ge \tau(v)$，从而 $k(e)\ge \lceil \log_2\tau(v)\rceil=k_{\min}(v)$。
- 长度：由定义
  $$
  \begin{aligned}
  k_{\text{high}}-k_{\min}
  &\le \left(\log_2(2T)+1\right)-\log_2(\tau)
  =\log_2\!\left(\frac{2T}{\tau}\right)+1\\
  &=\log_2\!\left(\frac{2T}{\varepsilon_{\mathrm{clip}}T/m}\right)+1
  =\log_2\!\left(\frac{2m}{\varepsilon_{\mathrm{clip}}}\right)+1,
  \end{aligned}
  $$
  再加上端点取整产生的常数项即可。证毕。

---

## 3.5 桶内确定性更新：预分配数组段与 swap-delete 句柄

本节给出桶内数据结构，使得每条 canonical membership 的插删为确定最坏 $O(1)$，且运行期不依赖哈希、平衡树或均摊扩容。

### 3.5.1 对象与内存布局

**(1) membership record**  
每条 record 对应“某个元素 $e$ 在某个桶 $v$ 的出现”（即 $v\in\mathcal C(e)$）。record 维护字段：

- `eid`：元素编号；
- `v`：桶编号；
- `k`：该元素的 dyadic level $k(e)$（由权重静态决定）；
- `segptr`：指向数组段 $\mathrm{Seg}(v,k)$ 的指针/索引；
- `pos`：若该 record 处于活跃态，则为其在数组段 `cell` 中的位置；否则为 $-1$。

**(2) 桶-层数组段**  
对每个桶 $v$ 与其可能出现的 level $k$，维护数组段 $\mathrm{Seg}(v,k)$：

- `cap(v,k)`：固定容量（预处理阶段确定）；
- `sz(v,k)`：当前活跃大小；
- `cell(v,k)[0..cap-1]`：存放活跃 record id（即 rid）。

**(3) 元素的固定 membership 列表**  
对每个元素 $e$，其规范覆盖 $\mathcal C(e)$ 固定，长度 $d(e)=|\mathcal C(e)|=O(\log N)$。预处理阶段为其分配定长数组 `MemList(e)[0..d(e)-1]`，存放所有对应 record id。运行期插入/删除元素时仅遍历这段数组，不做任何查找。

**(4) 桶聚合量**  
每个桶 $v$ 维护：
$$
m(v)=\sum_k \mathrm{sz}(v,k),\qquad
T(v)=\sum_{e\in B(v)}w(e).
$$
插删 record 时同步更新 $m(v)$ 与 $T(v)$。

### 3.5.2 Insert/Delete：swap-delete 伪代码

下面用 `rid` 表示 record id，`Rec[rid]` 为 record 对象。

~~~text
Procedure InsertRecord(rid):
  rec ← Rec[rid]
  if rec.pos ≠ -1: return            // idempotent
  seg ← rec.segptr
  assert seg.sz < seg.cap

  p ← seg.sz
  seg.cell[p] ← rid
  seg.sz ← seg.sz + 1
  rec.pos ← p

  v ← rec.v
  eid ← rec.eid
  m(v) ← m(v) + 1
  T(v) ← T(v) + w(eid)

Procedure DeleteRecord(rid):
  rec ← Rec[rid]
  if rec.pos = -1: return            // idempotent
  seg ← rec.segptr
  p ← rec.pos
  last ← seg.sz - 1
  rid2 ← seg.cell[last]

  seg.cell[p] ← rid2                 // swap
  seg.sz ← seg.sz - 1
  rec.pos ← -1

  if rid2 ≠ rid:
      Rec[rid2].pos ← p              // fix pointer

  v ← rec.v
  eid ← rec.eid
  m(v) ← m(v) - 1
  T(v) ← T(v) - w(eid)
~~~

### 定理 3.11（每条 membership 更新确定最坏 $O(1)$）

在运行期不发生数组扩容且每个数组段容量 `cap(v,k)` 预先正确分配的前提下，`InsertRecord` 与 `DeleteRecord` 对任意 record 的最坏时间均为 $O(1)$。

**证明**  
过程包含常数次数组访问与常数次标量加减，无循环、无递归、无查找；且 swap-delete 只涉及常数次写回。容量不溢出的前提由引理 3.14 给出。证毕。

---

## 3.6 预处理：record 生成、分组与容量分配

预处理的目标是在运行期开始前完成以下事项：

1. 为每个元素 $e$ 生成其 membership records，并填充 `MemList(e)`；
2. 按键 $(v,k)$ 分组 records，得到每个 $(v,k)$ 的容量 `cap(v,k)`；
3. 分配并初始化每个数组段 $\mathrm{Seg}(v,k)$，同时为每条 record 写入 `segptr`；
4. 为每个桶 $v$ 生成静态有序 level 列表 $\mathrm{KList}(v)$（包含该桶中可能出现的所有 $k$ 值）。

### 3.6.1 record 数量界

令元素总数为 $N$，段树叶子数为 $q=O(N)$。

### 引理 3.12（总 record 数）

令
$$
P=\sum_{e}|\mathcal C(e)|.
$$
则 $P=O(N\log N)$。

**证明**  
对任意元素 $e$，规范覆盖大小 $|\mathcal C(e)|\le 2\lceil\log_2 q\rceil=O(\log N)$，对 $N$ 个元素求和即得。证毕。

### 3.6.2 键编码：处理负 level

为分组，需要把 record 的二元键 $(v,k)$ 编码为非负整数。设段树节点总数为 $|V|=2q-1$，编号 $\mathrm{bucketID}(v)\in[0,|V|-1]$。

level $k$ 可能为负。预处理阶段先扫描所有元素得到
$$
k_{\mathrm{low}}=\min_e k(e),\qquad k_{\mathrm{high}}=\max_e k(e),
$$
并定义偏移后的非负 level
$$
k'(e)=k(e)-k_{\mathrm{low}}\in[0,\ k_{\mathrm{high}}-k_{\mathrm{low}}].
$$
令
$$
B_k=\left\lceil \log_2\big(k_{\mathrm{high}}-k_{\mathrm{low}}+1\big)\right\rceil,\qquad
B_v=\lceil \log_2|V|\rceil,
$$
则可将键编码为
$$
\mathrm{key}(v,k)=\big(\mathrm{bucketID}(v)\ll B_k\big)\ \vert\ (k-k_{\mathrm{low}}),
$$
这是一个 $B=B_v+B_k$ 位非负整数。

### 3.6.3 分组算法：线性版与通用版

对 $P$ 条 record 按 key 分组的需求仅是得到“同 key 的连续段”，不要求稳定排序之外的性质。这里给出两种确定性实现口径。

#### 定理 3.13（确定性分组）

给定 $P$ 条记录及其 $B$ 位非负整数 key：

1. **通用版**：可以用确定性的比较排序在 $O(P\log P)$ 时间内按 key 排序，从而得到连续分组；
2. **线性版**：若在 word-RAM 模型下 $B=O(\log P)$（等价地 key 宇宙规模为 $P^{O(1)}$，且 key 可装入常数个机器字），则存在确定性的 LSD radix 分组算法，在 $O(P)$ 时间与 $O(P)$ 额外空间内得到按 key 的连续分组。

**证明要点**  
(1) 使用任意确定性比较排序即可。  
(2) 取基数 $R=2^r$，其中 $r=\lceil\log_2 P\rceil$，则计数数组大小 $R\le 2P$。当 $B\le 2r$ 时两趟 radix 即可完成；更一般地，趟数为 $\lceil B/r\rceil=O(1)$（因为 $B=O(\log P)$）。每趟时间 $O(P+R)=O(P)$，总计 $O(P)$。证毕。

> 注：在本书后续实例中，$\mathrm{bucketID}$ 来自段树节点编号（$O(N)$ 范围），而 $k$ 来自有限精度权重表示的指数位，通常满足 $B=O(\log N)=O(\log P)$。

### 3.6.4 容量分配与指针写回

按 key 分组后，对每个组（对应某个 $(v,k)$）：

- 组大小即为 `cap(v,k)`；
- 在全局数组中为 `cell(v,k)` 分配一段连续切片；
- 初始化 `sz(v,k)=0`；
- 令组内每条 record 的 `segptr` 指向该段。

同时，对每个桶 $v$，将所有出现过的 level $k$ 以升序写入静态数组 $\mathrm{KList}(v)$。该数组在运行期只读。

### 引理 3.14（容量安全）

对任意 $(v,k)$，运行期任意时刻都有
$$
\mathrm{sz}(v,k)\le \mathrm{cap}(v,k).
$$

**证明**  
`cap(v,k)` 在预处理中计为“所有元素中，规范覆盖包含 $v$ 且 level 为 $k$ 的 record 总数”。运行期 `sz(v,k)` 仅统计其中当前活跃的 record 数，是其子集大小，故不超过 `cap(v,k)`。证毕。

---

## 3.7 桶内 TV-近似抽样：窗口扫描、alias 与有界拒绝

本节给出桶内抽样器 `BucketSample(v)`，其目标是在桶 $v$ 内近似采样 $\mathcal P_v(e)\propto w(e)$，并满足：

- **最坏时间有界**：不允许“直到接受”为止”的随机循环；
- **TV 上界显式**：裁剪引入 $\varepsilon_{\mathrm{clip}}$，有界拒绝引入 $2^{-M_{\mathrm{bucket}}}$。

### 3.7.1 提案分布与拒绝校正

对保留集 $K(v)$，定义提案分布
$$
\mathcal Q_v(e)=\frac{q(e)}{\sum_{u\in K(v)}q(u)}\quad (e\in K(v)).
$$
以接受概率
$$
a(e)=\frac{w(e)}{q(e)}\in\left(\tfrac12,1\right]
$$
做拒绝采样：提案被接受的非归一化权重为 $\mathcal Q_v(e)\cdot a(e)\propto w(e)$，因此无限次重试的极限分布为 $\mathcal P_{v,\mathrm{keep}}$。

实际算法用 *有界拒绝*：最多尝试 $M_{\mathrm{bucket}}$ 次，若全部拒绝则走固定失败分支，以保证确定上界。

### 3.7.2 查询窗口与小 alias

对桶 $v$，若 $m(v)=0$ 则返回 $\bot$。否则计算
$$
\tau=\varepsilon_{\mathrm{clip}}\cdot \frac{T(v)}{m(v)},\quad
k_{\min}=\lceil\log_2\tau\rceil,\quad
k_{\text{high}}=\lceil\log_2(2T(v))\rceil.
$$
由引理 3.10，保留元素仅可能出现在 level 窗口 $[k_{\min},k_{\text{high}}]$，其长度为
$$
r_{\max}=k_{\text{high}}-k_{\min}+1
=O\!\left(\log\frac{m(v)}{\varepsilon_{\mathrm{clip}}}\right).
$$

实现上，桶 $v$ 有静态递增列表 $\mathrm{KList}(v)$。查询时用二分定位窗口起点，然后顺序扫描 $\mathrm{KList}(v)$ 中落入窗口的条目；由于 level 为整数且严格递增，该扫描次数不超过窗口长度。

对每个扫描到且当前活跃计数为 $c_k=\mathrm{sz}(v,k)>0$ 的 level $k$，定义层质量
$$
M_k=c_k\cdot 2^k.
$$
为避免处理大指数浮点，统一除以 $2^{k_{\min}}$，定义整数质量
$$
M_k'=c_k\cdot 2^{k-k_{\min}}\in\mathbb N.
$$
按 $\{M_k'\}$ 构建规模 $r\le r_{\max}$ 的 alias 表，随后每次“抽 level”均为 $O(1)$。

### 3.7.3 桶内采样算法（伪代码）

下面的伪代码使用两类随机原语：

- `RandInt(M)`：返回 $\{0,1,\dots,M-1\}$ 上均匀整数；
- `Bernoulli(p)`：返回 $1$ 的概率为 $p$（$0<p\le 1$）。

这些原语只消耗随机比特，不包含不定循环。

~~~text
Procedure BucketSample(v):
  if m(v) = 0: return ⊥
  T ← T(v); m ← m(v)

  τ ← ε_clip * T / m
  kmin  ← ceil_log2(τ)
  khigh ← ceil_log2(2*T)

  // 1) collect active levels inside [kmin, khigh]
  i ← LowerBound(KList(v), kmin)      // deterministic binary search
  r ← 0
  Wsum ← 0
  while i < |KList(v)| and KList(v)[i] ≤ khigh:
      k ← KList(v)[i]
      seg ← Seg(v,k)
      c ← seg.sz
      if c > 0:
          Level[r] ← k
          Mass[r]  ← c * 2^(k-kmin)    // integer
          Wsum ← Wsum + Mass[r]
          r ← r + 1
      i ← i + 1
  assert r > 0                         // since T>0 and ε_clip<1

  // 2) build alias table for Mass[0..r-1]
  BuildAliasInt(Mass, r, Wsum, Prob[0..r-1], Alias[0..r-1])

  // 3) bounded rejection (deterministic loop bound)
  for t = 1..M_bucket:
      j ← AliasSampleInt(Prob, Alias, r, Wsum)   // O(1)
      k ← Level[j]
      seg ← Seg(v,k)
      p ← RandInt(seg.sz)
      rid ← seg.cell[p]
      eid ← Rec[rid].eid

      if Bernoulli( w(eid) / 2^k ) = 1:
          return eid
      last ← eid
  return last
~~~

其中 `BuildAliasInt / AliasSampleInt` 可采用整数版 alias（避免浮点）：

~~~text
Procedure BuildAliasInt(Mass[0..r-1], r, Wsum, Prob, Alias):
  // scale: S[i] = Mass[i] * r, compare to Wsum
  for i=0..r-1:
      S[i] ← Mass[i] * r

  init stacks Small, Large
  for i=0..r-1:
      if S[i] < Wsum: push Small(i) else push Large(i)

  while Small not empty and Large not empty:
      i ← pop Small
      j ← pop Large
      Prob[i] ← S[i]            // integer threshold in [0, Wsum)
      Alias[i] ← j
      S[j] ← S[j] + S[i] - Wsum
      if S[j] < Wsum: push Small(j) else push Large(j)

  while Large not empty:
      j ← pop Large
      Prob[j] ← Wsum
      Alias[j] ← j
  while Small not empty:
      i ← pop Small
      Prob[i] ← Wsum
      Alias[i] ← i

Procedure AliasSampleInt(Prob, Alias, r, Wsum):
  i ← RandInt(r)
  u ← RandInt(Wsum)          // uniform in {0,...,Wsum-1}
  if u < Prob[i]: return i else return Alias[i]
~~~

### 3.7.4 TV 误差界：有界拒绝项

记 `BucketSample(v)` 的输出分布为 $\widetilde{\mathcal P}_v$。

### 引理 3.15（有界拒绝的 TV 上界）

在保留集 $K(v)$ 上，以提案分布 $\mathcal Q_v$ 独立提案，并以接受概率 $a(e)=w(e)/q(e)$ 拒绝采样；最多尝试 $M_{\mathrm{bucket}}$ 次，若全部拒绝则输出第 $M_{\mathrm{bucket}}$ 次提案。则
$$
d_{\mathrm{TV}}\big(\widetilde{\mathcal P}_v,\ \mathcal P_{v,\mathrm{keep}}\big)
\le (1-a_{\min})^{M_{\mathrm{bucket}}}
\le 2^{-M_{\mathrm{bucket}}},
$$
其中 $a_{\min}=\inf_{e\in K(v)}a(e)\ge \tfrac12$（引理 3.5）。

**证明**  
令事件 $\mathsf{Fail}$ 表示前 $M_{\mathrm{bucket}}$ 次均拒绝。无限拒绝（直到接受）诱导分布为 $\mathcal P_{v,\mathrm{keep}}$。有界过程的输出分布可写为凸组合
$$
\widetilde{\mathcal P}_v=\Pr[\neg\mathsf{Fail}]\cdot \mathcal P_{v,\mathrm{keep}}
+\Pr[\mathsf{Fail}]\cdot \mathcal Q_{\mathsf{fail}},
$$
其中 $\mathcal Q_{\mathsf{fail}}$ 为失败分支输出分布。于是
$$
d_{\mathrm{TV}}(\widetilde{\mathcal P}_v,\mathcal P_{v,\mathrm{keep}})
\le \Pr[\mathsf{Fail}].
$$
每次尝试的接受概率至少为 $a_{\min}$，各次尝试独立，因此
$$
\Pr[\mathsf{Fail}]\le (1-a_{\min})^{M_{\mathrm{bucket}}}\le 2^{-M_{\mathrm{bucket}}}.
$$
证毕。

### 定理 3.16（桶内总误差）

当 $m(v)>0$ 时，
$$
d_{\mathrm{TV}}\big(\widetilde{\mathcal P}_v,\ \mathcal P_v\big)
\le \varepsilon_{\mathrm{clip}}+2^{-M_{\mathrm{bucket}}}
=: \varepsilon_{\mathrm{bucket}}.
$$

**证明**  
由三角不等式
$$
d_{\mathrm{TV}}(\widetilde{\mathcal P}_v,\mathcal P_v)
\le d_{\mathrm{TV}}(\widetilde{\mathcal P}_v,\mathcal P_{v,\mathrm{keep}})
+d_{\mathrm{TV}}(\mathcal P_{v,\mathrm{keep}},\mathcal P_v)
\le 2^{-M_{\mathrm{bucket}}}+\varepsilon_{\mathrm{clip}},
$$
其中第二项由推论 3.9 给出，第一项由引理 3.15 给出。证毕。

### 引理 3.17（桶内查询最坏时间）

`BucketSample(v)` 的最坏时间为
$$
O\!\left(\log|\mathrm{KList}(v)|
+\log\frac{m(v)}{\varepsilon_{\mathrm{clip}}}
+M_{\mathrm{bucket}}\right)
\subseteq
O\!\left(\log N+\log\frac{m(v)}{\varepsilon_{\mathrm{clip}}}+M_{\mathrm{bucket}}\right).
$$

**证明要点**  
二分定位窗口起点 $O(\log|\mathrm{KList}(v)|)$；窗口扫描访问的 level 数不超过 $r_{\max}=O(\log(m/\varepsilon_{\mathrm{clip}}))$；alias 构建 $O(r_{\max})$；有界拒绝循环次数固定为 $M_{\mathrm{bucket}}$，每次 $O(1)$。证毕。

---

## 3.8 回到 stabbing：Path 混合误差不放大与主定理

本节把桶内结论接入 Path 混合骨架，得到动态加权 stabbing 的整体复杂度与误差界。

### 引理 3.18（同权混合下 TV 不放大）

固定查询点 $y$，设路径上每个桶 $v\in\mathrm{Path}(y)$ 的精确桶内分布为 $\mathcal P_v$，近似桶内分布为 $\widetilde{\mathcal P}_v$，并满足
$$
d_{\mathrm{TV}}(\widetilde{\mathcal P}_v,\mathcal P_v)\le \varepsilon_{\mathrm{bucket}}.
$$
定义精确与近似的路径混合分布
$$
\mathcal P_y=\sum_{v\in\mathrm{Path}(y)}\frac{T(v)}{T_{\mathrm{path}}(y)}\cdot \mathcal P_v,
\qquad
\widetilde{\mathcal P}_y=\sum_{v\in\mathrm{Path}(y)}\frac{T(v)}{T_{\mathrm{path}}(y)}\cdot \widetilde{\mathcal P}_v,
$$
其中混合系数使用同一组精确权重 $T(v)$。则
$$
d_{\mathrm{TV}}(\widetilde{\mathcal P}_y,\mathcal P_y)\le \varepsilon_{\mathrm{bucket}}.
$$

**证明**  
用 TV 的凸性：
$$
d_{\mathrm{TV}}\!\left(\sum_v \lambda_v\widetilde{\mathcal P}_v,\sum_v \lambda_v \mathcal P_v\right)
\le \sum_v \lambda_v\, d_{\mathrm{TV}}(\widetilde{\mathcal P}_v,\mathcal P_v),
\quad
\lambda_v=\frac{T(v)}{T_{\mathrm{path}}(y)}.
$$
并用每项上界 $\varepsilon_{\mathrm{bucket}}$ 即得。证毕。

### 定理 3.19（动态加权 stabbing：确定性更新 + TV-近似查询）

给定端点集合 $Y$ 上的段树骨架 $\mathcal T_y$，以及一组静态带权区间集合（权重 $w(I)>0$）。存在一个 stabbing 采样器满足：

1. **预处理**：生成 $P=O(N\log N)$ 条 membership records，分配所有桶-层数组段并填充 `MemList`。  
   - 在通用版分组口径下，预处理时间 $O(P\log P)$，空间 $O(P)$；  
   - 在定理 3.13 的线性版口径下，预处理时间 $O(P)$，空间 $O(P)$。
2. **更新最坏时间**：一次 `Insert(I)` 或 `Delete(I)` 仅遍历 $\mathcal C(I)$ 的 $O(\log N)$ 条 membership record，并对每条 record 做 $O(1)$ 插删（定理 3.11），因此更新最坏时间为
   $$
   O(\log N).
   $$
3. **查询最坏时间**：`Sample(y)` 沿 $\mathrm{Path}(y)$ 计算并按 $\{T(v)\}$ 选桶，代价 $O(\log N)$；随后仅对被选中桶 $v^\*$ 调用 `BucketSample(v^\*)`，代价见引理 3.17。因此整体查询最坏时间为
   $$
   O\!\left(\log N+\log\frac{m(v^*)}{\varepsilon_{\mathrm{clip}}}+M_{\mathrm{bucket}}\right).
   $$
4. **TV 误差**：对任意 $y$，输出分布满足
   $$
   d_{\mathrm{TV}}(\widetilde{\mathcal P}_y,\mathcal P_y)
   \le \varepsilon_{\mathrm{bucket}}
   =\varepsilon_{\mathrm{clip}}+2^{-M_{\mathrm{bucket}}}.
   $$

**证明要点**  
更新与查询复杂度分别由段树规范覆盖规模、定理 3.11 与引理 3.17 直接给出；TV 上界由定理 3.16（桶内）与引理 3.18（混合不放大）组合得到。证毕。

---

## 3.9 参数接口（供后续章节直接调用）

本章输出的 stabbing 模块可用两参数描述其逐点 TV 误差上界：
$$
\varepsilon_{\mathrm{stab}}=\varepsilon_{\mathrm{clip}}+2^{-M_{\mathrm{bucket}}}.
$$
给定目标 $\varepsilon_{\mathrm{stab}}$，一种直接选择是：
$$
\varepsilon_{\mathrm{clip}}=\frac{\varepsilon_{\mathrm{stab}}}{2},
\qquad
M_{\mathrm{bucket}}=\left\lceil \log_2\frac{2}{\varepsilon_{\mathrm{stab}}}\right\rceil,
$$
从而 $\varepsilon_{\mathrm{clip}}+2^{-M_{\mathrm{bucket}}}\le \varepsilon_{\mathrm{stab}}$，且查询最坏时间中额外循环项仅为 $M_{\mathrm{bucket}}=\Theta(\log(1/\varepsilon_{\mathrm{stab}}))$。

本章的结构性特征是：更新端完全确定 $O(\log N)$；查询端仅对 *一个* 被选中桶做窗口扫描与 alias 构建；所有近似均以 TV 量化并可在更高层的混合结构中直接传递。

---

# 第 4 章 单项几何子问题的采样器：slab 扫描 + 乘积线段树 + stabbing

本章固定一个指数项编号 $\ell$，研究如下单项几何子权重（定义在跨类有序对空间上）：
$$
W_\ell(r,s)=I(r,s)\,w_\ell(r)\,v_\ell(s),
\qquad 
w_\ell(r)=e^{-\alpha_\ell A(r)},\ 
v_\ell(s)=\omega_\ell e^{-\alpha_\ell A(s)}.
$$
其中 $r\in R_c$（红集合），$s\in R_{\bar c}$（蓝集合），矩形面积 $A(\cdot)>0$，交面积 $I(r,s)=\mathrm{Area}(r\cap s)$。本章的目标是：在不枚举所有跨集合相交对的前提下，构造一个对分布
$$
\mathcal D_\ell(r,s)=\frac{W_\ell(r,s)}{Z_\ell},
\qquad 
Z_\ell=\sum_{r\in R_c}\sum_{s\in R_{\bar c}} W_\ell(r,s)
$$
进行 i.i.d.（有放回）采样的算法，输出 $t$ 个独立同分布样本，并给出严格的最坏时间复杂度与 TV（总变差）误差闭环。

本章的核心结构由三部分组成：

1. **slab 混合**：沿 $x$ 轴把平面分解为不交竖直 slab，使 $Z_\ell$ 写成 slab 质量之和 $Z_\ell=\sum_i M_i$，并将 $\mathcal D_\ell$ 表达为 slab 条件分布的混合；
2. **乘积线段树**：在固定 slab 内把 $M_i$ 写成一维积分
   $$
   M_i=\Delta x_i \int F_i(y)G_i(y)\,dy,
   $$
   并用线段树动态维护 $\int F_i$、$\int G_i$、$\int F_iG_i$，支持按密度 $\propto F_iG_i$ 采样 $y$；
3. **加权 stabbing**：给定 $y$，分别从覆盖 $y$ 的红/蓝矩形中按权重抽样，从而完成 slab 内条件采样。桶内采用确定最坏更新的数据结构与 TV-近似抽样，使扫描更新端避免 $O(\log^2 N)$ 回潮。

---

## 4.0 计算模型与随机性原语

本章的正确性与复杂度分析在如下口径下给出：

- 系统模型为 RAM（real-RAM 或 word-RAM 均可），并假设实现提供常数时间的基础算术与比较原语，用于处理输入坐标、面积、以及本章涉及的权重表示。
- 随机性以可分离的随机流（RNG streams）提供。每个采样过程消耗有限随机位；所有循环均具有确定上界，不允许“直到接受为止”的不定循环。
- 对均匀随机数原语，允许使用：
  - `RandInt(m)`：返回 $\{0,1,\dots,m-1\}$ 上均匀整数；
  - `Unif(a,b)`：在半开区间 $[a,b)$ 上均匀采样（实现上可由 `RandInt(2^b)` 与缩放得到）。

复杂度以这些原语调用次数为计数单位；其中所有数据结构更新与查询均给出确定最坏上界。

---

## 4.1 slab 混合：质量 $M_i$ 与条件分布 $\mathcal D_{\ell,i}$

### 4.1.1 输入几何与半开语义

每个矩形 $r$ 表示为二维轴对齐半开矩形
$$
r=[L_1(r),R_1(r))\times [L_2(r),R_2(r)),
\qquad L_1(r)<R_1(r),\ \ L_2(r)<R_2(r).
$$
半开语义的作用是：边界接触不会产生正面积交叠，从而扫描线分解与端点归属具有严格一致性。

定义一维半开区间长度
$$
\mathrm{len}([a,b))=\max\{0,b-a\}.
$$
对两矩形 $r,s$，定义 $y$ 方向交叠长度
$$
h(r,s)=\mathrm{len}\Big([\max\{L_2(r),L_2(s)\},\ \min\{R_2(r),R_2(s)\})\Big),
$$
并注意 $h(r,s)$ 与 $x$ 无关。

---

### 4.1.2 slab 定义与事件顺序

收集所有矩形的 $x$ 端点集合
$$
X=\{L_1(u),R_1(u)\mid u\in R_c\cup R_{\bar c}\},
$$
去重排序得到
$$
x_0<x_1<\cdots<x_m,\qquad m\le 2N-1.
$$
定义第 $i$ 个竖直 slab
$$
S_i=[x_i,x_{i+1}),\qquad \Delta x_i=x_{i+1}-x_i>0,\qquad i=0,1,\dots,m-1.
$$

扫描按 $x$ 递增处理事件，并规定在同一 $x$ 处先处理 `END`（删除满足 $R_1(\cdot)=x$ 的矩形）后处理 `START`（插入满足 $L_1(\cdot)=x$ 的矩形）。在这一顺序下，处理完 $x=x_i$ 的所有事件后，当前活跃集合恰好对应覆盖 slab $S_i=[x_i,x_{i+1})$ 的矩形集合。

定义 slab 覆盖指示
$$
\chi_r(i)=\mathbf 1\big(S_i\subseteq [L_1(r),R_1(r))\big).
$$
由于 slab 端点来自所有矩形端点集合，任意矩形要么覆盖整个 slab，要么与其不相交，因此在 slab 内 $\chi_r(i)$ 为常量。

---

### 4.1.3 slab 分量权重与总质量分解

两矩形在 slab $S_i$ 内的交面积贡献为
$$
I_i(r,s)=\Delta x_i\cdot \chi_r(i)\chi_s(i)\cdot h(r,s).
$$
进而定义 slab 分量权重
$$
W_{\ell,i}(r,s)=I_i(r,s)\,w_\ell(r)\,v_\ell(s).
$$

**引理 4.1（slab 分量可加性）**  
对任意 $r\in R_c$、$s\in R_{\bar c}$，有
$$
I(r,s)=\sum_{i=0}^{m-1} I_i(r,s),
\qquad
W_\ell(r,s)=\sum_{i=0}^{m-1} W_{\ell,i}(r,s).
$$

*证明*：在每个 slab 上，两矩形在 $x$ 方向的交叠长度要么为 $\Delta x_i$（两者均覆盖该 slab），要么为 $0$。因此 $x$ 方向交叠长度可写为 $\sum_i \Delta x_i\chi_r(i)\chi_s(i)$，乘以 $y$ 方向交叠长度 $h(r,s)$ 得 $I(r,s)=\sum_i I_i(r,s)$。再乘以 $w_\ell(r)v_\ell(s)$ 得第二式。$\square$

定义 slab 质量与总质量：
$$
M_i=\sum_{r\in R_c}\sum_{s\in R_{\bar c}} W_{\ell,i}(r,s),
\qquad
Z_\ell=\sum_{i=0}^{m-1}M_i.
$$
若 $Z_\ell=0$，表示该单项子权重在跨类对上处处为零，按约定采样器输出空序列。

---

### 4.1.4 slab 条件分布与混合形式

当 $M_i>0$ 时定义 slab 条件分布
$$
\mathcal D_{\ell,i}(r,s)=\frac{W_{\ell,i}(r,s)}{M_i}.
$$
由分量混合恒等式可得：若先按概率 $M_i/Z_\ell$ 抽取 slab $i$，再条件于 $i$ 从 $\mathcal D_{\ell,i}$ 采样 $(r,s)$，则边缘分布即为 $\mathcal D_\ell$。

因此本章的关键在于实现两类接口：

- `Mass(i)`：在扫描状态下得到 slab $S_i$ 的质量 $M_i$；
- `CondSample(i)`：在扫描状态下（已维护足够数据结构）从 $\mathcal D_{\ell,i}$ 采样。

---

## 4.2 乘积线段树：维护 $\int F_i$、$\int G_i$、$\int F_iG_i$ 并按密度采样 $y$

本节固定 slab $S_i$（扫描过程中 $i$ 随 $x$ 变化），把 $M_i$ 写成一维积分，并给出支持在线更新与采样的数据结构。

### 4.2.1 $y$ 端点压缩与基本段

收集所有矩形的 $y$ 端点集合
$$
Y=\{L_2(u),R_2(u)\mid u\in R_c\cup R_{\bar c}\},
$$
去重排序得到
$$
y_0<y_1<\cdots<y_q,\qquad q\le 2N-1.
$$
定义基本段
$$
I_j=[y_j,y_{j+1}),\qquad j=0,1,\dots,q-1.
$$
由于所有矩形的 $y$ 边界都来自 $Y$，在任意 slab 内下文定义的覆盖函数在每个基本段上为常数，且在区间 $[y_0,y_q)$ 外恒为 $0$。

---

### 4.2.2 slab 内加权覆盖函数

在 slab $S_i$ 内，定义红/蓝活跃集合：
$$
A_c(i)=\{r\in R_c\mid \chi_r(i)=1\},\qquad
A_{\bar c}(i)=\{s\in R_{\bar c}\mid \chi_s(i)=1\}.
$$

定义加权覆盖函数：
$$
F_i(y)=\sum_{r\in A_c(i)} w_\ell(r)\cdot \mathbf 1[y\in r^y],
\qquad
G_i(y)=\sum_{s\in A_{\bar c}(i)} v_\ell(s)\cdot \mathbf 1[y\in s^y],
$$
其中 $r^y=[L_2(r),R_2(r))$。

---

### 4.2.3 slab 质量的积分形式

**引理 4.2（质量积分表达）**  
对任意 slab $S_i$，有
$$
M_i=\Delta x_i\int_{y_0}^{y_q}F_i(y)G_i(y)\,dy.
$$

*证明*：展开乘积并交换有限求和与积分：
$$
\int_{y_0}^{y_q}F_i(y)G_i(y)\,dy
=\sum_{r\in A_c(i)}\sum_{s\in A_{\bar c}(i)} w_\ell(r)v_\ell(s)\int_{y_0}^{y_q} \mathbf 1[y\in r^y]\mathbf 1[y\in s^y]\,dy.
$$
积分项等于 $r^y$ 与 $s^y$ 的交叠长度 $h(r,s)$，因此
$$
\int_{y_0}^{y_q}F_iG_i
=\sum_{r\in A_c(i)}\sum_{s\in A_{\bar c}(i)} w_\ell(r)v_\ell(s)\,h(r,s).
$$
两边乘以 $\Delta x_i$，并注意 $W_{\ell,i}(r,s)=\Delta x_i\,h(r,s)\,w_\ell(r)v_\ell(s)$，得到 $M_i=\sum_{r,s}W_{\ell,i}(r,s)$。$\square$

---

### 4.2.4 乘积线段树：节点量与闭式 lazy 更新

建立线段树覆盖区间 $[y_0,y_q)$，叶子对应基本段 $I_j=[y_j,y_{j+1})$。

对每个节点 $u$，令 $\mathrm{seg}(u)=[a(u),b(u))$，$\mathrm{len}(u)=b(u)-a(u)$。节点维护三项积分量：
- $\mathrm{SF}(u)=\int_{\mathrm{seg}(u)}F(y)\,dy$；
- $\mathrm{SG}(u)=\int_{\mathrm{seg}(u)}G(y)\,dy$；
- $\mathrm{SP}(u)=\int_{\mathrm{seg}(u)}F(y)G(y)\,dy$；

并维护懒标记 $\mathrm{addF}(u),\mathrm{addG}(u)$，表示对整个区间执行 $F\leftarrow F+\mathrm{addF}$、$G\leftarrow G+\mathrm{addG}$ 但尚未下推。

由于 $F,G$ 在叶子段上恒为常数，区间加常数的积分更新可以用闭式公式完成：

- 若在节点区间上执行 $F\leftarrow F+\Delta$（$\Delta$ 为常数），则
  $$
  \mathrm{SF}\leftarrow \mathrm{SF}+\Delta\cdot \mathrm{len},\qquad
  \mathrm{SP}\leftarrow \mathrm{SP}+\Delta\cdot \mathrm{SG},\qquad
  \mathrm{addF}\leftarrow \mathrm{addF}+\Delta.
  $$
- 若执行 $G\leftarrow G+\Delta$，则
  $$
  \mathrm{SG}\leftarrow \mathrm{SG}+\Delta\cdot \mathrm{len},\qquad
  \mathrm{SP}\leftarrow \mathrm{SP}+\Delta\cdot \mathrm{SF},\qquad
  \mathrm{addG}\leftarrow \mathrm{addG}+\Delta.
  $$

使用标准线段树区间更新（递归覆盖并在必要时下推懒标记）即可在 $O(\log q)=O(\log N)$ 最坏时间内完成一次“对 $F$ 或 $G$ 的区间加常数”更新。

---

### 4.2.5 按密度 $\propto F_iG_i$ 采样 $y$

记当前 slab 的乘积积分
$$
T_i=\int_{y_0}^{y_q} F_i(y)G_i(y)\,dy,
$$
它等于根节点的 $\mathrm{SP}(\mathrm{root})$。当 $T_i>0$ 时，定义密度
$$
p_i(y)=\frac{F_i(y)G_i(y)}{T_i},\qquad y\in[y_0,y_q).
$$

下面给出 `SampleY()` 的实现方式：从根开始向下，每次在左右子节点中按子区间的 $\mathrm{SP}$ 比例选择，直到到达叶子段，再在该基本段内均匀采样。

```text
Procedure SampleY():
  u ← root
  while u is not a leaf:
      PushDown(u)                       // ensure children values are correct
      L ← left(u), R ← right(u)
      tL ← SP(L), tR ← SP(R)            // both ≥ 0
      x ← RandReal(0, tL + tR)          // or RandInt on scaled integer
      if x < tL: u ← L else u ← R
  // u is a leaf with segment [y_j, y_{j+1})
  return Unif(a(u), b(u))
```

**引理 4.3（按 $F_iG_i$ 密度采样的正确性）**
 当 $T_i>0$ 时，上述 `SampleY()` 输出随机变量 $Y$ 的密度为
$$
p_i(y)=\frac{F_i(y)G_i(y)}{\int_{y_0}^{y_q} F_i(u)G_i(u)\,du}.
$$
*证明*：叶子段上 $F_i,G_i$ 为常数，因此 $F_iG_i$ 在每个基本段 $I_j$ 上为常数。`SampleY()` 以概率
$$
\frac{\int_{I_j}F_iG_i}{\int_{y_0}^{y_q}F_iG_i}
$$
选择叶子段 $I_j$，并在该段内均匀采样，因而得到分段常数密度 $F_iG_i/\int F_iG_i$。$\square$

------

## 4.3 slab 内条件采样：先采 $y$，再两次 stabbing

本节构造从 slab 条件分布 $\mathcal D_{\ell,i}$ 采样的过程。关键点在于：利用一个辅助变量 $y$，把二维权重采样分解为一次连续采样与两次一维 stabbing 采样。

### 4.3.1 条件分布的分解形态

在 slab $S_i$ 中，$W_{\ell,i}(r,s)$ 可写为
$$
W_{\ell,i}(r,s)=\Delta x_i\cdot h(r,s)\cdot w_\ell(r)\cdot v_\ell(s)\cdot \chi_r(i)\chi_s(i).
$$
由引理 4.2，$M_i=\Delta x_i\int F_iG_i$。

令 $Y$ 的密度为 $p_i(y)\propto F_i(y)G_i(y)$（由引理 4.3），并定义给定 $y$ 的红侧与蓝侧条件分布：
$$
\mathcal P_r(r\mid y)=\frac{w_\ell(r)\,\mathbf 1[\chi_r(i)=1]\,\mathbf 1[y\in r^y]}{F_i(y)},
\qquad
\mathcal P_b(s\mid y)=\frac{v_\ell(s)\,\mathbf 1[\chi_s(i)=1]\,\mathbf 1[y\in s^y]}{G_i(y)}.
$$
当 $p_i(y)>0$ 时有 $F_i(y)>0$ 且 $G_i(y)>0$，因此上述条件分布均良定义，并且对应的覆盖集合均非空；算法不会在空覆盖集合上调用 stabbing。

------

### 4.3.2 slab 内条件采样过程

**算法 4.1（slab 内条件采样）**
 输入：当前 slab $S_i$ 的数据结构状态（乘积线段树 + 红 stabbing + 蓝 stabbing）。
 输出：一对 $(r,s)$。

1. $y\leftarrow \mathrm{SampleY}()$（按密度 $p_i(y)\propto F_i(y)G_i(y)$）；
2. $r\leftarrow \mathrm{RedStab.Sample}(y)$（从覆盖 $y$ 的红矩形中按权重 $w_\ell$ 抽样）；
3. $s\leftarrow \mathrm{BlueStab.Sample}(y)$（从覆盖 $y$ 的蓝矩形中按权重 $v_\ell$ 抽样）；
4. 输出 $(r,s)$。

上面两次 stabbing 使用独立随机性，因此在给定 $y$ 时 $r,s$ 条件独立。

------

### 4.3.3 精确 stabbing 下的分布正确性

**定理 4.4（精确 stabbing 时的 slab 条件分布）**
 若红/蓝两侧 stabbing 在任意查询点 $y$ 上都返回精确的加权抽样结果（分别为 $\mathcal P_r(\cdot\mid y)$ 与 $\mathcal P_b(\cdot\mid y)$），则算法 4.1 输出分布为 $\mathcal D_{\ell,i}$。

*证明*：对任意 $r\in R_c,s\in R_{\bar c}$，由全概率公式与条件独立性，
$$
\Pr[(r,s)]
=\int_{y_0}^{y_q} p_i(y)\,\mathcal P_r(r\mid y)\,\mathcal P_b(s\mid y)\,dy.
$$
代入
$$
p_i(y)=\frac{F_i(y)G_i(y)}{\int F_iG_i},\quad
\mathcal P_r(r\mid y)=\frac{w_\ell(r)\mathbf 1[\chi_r(i)=1]\mathbf 1[y\in r^y]}{F_i(y)},\quad
\mathcal P_b(s\mid y)=\frac{v_\ell(s)\mathbf 1[\chi_s(i)=1]\mathbf 1[y\in s^y]}{G_i(y)},
$$
可得分母 $F_i(y)$ 与 $G_i(y)$ 抵消，得到
$$
\Pr[(r,s)]
=\frac{w_\ell(r)v_\ell(s)\mathbf 1[\chi_r(i)=1]\mathbf 1[\chi_s(i)=1]}{\int F_iG_i}
\int_{y_0}^{y_q}\mathbf 1[y\in r^y]\mathbf 1[y\in s^y]\,dy.
$$
积分项为 $h(r,s)$。再用引理 4.2 的 $M_i=\Delta x_i\int F_iG_i$，得到
$$
\Pr[(r,s)]=\frac{\Delta x_i\,h(r,s)\,w_\ell(r)v_\ell(s)\,\chi_r(i)\chi_s(i)}{M_i}
=\frac{W_{\ell,i}(r,s)}{M_i}=\mathcal D_{\ell,i}(r,s).
$$
$\square$

------

### 4.3.4 近似 stabbing 下的误差合成

本节把一维 stabbing 的 TV-近似误差合并到 slab 条件分布上。

记对固定 $y$：

- 红侧真实分布为 $\mathcal P_r(\cdot\mid y)$，近似分布为 $\tilde{\mathcal P}_r(\cdot\mid y)$；
- 蓝侧真实分布为 $\mathcal P_b(\cdot\mid y)$，近似分布为 $\tilde{\mathcal P}_b(\cdot\mid y)$。

假设对所有满足 $p_i(y)>0$ 的 $y$，有
$$
d_{\mathrm{TV}}(\tilde{\mathcal P}_r(\cdot\mid y),\mathcal P_r(\cdot\mid y))\le \varepsilon_{\mathrm{stab}},
\qquad
d_{\mathrm{TV}}(\tilde{\mathcal P}_b(\cdot\mid y),\mathcal P_b(\cdot\mid y))\le \varepsilon_{\mathrm{stab}}.
$$
**引理 4.5（给定 $y$ 的红蓝误差加法）**
 若红/蓝采样使用独立随机性，则对联合分布有
$$
d_{\mathrm{TV}}\Big(\tilde{\mathcal P}_r(\cdot\mid y)\times \tilde{\mathcal P}_b(\cdot\mid y),\ 
\mathcal P_r(\cdot\mid y)\times \mathcal P_b(\cdot\mid y)\Big)
\le 2\varepsilon_{\mathrm{stab}}.
$$
*证明*：对任意有限域上分布 $P,P',Q,Q'$，有
$$
\|P\times Q - P'\times Q'\|_1
\le \|P-P'\|_1+\|Q-Q'\|_1.
$$
两边除以 $2$ 即得 TV 上界。$\square$

**引理 4.6（对 $y$ 的混合不放大 TV）**
 设 $p(y)$ 是任意密度（或离散混合权重），${\mathcal Q_y}$ 与 ${\tilde{\mathcal Q}_y}$ 是同一域上的分布族，且对所有 $y$ 有
$$
d_{\mathrm{TV}}(\tilde{\mathcal Q}_y,\mathcal Q_y)\le \eta,
$$
则
$$
d_{\mathrm{TV}}\Big(\int p(y)\tilde{\mathcal Q}_y\,dy,\ \int p(y)\mathcal Q_y\,dy\Big)\le \eta.
$$
*证明*：用 $\ell_1$ 形式与凸性：
$$
\left\|\int p(y)(\tilde{\mathcal Q}_y-\mathcal Q_y)\,dy\right\|_1
\le \int p(y)\|\tilde{\mathcal Q}_y-\mathcal Q_y\|_1\,dy
\le \int p(y)\cdot 2\eta\,dy
=2\eta.
$$
两边除以 $2$ 得结论。$\square$

由引理 4.5 与引理 4.6，得到 slab 条件分布的误差：

**定理 4.7（slab 条件分布的 TV 误差）**
 令 $\tilde{\mathcal D}*{\ell,i}$ 表示算法 4.1 在使用近似 stabbing 时的输出分布。若两侧 stabbing 的逐点误差均不超过 $\varepsilon*{\mathrm{stab}}$，则
$$
d_{\mathrm{TV}}(\tilde{\mathcal D}_{\ell,i},\mathcal D_{\ell,i})\le 2\varepsilon_{\mathrm{stab}}.
$$
*证明*：对每个 $y$ 的联合误差由引理 4.5 上界为 $2\varepsilon_{\mathrm{stab}}$，再以 $p_i(y)$ 混合并用引理 4.6 得到同样上界。$\square$

------

## 4.4 固定 $\ell$ 的两遍扫描回填：i.i.d. 采样器组织

本节将 slab 混合与 slab 内条件采样组织成一个固定 $\ell$ 的 i.i.d. 采样器。整体结构分为 Pass1（仅统计质量）、Planning（独立指派 slab 并分组）、Pass2（条件采样并回填）。

### 4.4.1 预处理：事件表、端点骨架与 stabbing 内存布局

固定 $\ell$ 后，以下对象在整个两遍扫描过程中视为静态：

- 所有矩形的几何端点（用于生成 slab 与 $y$ 基本段）；
- 红侧权重 ${w_\ell(r)}*{r\in R_c}$ 与蓝侧权重 ${v*\ell(s)}*{s\in R*{\bar c}}$；
- stabbing 结构所需的段树骨架与 canonical cover 归属关系；
- stabbing 桶内为确定最坏更新所需的 membership records、按 $(\text{bucket},\text{level})$ 的分组与数组段容量分配。

预处理阶段执行：

1. 生成并排序 $X$ 与 $Y$，建立 slab 列表与 $y$ 基本段；
2. 构建乘积线段树骨架（覆盖 $[y_0,y_q)$）；
3. 为红侧与蓝侧分别构建一套动态加权 stabbing 结构：
   - 为每个矩形生成其规范覆盖（canonical cover）对应的 membership records；
   - 根据权重计算 dyadic level（例如 $k=\lceil\log_2 w\rceil$）并按键 $(\text{bucket},k)$ 分组；
   - 为每个桶-层分配定长数组段，并为每条 record 写回段指针；
   - 为每个矩形写入定长的 membership 列表，以便运行期插删仅遍历该列表，不做查找。

该布局保证运行期每条 membership 的插入/删除为确定 $O(1)$，从而一次矩形的插入/删除（触达 $O(\log N)$ 个 canonical bucket）为确定 $O(\log N)$。

------

### 4.4.2 Pass1：仅统计 slab 质量

第一遍扫描仅维护乘积线段树（不维护 stabbing），以便得到每个 slab 的质量：
$$
M_i=\Delta x_i\cdot \mathrm{SP}_i,\qquad \mathrm{SP}_i=\int_{y_0}^{y_q}F_i(y)G_i(y)\,dy.
$$
扫描过程如下：

- 初始时 $F\equiv 0,G\equiv 0$，线段树的 $\mathrm{SF},\mathrm{SG},\mathrm{SP}$ 全为 $0$；
- 对每个 $x=x_i$，先处理 `END` 事件再处理 `START` 事件：
  - 若红矩形 $r$ 插入，则对区间 $r^y$ 执行 $F\leftarrow F+w_\ell(r)$；
  - 若红矩形 $r$ 删除，则执行 $F\leftarrow F-w_\ell(r)$；
  - 蓝矩形 $s$ 同理在 $G$ 上做区间加减 $v_\ell(s)$；
- 处理完 $x_i$ 处事件后，当前状态对应 slab $S_i=[x_i,x_{i+1})$，此时根节点 $\mathrm{SP}(\mathrm{root})$ 即为 $\int F_iG_i$，从而 $M_i=\Delta x_i\cdot \mathrm{SP}(\mathrm{root})$。

完成所有 slab 后得到 $Z_\ell=\sum_i M_i$。若 $Z_\ell=0$，返回空序列。

------

### 4.4.3 Planning：独立指派 slab 索引并分组

当 $Z_\ell>0$ 时，建立 slab 离散分布
$$
\Pr[I=i]=\frac{M_i}{Z_\ell},\qquad i=0,\dots,m-1.
$$
对每个输出位置 $j=1,\dots,t$ 独立抽取 $I_j$，并按 slab 分组为列表
$$
L[i]=\{j\in[t]: I_j=i\}.
$$
该阶段使用独立随机性，与 Pass2 的条件采样随机性相互独立。

------

### 4.4.4 Pass2：维护乘积线段树 + 两套 stabbing 并回填

第二遍扫描从空结构开始，同时维护三类对象：

- 乘积线段树（维护 $F_i,G_i,\int F_iG_i$ 并支持 `SampleY()`）；
- 红 stabbing（动态维护活跃红矩形的 $y$ 区间，权重为 $w_\ell(r)$）；
- 蓝 stabbing（动态维护活跃蓝矩形的 $y$ 区间，权重为 $v_\ell(s)$）。

扫描到 $x=x_i$ 时，仍按 `END` 再 `START` 更新活跃集合：

- 对乘积线段树做区间加减更新 $F,G$；
- 对两套 stabbing 分别做对应矩形的 `Insert/Delete`（每次更新触达 $O(\log N)$ 个 bucket，每个 bucket 的 membership 更新为确定 $O(1)$）。

处理完 $x_i$ 处事件后进入 slab $S_i$。若 $k=|L[i]|>0$，则独立重复算法 4.1 共 $k$ 次生成 $k$ 个样本，并按 $L[i]$ 指定位置回填到输出数组。

------

### 4.4.5 i.i.d. 组织与随机性分离

本章默认如下随机性分离：

- Planning 阶段生成 ${I_j}_{j=1}^t$ 的随机性相互独立；
- Pass2 中每次调用算法 4.1 的随机性相互独立；
- `SampleY()` 的随机性与两次 stabbing 的随机性相互独立；
- 红/蓝两次 stabbing 在同一次条件采样中使用独立随机性。

在这些条件下，回填只是对一族独立随机变量的确定性重标号，因此输出样本相互独立。

------

## 4.5 固定 $\ell$ 的误差与复杂度定理（账本结论）

本节给出固定 $\ell$ 子采样器的误差与最坏复杂度。误差仅来自 stabbing 的 TV-近似；乘积线段树的质量统计与 `SampleY()` 在本章模型下为精确实现。

### 4.5.1 stabbing 子模块参数与逐点误差

对每种颜色各维护一套动态加权 stabbing 采样器。其逐点 TV 误差由两参数刻画：

- 裁剪参数 $\varepsilon_{\mathrm{clip}}\in(0,1)$；
- 桶内有界拒绝次数 $M_{\mathrm{bucket}}\in\mathbb N$。

定义 stabbing 的逐点 TV 误差上界
$$
\varepsilon_{\mathrm{stab}}=\varepsilon_{\mathrm{clip}}+2^{-M_{\mathrm{bucket}}}.
$$
该上界对所有查询点 $y\in[y_0,y_q)$ 成立；在本章的调用中，实际查询的 $y$ 总满足 $p_i(y)>0$，因此不会触发空覆盖集合的情形。

------

### 4.5.2 固定 $\ell$ 的正确性与 TV 误差

记 Pass2 中 slab $i$ 的近似条件分布为 $\tilde{\mathcal D}_{\ell,i}$，并定义整体输出的单次边缘分布为
$$
\tilde{\mathcal D}_\ell
=\sum_{i=0}^{m-1}\frac{M_i}{Z_\ell}\tilde{\mathcal D}_{\ell,i}.
$$
真实目标为
$$
\mathcal D_\ell
=\sum_{i=0}^{m-1}\frac{M_i}{Z_\ell}\mathcal D_{\ell,i}.
$$
**定理 4.8（固定 $\ell$ 的分布与 TV 误差）**
 当 $Z_\ell>0$ 时，两遍扫描回填算法输出 $t$ 个 i.i.d. 样本，其单次边缘分布为 $\tilde{\mathcal D}_\ell$，并且
$$
d_{\mathrm{TV}}(\tilde{\mathcal D}_\ell,\mathcal D_\ell)\le 2\varepsilon_{\mathrm{stab}}
=2\big(\varepsilon_{\mathrm{clip}}+2^{-M_{\mathrm{bucket}}}\big).
$$
*证明*：对每个 slab $i$，由定理 4.7 得
 $d_{\mathrm{TV}}(\tilde{\mathcal D}*{\ell,i},\mathcal D*{\ell,i})\le 2\varepsilon_{\mathrm{stab}}$。
 混合系数 ${M_i/Z_\ell}$ 精确且相同，因此对 slab 混合使用 TV 的凸性可得
$$
d_{\mathrm{TV}}(\tilde{\mathcal D}_\ell,\mathcal D_\ell)
\le \sum_i \frac{M_i}{Z_\ell}\cdot 2\varepsilon_{\mathrm{stab}}
=2\varepsilon_{\mathrm{stab}}.
$$
i.i.d. 性由第 4.4.5 节的独立性组织得到。$\square$

------

### 4.5.3 固定 $\ell$ 的最坏时间复杂度

记 $N=|R_c|+|R_{\bar c}|$，$m\le 2N-1$ 为 slab 数，$q\le 2N-1$ 为 $y$ 基本段数，因此 $\log m,\log q=O(\log N)$。

为强调计费口径，将总时间分为“预处理、两遍扫描更新、采样生成”三部分。

**定理 4.9（固定 $\ell$ 的最坏时间账本）**
 固定 $\ell$ 且需要 $t$ 个样本。上述两遍扫描回填采样器的最坏时间满足：
$$
T_\ell(N,t)
=
T_{\mathrm{prep}}
+
T_{\mathrm{Pass1}}
+
T_{\mathrm{plan}}
+
T_{\mathrm{Pass2\_update}}
+
T_{\mathrm{sample}},
$$
其中：

1. **预处理**

- 构造 $X,Y$ 并排序：$O(N\log N)$；
- 建立 slab 与事件表：$O(N)$；
- 建立乘积线段树骨架：$O(N)$；
- 为红/蓝两侧各构建一套 stabbing 内存布局：membership records 总数为 $P=O(N\log N)$，分组与数组段分配可在
  - 通用比较分组口径下用 $O(P\log P)$；
  - 当键可装入 $O(1)$ 个机器字时，用确定的 radix 分组实现 $O(P)$；
     因此 $T_{\mathrm{prep}}=O(N\log N)$（在 $O(P)$ 分组口径下）或 $T_{\mathrm{prep}}=O(P\log P)$（通用口径）。

1. **Pass1（mass-only）**
    每个矩形触发一次 `START` 与一次 `END`，每次对乘积线段树做一次区间加法更新 $O(\log N)$，因此

$$
T_{\mathrm{Pass1}}=O(N\log N).
$$

1. **Planning（slab 指派）**
    建立 slab 离散采样器并独立抽取 $t$ 次，时间

$$
T_{\mathrm{plan}}=O(m+t)=O(N+t).
$$

1. **Pass2 更新**

- 乘积线段树更新：同 Pass1，总计 $O(N\log N)$；
- stabbing 更新：每次矩形插入/删除触达 $O(\log N)$ 个 canonical bucket，每个 bucket 的 membership 更新为确定 $O(1)$，因此每次事件 $O(\log N)$，总计

$$
T_{\mathrm{Pass2\_update}}=O(N\log N).
$$

1. **样本生成**
    每个样本包含：

- 一次 `SampleY()`：沿线段树下降，$O(\log N)$；
- 两次 stabbing 查询：每次由“Path 选桶 + 单桶窗口扫描 + 小 alias + 有界拒绝”组成，其确定最坏时间为

$$
O\Big(\log N+\log(N/\varepsilon_{\mathrm{clip}})+M_{\mathrm{bucket}}\Big).
$$

因此
$$
T_{\mathrm{sample}}
=
O\Big(
t\big(\log N+\log(N/\varepsilon_{\mathrm{clip}})+M_{\mathrm{bucket}}\big)
\Big).
$$
在采用 $T_{\mathrm{prep}}=O(N\log N)$ 的分组口径时，总体最坏时间可写为
$$
T_\ell(N,t)
=
O\Big(
N\log N
+t\big(\log N+\log(N/\varepsilon_{\mathrm{clip}})+M_{\mathrm{bucket}}\big)
\Big).
$$
该界中扫描更新部分不包含 $\log^2 N$ 项。

------

### 4.5.4 空间复杂度

固定 $\ell$ 的运行时空间由以下部分构成：

1. 乘积线段树：$O(q)=O(N)$ 节点，每节点常数域；
2. 两套 stabbing：membership records 总数 $P=O(N\log N)$，桶-层数组段按预处理精确分配，总空间 $O(N\log N)$；
3. slab 质量数组 ${M_i}$ 与分组列表 ${L[i]}$：$O(N+t)$；
4. 输出样本数组：$O(t)$。

因此空间上界为
$$
S_\ell(N,t)=O(N\log N+t).
$$

---

# 第 5 章 端到端装配：多层回填、独立性、TV 闭环与最坏复杂度

本章把前述构件组装成一个端到端的 IoU 采样器。整体结构是：

- **提案端**：以 $\rho(P)=I(P)/(A(P)+B(P))$ 为常数因子提案，并用正系数指数和把 $\rho$ 写成 $L$ 个非负可分离分量的混合；
- **分量端**：对每个分量 $\ell$，调用固定 $\ell$ 的几何采样器生成 i.i.d. 提案；
- **外层校正**：用接受概率 $a(P)=1/[2(1-\rho(P))]$ 将提案分布映射到 $\mathcal D_{\mathrm{IoU}}$；为获得确定最坏时间，使用有界拒绝（每个输出最多尝试 $M$ 次）。

本章给出两类结论：

1) **分布正确性与 TV 误差账本**：把端到端的总变差距离拆为四类可加项，并保持外层重加权算子的 TV-Lipschitz 常数为 $4$，从而误差闭环稳定。

2) **确定最坏复杂度**：给出端到端时间的逐项加和式，明确哪些代价来自共享几何骨架、哪些来自对所有 $\ell$ 的质量统计、哪些来自仅对被使用到的 $\ell$ 的条件采样扫描，以及外层拒绝的开销；同时给出两种等价的内存组织选项与对应的空间界。

---

## 5.0 记号与基本约定（端到端口径）

样本空间为跨类有序对集合 $\Omega=R_c\times R_{\bar c}$。对 $P=(r,s)\in\Omega$ 记
- $A(P)=A(r)$，$B(P)=A(s)$，
- $I(P)=I(r,s)$ 为交面积。

定义正权重支撑集
$$
J=\{P\in\Omega:\ I(P)>0\}.
$$
若 $J=\varnothing$（等价于所有跨类对均无正交面积），则 IoU 总权重为 $0$，采样任务按约定输出空序列。以下默认 $J\neq\varnothing$。

对任意 $P\in J$ 定义
$$
\rho(P)=\frac{I(P)}{A(P)+B(P)}\in\Big(0,\frac12\Big],
\qquad
W_{\mathrm{IoU}}(P)=\frac{I(P)}{A(P)+B(P)-I(P)}=\frac{\rho(P)}{1-\rho(P)}.
$$
所有分布的总变差距离（TV）均以有限域 $J$ 为支撑计算：
$$
d_{\mathrm{TV}}(\mathcal P,\mathcal Q)=\frac12\sum_{P\in J}|\mathcal P(P)-\mathcal Q(P)|.
$$

---

## 5.1 $\ell$ 上混合与提案序列的 i.i.d. 组织

### 5.1.1 近似提案的分量混合形式

提案目标分布为 $\mathcal D_\rho(P)\propto \rho(P)$。为实现可分离结构，在区间 $x\in[S_{\min},S_{\max}]$ 上采用正系数指数和相对逼近
$$
\left|\frac{1}{x}-\sum_{\ell=1}^L \omega_\ell e^{-\alpha_\ell x}\right|\le \frac{\delta}{x},
\qquad \omega_\ell>0,\ \alpha_\ell>0.
$$
于是定义近似提案的未归一化权重（在 $J$ 上）
$$
\tilde w_\rho(P)=I(P)\sum_{\ell=1}^L \omega_\ell e^{-\alpha_\ell(A(P)+B(P))}
=\sum_{\ell=1}^L W_\ell(P),
$$
其中分量权重为
$$
W_\ell(r,s)=I(r,s)\,w_\ell(r)\,v_\ell(s),
\qquad
w_\ell(r)=e^{-\alpha_\ell A(r)},\quad
v_\ell(s)=\omega_\ell e^{-\alpha_\ell A(s)}.
$$
对每个 $\ell$ 定义分量总质量与总质量
$$
Z_\ell=\sum_{P\in J}W_\ell(P),
\qquad
Z=\sum_{\ell=1}^L Z_\ell=\sum_{P\in J}\tilde w_\rho(P).
$$
由于 $J\neq\varnothing$ 且每项非负，可有 $Z>0$；同时允许存在某些 $Z_\ell=0$（该分量在 $J$ 上处处为 $0$）。

对 $Z_\ell>0$ 的分量定义条件分布与混合系数
$$
\mathcal D_\ell(P)=\frac{W_\ell(P)}{Z_\ell},
\qquad
\pi_\ell=\frac{Z_\ell}{Z},
\qquad
\sum_{\ell:Z_\ell>0}\pi_\ell=1.
$$
则近似提案分布可写为分量混合
$$
\tilde{\mathcal D}_\rho(P)=\frac{\tilde w_\rho(P)}{\sum_{Q\in J}\tilde w_\rho(Q)}
=\sum_{\ell:Z_\ell>0}\pi_\ell\,\mathcal D_\ell(P).
$$

---

### 5.1.2 多层回填：分量指派与分量内批量生成

外层有界拒绝中，每个 IoU 输出位置最多尝试 $M$ 次，因此总提案数固定为
$$
t_{\mathrm{prop}}=t\cdot M.
$$

提案序列用“两层独立指派 + 分组回填”的方式生成：

**层 1（分量指派）**  
对每个提案位置 $j\in[t_{\mathrm{prop}}]$，独立抽取分量索引
$$
L_j\sim \{\hat\pi_\ell\}_{\ell=1}^L,
$$
其中 $\hat\pi$ 是实现中使用的混合系数（可等于真实 $\pi$，也可为其近似）。令
$$
t_\ell=|\{j:\ L_j=\ell\}|,
\qquad
J_\ell=\{j\in[t_{\mathrm{prop}}]:\ L_j=\ell\}
$$
为计数与位置列表。

**层 2（分量内批量生成并回填）**  
对每个满足 $t_\ell>0$ 的分量 $\ell$，调用“固定 $\ell$ 分量采样器”一次性生成 $t_\ell$ 个相互独立样本
$$
P^{(\ell)}_1,\dots,P^{(\ell)}_{t_\ell}\ \sim\ \widetilde{\mathcal D}_\ell,
$$
再按列表 $J_\ell$ 的顺序把这些样本回填到提案数组 $P_1,\dots,P_{t_{\mathrm{prop}}}$ 的对应位置。

固定 $\ell$ 分量采样器自身再采用 slab 两遍扫描回填（先指派 slab、再在 slab 内条件采样并回填），因此整体提案生成过程是“分量—slab”的两级回填组合，但所有指派变量都在外层先生成且相互独立。

---

### 5.1.3 随机性流分离与 i.i.d. 充分条件

为避免实现层面引入隐式相关性，采用可分离的随机性流（RNG streams），并要求以下四类随机源两两独立：

- $\mathsf{RNG}_{\ell\text{-assign}}$：生成 $\{L_j\}_{j=1}^{t_{\mathrm{prop}}}$；
- $\mathsf{RNG}_{\text{slab-assign}}[\ell]$：固定 $\ell$ 的 slab 指派；
- $\mathsf{RNG}_{\text{cond}}[\ell]$：固定 $\ell$ 的 slab 内条件采样（包含 `SampleY()` 与红/蓝 stabbing 的随机位）；
- $\mathsf{RNG}_{\text{accept}}$：外层有界拒绝的均匀随机数与比较。

此外，所有循环均为确定上界（外层 $M$、桶内 $M_{\mathrm{bucket}}$ 等），因此每条随机流的调用序列长度对任意输入与任意随机结果都是有界且可预先上界的。

**引理 5.1（两层回填生成的提案序列 i.i.d.）**  
设对每个 $\ell$，分量采样器返回 $t_\ell$ 个相互独立、边缘分布为 $\widetilde{\mathcal D}_\ell$ 的样本；并且
1) $\{L_j\}$ 相互独立，且 $\Pr[L_j=\ell]=\hat\pi_\ell$；
2) 各 $\ell$ 的分量采样过程所用随机性相互独立，且与 $\{L_j\}$ 独立；  
则回填得到的提案序列 $P_1,\dots,P_{t_{\mathrm{prop}}}$ 相互独立同分布，且单次边缘分布为
$$
\widehat{\mathcal D}_\rho=\sum_{\ell: t_\ell>0}\hat\pi_\ell\,\widetilde{\mathcal D}_\ell.
$$

**证明**  
对任意位置 $j$ 与任意 $P\in J$，由全概率公式
$$
\Pr[P_j=P]=\sum_{\ell}\Pr[L_j=\ell]\Pr[P_j=P\mid L_j=\ell]
=\sum_{\ell}\hat\pi_\ell\,\widetilde{\mathcal D}_\ell(P).
$$
独立性：给定 $\{L_j\}$ 后，各分量内部样本独立、不同分量的样本族独立；回填是确定性重标号，不产生新的相关性，因此 $P_1,\dots,P_{t_{\mathrm{prop}}}$ 相互独立同分布。$\square$

---

## 5.2 外层有界拒绝：从提案到 IoU

### 5.2.1 接受概率与重加权算子

由恒等式 $W_{\mathrm{IoU}}(P)=\rho(P)/(1-\rho(P))$，以 $2\rho$ 作为常数包络，定义接受概率
$$
a(P)=\frac{W_{\mathrm{IoU}}(P)}{2\rho(P)}=\frac{1}{2(1-\rho(P))}.
$$
因为 $\rho(P)\le 1/2$，所以对所有 $P\in J$ 有
$$
a(P)\in\left[\frac12,1\right].
$$

对任意分布 $\mathcal P$（定义在 $J$ 上）定义重加权算子
$$
\mathcal T_a(\mathcal P)(P)=\frac{\mathcal P(P)\,a(P)}{\sum_{Q\in J}\mathcal P(Q)\,a(Q)}.
$$
当 $\mathcal P=\mathcal D_\rho$ 时，上式满足 $\mathcal T_a(\mathcal D_\rho)=\mathcal D_{\mathrm{IoU}}$，即“从 $\rho$ 校正到 IoU”。

---

### 5.2.2 有界拒绝算子与 TV 误差

经典拒绝采样“直到接受”为止能精确得到 $\mathcal T_a(\mathcal P)$，但运行时间随机。这里采用确定最坏时间的有界版本：每个输出最多尝试 $M$ 次，若全部拒绝则输出第 $M$ 次提案。

形式化地，给定提案分布 $\mathcal P$，定义 $\mathcal T_{a,M}(\mathcal P)$ 为如下过程的输出分布：

- 独立提案 $P_1,\dots,P_M\sim \mathcal P$；
- 独立均匀随机数 $U_1,\dots,U_M\sim\mathrm{Unif}(0,1)$；
- 输出 $P_\tau$，其中 $\tau$ 是最小满足 $U_\tau\le a(P_\tau)$ 的下标；若不存在，则输出 $P_M$。

**引理 5.2（外层有界拒绝的 TV 误差）**  
对任意提案分布 $\mathcal P$，有
$$
d_{\mathrm{TV}}\big(\mathcal T_{a,M}(\mathcal P),\mathcal T_a(\mathcal P)\big)\le 2^{-M}.
$$

**证明**  
由于 $a(P)\ge 1/2$，每次尝试的接受概率至少为 $1/2$。令 $\mathsf{Fail}$ 表示前 $M$ 次均拒绝，则
$$
\Pr[\mathsf{Fail}]\le (1-1/2)^M=2^{-M}.
$$
有界拒绝的输出分布可以写成“成功时的无界拒绝输出分布”的凸组合加上失败分支分布，因此与无界拒绝分布的 TV 距离至多为失败概率，得到结论。$\square$

---

### 5.2.3 外层 i.i.d. 性

外层对每个输出位置使用 $M$ 个独立提案与独立均匀随机数，并且不同位置之间使用独立随机性。于是：

**推论 5.2.1（外层有界拒绝保持 i.i.d.）**  
若提案序列按引理 5.1 是 i.i.d.，并且外层接受随机性对每次尝试均独立且与提案独立，则最终输出 $t$ 个样本相互独立同分布，其单次边缘分布为 $\mathcal T_{a,M}(\widehat{\mathcal D}_\rho)$。

---

## 5.3 端到端 TV 闭环（保持 Lipschitz 常数 4）

本节给出端到端的 TV 上界。核心分为两部分：
- **提案侧误差**：$\widehat{\mathcal D}_\rho$ 与理想提案 $\mathcal D_\rho$ 的 TV 距离；
- **外层传播**：有界拒绝引入 $2^{-M}$，重加权算子把提案侧 TV 误差放大至多常数 $4$。

---

### 5.3.1 提案侧误差的三项分解

理想提案分布为
$$
\mathcal D_\rho(P)=\frac{\rho(P)}{\sum_{Q\in J}\rho(Q)}.
$$

**(1) 核展开误差：$\mathcal D_\rho \to \tilde{\mathcal D}_\rho$**  
由指数和相对误差与“乘性权重误差 $\Rightarrow$ TV”的标准不等式，可得
$$
d_{\mathrm{TV}}(\mathcal D_\rho,\tilde{\mathcal D}_\rho)\le \frac{\delta}{1-\delta}.
$$

**(2) 混合系数误差：$\pi\to\hat\pi$（可选）**  
若实现中采用 $\hat\pi$ 代替真实 $\pi$，定义
$$
\varepsilon_{\mathrm{mix}}=\frac12\sum_{\ell=1}^L|\hat\pi_\ell-\pi_\ell|.
$$
则对同一组分量分布 $\{\mathcal D_\ell\}$，有
$$
d_{\mathrm{TV}}\Big(\sum_\ell \hat\pi_\ell\mathcal D_\ell,\ \sum_\ell \pi_\ell\mathcal D_\ell\Big)\le \varepsilon_{\mathrm{mix}}.
$$
若 $\pi$ 精确可得，则取 $\varepsilon_{\mathrm{mix}}=0$。

**(3) 分量几何误差：$\mathcal D_\ell \to \widetilde{\mathcal D}_\ell$**  
固定 $\ell$ 的几何采样器返回分布 $\widetilde{\mathcal D}_\ell$，并满足
$$
d_{\mathrm{TV}}(\widetilde{\mathcal D}_\ell,\mathcal D_\ell)\le \varepsilon_{\mathrm{geom}},
\qquad
\varepsilon_{\mathrm{geom}}:=2\,\varepsilon_{\mathrm{stab}},
\qquad
\varepsilon_{\mathrm{stab}}=\varepsilon_{\mathrm{clip}}+2^{-M_{\mathrm{bucket}}}.
$$

---

### 5.3.2 分量混合下的 TV 不放大

实际提案分布为
$$
\widehat{\mathcal D}_\rho=\sum_{\ell}\hat\pi_\ell\,\widetilde{\mathcal D}_\ell,
\qquad
\tilde{\mathcal D}_\rho=\sum_{\ell}\pi_\ell\,\mathcal D_\ell.
$$

**引理 5.3（分量混合下 TV 的可组合性）**  
若对所有 $\ell$ 有 $d_{\mathrm{TV}}(\widetilde{\mathcal D}_\ell,\mathcal D_\ell)\le \varepsilon_{\mathrm{geom}}$，则
$$
d_{\mathrm{TV}}\Big(\sum_\ell \pi_\ell\widetilde{\mathcal D}_\ell,\ \sum_\ell \pi_\ell\mathcal D_\ell\Big)\le \varepsilon_{\mathrm{geom}}.
$$
进一步，允许混合系数扰动时有
$$
d_{\mathrm{TV}}(\widehat{\mathcal D}_\rho,\tilde{\mathcal D}_\rho)
\le \varepsilon_{\mathrm{mix}}+\varepsilon_{\mathrm{geom}}.
$$

**证明**  
第一式由 TV 的凸性：
$$
d_{\mathrm{TV}}\Big(\sum_\ell \pi_\ell\widetilde{\mathcal D}_\ell,\ \sum_\ell \pi_\ell\mathcal D_\ell\Big)
\le \sum_\ell \pi_\ell\,d_{\mathrm{TV}}(\widetilde{\mathcal D}_\ell,\mathcal D_\ell)
\le \sum_\ell \pi_\ell\,\varepsilon_{\mathrm{geom}}=\varepsilon_{\mathrm{geom}}.
$$
第二式把差异分解为“系数差异”与“分量分布差异”，并对系数差异使用 $\ell_1/2$ 上界即可。$\square$

---

### 5.3.3 外层重加权的 Lipschitz 常数

**引理 5.4（重加权算子的 TV-Lipschitz 常数）**  
当 $a(P)\in[1/2,1]$ 时，对任意两分布 $\mathcal P,\mathcal Q$（定义在 $J$ 上），
$$
d_{\mathrm{TV}}\big(\mathcal T_a(\mathcal P),\mathcal T_a(\mathcal Q)\big)\le 4\,d_{\mathrm{TV}}(\mathcal P,\mathcal Q).
$$

**证明**  
令未归一化向量 $\mu_{\mathcal P}(P)=\mathcal P(P)a(P)$，$\mu_{\mathcal Q}(P)=\mathcal Q(P)a(P)$。归一化常数为
$$
\|\mu_{\mathcal P}\|_1=\mathbb E_{\mathcal P}[a]\ge 1/2,\qquad
\|\mu_{\mathcal Q}\|_1=\mathbb E_{\mathcal Q}[a]\ge 1/2.
$$
归一化的标准不等式给出
$$
\left\|\frac{\mu_{\mathcal P}}{\|\mu_{\mathcal P}\|_1}-\frac{\mu_{\mathcal Q}}{\|\mu_{\mathcal Q}\|_1}\right\|_1
\le \frac{2\|\mu_{\mathcal P}-\mu_{\mathcal Q}\|_1}{\min\{\|\mu_{\mathcal P}\|_1,\|\mu_{\mathcal Q}\|_1\}}
\le 4\|\mu_{\mathcal P}-\mu_{\mathcal Q}\|_1.
$$
并且
$$
\|\mu_{\mathcal P}-\mu_{\mathcal Q}\|_1
=\sum_{P\in J}a(P)|\mathcal P(P)-\mathcal Q(P)|
\le \sum_{P\in J}|\mathcal P(P)-\mathcal Q(P)|
=2d_{\mathrm{TV}}(\mathcal P,\mathcal Q).
$$
代入并除以 $2$ 得到常数 $4$。$\square$

---

### 5.3.4 端到端 TV 闭环定理

**定理 5.5（端到端 TV 闭环）**  
设端到端算法的单次输出边缘分布为 $\tilde{\mathcal D}$，理想 IoU 分布为 $\mathcal D_{\mathrm{IoU}}$。若提案端实现得到 $\widehat{\mathcal D}_\rho$，外层使用 $M$ 次有界拒绝，则
$$
d_{\mathrm{TV}}(\tilde{\mathcal D},\mathcal D_{\mathrm{IoU}})
\le 2^{-M}+4\left(\frac{\delta}{1-\delta}+\varepsilon_{\mathrm{mix}}+\varepsilon_{\mathrm{geom}}\right),
$$
其中 $\varepsilon_{\mathrm{geom}}=2(\varepsilon_{\mathrm{clip}}+2^{-M_{\mathrm{bucket}}})$。

**证明**  
记无界校正算子为 $\mathcal T_a$，有界校正为 $\mathcal T_{a,M}$。由三角不等式：
$$
d_{\mathrm{TV}}(\mathcal T_{a,M}(\widehat{\mathcal D}_\rho),\mathcal T_a(\mathcal D_\rho))
\le d_{\mathrm{TV}}(\mathcal T_{a,M}(\widehat{\mathcal D}_\rho),\mathcal T_a(\widehat{\mathcal D}_\rho))
+ d_{\mathrm{TV}}(\mathcal T_a(\widehat{\mathcal D}_\rho),\mathcal T_a(\mathcal D_\rho)).
$$
第一项由引理 5.2 得 $\le 2^{-M}$。第二项由引理 5.4 得
$$
d_{\mathrm{TV}}(\mathcal T_a(\widehat{\mathcal D}_\rho),\mathcal T_a(\mathcal D_\rho))
\le 4\,d_{\mathrm{TV}}(\widehat{\mathcal D}_\rho,\mathcal D_\rho).
$$
再用三角不等式分解
$$
d_{\mathrm{TV}}(\widehat{\mathcal D}_\rho,\mathcal D_\rho)
\le d_{\mathrm{TV}}(\widehat{\mathcal D}_\rho,\tilde{\mathcal D}_\rho)
+ d_{\mathrm{TV}}(\tilde{\mathcal D}_\rho,\mathcal D_\rho).
$$
其中第一项由引理 5.3 上界为 $\varepsilon_{\mathrm{mix}}+\varepsilon_{\mathrm{geom}}$，第二项由核展开误差上界为 $\delta/(1-\delta)$。合并即得。$\square$

---

### 5.3.5 多样本联合分布的 TV 上界（可选）

当输出 $t$ 个样本为 i.i.d. 时，联合分布为 $\tilde{\mathcal D}^{\otimes t}$，理想为 $\mathcal D_{\mathrm{IoU}}^{\otimes t}$，并有线性放大界
$$
d_{\mathrm{TV}}(\tilde{\mathcal D}^{\otimes t},\mathcal D_{\mathrm{IoU}}^{\otimes t})
\le t\cdot d_{\mathrm{TV}}(\tilde{\mathcal D},\mathcal D_{\mathrm{IoU}}).
$$
因此定理 5.5 的单次边缘界可直接转化为联合分布界。

---

### 5.3.6 参数分配示例（保持端到端对数结构）

给定目标误差 $\varepsilon\in(0,1)$，一种直接的预算分配为
$$
2^{-M}\le \frac{\varepsilon}{4},\qquad
\frac{\delta}{1-\delta}\le \frac{\varepsilon}{16},\qquad
\varepsilon_{\mathrm{mix}}\le \frac{\varepsilon}{16},\qquad
\varepsilon_{\mathrm{geom}}\le \frac{\varepsilon}{16}.
$$
例如可取
$$
M=\left\lceil \log_2\frac{4}{\varepsilon}\right\rceil,\qquad
\delta=\frac{\varepsilon}{32},
$$
并令
$$
\varepsilon_{\mathrm{clip}}=\frac{\varepsilon}{64},\qquad
M_{\mathrm{bucket}}=\left\lceil \log_2\frac{64}{\varepsilon}\right\rceil,
$$
则
$$
\varepsilon_{\mathrm{geom}}=2\left(\frac{\varepsilon}{64}+2^{-M_{\mathrm{bucket}}}\right)\le \frac{\varepsilon}{16}.
$$
该选择仅引入 $M=\Theta(\log(1/\varepsilon))$ 与 $M_{\mathrm{bucket}}=\Theta(\log(1/\varepsilon))$，并以加法形式进入桶内扫描窗口 $\log(N/\varepsilon_{\mathrm{clip}})=\log N+\Theta(\log(1/\varepsilon))$，不与扫描过程中的 $\log N$ 发生新的乘法叠加。

---

## 5.4 端到端最坏复杂度主定理与分解账本

本节给出端到端的最坏时间与空间界，并把“共享结构”和“按 $\ell$ 重复的结构”严格区分。

### 5.4.1 指数和项数与提案数

指数和项数记为 $L=L(\delta)$，满足
$$
L(\delta)=O\big((\log\kappa+\log(1/\delta))\log(1/\delta)\big),
\qquad
\kappa=\frac{S_{\max}}{S_{\min}}.
$$
外层有界拒绝次数为 $M$，因此提案数固定为 $t_{\mathrm{prop}}=tM$。

令
$$
U=\{\ell\in[L]: t_\ell>0\},
\qquad
L'=|U|\le \min\{L,tM\}
$$
为被实际使用到的分量集合及其大小。

---

### 5.4.2 共享几何骨架与按 $\ell$ 的可复用工作区

端到端实现中，有两类几何对象：

**(A) 共享几何骨架（一次性构建）**
- 排序得到 $x$ 端点与 slab 列表、事件表；
- 排序得到 $y$ 端点与基本段；
- 构建乘积线段树的拓扑骨架；
- 构建 stabbing 的段树骨架与每个矩形 $y$ 区间的 canonical cover 列表（这部分与权重无关）。

这些对象只依赖输入几何与端点压缩，因此可在端到端只构建一次。

**(B) 按 $\ell$ 的工作区（可复用）**
- 乘积线段树的节点值（$\mathrm{SF},\mathrm{SG},\mathrm{SP}$ 与懒标记）依赖当前 $\ell$ 的权重，按 $\ell$ 迭代时可在同一内存上清零复用；
- stabbing 的桶内数组段与 record 的 level 分组依赖当前 $\ell$（因为 dyadic level 由 $w_\ell,v_\ell$ 决定），按 $\ell$ 串行执行时同样可复用同一块工作内存。

本章默认按 $\ell$ 串行执行（不并行维护多个 $\ell$ 的运行态），因此空间界不随 $L'$ 乘法增长。

---

### 5.4.3 固定 $\ell$ 子程序代价（口径化定义）

为了端到端账本清晰，将固定 $\ell$ 的工作拆成两类子程序：

**(1) 质量统计（mass-only）子程序**  
仅维护乘积线段树，扫描一遍 $x$ 事件以得到
$$
Z_\ell=\sum_{i}M_i^{(\ell)},
\qquad
M_i^{(\ell)}=\Delta x_i\int F_i^{(\ell)}(y)G_i^{(\ell)}(y)\,dy.
$$
该过程不维护 stabbing。其最坏时间为
$$
T_{\ell,\mathrm{mass}}(N)=O(N\log N).
$$

**(2) 条件采样（cond）子程序**  
需要生成 $t_\ell$ 个来自分量分布（或其 TV-近似实现）的样本。典型组织为：
- （可选）一遍扫描得到 slab 质量并完成 slab 指派；
- 构建当前 $\ell$ 的两套 stabbing 运行态（包括 dyadic level 计算、按 $(\text{bucket},\text{level})$ 分组、数组段分配与初始化）；
- 第二遍扫描维护乘积线段树与两套 stabbing，并在每个被指派的 slab 内重复进行 $t_\ell$ 次 slab 条件采样。

其最坏时间写作
$$
T_{\ell,\mathrm{cond}}(N,t_\ell)=O\Big(
N\log N
+
t_\ell\cdot(\log N+\log(N/\varepsilon_{\mathrm{clip}})+M_{\mathrm{bucket}})
\Big).
$$
其中 $N\log N$ 项包含：当前 $\ell$ 的 stabbing 运行态构建（在本章的计算模型下可用确定的分组口径做到 $O(N\log N)$）以及 Pass2 扫描更新；样本项包含一次 `SampleY()` 与两次 stabbing 查询的确定最坏代价。

---

### 5.4.4 端到端时间账本（逐项加和）

端到端时间分解为：

1) **共享几何骨架构建**
$$
T_{\mathrm{skel}}=O(N\log N).
$$

2) **对所有 $\ell$ 的质量统计（用于得到 $\{Z_\ell\}$ 与混合系数）**
$$
T_{\mathrm{all-mass}}=\sum_{\ell=1}^L T_{\ell,\mathrm{mass}}(N)=O(L\cdot N\log N).
$$

3) **分量指派（生成 $\{L_j\}_{j=1}^{tM}$ 并分组）**  
构建 $\ell$ 上离散采样器并抽样 $tM$ 次：
$$
T_{\ell\text{-assign}}=O(L+tM).
$$

4) **仅对被使用到的 $\ell$ 执行条件采样并回填提案**  
对每个 $\ell\in U$，执行 $T_{\ell,\mathrm{cond}}(N,t_\ell)$：
$$
T_{\mathrm{cond}}
=\sum_{\ell\in U}T_{\ell,\mathrm{cond}}(N,t_\ell)
=O\Big(
L'\cdot N\log N
+
tM(\log N+\log(N/\varepsilon_{\mathrm{clip}})+M_{\mathrm{bucket}})
\Big),
$$
其中用到 $\sum_{\ell\in U}t_\ell=tM$。

5) **外层有界拒绝校正**  
对每个输出位置最多尝试 $M$ 次，每次做常数次算术与一次比较，因此
$$
T_{\mathrm{accept}}=O(tM).
$$

合并得：

**定理 5.6（端到端最坏时间复杂度）**  
端到端算法的最坏时间满足
$$
\begin{aligned}
T(N,t,\varepsilon)
&=
T_{\mathrm{skel}}
+T_{\mathrm{all-mass}}
+T_{\ell\text{-assign}}
+T_{\mathrm{cond}}
+T_{\mathrm{accept}}\\
&=
O\Big(
N\log N
+L\cdot N\log N
+(L+tM)
+L'\cdot N\log N\\
&\hspace{3.5em}
+tM(\log N+\log(N/\varepsilon_{\mathrm{clip}})+M_{\mathrm{bucket}})
+tM
\Big).
\end{aligned}
$$
在最坏情形 $L'=L$ 下，上界化简为
$$
T(N,t,\varepsilon)
=
O\Big(
L(\delta)\cdot N\log N
+
t\log(1/\varepsilon)\cdot(\log N+\log(N/\varepsilon_{\mathrm{clip}})+M_{\mathrm{bucket}})
\Big).
$$
其中扫描更新部分不出现 $L'\cdot N\log^2 N$ 项；原因是 stabbing 的 membership 更新为确定 $O(1)$，每次插删触达的 canonical bucket 数为 $O(\log N)$，从而 Pass2 的更新扫描为 $O(N\log N)$（逐 $\ell$）。

---

### 5.4.5 空间复杂度与两种内存组织

端到端实现需要存储以下对象：

- 共享几何骨架：端点数组、事件表、slab 列表、乘积线段树拓扑、stabbing 段树拓扑与 canonical cover 列表，整体 $O(N\log N)$；
- 混合相关数组：$\{(\omega_\ell,\alpha_\ell)\}$、$\{Z_\ell\}$、$\{\hat\pi_\ell\}$，整体 $O(L)$；
- 提案缓冲：长度 $tM$ 的提案数组，整体 $O(tM)$；
- 按 $\ell$ 复用的工作区：乘积线段树节点值与 stabbing 的桶内数组段/record，整体 $O(N\log N)$。

因此在“按 $\ell$ 串行执行并复用工作区”的默认口径下，空间界为
$$
S(N,t,\varepsilon)=O(N\log N + L + tM).
$$

此外，slab 质量数组的保存方式有两种等价选项：

**选项 A（低空间）**  
只保存 $\{Z_\ell\}$，不保存所有 $M_i^{(\ell)}$。当某个 $\ell$ 被使用到时，在执行 $T_{\ell,\mathrm{cond}}$ 前为该 $\ell$ 再计算一次 slab 质量用于 slab 指派。该选项维持上述空间界。

**选项 B（低重复扫描）**  
在执行 $T_{\mathrm{all-mass}}$ 时同时保存每个 $\ell$ 的 slab 质量数组 $\{M_i^{(\ell)}\}_{i=0}^{m-1}$，从而对 $\ell\in U$ 的条件采样阶段可直接进行 slab 指派而无需再次计算 slab 质量。该选项额外占用
$$
O(L\cdot m)=O(LN)
$$
空间，用于换取更少的重复扫描。

---

## 5.5 计算模型与结构性取舍（可审查口径）

本节给出两套自洽的计算口径，用于支撑以下操作的“确定最坏时间”叙述：
- dyadic level 的计算与比较；
- 桶内阈值 $\tau=\varepsilon_{\mathrm{clip}}T/m$ 的计算；
- Bernoulli$(w/2^k)$ 与外层 Bernoulli$(a(P))$ 的实现；
- 预处理中的按键分组（radix）与数组段分配。

---

### 5.5.1 模型 A：word-RAM + 对数域权重（工程友好）

- 字长 $W=\Theta(\log N+\log(1/\varepsilon))$；
- 权重用对数域表示（例如存 $\log w_\ell(r)=-\alpha_\ell A(r)$ 与 $\log v_\ell(s)=\log\omega_\ell-\alpha_\ell A(s)$）；
- 支持 $O(1)$ 的 `msb/ceil_log2`、整数比较、以及 `RandInt(m)`；
- Bernoulli$(p)$ 通过 `RandInt(2^b)` 与阈值比较实现（$b=O(W)$）。

在该口径下，桶内扫描窗口长度为 $O(\log(N/\varepsilon_{\mathrm{clip}}))=O(W)$，因此每次桶内构建 alias 的整数规模与随机位消耗均为确定上界。

---

### 5.5.2 模型 B：纯整数/定点位宽（更保守、可审查）

- 坐标为 $b$ 位整数；面积与交面积为整数；
- 将所有权重与概率用 $b_w=O(W)$ 位定点（或统一分母的有理数）表示；
- dyadic level 通过最高位位置 `msb` 得到；
- Bernoulli$(w/2^k)$ 用“生成 $k$ 位随机整数并比较”实现，且 $k$ 被窗口上界控制在 $O(W)$；
- 外层接受概率 $a(P)=1/[2(1-\rho(P))]$ 可用定点形式计算并以同样的比较方式完成一次 Bernoulli。

该口径把位宽相关的代价吸收进常数次 word 操作，不引入额外的 $\log\log$ 因子，也不改变第 5.4 节的最坏复杂度账本。

---

### 5.5.3 结构性取舍：为何“单桶扫描 + 显式 TV 误差”是自然的代价

动态按权重抽样若同时追求三点：
1) **精确分布**（零误差）；
2) **更新确定 $O(1)$**（对任意插删序列）；
3) **查询确定 $O(1)$**（对任意时刻与任意查询点）；  
通常会迫使数据结构在查询端隐式完成“动态前缀质量定位 / 动态选择”一类操作。对一般权重与一般更新序列，这类操作往往要求非平凡的全局信息维护。

本书采用的路线将困难显式移位为两类可控代价：

- **更新端保持确定性**：每条 canonical membership 的桶内插删用预分配数组段与 swap-delete 句柄完成，避免动态映射与均摊扩容；
- **查询端仅做局部扫描**：在整条 Path 上只选中一个桶进行窗口扫描与小 alias 构建；
- **误差端用 TV 量化**：用平均阈值整层裁剪产生 $\varepsilon_{\mathrm{clip}}$，用有界拒绝产生 $2^{-M_{\mathrm{bucket}}}$，并通过混合凸性与外层 Lipschitz 常数把误差闭环为可加账本。

因此，“单桶扫描 + 显式 TV 误差”并不是额外的技巧负担，而是让更新与查询同时保持确定上界时，一种结构上自然且可审查的代价交换方式。
