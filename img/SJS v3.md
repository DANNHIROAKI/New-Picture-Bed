# 第 1 章 目标分布与可采样分解

本章只做三件事：

1. 精确定义要在什么集合上采样、采什么分布；
2. 给出一个与半开语义严格一致的扫描事件序；
3. 证明一个**无重分解**：把全体相交对 $J$ 分解成若干互不重叠的“事件块” $J_e$，从而把“在 $J$ 上均匀采样”化为“按块权重采块 + 块内均匀”。

------

## 1.1 输入模型与半开语义下的相交

给定两类 $d$ 维轴对齐半开盒集合：
$$
R_c=\{r_{c1},\dots,r_{cn}\},\qquad 
R_{\bar c}=\{r_{\bar c1},\dots,r_{\bar cm}\}.
$$
每个盒
$$
r=\prod_{i=1}^{d}[L_i(r),R_i(r)),\qquad L_i(r)<R_i(r).
$$
半开语义下，两盒 $r,s$ 的交集非空当且仅当对所有维度 $i$：
$$
r\cap s\neq\varnothing
\iff
\big(L_i(r)<R_i(s)\big)\ \wedge\ \big(L_i(s)<R_i(r)\big).
$$
该判定使用严格不等号，保证“贴边”（例如 $R_i(r)=L_i(s)$）不计为相交。

------

## 1.2 目标集合与采样目标分布

仅考虑跨集合、有序方向固定的相交对集合：
$$
J=\{(r_c,r_{\bar c})\mid r_c\in R_c,\ r_{\bar c}\in R_{\bar c},\ r_c\cap r_{\bar c}\neq\varnothing\}.
$$
给定 $t\ge 1$，输出序列
$$
Z_1,\dots,Z_t\in J
$$
满足 **i.i.d. 均匀（有放回）**：
$$
\Pr\{Z_j=P\}=\frac{1}{|J|}\quad(\forall P\in J),\qquad Z_1,\dots,Z_t\ \text{相互独立}.
$$
若 $|J|=0$，输出空序列。

------

## 1.3 扫描轴事件序与相交对的无重块分解

本节在第一维上做一次 plane sweep，并将每个相交对唯一归属到某个“START 事件块”。

### 1.3.1 扫描维度、事件与全序

固定扫描维度为第 1 维。对每个盒 $r\in R_c\cup R_{\bar c}$ 定义两类事件：

- $\mathrm{START}(r)$：坐标 $x=L_1(r)$；
- $\mathrm{END}(r)$：坐标 $x=R_1(r)$。

令 $E$ 为所有事件的集合，令 $E^+\subset E$ 为所有 START 事件集合。

为了与半开语义严格一致，并消除同坐标歧义，对 $E$ 定义一个总序 $\prec$：

1. 若 $x(e)<x(e')$，则 $e\prec e'$；
2. 若 $x(e)=x(e')$，则所有 END 事件严格先于所有 START 事件：
    $\mathrm{END}(\cdot)\prec \mathrm{START}(\cdot)$；
3. 若同为 START 且坐标相同，则按盒子的一个固定全序（例如唯一 id 的升序）打破平局；
    若同为 END 且坐标相同，也用固定规则打破平局（仅为确定性，不影响半开语义）。

> 这一序使得当 $R_1(r)=L_1(s)$ 时，$\mathrm{END}(r)\prec \mathrm{START}(s)$，从而在处理 $\mathrm{START}(s)$ 之前 $r$ 已被移出活跃集合，符合“贴边不相交”。

### 1.3.2 活跃集合

在按 $\prec$ 顺序处理事件的过程中，维护两类活跃集合：
$$
A_c,\ A_{\bar c}\subseteq R_c,\ R_{\bar c}.
$$
处理规则为：

- 遇到 $\mathrm{END}(r)$：将 $r$ 从其所属类的活跃集合中删除；
- 遇到 $\mathrm{START}(r)$：先定义该事件的“事件块”（见下节），然后将 $r$ 插入其所属类的活跃集合。

因此，在处理某个 START 事件 $e=\mathrm{START}(q)$ 的时刻，活跃集合精确包含了所有满足
$$
L_1(s)<L_1(q)\quad\text{或}\quad \big(L_1(s)=L_1(q)\ \wedge\ \mathrm{START}(s)\prec \mathrm{START}(q)\big),
$$
且同时
$$
L_1(q)<R_1(s)
$$
的盒 $s$（并按类别分在 $A_c$ 或 $A_{\bar c}$ 中）。后一条件来自“尚未 END”，与半开语义严格一致。

### 1.3.3 剩余维度投影与事件块定义

记剩余维度数 $m=d-1$。对任意盒 $r$ 定义其在维度 $2..d$ 上的投影：
$$
r^\downarrow=\prod_{i=2}^{d}[L_i(r),R_i(r)).
$$
并用同样的严格判定定义 $r^\downarrow\cap s^\downarrow\neq\varnothing$。

对任意 START 事件 $e=\mathrm{START}(q)\in E^+$，定义其对应的 partner 集合 $K_e$ 与事件块 $J_e$：

- 若 $q\in R_c$（事件盒来自左类）：
  $$
  K_e=\{s\in A_{\bar c}\mid q^\downarrow\cap s^\downarrow\neq\varnothing\},
  \qquad
  J_e=\{(q,s)\mid s\in K_e\}.
  $$

- 若 $q\in R_{\bar c}$（事件盒来自右类）：
  $$
  K_e=\{s\in A_{c}\mid s^\downarrow\cap q^\downarrow\neq\varnothing\},
  \qquad
  J_e=\{(s,q)\mid s\in K_e\}.
  $$

注意：$J$ 的方向固定为 $(r_c,r_{\bar c})$，因此当事件盒 $q\in R_{\bar c}$ 时，块内元素写作 $(s,q)$ 以保持方向一致。

定义块权重（块大小）：
$$
w_e:=|J_e|=|K_e|.
$$

### 1.3.4 核心定理：全体相交对的无重分解

**定理 1.1（按“较晚 START”归属的块分解）**
 在上述事件序 $\prec$ 与事件块定义下，有
$$
J=\biguplus_{e\in E^+} J_e,
\qquad
|J|=\sum_{e\in E^+} w_e.
$$
其中 $\biguplus$ 表示不交并（互不重叠的并）。

**证明** 分两步：覆盖性与互斥性。

**(1) 覆盖性：任意 $(r_c,r_{\bar c})\in J$ 必落入某个 $J_e$。**
 令
$$
e_c=\mathrm{START}(r_c),\qquad e_{\bar c}=\mathrm{START}(r_{\bar c}),
$$
并令
$$
e^\star=\max_{\prec}\{e_c,e_{\bar c}\}
$$
为二者中在总序 $\prec$ 下更晚的 START 事件。设 $e^\star=\mathrm{START}(q)$，则 $q$ 是 $\{r_c,r_{\bar c}\}$ 中“较晚开始”的那个（若同坐标则由 tie-break 决定）。

不妨设 $q=r_c$（另一种情形完全对称）。因为 $(r_c,r_{\bar c})\in J$，在第 1 维必有
$$
L_1(r_c)<R_1(r_{\bar c}).
$$
而 $e^\star$ 发生在 $x=L_1(r_c)$，且所有同坐标 END 事件已先处理，因此 $r_{\bar c}$ 在该时刻尚未被删除；又因为 $e_{\bar c}\prec e^\star$，$r_{\bar c}$ 的 START 已先处理，因此 $r_{\bar c}\in A_{\bar c}$ 在 $e^\star$ 时刻成立。再结合剩余维度相交 $r_c^\downarrow\cap r_{\bar c}^\downarrow\neq\varnothing$，得到 $r_{\bar c}\in K_{e^\star}$，从而
$$
(r_c,r_{\bar c})\in J_{e^\star}.
$$
因此每个 $(r_c,r_{\bar c})\in J$ 都属于某个事件块。

**(2) 互斥性：任意 $(r_c,r_{\bar c})\in J$ 至多属于一个 $J_e$。**
 由构造，若 $(r_c,r_{\bar c})\in J_e$，则 $e$ 必为 $e_c,e_{\bar c}$ 中较晚的 START：因为在处理较早 START 时，较晚 START 的盒尚未进入活跃集合，不可能成为 partner；而在处理较晚 START 时，较早者若与之相交则必仍活跃并被纳入 partner 集合。由于 $\prec$ 是总序，较晚 START 唯一，因此归属的事件块唯一。

由 (1)(2) 得 $J=\biguplus_{e\in E^+}J_e$。权重求和式 $|J|=\sum_e |J_e|=\sum_e w_e$ 由不交并立即得到。证毕。

------

## 1.4 从全局均匀采样到“按块加权 + 块内均匀”的统一原理

有了定理 1.1，均匀采样 $J$ 可以完全等价为对事件块的两阶段抽样。

令
$$
W:=|J|=\sum_{e\in E^+} w_e.
$$
若 $W=0$，则 $J=\varnothing$ 且应输出空序列。

否则定义如下抽样过程（单次输出）：

1. 以概率
   $$
   \Pr\{E=e\}=\frac{w_e}{W}
   $$
   从 $E^+$ 中抽取一个 START 事件 $e$；

2. 在该事件块 $J_e$ 内做一次均匀抽样，得到 $Z\in J_e$。

**命题 1.2（单次输出的均匀性）**
 上述过程输出的 $Z$ 在 $J$ 上均匀：
$$
\Pr\{Z=P\}=\frac{1}{|J|},\quad\forall P\in J.
$$
**证明** 任取 $P\in J$。由定理 1.1，存在唯一事件 $e(P)$ 使得 $P\in J_{e(P)}$。因此
$$
\Pr\{Z=P\}
=\Pr\{E=e(P)\}\cdot \Pr\{Z=P\mid E=e(P)\}
=\frac{w_{e(P)}}{W}\cdot \frac{1}{w_{e(P)}}
=\frac{1}{W}
=\frac{1}{|J|}.
$$
证毕。

**命题 1.3（多次输出的 i.i.d.）**
 若对每次输出都使用相互独立的随机性重复执行该两阶段过程，并且块内抽样为有放回，则 $Z_1,\dots,Z_t$ 相互独立且同分布为 $J$ 上的均匀分布。

------

### 1.4.1 事件块原语

定理 1.1 与命题 1.2–1.3 将整个问题规约为：对每个 START 事件 $e$，能够有效获得/生成事件块 $J_e$ 的信息。等价地，只需对 partner 集 $K_e$ 支持三类原语：

- **计数**：返回 $w_e=|K_e|$；
- **枚举**：输出 $K_e$（从而得到 $J_e$）；
- **均匀采样**：从 $K_e$ 上 i.i.d. 均匀（有放回）抽取 $k$ 个 partner。

由于活跃集合已经确保第 1 维与事件位置严格重叠，$K_e$ 的判定完全发生在投影空间 $2..d$ 上：
$$
q^\downarrow\cap s^\downarrow\neq\varnothing.
$$
因此后续章节的核心任务，就是在动态活跃集上实现对投影盒的 `COUNT/REPORT/SAMPLE`，并在第 3 章把这些原语组合成对 $J$ 的 i.i.d. 均匀采样。

# 第 2 章 事件块查询的数据结构：模式化分解与段树递归

本章在扫描线过程中维护“对侧类”的动态活跃集，并在每个 $\mathrm{START}(q)$ 事件处实现事件块原语：

- `COUNT(e)`：返回事件块 $e$ 的 partner 集大小；
- `REPORT(e)`：枚举 partner 集；
- `SAMPLE(e,k)`：在 partner 集上返回 $k$ 个**独立同分布**（i.i.d.）、**均匀**、**有放回**样本。

扫描线沿第 $1$ 维推进。对任意盒 $r\subset\mathbb{R}^d$，写成半开轴对齐盒：
$$
r=\prod_{i=1}^{d}[L_i(r),R_i(r)),\qquad L_i(r)<R_i(r).
$$
当处理某个 $\mathrm{START}(q)$ 时，只需检查其余 $m=d-1$ 维的投影是否相交。记投影盒
$$
r^\downarrow=\prod_{i=2}^{d}[L_i(r),R_i(r)).
$$
在本章中，所有对象均带唯一标识 $\mathrm{id}(r)$；即使几何参数完全相同，也视为不同元素。

本章核心路线：

1. 将“每一维相交”拆成两类互斥形态（stabbing / range），从而把 $m$ 维相交问题划分为 $2^m$ 个互不重叠的“模式”子问题；
2. 为每个固定模式 $g\in\{A,B\}^m$ 构造递归段树结构 $\mathcal{D}_{m,g}$，支持动态更新与 `COUNT/REPORT/SAMPLE`；
3. 在模式间按权重混合，得到完整的事件块原语。

---

## 2.1 一维相交的二分结构与 $m$ 维模式分解

### 2.1.1 一维二分：stabbing 形态与 range 形态

考虑一维半开区间
$$
I=[L(I),R(I)),\qquad J=[L(J),R(J)),\qquad L(\cdot)<R(\cdot).
$$
半开语义下，相交当且仅当
$$
I\cap J\neq\varnothing
\iff
L(I)<R(J)\ \wedge\ L(J)<R(I).
$$
令查询点为 $x:=L(J)$（即查询区间的左端点），定义两类条件：

- **A（stabbing 条件）**
  $$
  \mathbf{A}(I\mid J):\quad L(I)\le x < R(I).
  $$
  直观上表示 $I$ “刺中”查询点 $x$。

- **B（range 条件）**
  $$
  \mathbf{B}(I\mid J):\quad L(J)<L(I)<R(J).
  $$
  直观上表示 $I$ 的左端点严格落在 $J$ 的内部。

**引理 2.1（一维相交的互斥完备划分）**  
对任意区间 $I,J$，有
$$
I\cap J\neq\varnothing
\iff
\mathbf{A}(I\mid J)\ \vee\ \mathbf{B}(I\mid J),
$$
并且 $\mathbf{A}$ 与 $\mathbf{B}$ 互斥。

**证明**  

- （充分性）  
  若 $\mathbf{A}(I\mid J)$ 成立，则 $x\in I$。又 $x=L(J)\in J$（因为 $L(J)<R(J)$），故 $x\in I\cap J$，从而相交。  
  若 $\mathbf{B}(I\mid J)$ 成立，则 $L(I)\in (L(J),R(J))\subseteq J$，且 $L(I)\in I$（半开区间包含左端点），故 $L(I)\in I\cap J$，相交成立。

- （必要性）  
  若 $I\cap J\neq\varnothing$，则 $L(I)<R(J)$ 且 $L(J)<R(I)$。分两种情况：  
  1) 若 $L(I)\le L(J)=x$，结合 $x<L(J)<R(I)$ 得到 $L(I)\le x<R(I)$，即 $\mathbf{A}(I\mid J)$；  
  2) 若 $L(I)>L(J)$，由 $L(I)<R(J)$ 得 $L(J)<L(I)<R(J)$，即 $\mathbf{B}(I\mid J)$。

- （互斥性）  
  $\mathbf{A}$ 要求 $L(I)\le L(J)$，而 $\mathbf{B}$ 要求 $L(I)>L(J)$，不可同时成立。证毕。

---

### 2.1.2 $m$ 维投影盒的模式向量与无重分解

对扫描维以外的 $m=d-1$ 个维度，将其重编号为 $j=1,\dots,m$ 对应原维 $i=j+1$，定义
$$
L'_j(r):=L_{j+1}(r),\qquad R'_j(r):=R_{j+1}(r),
$$
从而
$$
r^\downarrow=\prod_{j=1}^{m}[L'_j(r),R'_j(r)).
$$

令对侧活跃集合（投影对象集合）为 $S$。给定查询盒 $q^\downarrow$，定义其余维相交集合：
$$
S_\cap(q):=\{\,s\in S\mid s^\downarrow\cap q^\downarrow\neq\varnothing\,\}.
$$

对每一维 $j$，应用引理 2.1 得到两种互斥形态。对一个模式向量
$$
g=(g_1,\dots,g_m)\in\{A,B\}^m,
$$
定义模式集合
$$
S_g(q)
=
\big\{s\in S\ \big|\ \forall j\in[m],\
\begin{cases}
L'_j(s)\le L'_j(q)<R'_j(s), & g_j=A,\\
L'_j(q)<L'_j(s)<R'_j(q), & g_j=B
\end{cases}
\big\}.
$$

**定理 2.2（模式分解无重且完备）**  
对任意查询盒 $q^\downarrow$，
$$
S_\cap(q)=\biguplus_{g\in\{A,B\}^m} S_g(q),
$$
即所有 $S_g(q)$ 两两不交，且并为全部相交对象。

**证明**  

- 覆盖性：若 $s\in S_\cap(q)$，则每一维区间 $[L'_j(s),R'_j(s))$ 与 $[L'_j(q),R'_j(q))$ 相交。由引理 2.1，该维必满足 A 或 B 且二者互斥。将每维的选择组成 $g$，则 $s\in S_g(q)$。  
- 不交性：若 $g\neq g'$，存在某维 $j$ 使 $g_j\neq g'_j$。该维 A/B 条件互斥，故不存在 $s$ 同时属于两个模式集合。证毕。

---

## 2.2 一维原语 I：stabbing 型段树结构

本节处理 A 形态：给定点 $x$，维护所有包含 $x$ 的活跃区间，并支持 `COUNT/REPORT/SAMPLE`。

### 2.2.1 坐标压缩与段树

固定一维。设该维所有区间端点来自一个静态全集 $\mathcal{U}$（例如某一类盒的所有投影区间）。取端点集合
$$
X=\{L(I),R(I)\mid I\in\mathcal{U}\},
$$
去重排序得
$$
x_0<x_1<\cdots<x_{M}.
$$
建立段树 $T$，每个结点 $v$ 对应半开段
$$
\mathrm{seg}(v)=[x_\ell,x_r),
$$
叶子对应基本段 $[x_i,x_{i+1})$。段树高度为 $O(\log M)$。

对任意区间 $I=[a,b)$，其**规范覆盖**（canonical cover）$\mathcal{C}(I)$ 是一组结点集合，满足：

1. $\bigcup_{v\in\mathcal{C}(I)}\mathrm{seg}(v)=[a,b)$；
2. 对所有 $v\in\mathcal{C}(I)$，$\mathrm{seg}(v)\subseteq [a,b)$；
3. 结点段两两不交；
4. 覆盖是“最粗”的：对任意 $v\in\mathcal{C}(I)$，其父结点段不完全包含于 $[a,b)$。

标准段树性质给出 $|\mathcal{C}(I)|=O(\log M)$，且可在 $O(\log M)$ 时间求得。

### 2.2.2 桶、句柄与动态更新

对每个结点 $v$ 维护一个桶 $B(v)$，用于存放当前活跃区间中满足 $v\in\mathcal{C}(I)$ 的那些区间 $I$：
$$
B(v):=\{\,I\ \text{活跃}\mid v\in\mathcal{C}(I)\,\}.
$$

桶实现采用“致密数组 + 位置句柄”的结构：

- 桶内部用数组存储元素；
- 对每个被插入的元素，在每个桶里保存其数组下标句柄；
- 删除时用 swap-with-last 维持致密性，从而每次删除 $O(1)$。

于是：

- `INSERT(I)`：对所有 $v\in\mathcal{C}(I)$ 执行 `bucket_insert(B(v),I)`；
- `DELETE(I)`：用句柄从每个相关桶中 $O(1)$ 删除。

单次插入/删除触及 $O(\log M)$ 个桶。

### 2.2.3 路径分解与无重性

给定查询点 $x$，令 $\mathsf{Path}(x)$ 为从根到包含 $x$ 的叶子的路径结点集合。由于结点段为半开段，$x$ 属于唯一叶子段，从而路径唯一。

**引理 2.3（路径唯一归属）**  
若区间 $I$ 满足 $x\in I$，则存在唯一结点 $v^\star\in\mathcal{C}(I)$ 使得 $v^\star\in\mathsf{Path}(x)$。

**证明**  
$\mathcal{C}(I)$ 的结点段两两不交且并为 $I$。点 $x\in I$ 必落入其中唯一一段 $\mathrm{seg}(v^\star)$，故 $v^\star$ 唯一。又 $\mathrm{seg}(v^\star)$ 包含 $x$，因此 $v^\star$ 位于 $\mathsf{Path}(x)$ 上。证毕。

由此得到答案集合的无重分解：

**推论 2.4（stabbing 答案的无重并）**  
令
$$
\mathcal{V}(x):=\mathsf{Path}(x),
$$
则
$$
\{\,I\ \text{活跃}\mid x\in I\,\}
=
\biguplus_{v\in\mathcal{V}(x)} B(v).
$$

### 2.2.4 三原语与均匀采样

记
$$
C(x):=\sum_{v\in\mathcal{V}(x)} |B(v)|.
$$

- `COUNT(x)`：返回 $C(x)$；
- `REPORT(x)`：依次枚举所有 $v\in\mathcal{V}(x)$ 的桶 $B(v)$ 并输出（由推论 2.4 不重）；
- `SAMPLE(x,k)`：若 $C(x)=0$ 返回空序列；否则对每个样本位置 $t=1,\dots,k$ 独立执行：
  1) 以权重 $|B(v)|$ 在 $v\in\mathcal{V}(x)$ 上抽取结点 $V_t$；  
  2) 在桶 $B(V_t)$ 内均匀抽取一个区间作为第 $t$ 个样本。

为了将“按权重抽结点”实现为 $O(1)$ 每次抽样，引入子程序：

**WeightedChoice 子程序（alias 版本）**  
给定正权重向量 $(w_1,\dots,w_t)$，可在 $O(t)$ 预处理后，用 alias 方法在 $O(1)$ 时间抽取一个索引 $i$ 满足 $\Pr[i]=w_i/\sum_j w_j$。  
当 $t=O(\log M)$ 时，预处理代价为 $O(\log M)$。

**定理 2.5（均匀性与独立性）**  
`SAMPLE(x,k)` 返回的 $k$ 个样本在集合 $\{I\mid x\in I\}$ 上均匀、独立同分布（有放回）。

**证明**  
任取答案区间 $I$。由引理 2.3，$I$ 唯一归属到某个路径结点 $v^\star\in\mathcal{V}(x)$，且 $I\in B(v^\star)$。一次采样命中 $I$ 的概率为
$$
\Pr[I]
=\Pr[V=v^\star]\cdot \Pr[I\mid V=v^\star]
=\frac{|B(v^\star)|}{C(x)}\cdot\frac{1}{|B(v^\star)|}
=\frac{1}{C(x)},
$$
与 $I$ 无关，故均匀。每个样本位置使用独立随机性且有放回，故样本 i.i.d.。证毕。

### 2.2.5 复杂度

令段树高度为 $O(\log M)$。

- `INSERT/DELETE`：$O(\log M)$ 个桶更新；
- `COUNT`：$O(\log M)$；
- `REPORT`：$O(\log M + K)$，其中 $K$ 为输出规模；
- `SAMPLE(x,k)`：对路径桶构建一次 alias 表 $O(\log M)$，随后每次抽样 $O(1)$，总计 $O(\log M + k)$。

（若不使用 alias，而以前缀和 + 二分抽桶，则总复杂度为 $O(\log M + k\log\log M)$ 或 $O((k+1)\log M)$，取决于具体实现。）

---

## 2.3 一维原语 II：秩化的 range 型段树结构

本节处理 B 形态：给定开区间 $(\ell,r)$，返回所有活跃对象中满足 $\ell < L(\cdot) < r$ 的元素（这里 $L(\cdot)$ 表示该维的左端点），并支持均匀采样。

### 2.3.1 全局秩域与开区间映射

固定一个静态全集 $\mathcal{U}$，每个对象 $u\in\mathcal{U}$ 带键值 $a(u)\in\mathbb{R}$（本章中取 $a(u)=L(u)$）。为处理重复值与严格不等号，引入唯一键
$$
\alpha(u):=(a(u),\mathrm{id}(u)),
$$
按字典序排序得到全局秩
$$
\mathrm{rank}(u)\in\{1,2,\dots,|\mathcal{U}|\}.
$$
记 $N:=|\mathcal{U}|$。

对查询开区间 $(\ell,r)$，定义两个边界秩：

- $L:=\mathrm{upper\_bound}(\ell)$：满足 $a(\text{rank}=L)>\ell$ 的最小秩；若不存在，则 $L:=N+1$；
- $R:=\mathrm{lower\_bound}(r)$：满足 $a(\text{rank}=R)\ge r$ 的最小秩；若不存在，则 $R:=N+1$。

于是
$$
\{\,u\in\mathcal{U}\mid \ell<a(u)<r\,\}
=
\{\,u\in\mathcal{U}\mid \mathrm{rank}(u)\in [L,R)\,\}.
$$
当 $L\ge R$ 时答案为空。

### 2.3.2 段树、路径累积桶与规范分解

在秩域 $\{1,\dots,N\}$ 上建立段树。每个结点 $v$ 对应秩段 $\mathrm{seg}(v)=[p_\ell,p_r)$，维护桶
$$
B(v):=\{\,u\ \text{活跃}\mid \mathrm{rank}(u)\in \mathrm{seg}(v)\ \text{且}\ v\ \text{在}\ \mathrm{rank}(u)\ \text{到根的路径上}\,\}.
$$

动态维护采用“叶到根路径累积”：

- `INSERT(u)`：从叶结点 $\mathrm{rank}(u)$ 沿父指针到根，对路径上所有结点 $v$ 将 $u$ 插入桶 $B(v)$；
- `DELETE(u)`：用句柄从同一路径上的桶中删除。

每次更新触及 $O(\log N)$ 个桶；桶同样使用“致密数组 + 句柄”，支持桶内均匀抽样与 $O(1)$ 删除。

对查询秩区间 $[L,R)$，令 $\mathcal{U}(L,R)$ 为该区间在段树上的规范分解结点集合（结点段两两不交，且并为 $[L,R)$），满足 $|\mathcal{U}(L,R)|=O(\log N)$。

**引理 2.6（range 答案的无重并）**  
对任意秩区间 $[L,R)$，
$$
\{\,u\ \text{活跃}\mid \mathrm{rank}(u)\in[L,R)\,\}
=
\biguplus_{v\in\mathcal{U}(L,R)} B(v).
$$

**证明**  
取任意活跃对象 $u$ 且 $\mathrm{rank}(u)\in[L,R)$。由于 $\mathcal{U}(L,R)$ 的结点段不交且覆盖 $[L,R)$，存在唯一结点 $v^\star\in\mathcal{U}(L,R)$ 满足 $\mathrm{rank}(u)\in\mathrm{seg}(v^\star)$。  
另一方面，$v^\star$ 是秩叶 $\mathrm{rank}(u)$ 的祖先结点，且插入时沿叶到根路径把 $u$ 加入到所有祖先桶中，因此 $u\in B(v^\star)$。唯一性由 $v^\star$ 的唯一性保证。证毕。

### 2.3.3 三原语与均匀采样

给定开区间 $(\ell,r)$，计算秩区间 $[L,R)$ 与分解集合 $\mathcal{U}(L,R)$，记
$$
C(\ell,r):=\sum_{v\in\mathcal{U}(L,R)} |B(v)|.
$$

- `COUNT(ℓ,r)`：返回 $C(\ell,r)$；
- `REPORT(ℓ,r)`：枚举所有 $v\in\mathcal{U}(L,R)$ 的桶 $B(v)$ 并输出（由引理 2.6 不重）；
- `SAMPLE(ℓ,r,k)`：若 $C(\ell,r)=0$ 返回空序列；否则对每个样本位置独立执行：
  1) 以权重 $|B(v)|$ 在 $v\in\mathcal{U}(L,R)$ 上抽取结点；  
  2) 在选中桶内均匀抽取一个对象。

均匀性与独立性证明与定理 2.5 完全同构。

### 2.3.4 复杂度

- `INSERT/DELETE`：$O(\log N)$；
- `COUNT`：$O(\log N)$；
- `REPORT`：$O(\log N + K)$；
- `SAMPLE(ℓ,r,k)`：构建 alias 表 $O(\log N)$，随后 $k$ 次抽样 $O(k)$，总计 $O(\log N + k)$。

---

## 2.4 固定模式的多维递归结构 $\mathcal{D}_{m,g}$

本节把 2.2 与 2.3 的一维结构递归组合，支持对固定模式 $g\in\{A,B\}^m$ 的动态 `COUNT/REPORT/SAMPLE`。

### 2.4.1 固定模式查询问题

维护动态集合 $S$（对侧活跃集的投影盒集合），每个元素为
$$
s^\downarrow=\prod_{j=1}^{m}[L'_j(s),R'_j(s)).
$$
给定查询盒 $q^\downarrow$ 与模式 $g=(g_1,\dots,g_m)\in\{A,B\}^m$，定义模式答案集 $S_g(q)$（见 2.1.2）。目标是在动态更新下支持：

- `INSERT(s)`, `DELETE(s)`；
- `COUNT_g(q)=|S_g(q)|`；
- `REPORT_g(q)`：枚举 $S_g(q)$；
- `SAMPLE_g(q,k)`：在 $S_g(q)$ 上返回 $k$ 个 i.i.d. 均匀样本（有放回）。

### 2.4.2 递归定义

记 $g_{2..m}=(g_2,\dots,g_m)$，并定义剩余维投影
$$
s^{\downarrow\downarrow}:=\prod_{j=2}^{m}[L'_j(s),R'_j(s)),\qquad
q^{\downarrow\downarrow}:=\prod_{j=2}^{m}[L'_j(q),R'_j(q)).
$$

**基例 $m=1$**

- 若 $g_1=A$：使用 2.2 的 stabbing 结构，查询点为 $x=L'_1(q)$，区间为 $[L'_1(s),R'_1(s))$；
- 若 $g_1=B$：使用 2.3 的 range 结构，键为 $a(s)=L'_1(s)$，查询开区间为 $(L'_1(q),R'_1(q))$。

**递归步 $m>1$**

在第 1 个投影维度 $j=1$ 上建立外层结构，并在外层结点上挂接内层结构处理维度 $2..m$：

- 若 $g_1=A$：外层为 stabbing 段树（2.2），其中每个外层结点 $v$ 维护一个内层结构
  $$
  \mathcal{D}_{m-1,g_{2..m}}(v),
  $$
  用于索引所有被分配到 $v$ 的对象集合的 $s^{\downarrow\downarrow}$。

- 若 $g_1=B$：外层为秩化 range 段树（2.3），同样每个外层结点 $v$ 维护一个内层结构 $\mathcal{D}_{m-1,g_{2..m}}(v)$。

外层结构如何将对象“分配到结点”，遵循其一维更新规则：

- $g_1=A$：对象 $s$ 被插入到其区间 $[L'_1(s),R'_1(s))$ 的规范覆盖结点集合 $\mathcal{C}_1(s)$（大小 $O(\log N)$）；
- $g_1=B$：对象 $s$ 被插入到其左端点秩叶到根的路径结点集合 $\mathsf{Path}_1(s)$（大小 $O(\log N)$）。

据此定义统一的更新过程：

- `INSERT(s)`：对所有外层结点 $v\in\mathcal{N}_1(s)$（对应 $\mathcal{C}_1(s)$ 或 $\mathsf{Path}_1(s)$），调用内层 `INSERT(s^{\downarrow\downarrow})`；
- `DELETE(s)`：对所有 $v\in\mathcal{N}_1(s)$，调用内层 `DELETE(s^{\downarrow\downarrow})`。

查询时，外层结构给出一个结点集合 $\mathcal{V}_1(q)$：

- $g_1=A$：$\mathcal{V}_1(q)=\mathsf{Path}(L'_1(q))$；
- $g_1=B$：$\mathcal{V}_1(q)=\mathcal{U}(L,R)$，其中 $[L,R)$ 来自 2.3.1 对 $(L'_1(q),R'_1(q))$ 的秩映射。

对每个 $v\in\mathcal{V}_1(q)$，在内层结构上查询 $q^{\downarrow\downarrow}$，得到结点内答案集 $S_{g,v}(q)$。于是全体答案可写成
$$
S_g(q)=\bigcup_{v\in\mathcal{V}_1(q)} S_{g,v}(q).
$$

### 2.4.3 正确性：无重分解

关键事实是：外层结点集合 $\mathcal{V}_1(q)$ 将所有满足第 1 维模式条件的对象划分为互不重叠的“块”，每个对象在外层唯一归属到其中某个结点，然后由内层结构筛选剩余维。

**定理 2.7（固定模式答案的无重分解）**  
对任意 $m\ge 1$、模式 $g\in\{A,B\}^m$ 与查询盒 $q^\downarrow$，
$$
S_g(q)=\biguplus_{v\in\mathcal{V}_1(q)} S_{g,v}(q),
$$
并且递归结构 $\mathcal{D}_{m,g}$ 返回的结果恰为 $S_g(q)$。

**证明（归纳）**  

- 基例 $m=1$：  
  - $g_1=A$ 时，外层即 stabbing 结构，推论 2.4 给出无重并；  
  - $g_1=B$ 时，外层即 range 结构，引理 2.6 给出无重并。  
  因而成立。

- 归纳步 $m>1$：  
  外层第 1 维按 $g_1$ 使用 stabbing 或 range。对任意满足第 1 维模式条件的对象 $s$：  
  - 若 $g_1=A$，由引理 2.3，$s$ 在查询路径 $\mathcal{V}_1(q)$ 上唯一归属到某个结点桶；  
  - 若 $g_1=B$，由引理 2.6，$s$ 在区间分解结点集合 $\mathcal{V}_1(q)$ 中唯一归属到某个桶。  
  因此外层结点间对象集合互不重叠。固定某个外层结点 $v$，其内层结构索引的是恰好分配到 $v$ 的对象，并按归纳假设在剩余 $m-1$ 维上正确筛选得到 $S_{g,v}(q)$。将所有结点 $v\in\mathcal{V}_1(q)$ 的结果并起来，得到 $S_g(q)$，且不重。证毕。

### 2.4.4 三原语与批量采样的实现方式

由定理 2.7，无重分解直接给出三原语：

- `COUNT_g(q)`：对每个 $v\in\mathcal{V}_1(q)$ 调用内层计数并求和：
  $$
  \mathrm{COUNT}_g(q)=\sum_{v\in\mathcal{V}_1(q)} \mathrm{COUNT}_{m-1,g_{2..m}}(v,q^{\downarrow\downarrow}).
  $$

- `REPORT_g(q)`：对每个 $v\in\mathcal{V}_1(q)$ 调用内层 `REPORT` 并拼接输出（无重复）。

- `SAMPLE_g(q,k)`：令
  $$
  c_v:=|S_{g,v}(q)|,\qquad C:=\sum_{v\in\mathcal{V}_1(q)} c_v.
  $$
  若 $C=0$ 返回空序列。否则定义每个样本位置的抽样过程为：
  $$
  \Pr[V_t=v]=\frac{c_v}{C},
  \qquad
  \Pr[\text{sample}=s\mid V_t=v]=\frac{1}{c_v}\ \ (s\in S_{g,v}(q)).
  $$

为了将批量采样实现为 $O(1)$ 每次抽样的代价，并保持输出序列的 i.i.d. 分布，采用“先指派结点、再分组采样、最后按位置写回”的流程：

1. 计算所有 $c_v$ 并构建对 $\mathcal{V}_1(q)$ 的 alias 表；
2. 对每个样本位置 $t=1,\dots,k$ 独立抽取结点 $V_t$，并记录该位置属于哪个结点；
3. 对每个结点 $v$，设其被指派的位置集合为 $P_v$，大小为 $k_v=|P_v|$，调用内层结构一次 `SAMPLE_{m-1,g_{2..m}}(v,q^{\downarrow\downarrow},k_v)` 得到 $k_v$ 个样本；
4. 将这 $k_v$ 个样本依次写回到位置集合 $P_v$ 中。

该实现与逐次独立采样在分布上完全一致。

**定理 2.8（固定模式采样的均匀性与独立性）**  
`SAMPLE_g(q,k)` 产生的输出序列 $(X_1,\dots,X_k)$ 在 $S_g(q)$ 上独立同分布且均匀：
$$
\Pr[X_t=s]=\frac{1}{|S_g(q)|}\quad(\forall s\in S_g(q),\ \forall t),
$$
并且 $X_1,\dots,X_k$ 相互独立。

**证明**  
由定理 2.7，$S_g(q)$ 为各块 $S_{g,v}(q)$ 的无重并。任取 $s\in S_g(q)$，其唯一属于某个块 $v^\star$。单次抽样命中 $s$ 的概率为
$$
\Pr[s]
=\Pr[V=v^\star]\cdot\Pr[s\mid V=v^\star]
=\frac{c_{v^\star}}{C}\cdot\frac{1}{c_{v^\star}}
=\frac{1}{C}
=\frac{1}{|S_g(q)|}.
$$
对不同样本位置，结点指派与块内抽样均使用独立随机性，并且有放回，故样本独立同分布。证毕。

### 2.4.5 复杂度与空间

令 $N$ 为该类对象的静态全集大小（活跃集是其动态子集）。各层段树高度均为 $O(\log N)$（坐标压缩后端点规模与 $N$ 同阶时成立；否则将 $\log N$ 替换为对应维度端点规模的对数即可）。

**定理 2.9（固定模式时间复杂度）**  
固定模式 $g\in\{A,B\}^m$ 下，结构 $\mathcal{D}_{m,g}$ 支持：

- `INSERT/DELETE`：$O((\log N)^m)$；
- `COUNT`：$O((\log N)^m)$；
- `REPORT`：$O((\log N)^m + K)$；
- `SAMPLE(k)`：$O((\log N)^m + k)$（使用 alias 抽样）。

**证明（递推）**  
外层一次操作触及 $O(\log N)$ 个结点，每个结点上调用一次 $m-1$ 维内层操作，因此
$$
T(m)=O(\log N)\cdot T(m-1),\qquad T(1)=O(\log N),
$$
解得 $T(m)=O((\log N)^m)$。  
对 `SAMPLE(k)`：计算所有 $c_v$ 的代价已包含在 $O((\log N)^m)$ 中；对结点指派构建 alias 表与抽取 $k$ 次分别是 $O(\log N)$ 与 $O(k)$，并对各结点内层批量采样的总代价为 $O((\log N)^m + k)$。证毕。

**定理 2.10（固定模式空间复杂度）**  
固定模式 $g$ 的总存储为
$$
O\!\left(N(\log N)^m\right).
$$

**证明（按元素出现次数计数）**  
在一维外层，元素被复制到 $O(\log N)$ 个桶结点；递归到下一维时，每一份复制再进入对应结点的内层结构并继续产生 $O(\log N)$ 倍的复制。由此每个元素在 $m$ 层递归中总出现次数为 $O((\log N)^m)$，从而桶中存储的元素副本总数为 $O(N(\log N)^m)$。  
实现上，内层结构与桶均按需创建：仅当某结点首次接收元素时才分配其桶与（必要的）子结构节点，因此结构节点开销同阶受元素副本数上界控制。证毕。

---

## 2.5 事件块原语：跨模式聚合与扫描对接

本节把固定模式结构汇总为“其余维相交”的完整 `COUNT/REPORT/SAMPLE`，并说明如何嵌入扫描线过程。

### 2.5.1 partner 集的模式无重并

在某个 $\mathrm{START}(q)$ 事件处，对侧活跃集投影为 $S$，事件块的 partner 集为
$$
K(q):=\{\,s\in S\mid s^\downarrow\cap q^\downarrow\neq\varnothing\,\}=S_\cap(q).
$$
由定理 2.2，
$$
K(q)=\biguplus_{g\in\{A,B\}^m} S_g(q).
$$
定义模式权重
$$
w^{(g)}(q):=|S_g(q)|,\qquad w(q):=|K(q)|=\sum_{g} w^{(g)}(q).
$$

### 2.5.2 跨模式的 COUNT / REPORT / SAMPLE

为每个模式 $g$ 维护一套结构 $\mathcal{D}_{m,g}$（对象集合为对侧活跃集的投影盒）。则：

- `COUNT(q)`：
  $$
  \mathrm{COUNT}(q)=\sum_{g} \mathrm{COUNT}_g(q).
  $$

- `REPORT(q)`：输出所有 $\mathrm{REPORT}_g(q)$ 的拼接即可；由于模式间无重，不会重复。

- `SAMPLE(q,k)`：  
  先计算所有模式计数 $w^{(g)}(q)$ 并令 $W=\sum_g w^{(g)}(q)$。若 $W=0$ 返回空序列。  
  对每个样本位置 $t=1,\dots,k$ 独立抽取模式 $G_t$ 满足
  $$
  \Pr[G_t=g]=\frac{w^{(g)}(q)}{W}.
  $$
  条件于 $G_t=g$，调用 `SAMPLE_g(q,1)` 返回一个在 $S_g(q)$ 上均匀的样本。  

  为提高批量效率，使用与 2.4.4 相同的“按位置指派、再分组采样、最后写回”的流程：  
  1) 在模式集合上用权重 $w^{(g)}(q)$ 构建 alias 表；  
  2) 对每个位置 $t$ 抽取 $G_t$ 并记录位置集合；  
  3) 对每个模式 $g$，令其被指派次数为 $k_g$，调用一次 `SAMPLE_g(q,k_g)`；  
  4) 将结果写回对应位置。

**命题 2.11（跨模式采样的均匀性与独立性）**  
`SAMPLE(q,k)` 返回的输出序列在 $K(q)$ 上独立同分布且均匀。

**证明**  
由定理 2.2 的无重分解，任取 $s\in K(q)$，其唯一属于某个模式 $g^\star$。单次采样命中 $s$ 的概率为
$$
\Pr[s]
=\Pr[G=g^\star]\cdot\Pr[s\mid G=g^\star]
=\frac{w^{(g^\star)}(q)}{W}\cdot\frac{1}{w^{(g^\star)}(q)}
=\frac{1}{W}
=\frac{1}{|K(q)|}.
$$
每个位置的模式抽取与模式内抽样均独立，且有放回，因此得到 i.i.d. 均匀样本序列。证毕。

### 2.5.3 与扫描线过程的对接

扫描过程中分别对两类对象维护对侧索引。对于任意时刻的活跃集：

- 以结构族 $\{\mathcal{D}^{(c)}_{m,g}\}_{g}$ 维护当前活跃的 $R_c$ 投影集合；
- 以结构族 $\{\mathcal{D}^{(\bar c)}_{m,g}\}_{g}$ 维护当前活跃的 $R_{\bar c}$ 投影集合。

处理事件时遵循如下顺序：

- 若事件为 $\mathrm{END}(r)$：将 $r$ 从其所属类的所有模式结构中 `DELETE`；
- 若事件为 $\mathrm{START}(q)$：先用对侧类的结构族对 $q^\downarrow$ 执行 `COUNT/REPORT/SAMPLE` 得到事件块 partner 信息，再将 $q$ 插入其所属类结构族（`INSERT`）。

该顺序保证查询仅依赖于“先于当前 START 已处于活跃状态”的对侧对象，与半开区间语义一致。

### 2.5.4 总复杂度（常数维）

模式数为 $2^m$。当 $d$ 视为常数时，$m=d-1$ 亦为常数，因此跨模式聚合仅引入常数因子。

对单次事件块查询（对侧全集规模为 $N$）：

- `COUNT`：$O\!\left(2^m(\log N)^m\right)$；
- `REPORT`：$O\!\left(2^m(\log N)^m + K\right)$；
- `SAMPLE(k)`：$O\!\left(2^m(\log N)^m + k\right)$（使用 alias 抽样）。

其中 $K$ 为输出规模。

# 第 3 章 从事件块原语到 i.i.d. 均匀 Join 采样的三种采样框架

本章讨论如下问题：给定两类 $d$ 维轴对齐半开盒集合 $R_c$ 与 $R_{\bar c}$，令
$$
J=\{(r_c,r_{\bar c})\mid r_c\in R_c,\ r_{\bar c}\in R_{\bar c},\ r_c\cap r_{\bar c}\neq\varnothing\}
$$
为跨集合、有序方向固定的相交对集合。目标是在 $J$ 上生成 $t$ 个样本 $Z_1,\dots,Z_t$，满足 **独立同分布（i.i.d.）且均匀（有放回）**。

本章给出三套从“事件块原语”出发的采样框架，分别在扫描次数、内存与原语调用代价之间提供不同取舍：

| 框架 | 扫描次数 | 额外空间 | 主要优点 | 主要代价 |
|---|---:|---:|---|---|
| 框架 I：显式枚举 | 1 | $\Theta(|J|)$ | 逻辑最直接 | $|J|$ 可能极大 |
| 框架 II：两遍回填 | 2 | $O(M+t)$ | 无需存 $J$，分布精确 | 需要可重放第二遍扫描 |
| 框架 III：预算自适应 | $\le 2$ | $\le B + O(M+t)$ | 在预算 $B$ 下尽量减少第二遍工作，可一遍结束 | 需要缓存策略与代价模型（仅影响性能，不影响分布） |

其中 $M$ 表示 START 事件数量（等于 $|R_c|+|R_{\bar c}|$）。

---

## 3.1 采样目标、扫描事件与事件块

### 3.1.1 半开盒、相交与 Join 目标分布

每个盒 $r$ 表示为
$$
r=\prod_{k=1}^d [L_k(r),R_k(r)),\qquad L_k(r)<R_k(r).
$$
半开语义下，两盒 $r,s$ 相交当且仅当对所有维度 $k$ 同时满足
$$
L_k(r) < R_k(s)\quad\wedge\quad L_k(s) < R_k(r).
$$
严格不等号确保“贴边”（例如 $R_k(r)=L_k(s)$）不算相交。

采样目标：若 $|J|>0$，输出 $Z_1,\dots,Z_t\in J$，满足
$$
\Pr\{Z_j=P\}=\frac{1}{|J|}\quad(\forall P\in J,\ \forall j),
\qquad Z_1,\dots,Z_t\ \text{相互独立}.
$$
若 $|J|=0$，输出空序列。

---

### 3.1.2 扫描维、事件总序与活跃集

固定第 1 维为扫描维。对任意盒 $r\in R_c\cup R_{\bar c}$ 定义两个事件：
- $\mathrm{START}(r)$，坐标 $x=L_1(r)$；
- $\mathrm{END}(r)$，坐标 $x=R_1(r)$。

对所有事件定义一个可复现的总序 $\prec$：
1. 坐标小者先；
2. 坐标相同则 **END 严格先于 START**；
3. 若同为 START（或同为 END）且坐标相同，则按唯一 id（或任意固定全序）打破平局。

扫描过程中维护两侧活跃集 $A_c,A_{\bar c}$：
- 处理 $\mathrm{END}(r)$：将 $r$ 从其所属活跃集中删除；
- 处理 $\mathrm{START}(r)$：先执行事件块查询（定义见下节），再将 $r$ 插入其所属活跃集。

上述规则与半开语义一致：当 $R_1(r)=L_1(s)$ 时，$\mathrm{END}(r)\prec \mathrm{START}(s)$，因此 $r$ 不会被视为 $s$ 的活跃候选，从而贴边不产生相交对。

---

### 3.1.3 partner 集与事件块

记剩余维度数 $m=d-1$。对盒 $r$ 定义其余维投影
$$
r^\downarrow=\prod_{k=2}^d [L_k(r),R_k(r)).
$$

令 $E^+$ 为所有 START 事件集合，并按 $\prec$ 顺序编号
$$
E^+=\{e_1,e_2,\dots,e_M\},\qquad e_i=\mathrm{START}(q_i).
$$

**引理 3.0（扫描维自动满足相交的一半）**  
当处理 START 事件 $e=\mathrm{START}(q)$（坐标 $x=L_1(q)$）时，任意对侧活跃盒 $s$ 都满足
$$
L_1(s)<R_1(q)\quad\text{且}\quad L_1(q)<R_1(s).
$$

**证明**  
由于 $s$ 在对侧活跃集中，其 START 已发生且 END 尚未发生，故 $L_1(s)\le x$ 且 $x<R_1(s)$。又 $R_1(q)>L_1(q)=x$（盒非空），于是 $L_1(s)\le x<R_1(q)$，从而 $L_1(s)<R_1(q)$；同时 $L_1(q)=x<R_1(s)$ 亦成立。证毕。

因此，在 START 事件处判断相交只需检查其余维投影是否相交。

**定义（partner 集合）**  
当处理 START 事件 $e=\mathrm{START}(q)$ 时：
- 若 $q\in R_c$，定义
  $$
  K_e=\{s\in A_{\bar c}\mid q^\downarrow\cap s^\downarrow\neq\varnothing\};
  $$
- 若 $q\in R_{\bar c}$，定义
  $$
  K_e=\{s\in A_{c}\mid q^\downarrow\cap s^\downarrow\neq\varnothing\}.
  $$

为统一输出方向恒为 $(r_c,r_{\bar c})$，定义映射 $\Phi_e$：
- 若 $q\in R_c$，则对 $s\in K_e$，$\Phi_e(s)=(q,s)$；
- 若 $q\in R_{\bar c}$，则对 $s\in K_e$，$\Phi_e(s)=(s,q)$。

**定义（事件块）**  
对每个 START 事件 $e\in E^+$ 定义事件块与权重
$$
J_e=\{\Phi_e(s)\mid s\in K_e\}\subseteq J,\qquad w_e:=|J_e|=|K_e|.
$$

下面给出 $J$ 的无重分解。

**定理 3.1（事件块无重分解）**  
有
$$
J=\biguplus_{e\in E^+} J_e,\qquad |J|=\sum_{e\in E^+} w_e,
$$
其中 $\biguplus$ 表示不交并。

**证明**  
*覆盖性*：任取 $(r_c,r_{\bar c})\in J$。考虑二者的 START 事件 $e_c=\mathrm{START}(r_c)$ 与 $e_{\bar c}=\mathrm{START}(r_{\bar c})$，令
$$
e^\star=\max_{\prec}\{e_c,e_{\bar c}\}.
$$
设 $e^\star=\mathrm{START}(q)$，另一侧盒为 $s$（即 $\{q,s\}=\{r_c,r_{\bar c}\}$）。由于 $e^\star$ 为较晚 START，处理 $e^\star$ 时 $s$ 的 START 已发生。又两盒在第 1 维相交，结合 END-before-START 的事件序，$s$ 在时刻 $L_1(q)$ 尚未 END，必在对侧活跃集中。并且两盒其余维投影相交，故 $s\in K_{e^\star}$，于是 $\Phi_{e^\star}(s)=(r_c,r_{\bar c})\in J_{e^\star}$。

*互斥性*：对任意相交对 $(r_c,r_{\bar c})$，它只可能在两者 START 中较晚的那一次被“看见”。在较早 START 时，较晚者尚未插入活跃集，不可能成为 partner；因此该对不可能同时落入两个不同事件块。证毕。

---

### 3.1.4 事件块原语与随机性约定

本章三套框架仅依赖如下原语语义。对每个 START 事件 $e$（即对每个 $K_e$）假定可调用：

- `COUNT(e)`：返回 $w_e=|K_e|$；
- `REPORT(e)`：枚举 $K_e$ 中每个元素一次（顺序任意）；
- `SAMPLE(e,k)`：返回 $k$ 个在 $K_e$ 上 **i.i.d. 均匀（有放回）** 的样本。

为严格保证全局 i.i.d.，需要不同阶段与不同调用使用相互独立的随机性。一个可复现的实现方式是：给每个事件 $e_i$ 分配整数编号 $i$，并为每个阶段（如“位置指派”“预取采样”“残差采样”等）分配固定的 `phase_id`，用哈希派生子种子：
$$
\texttt{seed}(i,\texttt{phase\_id},\texttt{ctr})=\mathrm{Hash}(\texttt{master\_seed},i,\texttt{phase\_id},\texttt{ctr}),
$$
其中 `ctr` 为该阶段的调用计数器。这样可在不共享随机流的前提下保证可重现实验。

---

## 3.2 事件块混合采样：全局均匀与 i.i.d.

令
$$
W:=|J|=\sum_{e\in E^+} w_e.
$$
若 $W=0$，则 $J=\varnothing$。

否则考虑以下单次输出过程：
1. 按分布 $\Pr\{E=e\}=w_e/W$ 抽取一个 START 事件 $E$；
2. 在 $K_E$ 上均匀抽取一个 partner $S$（有放回），输出 $Z=\Phi_E(S)$。

**引理 3.2（单次输出均匀）**  
上述过程输出的 $Z$ 在 $J$ 上均匀：
$$
\Pr\{Z=P\}=\frac{1}{|J|},\qquad \forall P\in J.
$$

**证明**  
任取 $P\in J$。由定理 3.1，存在唯一事件 $e(P)$ 使 $P\in J_{e(P)}$。因此
$$
\Pr\{Z=P\}
=\Pr\{E=e(P)\}\cdot\Pr\{Z=P\mid E=e(P)\}
=\frac{w_{e(P)}}{W}\cdot\frac{1}{w_{e(P)}}
=\frac{1}{W}
=\frac{1}{|J|}.
$$
证毕。

**引理 3.3（多次输出 i.i.d.）**  
若独立重复上述两步 $t$ 次（每次事件选择与块内抽样使用独立随机性，且为有放回），得到 $Z_1,\dots,Z_t$，则它们相互独立且同分布为 $J$ 上的均匀分布。

---

## 3.3 框架 I：显式枚举 Join 后均匀抽样

框架 I 在一次扫描中显式构造 $J$ 的全部元素，再在数组上独立均匀抽样。

### 3.3.1 算法

```text
Algorithm 3.1  Materialize-and-Sample
Input: R_c, R_{\bar c}, sample size t
Output: Z_1..Z_t

1: Pairs ← empty dynamic array
2: Scan all events in order ≺, maintaining active sets
3: For each START event e:
4:     For each s in REPORT(e):            // enumerates K_e
5:         Pairs.append( Φ_e(s) )          // append one element of J_e
6: If |Pairs| = 0: return empty sequence
7: For j = 1..t:
8:     Draw I_j uniformly from {1..|Pairs|} (with replacement)
9:     Z_j ← Pairs[I_j]
10: return Z_1..Z_t
```

### 3.3.2 正确性

**定理 3.4（框架 I 的正确性）**  
算法 3.1 输出 $Z_1,\dots,Z_t$ 为 $J$ 上 i.i.d. 均匀（有放回）样本。

**证明**  
扫描阶段对每个事件 $e$ 追加 $J_e$ 的全部元素。由定理 3.1，`Pairs` 恰为 $J$ 的无重枚举。随后对 `Pairs` 做独立均匀索引抽样即在 $J$ 上独立均匀抽样。证毕。

### 3.3.3 代价

时间与空间均由 $|J|$ 主导：时间 $\Theta(|J|)+O(t)$，空间 $\Theta(|J|)$（加上扫描维护开销）。

---

## 3.4 框架 II：两遍扫描的精确计数—位置指派—回填生成

框架 II 不显式构造 $J$，而是用第一遍扫描精确获得每个事件块权重 $w_i$，再按 $w_i/W$ 为每个输出位置指派事件索引，第二遍仅对被指派到的事件块生成所需样本并回填。

### 3.4.1 算法

将 START 事件按 $\prec$ 顺序记为 $e_1,\dots,e_M$，令 $w_i:=w_{e_i}$。

```text
Algorithm 3.2  Two-Pass Count-and-Backfill
Input: R_c, R_{\bar c}, sample size t
Output: Z_1..Z_t

// Pass 1: exact block weights
1: Scan events in order ≺
2: For each START event e_i:
3:     w_i ← COUNT(e_i)
4: W ← sum_{i=1}^M w_i
5: If W = 0: return empty sequence

// Planning: assign each output position to an event index
6: Build a discrete sampler over {1..M} with Pr{I=i} = w_i / W
7: For i = 1..M: L_i ← empty list
8: For j = 1..t:
9:     Draw I_j independently from Pr{I=i} = w_i / W
10:    Append position j to L_{I_j}

// Pass 2: backfill samples by blocks
11: Initialize output array Z[1..t]
12: Scan events in order ≺ again (same tie-break rules)
13: For each START event e_i:
14:     k_i ← |L_i|
15:     If k_i = 0: continue
16:     Draw partners s_1..s_{k_i} ← SAMPLE(e_i, k_i)   // i.i.d. on K_{e_i}
17:     For ℓ = 1..k_i:
18:         j ← L_i[ℓ]
19:         Z[j] ← Φ_{e_i}( s_ℓ )
20: return Z[1..t]
```

**可重放约束**  
第二遍扫描必须使用与第一遍完全相同的事件总序 $\prec$（包括坐标相等时 END-before-START 与 tie-break 规则），以保证 $K_{e_i}$ 的定义在两遍中一致。

---

### 3.4.2 正确性

**定理 3.5（框架 II 的 i.i.d. 均匀性）**  
若 `COUNT` 返回精确 $w_i$，`SAMPLE(e_i,k_i)` 返回 $K_{e_i}$ 上 i.i.d. 均匀（有放回）样本，且各阶段随机性独立，则算法 3.2 输出 $Z_1,\dots,Z_t$ 在 $J$ 上 i.i.d. 均匀（有放回）。

**证明**  
规划阶段对每个位置 $j$ 独立抽取事件索引 $I_j$，满足 $\Pr\{I_j=i\}=w_i/W$。对任意 $j$，条件于 $I_j=i$，回填阶段令 $Z_j=\Phi_{e_i}(S)$，其中 $S$ 在 $K_{e_i}$ 上均匀（有放回）。因此由引理 3.2，$Z_j$ 的边缘分布在 $J$ 上均匀。

独立性来自两点：$\{I_j\}$ 相互独立；对每个事件 $i$，`SAMPLE(e_i,k_i)` 输出 $k_i$ 个独立样本，并且不同事件的采样调用使用独立随机性。位置回填是确定性写入，不引入相关性。因此 $Z_1,\dots,Z_t$ 相互独立且同分布，证毕。

---

### 3.4.3 原语调用代价（以次数计）

- Pass 1：对每个 START 事件一次 `COUNT`，共 $M$ 次；
- 规划阶段：$t$ 次离散抽样；
- Pass 2：对所有 $k_i>0$ 的事件各调用一次 `SAMPLE(e_i,k_i)`；总样本数 $\sum_i k_i=t$。

额外空间主要为数组 $w[1..M]$ 与位置列表总大小 $t$。

---

## 3.5 框架 III：预算约束下的自适应生成（预取 + 缓存 + 流式回填）

框架 III 在不改变目标分布的前提下，引入预算 $B$（以“可缓存的 partner 记录条数”计量），通过第一遍扫描的缓存与预取来尽量减少第二遍仍需执行的原语工作量，并在缓存足够时实现“一遍结束”。

### 3.5.1 预算与两类缓存

对每个事件 $e_i$，允许保存两种对象（可同时存在）：

1. **全量缓存** $\mathcal C_i$：保存完整 partner 列表 $K_{e_i}$（由一次 `REPORT(e_i)` 得到），占用 $|\mathcal C_i|=w_i$ 单位预算。此后可以通过均匀随机访问在 $K_{e_i}$ 上进行任意次数的均匀抽样。

2. **预取样本缓存** $\mathcal S_i$：保存 $s_i$ 个由 `SAMPLE(e_i,s_i)` 得到的 i.i.d. 样本，占用 $|\mathcal S_i|=s_i$ 单位预算。

总预算约束为
$$
\sum_i |\mathcal C_i|+\sum_i |\mathcal S_i|\le B.
$$

> 重要事实：缓存策略只影响生成过程的代价，不影响分布正确性。分布正确性只依赖于：位置指派按 $w_i/W$；每个被填入的位置对应的 partner 是从 $K_{e_i}$ 上独立均匀（有放回）生成。

---

### 3.5.2 预算分配的离线基准：边际收益与最优槽位选择

规划阶段把 $t$ 个输出位置按分布 $p_i=w_i/W$ 指派给事件块。令
$$
k_i := |\mathcal L_i|
$$
为事件 $i$ 被指派的位置数。向量 $(k_1,\dots,k_M)$ 服从多项分布 $\mathrm{Multinomial}(t;p_1,\dots,p_M)$，从而每个边缘 $k_i\sim\mathrm{Binomial}(t,p_i)$。

若第一遍为事件 $i$ 预取 $s_i$ 个样本，则需在第二遍补齐的残差为
$$
r_i = (k_i-s_i)_+.
$$
一个自然目标是使 $\mathbb E[\sum_i r_i]$ 尽可能小。

**引理 3.6（边际收益）**  
对固定事件 $i$，将预取数从 $s$ 增加到 $s+1$ 带来的期望残差减少量为
$$
\Delta_i(s+1)
=\mathbb E[(k_i-s)_+]-\mathbb E[(k_i-(s+1))_+]
=\Pr(k_i\ge s+1).
$$

**证明**  
利用恒等式 $(k-s)_+=\sum_{r=s+1}^{t}\mathbf 1[k\ge r]$，取期望得
$$
\mathbb E[(k_i-s)_+]=\sum_{r=s+1}^{t}\Pr(k_i\ge r),
$$
相邻相减即得结论。证毕。

由于 $\Pr(k_i\ge r)$ 随 $r$ 单调递减，每个事件的预取具有“边际收益递减”。

**定理 3.7（离线最优：选择最大价值的 $B$ 个槽位）**  
若已知所有 $p_i$，且预算全部用于预取样本（每单位预算对应 1 个预取样本），则使 $\mathbb E[\sum_i (k_i-s_i)_+]$ 最小的整数解可表述为：从槽位集合
$$
\{(i,r)\mid i\in\{1,\dots,M\},\ r\in\{1,\dots,t\}\}
$$
中选择恰好 $B$ 个槽位，使总价值最大，其中槽位 $(i,r)$ 的价值为
$$
v_{i,r}=\Pr(k_i\ge r).
$$
并且最优解对每个事件 $i$ 取一个前缀：若选择了 $(i,r)$，则必选择所有 $(i,1),\dots,(i,r-1)$。

**证明要点**  
由引理 3.6，每增加 1 个预取样本相当于选择一个槽位，其价值等于期望残差的减少量。预算为 $B$ 时，总减少量是所选槽位价值之和，因此最优策略为选择价值最大的 $B$ 个槽位。由于 $v_{i,r}$ 对 $r$ 递减，最优解对每个事件必为前缀。证毕。

---

### 3.5.3 在线预取与可撤销的槽位堆

第一遍扫描时尚未知 $W$，从而 $p_i=w_i/W$ 不可用。本节给出一个在线策略：用估计值构造每个槽位的“价值评分”，并用一个全局堆维护在预算内优先保留的槽位。该策略仅用于性能优化；分布正确性不依赖评分是否准确。

**在线估计与评分**  
扫描到第 $i$ 个 START 事件时已知 $w_i$ 与前缀和 $W_{\text{sofar}}=\sum_{\ell\le i} w_\ell$。可用任意外推得到 $\widehat W_i$，例如
$$
\widehat W_i = W_{\text{sofar}} + (M-i)\cdot \frac{W_{\text{sofar}}}{i}.
$$
并令
$$
\widehat p_i = \frac{w_i}{\widehat W_i},\qquad \widehat\mu_i=t\widehat p_i.
$$
对槽位 $(i,r)$ 定义一个单调递减的评分序列 $\widehat v_{i,1}\ge \widehat v_{i,2}\ge\cdots$ 用于近似 $v_{i,r}=\Pr(k_i\ge r)$。常用做法是用 Poisson 近似：
$$
k_i\approx \mathrm{Poisson}(\widehat\mu_i),\qquad \widehat v_{i,r}\approx \Pr(\mathrm{Poisson}(\widehat\mu_i)\ge r).
$$
也可以直接用二项尾概率（以 $\widehat p_i$ 代替 $p_i$）或其他更保守的上/下界；选择不影响分布，仅影响预取效果与计算开销。

**可撤销堆维护**  
令 $B_{\text{samp}}=B-\sum_i|\mathcal C_i|$ 为分配给预取样本缓存的剩余预算。维护一个最小堆 $H$，其中每个元素表示一个被保留的预取槽位，记录二元组 $(\widehat v,i)$。并维护计数 $s_i$ 表示事件 $i$ 当前保留的预取槽位数（因此最终应缓存 $s_i$ 个样本）。

堆维护不变式：
- $\sum_i s_i = |H|\le B_{\text{samp}}$；
- 对每个 $i$，被保留的槽位对应前缀 $\{(i,1),\dots,(i,s_i)\}$。

当处理事件 $i$ 时，按 $r=1,2,\dots$ 依次尝试加入槽位 $(i,r)$：
- 若 $|H|<B_{\text{samp}}$，直接加入；
- 若 $|H|=B_{\text{samp}}$ 且 $\widehat v_{i,r}>\min(H)$，则加入并弹出堆顶（全局最小评分槽位）；
- 否则停止（由于 $\widehat v_{i,r}$ 随 $r$ 递减，后续槽位更不可能进入堆）。

被弹出的槽位一定来自某个事件 $j$，则令 $s_j\leftarrow s_j-1$。若事件 $j$ 已有预取样本缓存 $\mathcal S_j$，则从其末尾丢弃一个样本以保持 $|\mathcal S_j|=s_j$。这一“丢弃”仅改变缓存内容，不影响分布正确性（丢弃的规则与样本取值无关，保留样本仍为 i.i.d. 子序列）。

在事件 $i$ 的槽位决策结束后，调用一次
$$
\mathcal S_i\leftarrow \texttt{SAMPLE}(e_i,s_i)
$$
并保存返回的 $s_i$ 个样本。将来若槽位被挤出，同样按上述规则从尾部丢弃即可。

**小块全量缓存**  
当 $w_i$ 小于阈值 $w_{\text{small}}$ 且预算允许时，可设置 $\mathcal C_i=K_{e_i}$（一次 `REPORT` 获取并保存）。当新增全量缓存使 $B_{\text{samp}}$ 下降时，需要反复弹出堆顶直至 $|H|\le B_{\text{samp}}$，并同步减少对应事件的 $s_j$ 与丢弃样本，保持预算不变式。

---

### 3.5.4 仅一次 REPORT 的流式有放回采样

当某事件块在规划后需要补齐 $r$ 个样本时，除了 `SAMPLE(e,r)` 之外，还可以通过一次 `REPORT(e)` 实现 i.i.d. 有放回采样。

```text
Algorithm 3.3  STREAM-REPORT-SAMPLE(e, r, w)
Input: START event e, required samples r, block size w=|K_e|
Output: r i.i.d. samples from K_e (with replacement)

1: Draw u_1..u_r i.i.d. uniformly from {1..w}
2: Sort pairs (u_j, j) by u_j ascending
3: p ← 0
4: For each x in REPORT(e):                 // x is the next element of K_e
5:     p ← p + 1                            // x is the p-th element in this enumeration
6:     While next pair (u_j, j) has u_j = p:
7:         out[j] ← x
8:         advance to next (u_j, j)
9: return out[1..r]
```

**引理 3.8（流式回填正确性）**  
算法 3.3 输出在 $K_e$ 上 i.i.d. 均匀（有放回）。

**证明**  
`REPORT(e)` 给出 $K_e$ 的某个排列 $(x_1,\dots,x_w)$。算法输出 $out[j]=x_{u_j}$，其中 $u_j$ 独立且在 $\{1,\dots,w\}$ 上均匀。因此每个 $out[j]$ 在 $K_e$ 上均匀，且不同 $j$ 独立；重复的 $u_j$ 产生重复输出，对应有放回抽样。证毕。

---

### 3.5.5 完整算法：预算自适应 Join 采样

算法由“第一遍计数与缓存”“位置指派”“缓存回填”“必要时第二遍补齐”组成。

```text
Algorithm 3.4  Adaptive Budgeted Join Sampling (at most two passes)
Input: R_c, R_{\bar c}, sample size t, budget B
Output: Z_1..Z_t

// Pass 1: exact weights + caching/prefetch under budget
1: Initialize arrays w[1..M], C_i ← empty, S_i ← empty, s_i ← 0
2: mem_full ← 0, heap H ← empty
3: Scan events in order ≺
4: For each START event e_i:
5:     w_i ← COUNT(e_i)
6:     // optional: full cache for small blocks
7:     If (w_i ≤ w_small) and (mem_full + w_i ≤ B):
8:         C_i ← REPORT(e_i)                // store full partner list
9:         mem_full ← mem_full + w_i
10:        // shrink prefetch heap capacity if needed
11:        B_samp ← B - mem_full
12:        While |H| > B_samp:
13:            pop (v, j) from H; s_j ← s_j - 1; discard last sample from S_j
14:        continue
15:    // prefetch slots by global heap under remaining capacity
16:    B_samp ← B - mem_full
17:    While true:
18:        r ← s_i + 1
19:        compute score v_hat ← score(i, r)   // monotone decreasing in r
20:        If |H| < B_samp:
21:            push (v_hat, i) into H; s_i ← s_i + 1
22:        Else if v_hat > min(H):
23:            push (v_hat, i) into H; s_i ← s_i + 1
24:            pop (v_min, j) from H; s_j ← s_j - 1; discard last sample from S_j
25:        Else:
26:            break
27:    If s_i > 0:
28:        S_i ← SAMPLE(e_i, s_i)            // store i.i.d. samples (with replacement)

// Planning: assign t output positions to events by exact weights
29: W ← sum_{i=1}^M w_i
30: If W = 0: return empty sequence
31: Build discrete sampler over {1..M} with Pr{I=i} = w_i / W
32: For i = 1..M: L_i ← empty list
33: For j = 1..t:
34:     Draw I_j independently and append j to L_{I_j}

// Fill from caches/prefetch, record residual positions
35: Initialize output array Z[1..t], residual list R_i ← empty for all i
36: For i = 1..M:
37:     If |L_i| = 0: continue
38:     If C_i nonempty:
39:         For each j in L_i:
40:             draw s uniformly from C_i (with replacement)
41:             Z[j] ← Φ_{e_i}(s)
42:     Else:
43:         let a ← min(|L_i|, |S_i|)
44:         For ℓ = 1..a:
45:             j ← L_i[ℓ]
46:             Z[j] ← Φ_{e_i}( S_i[ℓ] )
47:         For ℓ = a+1..|L_i|:
48:             append L_i[ℓ] into R_i

// If no residual, finish in one pass
49: If R_i empty for all i: return Z[1..t]

// Pass 2: generate residual samples only where needed
50: Scan events in order ≺ again
51: For each START event e_i:
52:     If R_i empty or C_i nonempty: continue
53:     r_i ← |R_i|
54:     Choose one:
55:         partners ← SAMPLE(e_i, r_i)
56:         partners ← STREAM-REPORT-SAMPLE(e_i, r_i, w_i)
57:     For ℓ = 1..r_i:
58:         j ← R_i[ℓ]
59:         Z[j] ← Φ_{e_i}( partners[ℓ] )
60: return Z[1..t]
```

**选择 `SAMPLE` 还是流式 `REPORT`**  
可用简单代价比较决定第 55–56 行：设 `SAMPLE` 单样本代价近似为 $c_s$，`REPORT` 枚举单元素代价近似为 $c_r$，则当
$$
r_i\cdot c_s \gtrsim w_i\cdot c_r
$$
更倾向用一次 `REPORT` 的流式生成，否则更倾向直接 `SAMPLE(e_i,r_i)$。该选择只影响代价，不改变分布。

---

### 3.5.6 正确性与预算满足

**定理 3.9（框架 III 的 i.i.d. 均匀性与预算约束）**  
在以下前提下：
1) `COUNT/REPORT/SAMPLE` 满足其语义：`COUNT` 精确、`REPORT` 枚举一次、`SAMPLE` 在 $K_e$ 上返回 i.i.d. 均匀（有放回）样本；  
2) 位置指派与各次采样调用使用相互独立的随机性；  
3) 两遍扫描使用相同事件序 $\prec$；  

则算法 3.4 的输出 $Z_1,\dots,Z_t$ 在 $J$ 上 i.i.d. 均匀（有放回），且缓存占用始终不超过预算 $B$。

**证明**  

（预算）算法只在满足 `mem_full + w_i ≤ B` 时写入全量缓存；预取样本槽位通过堆维护始终满足 $|H|\le B_{\text{samp}}=B-mem_full$，并且每次堆弹出都会同步减少对应事件的样本缓存大小，因此
$$
\sum_i |\mathcal C_i| + \sum_i |\mathcal S_i|
= mem_full + |H| \le mem_full + (B-mem_full)=B.
$$

（分布）规划阶段对每个位置 $j$ 独立抽取事件索引 $I_j$，满足 $\Pr\{I_j=i\}=w_i/W$。因此只需证明：对任意事件 $i$，被填入到 $L_i$ 与 $R_i$ 的 partner 记录都是在 $K_{e_i}$ 上独立均匀（有放回）生成，且与 $\{I_j\}$ 独立。

- 若使用全量缓存 $\mathcal C_i=K_{e_i}$：每次通过均匀随机访问从 $\mathcal C_i$ 取出一个元素（有放回），得到 i.i.d. 均匀样本。
- 若使用预取样本缓存 $\mathcal S_i$：$\mathcal S_i$ 由一次 `SAMPLE(e_i,s_i)` 得到，内部为 i.i.d. 均匀。丢弃尾部样本不改变剩余样本的 i.i.d. 结构。将这些样本按固定顺序写入若干位置不会引入相关性。
- 若存在残差：第二遍对事件 $i$ 生成的 `partners` 来自 `SAMPLE(e_i,r_i)` 或算法 3.3（由引理 3.8），两者都在 $K_{e_i}$ 上输出 $r_i$ 个 i.i.d. 均匀（有放回）样本。该阶段随机性与预取与位置指派独立。

综上，对每个位置 $j$，条件于 $I_j=i$，$Z_j$ 等价于从 $K_{e_i}$ 上独立均匀抽取一个 partner 并经 $\Phi_{e_i}$ 映射得到的元素。由引理 3.2，$Z_j$ 在 $J$ 上均匀；由随机性独立性约定，$Z_1,\dots,Z_t$ 相互独立，从而为 i.i.d. 均匀（有放回）。证毕。

---

### 3.5.7 代价刻画与实践建议（不影响分布）

- 框架 III 的第一遍必做 $M$ 次 `COUNT`；同时可能做部分 `REPORT`（小块全量缓存）与部分 `SAMPLE`（预取样本）。
- 若缓存覆盖了所有被指派的需求（即所有 $R_i$ 为空），则无需第二遍扫描，实现“一遍结束”。
- 第二遍只对 $R_i\neq\varnothing$ 的事件执行补齐操作；对每个事件最多进行一次补齐（选择 `SAMPLE` 或一次 `REPORT` 流式生成）。
- 当 $B$ 很小或代价模型难以估计时，取 $w_{\text{small}}=0$ 且不做预取，框架 III 退化为框架 II；当 $B$ 足够大时，可通过全量缓存或预取使第二遍显著缩短。

---

# 第 4 章 对照方法一：正交范围树实现事件块原语

本章给出一种“直接型”的事件块原语实现路线：不做第 2 章的 $2^m$ 模式化分解，而是将“其余 $m=d-1$ 维的盒相交”整体提升为 $D=2m$ 维点集上的正交范围查询，并用 $D$ 维正交范围树（range tree）在扫描过程中动态维护对侧活跃集，从而在每个 $\mathrm{START}(q)$ 事件处实现：

- `COUNT(e)`：返回事件块 $e$ 的 partner 集大小 $w_e=|K_e|$；
- `REPORT(e)`：枚举 $K_e$ 中每个元素一次（顺序任意）；
- `SAMPLE(e,k)`：在 $K_e$ 上返回 $k$ 个 i.i.d.、均匀、有放回样本。

与第 1 章和第 3 章一致，扫描维固定为第 1 维，事件序采用 END-before-START 以严格匹配半开语义；因此在处理某个 START 事件时，对侧活跃集已经自动满足第 1 维相交的两条严格不等式，事件块判定只需要检查其余 $m$ 维投影是否相交。

---

## 4.1 其余维相交的 $2m$ 维提升

### 4.1.1 投影维度、重编号与 partner 集

固定扫描维度为第 1 维，其余维度数
$$
m=d-1.
$$
对任意盒 $r$，定义其在维度 $2..d$ 上的投影为
$$
r^\downarrow=\prod_{i=2}^{d}[L_i(r),R_i(r)).
$$
为统一记号，将投影维度重编号为 $j=1,\dots,m$（对应原维度 $i=j+1$）：
$$
L'_j(r):=L_{j+1}(r),\qquad R'_j(r):=R_{j+1}(r).
$$

在处理某个 START 事件 $e=\mathrm{START}(q)$ 时，对侧活跃集记为 $S$（其来源由 $q$ 属于哪一侧决定），partner 集为
$$
K_e=\{\,s\in S\mid q^\downarrow\cap s^\downarrow\neq\varnothing\,\}.
$$

---

### 4.1.2 投影盒 $\rightarrow$ $2m$ 维点

对对侧活跃集中的每个投影盒 $s^\downarrow$，构造 $D$ 维点
$$
p(s)=\big(L'_1(s),\dots,L'_m(s),\ R'_1(s),\dots,R'_m(s)\big)\in\mathbb{R}^{2m},
$$
并令
$$
D:=2m=2(d-1).
$$
点 $p(s)$ 携带对象的唯一标识 $\mathrm{id}(s)$，用于在查询后还原 partner 身份。

---

### 4.1.3 事件盒 $\rightarrow$ $2m$ 维正交范围

对事件盒 $q$（投影为 $q^\downarrow$），定义 $D$ 维正交范围
$$
Q(q)=
\Big(\prod_{j=1}^{m}(-\infty,\ R'_j(q))\Big)\times
\Big(\prod_{j=1}^{m}(L'_j(q),\ +\infty)\Big).
$$
注意：这里的边界均为**严格不等号**，与半开语义一致（贴边不算相交）。

---

### 4.1.4 等价性：投影相交 $\Longleftrightarrow$ 点落在范围内

**引理 4.1（$2m$ 维提升等价性）**  
对任意投影盒 $q^\downarrow$ 与 $s^\downarrow$，有
$$
q^\downarrow\cap s^\downarrow\neq\varnothing
\iff
p(s)\in Q(q).
$$

**证明**  
半开语义下，投影相交当且仅当对所有 $j\in[m]$：
$$
L'_j(q)<R'_j(s)\quad\wedge\quad L'_j(s)<R'_j(q).
$$
其中 $L'_j(s)<R'_j(q)$ 等价于点 $p(s)$ 的第 $j$ 个坐标落入 $(-\infty,R'_j(q))$；  
而 $L'_j(q)<R'_j(s)$ 等价于点 $p(s)$ 的第 $m+j$ 个坐标落入 $(L'_j(q),+\infty)$。  
逐维合并即得 $p(s)\in Q(q)$；反向同理。证毕。

由此，在 START 事件 $e$ 时刻，若对侧活跃点集为
$$
P_S:=\{p(s)\mid s\in S\}\subset\mathbb{R}^{D},
$$
则
$$
K_e=\{s\in S\mid p(s)\in Q(q(e))\},
\qquad
w_e=|K_e|=|P_S\cap Q(q(e))|.
$$

---

### 4.1.5 严格边界与重复坐标的实现约定

由于 $Q(q)$ 使用严格不等号，实现中必须确保“贴边点”不会被误计入。为此，本章统一采用以下**总序键**来消除坐标重复导致的歧义。

对任意点 $p$，其唯一 id 为 $\mathrm{id}(p)\in\{1,2,\dots\}$。在任何一维排序/比较中，不直接用实数坐标 $x$，而是用键
$$
\kappa(x,p):=(x,\mathrm{id}(p))
$$
并按字典序比较。这样即便 $x$ 相同，不同点也可严格区分。

当需要表达“严格阈值”时，引入两个概念哨兵：
- $\mathrm{id}_{\min}:=0$（小于任何真实 id）；
- $\mathrm{id}_{\max}:=+\infty$（大于任何真实 id）。

于是，对任意一维严格约束，都可转为键比较：

- 约束 $x<b$ 等价于 $\kappa(x,p)<(b,\mathrm{id}_{\min})$；
- 约束 $x>a$ 等价于 $\kappa(x,p)>(a,\mathrm{id}_{\max})$；
- 开区间 $a<x<b$ 等价于
  $$
  (a,\mathrm{id}_{\max})<\kappa(x,p)<(b,\mathrm{id}_{\min}).
  $$

后续所有范围查询均按上述严格语义执行。

---

## 4.2 $D$ 维正交范围树：结构定义与三原语实现

本节构造一个维度为 $D=2m$ 的数据结构，在扫描过程中动态维护对侧活跃点集 $P_S$，并在任意正交范围 $Q$ 上支持：

- `INSERT(p)` / `DELETE(p)`：活跃集变化；
- `COUNT(Q)`：返回 $|P_S\cap Q|$；
- `REPORT(Q)`：枚举 $P_S\cap Q$ 中的点（进而得到 partner id）；
- `SAMPLE(Q,k)`：从 $P_S\cap Q$ 上输出 $k$ 个 i.i.d. 均匀（有放回）样本。

核心思想分两层：
1. 用经典 range tree 把 $D$ 维正交范围查询递归降维；
2. 在最后一维用一个“可分裂/可合并的顺序统计树”显式实现区间计数、枚举与均匀采样，从而避免把采样/枚举托付给黑盒假设。

---

### 4.2.1 一维引擎：带顺序统计的 Treap（显式支持 COUNT/REPORT/SAMPLE）

先定义一个一维动态集合结构，维护当前活跃点的键集合 $\kappa$（本章约定为 $(x,\mathrm{id})$），支持如下操作：

- `insert(key)`：插入一个键（对应一个点）；
- `erase(key)`：删除一个键；
- `count((a,b))`：返回开区间 $(a,b)$ 内键的数量；
- `report((a,b))`：枚举开区间 $(a,b)$ 内所有键；
- `sample((a,b))`：在开区间 $(a,b)$ 内的键上均匀采样一个键。

为完全显式，本章采用随机化 Treap（笛卡尔树）作为顺序统计树。Treap 的每个节点包含：
- `key`：键 $(x,\mathrm{id})$；
- `prio`：独立随机优先级；
- 左右子树指针；
- `sz`：子树大小（用于顺序统计）。

Treap 的关键原语是 `split` 与 `merge`：

- `split(T, key0)`：把 Treap $T$ 分裂成 $(A,B)$，满足
  $$
  \text{keys}(A) < key0,\qquad \text{keys}(B)\ge key0,
  $$
  且保持两棵树各自仍是合法 Treap；
- `merge(A,B)`：在满足 $\max(\text{keys}(A)) < \min(\text{keys}(B))$ 时，将两棵 Treap 合并为一棵合法 Treap。

用 `split/merge` 可实现插入删除与区间裁切。下面给出区间裁切模板：对任意开区间 $(a,b)$（允许 $a=-\infty$ 或 $b=+\infty$），定义两条分裂边界键：
$$
key_L:=(a,\mathrm{id}_{\max}),\qquad key_R:=(b,\mathrm{id}_{\min}).
$$

- 若 $a=-\infty$，则不需要按 $key_L$ 分裂；
- 若 $b=+\infty$，则不需要按 $key_R$ 分裂。

对一般有限 $(a,b)$，将 $T$ 分裂为三段：
1) $(T_{\le a}, T_{>a}) := \mathrm{split}(T, key_L)$  
2) $(T_{(a,b)}, T_{\ge b}) := \mathrm{split}(T_{>a}, key_R)$

则中段 $T_{(a,b)}$ 的键集合恰为开区间内的所有活跃点键。于是：

- 计数：$\mathrm{count}((a,b)) = \mathrm{sz}(T_{(a,b)})$；
- 枚举：对 $T_{(a,b)}$ 做中序遍历即可输出全部键（输出规模为 $K$ 时为 $O(K)$）；
- 均匀采样：若 $K=\mathrm{sz}(T_{(a,b)})>0$，取 $r\sim\mathrm{Unif}\{1,\dots,K\}$，用顺序统计 `kth(T_{(a,b)}, r)` 得到第 $r$ 小键；它在该集合上均匀。

操作完成后用
$$
T := \mathrm{merge}\big(T_{\le a}, \mathrm{merge}(T_{(a,b)}, T_{\ge b})\big)
$$
恢复原 Treap 结构（集合不变）。

> 重要：以上过程只使用局部结构变换与子树大小维护，不依赖任何外部“区间采样黑盒”。

---

### 4.2.2 $k$ 维范围树的递归结构（$k=1$ 时落到 Treap）

下面用 $T^{(k)}$ 表示一个 $k$ 维范围树，维护一个动态活跃点集 $P\subset\mathbb{R}^k$（点携带 id）。当 $k=D$ 时即为本章需要的结构。

**结构定义**

- **基例 $k=1$**：  
  $T^{(1)}$ 仅包含一个 Treap，存储当前活跃点的键 $\kappa(x_1,p)=(x_1(p),\mathrm{id}(p))$。  
  它直接支持对任意开区间/半无限区间的 `COUNT/REPORT/SAMPLE`（见 4.2.1）。

- **递归步 $k>1$**：  
  $T^{(k)}$ 由一棵主树 $T_x$ 和若干关联结构（associated structures）组成。
  
  1) 主树 $T_x$：按第 1 坐标 $x_1$ 的键 $\kappa(x_1,p)$ 对点的**全集**建立一棵平衡二叉搜索树（可用“按排序数组取中位数递归建树”的静态完美平衡树），每个节点 $v$ 对应一个子树点集 $P(v)$（全集固定，活跃状态动态变化）。

  2) 关联结构：对主树每个节点 $v$，建立一个 $(k-1)$ 维范围树
     $$
     \mathrm{assoc}(v) := T^{(k-1)}\big(\pi(P(v))\big),
     $$
     其中投影 $\pi$ 丢弃第 1 坐标，保留 $(x_2,\dots,x_k)$，并保留同一个点 id 作为身份信息。

  直观上，$T_x$ 用于处理第一维的“落在某个区间内”，而每个 $\mathrm{assoc}(v)$ 用于在固定 $P(v)$ 的前提下处理其余 $k-1$ 维的正交范围约束。

---

### 4.2.3 动态更新：INSERT/DELETE 的显式递归过程

范围树维护的是“活跃点集”。在扫描过程中，点只会随盒进入/离开活跃集而被插入/删除一次，因此可直接使用 `INSERT/DELETE` 语义。

对 $T^{(k)}$ 插入/删除点 $p=(x_1,\dots,x_k)$ 的递归过程如下。

- **基例 $k=1$**：  
  `INSERT(p)` 等价于在 Treap 中插入键 $(x_1,\mathrm{id}(p))$；  
  `DELETE(p)` 等价于删除该键。

- **递归步 $k>1$**：  
  在主树 $T_x$ 中按键 $(x_1,\mathrm{id}(p))$ 搜索到对应叶子（静态树，高度为 $O(\log N)$），沿根到叶路径上的每个节点 $v$ 同步更新其关联结构：
  $$
  \mathrm{assoc}(v).\mathrm{INSERT}(\pi(p))
  \quad\text{或}\quad
  \mathrm{assoc}(v).\mathrm{DELETE}(\pi(p)).
  $$

由于每次更新访问主树路径上的 $O(\log N)$ 个节点，并在每个节点调用一次 $(k-1)$ 维更新，因此更新时间满足递归
$$
T_{\mathrm{upd}}(k)=O(\log N)\cdot T_{\mathrm{upd}}(k-1),\qquad T_{\mathrm{upd}}(1)=O(\log N),
$$
从而
$$
T_{\mathrm{upd}}(k)=O(\log^{k}N),
$$
其中 $N$ 为全集规模上界（活跃点数不超过全集）。

> 若一维 Treap 使用随机优先级，则上述界为**期望**界；若需要确定性界，可将一维结构替换为确定性平衡 BST + 顺序统计（实现方式类似，只是旋转规则不同）。

---

### 4.2.4 查询分解：$k$ 维正交范围 $\rightarrow$ 不交终端块

令
$$
Q=\prod_{i=1}^k I_i
$$
为 $k$ 维正交范围，各 $I_i$ 可为开区间或半无限开区间（严格语义见 4.1.5）。

为了实现 `COUNT/REPORT/SAMPLE`，本章把 $P\cap Q$ 表达成若干**不相交终端块**的并：每个终端块都落在某个一维 Treap 的“一个区间查询”上，从而能显式完成计数、枚举与均匀采样。

**定义（终端块）**  
对任意 $k$ 维范围树 $T^{(k)}$ 与范围 $Q$，定义过程 `BLOCKS_k(T^{(k)},Q)`，返回一个块列表 $\mathcal{B}(Q)$。每个块 $b\in\mathcal{B}(Q)$ 形如
$$
b=(\mathcal{T}, I),
$$
其中 $\mathcal{T}$ 是某个 $T^{(1)}$ 的 Treap 句柄，$I$ 是一维开区间/半无限开区间；该块对应的点集定义为
$$
P(b):=\{\,p\in P \mid p\ \text{处于活跃状态且其一维键落在 } I\,\}.
$$

**构造（递归）**

- **基例 $k=1$**：  
  $Q=I_1$，返回单块
  $$
  \mathrm{BLOCKS}_1(T^{(1)},I_1)=\{(T^{(1)}.\text{Treap},\ I_1)\}.
  $$

- **递归步 $k>1$**：  
  在主树 $T_x$ 上对第一维区间 $I_1$ 做规范分解，得到一组节点
  $$
  \mathcal{C}(I_1)=\{v_1,\dots,v_t\},\qquad t=O(\log N),
  $$
  其性质是：这些节点的子树点集两两不交，且恰好覆盖所有第一维落在 $I_1$ 内的点。  
  然后对每个 $v_j$ 在其关联结构上递归取块：
  $$
  \mathrm{BLOCKS}_{k-1}\big(\mathrm{assoc}(v_j),\ Q'\big),\qquad Q'=\prod_{i=2}^k I_i.
  $$
  将所有递归结果并在一起，即
  $$
  \mathrm{BLOCKS}_{k}(T^{(k)},Q)
  :=
  \bigcup_{j=1}^{t}\mathrm{BLOCKS}_{k-1}(\mathrm{assoc}(v_j),Q').
  $$

> 关键点：本章在 $k=1$ 时**不再**把一维区间继续拆成 $O(\log N)$ 个结点块，而是直接把“一维区间查询”作为终端块。这能将终端块数量稳定在 $O(\log^{k-1}N)$，并将区间内的计数/枚举/采样全部交给 4.2.1 的 Treap 显式完成。

---

### 4.2.5 分解正确性：覆盖与不交

**定理 4.2（终端块分解的覆盖与不交）**  
对任意 $k\ge 1$、活跃点集 $P$ 与正交范围 $Q=\prod_{i=1}^k I_i$，令
$$
\mathcal{B}(Q):=\mathrm{BLOCKS}_k(T^{(k)},Q).
$$
则：

1)（覆盖性）任意点 $p\in P\cap Q$ 必属于某个块 $b\in\mathcal{B}(Q)$ 的点集 $P(b)$；  
2)（不交性）任意两个不同块 $b\neq b'$，有 $P(b)\cap P(b')=\varnothing$。

因此
$$
P\cap Q=\biguplus_{b\in\mathcal{B}(Q)} P(b),
$$
其中 $\biguplus$ 为不交并。

**证明（归纳）**  
对 $k$ 归纳。

- $k=1$：只有一个块 $(\mathcal{T},I_1)$，其点集正是 $P\cap I_1$，覆盖与不交显然成立。

- $k>1$：考虑主树第一维区间分解 $\mathcal{C}(I_1)=\{v_1,\dots,v_t\}$。按规范分解性质，$\{P(v_j)\}$ 两两不交，且
  $$
  \{p\in P\mid x_1(p)\in I_1\}=\biguplus_{j=1}^t P(v_j).
  $$
  取任意 $p\in P\cap Q$，它必在某个唯一的 $P(v_j)$ 中；又 $p$ 的其余坐标满足 $Q'$，因此其投影 $\pi(p)$ 属于 $\mathrm{assoc}(v_j)$ 的查询结果。由归纳假设，$\pi(p)$ 必属于 $\mathrm{BLOCKS}_{k-1}(\mathrm{assoc}(v_j),Q')$ 中某个块的点集，从而 $p$ 属于 $\mathcal{B}(Q)$ 中某个块，覆盖性得证。

  不交性方面：来自不同 $v_j$ 的块分别只含 $P(v_j)$ 内的点，而 $P(v_j)$ 两两不交；来自同一 $v_j$ 的块的不交性由归纳假设保证。故整体不交。证毕。

---

### 4.2.6 在终端块之上实现 COUNT / REPORT / SAMPLE

记查询范围为 $Q$，终端块集合为 $\mathcal{B}(Q)=\{b_1,\dots,b_B\}$，其中
$$
b_i=(\mathcal{T}_i, I_i).
$$
定义每个块的权重（即该块内活跃点数）：
$$
w_i := |P(b_i)| = \mathcal{T}_i.\mathrm{count}(I_i),
\qquad
W:=\sum_{i=1}^B w_i = |P\cap Q|.
$$

#### `COUNT(Q)`
返回
$$
\mathrm{COUNT}(Q)=\sum_{i=1}^B w_i.
$$

#### `REPORT(Q)`
对 $i=1..B$，调用一维引擎执行 `report`：
- 输出 $\mathcal{T}_i.\mathrm{report}(I_i)$ 枚举到的每个点 id；
- 因为定理 4.2 不交，跨块不会重复。

#### `SAMPLE(Q,k)`
若 $W=0$ 返回空列表。否则要在 $P\cap Q$ 上有放回生成 $k$ 个 i.i.d. 均匀样本。

本章采用“按块权重混合 + 块内均匀”的显式两阶段抽样：

1) **按权重选择块**  
   构造前缀和数组
   $$
   S_0:=0,\qquad S_i:=\sum_{j=1}^i w_j\ (1\le i\le B),
   $$
   则 $S_B=W$。取随机整数
   $$
   R\sim\mathrm{Unif}\{1,2,\dots,W\},
   $$
   用二分查找最小的 $i$ 使 $S_i\ge R$，选中块 $b_i$。  
   这一步对每个样本耗时 $O(\log B)$，而 $B=O(\log^{D-1}N)$，故 $\log B=O(\log\log N)$ 为很小的额外因子。

2) **块内均匀采样**  
   在选中块 $b_i=(\mathcal{T}_i,I_i)$ 内调用
   $$
   p\leftarrow \mathcal{T}_i.\mathrm{sample}(I_i),
   $$
   由 4.2.1，`sample(I_i)` 在该区间内的活跃点上均匀。

重复上述两步 $k$ 次，即得到 $k$ 个样本。

---

### 4.2.7 采样正确性：均匀与独立（有放回）

**定理 4.3（`SAMPLE(Q,k)` 的均匀性与独立性）**  
设 $W=|P\cap Q|>0$。上述 `SAMPLE(Q,k)` 输出的每个样本在 $P\cap Q$ 上均匀；且 $k$ 次输出相互独立（有放回）。

**证明**  
由定理 4.2，$P\cap Q$ 可写为不交并
$$
P\cap Q=\biguplus_{i=1}^B P(b_i),\qquad |P(b_i)|=w_i,\quad \sum_i w_i=W.
$$
任取一个具体点 $p\in P\cap Q$，它属于唯一块 $b_{i^\star}$。一次采样命中 $p$ 的概率为
$$
\Pr[\text{命中 }p]
=
\Pr[\text{选中块 }i^\star]\cdot \Pr[p\mid \text{选中块 }i^\star]
=
\frac{w_{i^\star}}{W}\cdot \frac{1}{w_{i^\star}}
=
\frac{1}{W},
$$
与 $p$ 无关，故均匀。

独立性：每次采样均重新独立地生成随机数 $R$（选块）以及块内 `sample` 的随机选择；并且采样为有放回（不会改变集合语义，只是读查询）。因此各次输出相互独立。证毕。

---

### 4.2.8 复杂度与空间

设对侧全集规模上界为 $N$，维度为 $D=2(d-1)$。以下给出单次操作的渐近代价；其中一维 Treap 带来的对数因子为期望界。

- **动态更新**（`INSERT/DELETE`）：
  $$
  T_{\mathrm{upd}}=O(\log^{D}N).
  $$

- **计数**（`COUNT(Q)`）：
  终端块数 $B=O(\log^{D-1}N)$；每块一维计数 $O(\log N)$，总计
  $$
  T_{\mathrm{count}}=O(\log^{D}N).
  $$

- **枚举**（`REPORT(Q)`）：
  块构造与边界开销 $O(\log^{D}N)$，输出规模 $K=|P\cap Q|$，故
  $$
  T_{\mathrm{report}}=O(\log^{D}N+K).
  $$

- **采样**（`SAMPLE(Q,k)`）：
  先构造块并计算权重：$O(\log^{D}N)$；每个样本包含“二分选块” $O(\log B)=O(\log\log N)$ 与块内采样 $O(\log N)$，因此
  $$
  T_{\mathrm{sample}}
  =
  O(\log^{D}N + k(\log N+\log\log N)).
  $$
  若将“按权重选块”替换为别名表（alias table）等 $O(1)$ 选块方案，则可写成
  $$
  T_{\mathrm{sample}}=O(\log^{D}N + k\log N).
  $$

- **空间**：范围树的标准空间界为
  $$
  S=O(N\log^{D-1}N),
  $$
  且常数因子随维度 $D$ 增大而快速上升。对事件块问题而言 $D=2(d-1)$，因此该对照方法在维度较大时开销会非常显著。

---

## 4.3 与扫描事件块原语的对接

本节说明如何在第 1 章的扫描过程中使用本章范围树实现事件块 $K_e$ 的三原语，并直接供第 3 章三种采样框架调用。

### 4.3.1 两侧索引结构

分别为两类集合维护对侧活跃集的索引：

- 维护 $T_c$：存储当前活跃的 $R_c$ 盒对应的点 $p(r_c)$；
- 维护 $T_{\bar c}$：存储当前活跃的 $R_{\bar c}$ 盒对应的点 $p(r_{\bar c})$。

两棵树的点维度均为 $D=2(d-1)$，点构造与查询范围使用 4.1 的定义。

---

### 4.3.2 扫描事件处理

按第 1 章定义的事件总序 $\prec$ 处理事件：

- 若事件为 $\mathrm{END}(r)$：  
  从其所属集合对应的树中执行 `DELETE(p(r))`。

- 若事件为 $\mathrm{START}(q)$：  
  先对对侧树执行查询：
  - 若 $q\in R_c$，对 $T_{\bar c}$ 查询 $Q(q)$；
  - 若 $q\in R_{\bar c}$，对 $T_c$ 查询 $Q(q)$。

  于是事件块原语直接实现为：
  $$
  \texttt{COUNT}(e)=\texttt{COUNT}(Q(q(e))),
  $$
  $$
  \texttt{REPORT}(e)=\texttt{REPORT}(Q(q(e))),
  $$
  $$
  \texttt{SAMPLE}(e,k)=\texttt{SAMPLE}(Q(q(e)),k).
  $$

  对 `REPORT/SAMPLE` 返回的每个 partner id，按第 1 章/第 3 章定义的映射 $\Phi_e$ 还原为有序对并输出。

  最后将 $q$ 插入其所属集合对应的树中：`INSERT(p(q))`。

---

### 4.3.3 与三种采样框架的组合

第 3 章的三种采样框架只依赖事件块原语的语义正确性以及各次随机调用之间的独立性约定。本章已经在 4.2.7 中给出 `SAMPLE` 的均匀性与独立性证明，并且 `COUNT/REPORT` 均严格返回 $P_S\cap Q(q)$ 上的计数与无重枚举。因此，将本章原语实现替换进三种框架后，整体输出仍严格为 $J$ 上 i.i.d. 均匀（有放回）。

从代价角度，令 $M=|E^+|$ 为 START 事件数、$N$ 为对侧全集规模上界、$D=2(d-1)$，则单次事件：
- 更新：$O(\log^{D}N)$；
- 计数：$O(\log^{D}N)$；
- 枚举：$O(\log^{D}N+w_e)$；
- 采样（$k$ 个）：$O(\log^{D}N + k\log N)$（或带 $k\log\log N$ 的显式前缀和选块实现）。

这些代价可直接代入第 3 章三种框架的“原语调用计数”，得到总复杂度表达式。

# 第 5 章 对照方法二：kd-tree 上的 $2d$ 维范围计数与范围采样 Join

本章给出一套独立的 Join 采样对照体系：将“跨集合的 $d$ 维半开盒相交”提升为 $2d$ 维点集上的正交范围问题，在右侧集合上建立静态 kd-tree，并通过“按度加权选择左端 + 条件范围内均匀采样选择右端”实现对 Join 集合的 **i.i.d. 均匀（有放回）**采样。该方法结构直接、证明链条清晰，但由于查询维度为 $D=2d$，会呈现显著的维度效应，适合作为对照基线。

---

## 5.1 问题定义与 $2d$ 维提升

### 5.1.1 输入、半开语义与 Join

给定两类 $d$ 维轴对齐半开盒集合
$$
R_c=\{r_1,\dots,r_n\},\qquad R_{\bar c}=\{s_1,\dots,s_m\}.
$$
每个盒写作
$$
r=\prod_{i=1}^{d}[L_i(r),R_i(r)),\qquad L_i(r)<R_i(r),
$$
并假定每个对象携带唯一标识 $\mathrm{id}(\cdot)$（即便几何参数完全相同，也视为不同元素）。

半开语义下，两盒相交当且仅当对所有维度 $i$ 同时满足严格不等式
$$
r\cap s\neq\varnothing
\iff
\big(L_i(r)<R_i(s)\big)\ \wedge\ \big(L_i(s)<R_i(r)\big).
$$
严格不等号保证贴边（例如 $R_i(r)=L_i(s)$）不计为相交。

本章研究方向固定的 Join 集合
$$
J=\{(r,s)\mid r\in R_c,\ s\in R_{\bar c},\ r\cap s\neq\varnothing\},
$$
目标是在 $J$ 上输出 $t$ 个样本 $Z_1,\dots,Z_t$，满足
$$
\Pr\{Z_j=P\}=\frac{1}{|J|}\ (\forall P\in J),\qquad Z_1,\dots,Z_t\ \text{相互独立},
$$
即 **i.i.d. 且均匀（有放回）**。若 $|J|=0$，输出空序列。

---

### 5.1.2 右侧盒 $\rightarrow$ $2d$ 维点

对每个右侧盒 $s\in R_{\bar c}$，构造点
$$
p(s)=\big(L_1(s),\dots,L_d(s),\ R_1(s),\dots,R_d(s)\big)\in\mathbb{R}^{2d},
$$
得到点集
$$
P=\{p(s)\mid s\in R_{\bar c}\}.
$$
点 $p(s)$ 携带 $\mathrm{id}(s)$，用于在查询后恢复对应的右侧盒。

记提升维度
$$
D:=2d.
$$

---

### 5.1.3 左侧盒 $\rightarrow$ $2d$ 维正交范围

对每个左侧盒 $r\in R_c$，定义正交范围
$$
Q(r)=
\Big(\prod_{i=1}^{d}(-\infty,\ R_i(r))\Big)\times
\Big(\prod_{i=1}^{d}(L_i(r),\ +\infty)\Big).
$$
所有边界均为开边界（严格不等号），与半开盒相交判定完全一致。

---

### 5.1.4 等价性与 Join 度数

**引理 5.1（$2d$ 维提升等价性）**  
对任意 $r\in R_c$ 与 $s\in R_{\bar c}$，有
$$
p(s)\in Q(r)\iff r\cap s\neq\varnothing.
$$

**证明**  
$p(s)\in Q(r)$ 等价于两组逐维约束：
- 前 $d$ 个坐标要求 $L_i(s)<R_i(r)$；
- 后 $d$ 个坐标要求 $R_i(s)>L_i(r)$，等价于 $L_i(r)<R_i(s)$。

合并即得到半开语义下相交的充要条件。证毕。

据此定义左侧盒的 Join 度数（也称连接度）
$$
w(r)=|\{s\in R_{\bar c}\mid r\cap s\neq\varnothing\}|=|P\cap Q(r)|.
$$
并有
$$
|J|=\sum_{r\in R_c}w(r).
$$

---

## 5.2 静态 kd-tree：结构与构建

本节在静态点集 $P\subset\mathbb{R}^D$ 上建立平衡 kd-tree，用于支持对任意正交范围 $Q$ 的精确计数 `COUNT(Q)` 与范围内均匀采样 `SAMPLE(Q,k)`。

### 5.2.1 节点信息

每个节点 $v$ 存储：

- `left(v), right(v)`：左右子节点指针；
- `split_dim(v)\in\{0,\dots,D-1\}` 与 `split_val(v)`：分割维度与阈值（用于构建与剪枝）；
- `size(v)`：子树点数；
- 子树外包框（闭盒）
  $$
  \mathrm{MBR}(v)=\prod_{j=0}^{D-1}[\min_j(v),\max_j(v)],
  $$
  其中 $\min_j(v),\max_j(v)$ 分别是子树点在第 $j$ 维的最小/最大坐标。

此外，为支持“当 $\mathrm{MBR}(v)$ 被查询范围完全包含时，以 $O(1)$ 块内均匀抽样一个子树点”，采用如下静态构建约定：

- 构建过程对点数组做 in-place partition；
- 每个节点子树对应数组的一个连续切片 $[l_v,r_v)$；
- 因而 `size(v)=r_v-l_v`，并可在该切片上直接均匀抽取下标实现块内均匀采样。

该约定仅依赖 kd-tree 为静态索引（构建后不再重排点数组），不要求额外存储子树点列表。

---

### 5.2.2 构建方式与代价

设 $m=|P|$。递归构建平衡 kd-tree：

1. 在深度 $h$ 的节点选择 `split_dim = h mod D`；
2. 在该维以中位数划分，使左右子树规模尽量均衡；
3. 递归直到子树点数 $\le B$（叶子阈值，例如 32 或 64）；
4. 自底向上计算每个节点的 `MBR` 与 `size`。

构建代价取决于“取中位数并 partition”的实现：

- 若每层使用线性期望的选择算法（例如 `nth_element` 风格）完成中位数划分，则构建期望时间为
  $$
  T_{\mathrm{build}}=\mathbb{E}[O(m\log m)].
  $$
- 若每层对当前子数组排序后取中位数，则构建时间为
  $$
  O(m\log^2 m).
  $$

空间方面，kd-tree 节点数为 $O(m)$，每个节点存长度为 $D$ 的 $\min/\max$ 向量与常数指针字段，因此
$$
S_{\mathrm{tree}}=O(mD).
$$

---

## 5.3 开边界范围与闭盒 MBR 的关系判定

本章查询范围均来自 $Q(r)$，其每一维区间都是开边界的半无限区间：$(-\infty,b)$ 或 $(a,+\infty)$。为便于统一描述，本节允许更一般的开区间 $(a,b)$。

令查询范围
$$
Q=\prod_{j=0}^{D-1} I_j,
$$
其中 $I_j$ 为以下三类之一：

- $I_j=(-\infty,b)$（严格 $x<b$）；
- $I_j=(a,+\infty)$（严格 $x>a$）；
- $I_j=(a,b)$（严格 $a<x<b$）。

对节点 $v$ 的闭盒外包框
$$
\mathrm{MBR}(v)=\prod_{j=0}^{D-1}[\min_j(v),\max_j(v)].
$$

定义三类关系：

- `disjoint(v,Q)`：$\mathrm{MBR}(v)\cap Q=\varnothing$；
- `contained(v,Q)`：$\mathrm{MBR}(v)\subset Q$（注意 $Q$ 是开边界）；
- 其余为 `partial(v,Q)`。

逐维判定规则（严格语义）：

1. 若 $I_j=(-\infty,b)$：
   - 不相交：$\min_j(v)\ge b$；
   - 完全包含：$\max_j(v)<b$。
2. 若 $I_j=(a,+\infty)$：
   - 不相交：$\max_j(v)\le a$；
   - 完全包含：$\min_j(v)>a$。
3. 若 $I_j=(a,b)$：
   - 不相交：$\min_j(v)\ge b$ 或 $\max_j(v)\le a$；
   - 完全包含：$\min_j(v)>a$ 且 $\max_j(v)<b$。

这些不等式确保：当 $\mathrm{MBR}(v)$ 边界恰好等于查询开边界时，不会被误判为完全包含，从而不引入“贴边点”的错误计入。

---

## 5.4 `COUNT(Q)`：精确范围计数

### 5.4.1 递归计数算法

定义递归函数 `COUNT(v,Q)`：

- 若 `disjoint(v,Q)`：返回 $0$；
- 若 `contained(v,Q)`：返回 `size(v)`；
- 若 $v$ 为叶子：逐点检查叶内点 $p$ 是否满足 $p\in Q$（逐维严格比较），返回命中点数；
- 否则返回
  $$
  \texttt{COUNT}(\texttt{left}(v),Q)+\texttt{COUNT}(\texttt{right}(v),Q).
  $$

全树计数为 `COUNT(root,Q)`。

---

### 5.4.2 正确性

**定理 5.2（`COUNT` 的精确性）**  
对任意开边界正交范围 $Q$，`COUNT(root,Q)` 的返回值等于 $|P\cap Q|$。

**证明**  
对任意节点 $v$：

- 若不相交，则 $P(v)\cap Q=\varnothing$，返回 0 正确；
- 若完全包含，则 $P(v)\subset Q$，返回 `size(v)=|P(v)|` 正确；
- 否则递归分解到子节点，或在叶子逐点严格判定并计数。

三分情况覆盖且互斥，递归累加恰好统计所有命中点数，因此总计数为 $|P\cap Q|$。证毕。

---

### 5.4.3 访问代价的常见刻画

令 $m=|P|$，维度为 $D$。在平衡 kd-tree、$D$ 为常数且数据/查询不对抗的常见平均模型下，范围计数访问的节点数经常用
$$
T_{\mathrm{count}}(m,D)=\tilde O\!\left(m^{1-\frac{1}{D}}\right)
$$
来刻画，其中 $\tilde O(\cdot)$ 隐含多对数因子与叶阈值常数。  
在最坏情况下（例如高维或对抗分布/查询），访问代价仍可能退化到 $O(m)$。

---

## 5.5 `SAMPLE(Q,k)`：范围内 i.i.d. 均匀采样

本节给出在 $P\cap Q$ 上生成 $k$ 个 **i.i.d. 且均匀（有放回）**点样本的过程，并给出完整正确性论证。

### 5.5.1 不交块分解：contained 块与边界叶块

对 kd-tree 做一次遍历，按 5.3 的三类关系处理节点：

- 若某节点 $v$ 满足 `contained(v,Q)`：将 $v$ 作为一个 **contained 块**加入集合 $\mathcal{I}(Q)$，其块权重为
  $$
  w_v:=\texttt{size}(v).
  $$
  该块代表子树点集 $P(v)$ 整体属于 $P\cap Q$。

- 若节点与 $Q$ 不相交：剪枝跳过。

- 若节点与 $Q$ 部分相交：
  - 若非叶子，继续下探；
  - 若为叶子 $\ell$，逐点检查叶内点是否属于 $Q$，将命中点的数组位置收集为索引表 $\mathrm{hit}(\ell)$，并记
    $$
    w_\ell := |\mathrm{hit}(\ell)|.
    $$
    当 $w_\ell>0$ 时，将该叶子作为 **边界叶块**加入集合 $\mathcal{B}(Q)$。

定义总命中数
$$
C(Q)=\sum_{v\in\mathcal{I}(Q)} w_v+\sum_{\ell\in\mathcal{B}(Q)} w_\ell.
$$
若 $C(Q)=0$，返回空序列。

---

### 5.5.2 分解正确性

**引理 5.3（块分解覆盖且无重）**  
由上述遍历得到的块集合满足不交并分解
$$
P\cap Q
=
\biguplus_{v\in\mathcal{I}(Q)} P(v)
\ \uplus\
\biguplus_{\ell\in\mathcal{B}(Q)} \big(P(\ell)\cap Q\big),
$$
其中 $P(v)$ 为节点 $v$ 子树点集，$P(\ell)$ 为叶子 $\ell$ 内点集。

**证明**  

- 不交性：遍历遇到 `contained` 节点即停止下探，因此 $\mathcal{I}(Q)$ 中不会出现祖先与后代同时入块，不同 contained 块对应的子树点集两两不交。叶子之间点集天然不交。边界叶块来自“未被任何 contained 祖先截断”的路径末端，因此边界叶不可能位于任一 contained 块的子树内，从而两类块之间也不重叠。

- 覆盖性：任取 $p\in P\cap Q$。沿根到叶的路径下行：若途中出现第一个满足 `contained` 的节点 $v$，则 $p\in P(v)$ 并被该块覆盖；否则路径最终到达某个部分相交叶子 $\ell$，且逐点筛选会把 $p$ 纳入 $P(\ell)\cap Q$。因此所有命中点均被覆盖。

综上得到不交并分解。证毕。

由引理 5.3 立刻得到 $C(Q)=|P\cap Q|$。

---

### 5.5.3 两阶段采样：按块权重 + 块内均匀

对块集合 $\mathcal{I}(Q)\cup\mathcal{B}(Q)$ 建立按权重的离散分布，使一次选块满足
$$
\Pr\{\text{选中块 }b\}=\frac{w_b}{C(Q)}.
$$
块内均匀抽样：

- 若选中 contained 块 $v\in\mathcal{I}(Q)$：在其数组切片 $[l_v,r_v)$ 上均匀抽取
  $$
  u\sim \mathrm{Unif}\{l_v,\dots,r_v-1\},
  $$
  返回点数组第 $u$ 个点的 $\mathrm{id}$。

- 若选中边界叶块 $\ell\in\mathcal{B}(Q)$：在命中索引表 $\mathrm{hit}(\ell)$ 上均匀抽取一个位置并返回相应点的 $\mathrm{id}$。

重复 $k$ 次，且每次使用独立随机性，即得到 $k$ 个样本。

---

### 5.5.4 正确性：均匀与 i.i.d.

**定理 5.4（`SAMPLE(Q,k)` 的正确性）**  
当 $C(Q)>0$ 时，上述 `SAMPLE(Q,k)` 生成的 $k$ 个点样本在 $P\cap Q$ 上 **均匀且相互独立（有放回）**。

**证明**  
由引理 5.3，$P\cap Q$ 为块集合的不交并。任取 $p\in P\cap Q$，它属于唯一块 $b(p)$，且该块大小为 $w_{b(p)}$。一次采样命中 $p$ 的概率为
$$
\Pr[p]
=
\Pr[b(p)]\cdot \Pr[p\mid b(p)]
=
\frac{w_{b(p)}}{C(Q)}\cdot \frac{1}{w_{b(p)}}
=
\frac{1}{C(Q)},
$$
与 $p$ 无关，故单次输出在 $P\cap Q$ 上均匀。每次采样独立重做“选块 + 块内抽样”，且为有放回抽样，因此 $k$ 次输出相互独立且同分布。证毕。

---

### 5.5.5 代价：遍历 + 选块 + 块内抽样

记本次查询产生的块数为 $B_Q:=|\mathcal{I}(Q)|+|\mathcal{B}(Q)|$。一次 `SAMPLE(Q,k)` 的代价由两部分构成：

1. 遍历 kd-tree，形成块集合并计算权重（与 `COUNT(Q)` 同阶，并在边界叶做逐点筛选）；
2. 生成 $k$ 个样本：每个样本包含一次“按权重选块”与一次块内均匀抽样。

按权重选块可用两种常见实现：

- 前缀和 + 二分：预处理 $O(B_Q)$，每次选块 $O(\log B_Q)$；
- alias：预处理 $O(B_Q)$，每次选块期望 $O(1)$。

因此可写成更精确的两种表达：

- 前缀和版本：
  $$
  T_{\mathrm{sample}}=\tilde O\!\left(T_{\mathrm{visit}} + B_Q + k\log B_Q\right),
  $$
- alias 版本：
  $$
  T_{\mathrm{sample}}=\tilde O\!\left(T_{\mathrm{visit}} + B_Q + k\right),
  $$
  其中 $T_{\mathrm{visit}}$ 为遍历访问节点的代价。  
  在固定维度的常见平均模型下，$T_{\mathrm{visit}}$ 常用 $\tilde O(m^{1-1/D})$ 刻画；最坏情况下仍可能退化到 $O(m)$。

---

## 5.6 基于 kd-tree 的 Join 采样器

本节将 5.1 的提升与 5.4–5.5 的 `COUNT/SAMPLE` 组合，得到 $J$ 上的 i.i.d. 均匀（有放回）采样器。

### 5.6.1 预处理：构建 kd-tree 与计算度数

令
$$
n:=|R_c|,\qquad m:=|R_{\bar c}|,\qquad D:=2d.
$$

步骤如下：

1. 构造点集 $P=\{p(s)\mid s\in R_{\bar c}\}\subset\mathbb{R}^{D}$；
2. 在 $P$ 上构建 kd-tree 索引 $T$；
3. 对每个 $r\in R_c$，构造范围 $Q(r)$ 并计算
   $$
   w(r)\leftarrow \texttt{COUNT}(T,Q(r))=|P\cap Q(r)|.
   $$
4. 计算
   $$
   W:=\sum_{r\in R_c} w(r)=|J|.
   $$
   若 $W=0$，输出空序列并终止；
5. 在集合 $\{r\in R_c\mid w(r)>0\}$ 上建立离散采样结构，使一次抽取满足
   $$
   \Pr\{r\}=\frac{w(r)}{W}.
   $$

---

### 5.6.2 采样：生成 $t$ 个 Join 样本

对 $j=1,\dots,t$：

1. 按 $\Pr\{r\}=w(r)/W$ 抽取左端盒 $r\in R_c$；
2. 调用 `SAMPLE(T,Q(r),1)` 得到一个在 $P\cap Q(r)$ 上均匀的点 id，并据此恢复对应右侧盒 $s\in R_{\bar c}$；
3. 输出 $Z_j=(r,s)$。

为降低重复查询开销，可使用批量化生成：先独立抽取 $t$ 次左端盒得到序列 $r^{(1)},\dots,r^{(t)}$，统计每个左端 $r$ 被抽到的次数 $k_r$，然后对每个 $k_r>0$ 的 $r$ 只调用一次 `SAMPLE(T,Q(r),k_r)`，再按位置写回输出。该批量化操作不改变分布。

---

### 5.6.3 伪代码

~~~text
Algorithm 5.1  KD-Tree Join Sampler
Input : R_c, R_{\bar c}, sample size t
Output: Z_1..Z_t (i.i.d. uniform on J, with replacement)

1:  P ← { p(s) : s ∈ R_{\bar c} }, where p(s) = (L_1(s)..L_d(s), R_1(s)..R_d(s))
2:  Build a balanced kd-tree T on P   // static, with subtree array slices
3:  For each r ∈ R_c:
4:      Q(r) ← (∏_{i=1}^d (-∞, R_i(r))) × (∏_{i=1}^d (L_i(r), +∞))
5:      w(r) ← COUNT(T, Q(r))
6:  W ← Σ_{r∈R_c} w(r)
7:  If W = 0: return empty sequence
8:  Build a discrete sampler over {r ∈ R_c : w(r) > 0} with Pr[r] = w(r)/W

9:  For j = 1..t:
10:     Draw r ~ Pr[r] = w(r)/W
11:     Draw s_id ← SAMPLE(T, Q(r), 1)   // uniform in P ∩ Q(r)
12:     Let s be the box in R_{\bar c} with id = s_id
13:     Z_j ← (r, s)
14: return Z_1..Z_t
~~~

---

## 5.7 正确性：$J$ 上的均匀性与 i.i.d.

### 5.7.1 单次均匀性

**定理 5.5（Join 采样的均匀性）**  
当 $W=|J|>0$ 时，算法 5.1 的单次输出 $Z$ 在 $J$ 上均匀：
$$
\Pr\{Z=(r,s)\}=\frac{1}{|J|},\qquad \forall (r,s)\in J.
$$

**证明**  
任取 $(r,s)\in J$。由引理 5.1，$p(s)\in Q(r)$，故 $w(r)=|P\cap Q(r)|>0$。算法输出 $(r,s)$ 的概率为
$$
\Pr\{(r,s)\}
=
\Pr\{r\}\cdot \Pr\{s\mid r\}
=
\frac{w(r)}{W}\cdot \frac{1}{w(r)}
=
\frac{1}{W}
=
\frac{1}{|J|},
$$
其中 $\Pr\{s\mid r\}=1/w(r)$ 由 `SAMPLE(T,Q(r),1)` 在 $P\cap Q(r)$ 上均匀得到。证毕。

---

### 5.7.2 多次输出的独立同分布

**定理 5.6（Join 采样的 i.i.d.）**  
若每次左端抽取与每次 `SAMPLE` 调用都使用独立随机性，且 `SAMPLE` 为有放回采样，则算法 5.1 输出的 $Z_1,\dots,Z_t$ 相互独立且同分布为 $J$ 上的均匀分布。

**证明**  
由定理 5.5，每个 $Z_j$ 的边缘分布为 $J$ 上均匀。不同轮次的随机选择相互独立，且每轮采样均不改变后续轮次的候选集合（有放回）。因此 $Z_1,\dots,Z_t$ 相互独立且同分布。证毕。

---

## 5.8 复杂度与维度效应

令查询维度为 $D=2d$。

### 5.8.1 预处理

- kd-tree 构建：期望 $O(m\log m)$（使用线性期望选择算法取中位数），或 $O(m\log^2 m)$（每层排序取中位数）。
- 度数计算：对每个 $r\in R_c$ 调用一次 `COUNT(Q(r))`。在固定 $D$ 的常见平均模型下，单次计数常用
  $$
  \tilde O\!\left(m^{1-\frac{1}{D}}\right)=\tilde O\!\left(m^{1-\frac{1}{2d}}\right)
  $$
  刻画，因此预处理总计常写作
  $$
  \tilde O\!\left(n\cdot m^{1-\frac{1}{2d}}\right).
  $$
  在最坏情况下该部分可能退化到 $O(nm)$。
- 左端加权采样结构：$O(n)$ 级别（alias 或前缀和）。

---

### 5.8.2 采样阶段

朴素逐样本方式：每个样本调用一次 `SAMPLE(Q(r),1)`，其主要代价是一次范围遍历与块构造。在固定 $D$ 的常见平均模型下常写作
$$
\tilde O\!\left(m^{1-\frac{1}{2d}}\right),
$$
因此总采样代价近似
$$
\tilde O\!\left(t\cdot m^{1-\frac{1}{2d}}\right),
$$
最坏情况可到 $O(tm)$。

批量化方式：若将 $t$ 次左端抽样后分组，对每个出现过的左端 $r$ 只调用一次 `SAMPLE(Q(r),k_r)`，则总遍历次数等于不同左端的个数。设其为 $n_{\mathrm{dist}}$，则代价形态可写为
$$
\tilde O\!\left(n_{\mathrm{dist}}\cdot m^{1-\frac{1}{2d}} + t\right),
$$
其中 $t$ 来自回填写入与块内采样的线性部分。

---

### 5.8.3 空间

kd-tree 节点存 MBR（长度 $D$ 的 $\min/\max$）与常数结构字段，空间为
$$
S=O(mD).
$$
此外需要存储左侧度数数组与加权采样结构 $O(n)$。总空间量级
$$
S_{\mathrm{total}}=O(mD)+O(n).
$$

---

### 5.8.4 维度效应

本方法的查询维度为 $D=2d$。当 $d$ 增大时，
$$
1-\frac{1}{2d}\to 1,
$$
常见平均模型下的典型复杂度 $\tilde O(m^{1-1/(2d)})$ 会迅速趋近线性，导致 kd-tree 范围计数/采样在中高维上性能明显退化。这一维度依赖是该对照体系的主要结构性特征，也为后续更结构化的分解方法提供了对照基准。

---

# 第6章：实验详细方案

本章给出一套完整的实验方案，用于系统评估本文提出的 **SJS** 及两类对照基线在 **Join 采样（i.i.d. 且均匀、有放回）** 任务上的正确性、鲁棒性与性能表现。全文统一采用 **半开盒（half-open box）相交语义**：对任意两盒 $r,s$，仅当对所有维度 $k$ 同时满足
$$
L_k(r) < R_k(s)\ \wedge\ L_k(s) < R_k(r)
$$
才判定 $r\cap s\neq\varnothing$。贴边（如 $R_k(r)=L_k(s)$）不计为相交。扫描线事件序固定采用 **END-before-START**，以严格匹配半开语义。

---

## 6.1 实验目标与评估指标

### 6.1.1 任务定义：Join 上的 i.i.d. 均匀采样

给定两类 $d$ 维轴对齐半开盒集合 $R_c$ 与 $R_{\bar c}$，定义有序 Join 结果集合
$$
J=\{(r,s)\mid r\in R_c,\ s\in R_{\bar c},\ r\cap s\neq\varnothing\}.
$$
目标生成 $t$ 个样本 $Z_1,\dots,Z_t\in J$，满足 **相互独立且在 $J$ 上均匀（有放回）**。当 $|J|=0$ 时输出空序列。

### 6.1.2 评估指标

所有实验统一报告以下指标：

1. **正确性**
   - **几何正确性**：输出对 $(r,s)$ 必须满足半开语义下 $r\cap s\neq\varnothing$。
   - **分布正确性**：输出样本满足 $J$ 上 i.i.d. 且均匀（有放回）。小规模场景构造精确真值 $J$ 进行统计检验；中大规模场景用事件块边缘分布、块内均匀性与一致性约束做审计（详见 6.5）。

2. **性能**
   - **端到端耗时（Latency）**：生成 $t$ 个样本的总 wall-clock 时间。
   - **吞吐（Throughput）**：samples/sec。
   - **分阶段耗时**：将每个模型拆分为“事件准备/扫描更新/计数/枚举/采样/位置指派/回填/第二遍补齐”等阶段（详见 6.6、6.7）。
   - **峰值内存（Peak RSS）**：包含索引结构、缓存、输出缓冲区及临时结构。

3. **可复现性**
   - 固定 `master_seed`，并按“数据集标识 + 实验编号 + 运行编号 + 阶段编号 + 事件编号 + 调用计数”派生子种子，保证随机性独立且结果可复现。

---

## 6.2 实验平台与运行规约

### 6.2.1 硬件与软件环境（单机）

实验在单机服务器完成，所有模型在同一环境运行：

- CPU：x86-64 服务器级处理器（不少于 32 物理核）
- 内存：≥ 256GB（实际上我们的主机为1TB）
- 存储：NVMe SSD（用于数据落盘与加载）
- OS：Ubuntu 22.04 LTS
- 编译：C++17，`-O3 -march=native`
- Python：3.12（仅用于数据集下载与预处理；核心算法与计时点在编译态实现）

**运行约束**
- 统一采用 **单线程**（线程数固定为 1）。
- 固定 CPU 亲和性到同一物理核，避免调度抖动。
- 计时不包含数据集下载；数据加载与预处理阶段单独记录（不计入算法端到端采样耗时），算法计时从“输入已在内存并完成统一化表示”开始。

### 6.2.2 计时、热身与重复次数

- 每组配置先执行 **1 次热身运行**（不记录），用于稳定页缓存与内存分配路径。
- 正式记录 **5 次重复运行**，报告 **中位数**，并同时记录最小/最大值用于稳定性分析。
- 时间测量采用高精度单调时钟；内存采用 OS 级 RSS 统计。

**超限处理**
- 单次运行设置时间上限 2 小时；超过则标记为 `TIMEOUT`。
- 内存达到物理上限或触发 OOM 则标记为 `OOM`。
- 对 `TIMEOUT/OOM` 保留阶段日志与峰值资源占用，纳入结果汇总表。

---

## 6.3 数据集与预处理

本研究同时覆盖 **可控合成数据** 与 **三类真实数据**（建筑 / 地理 / CV），并统一到半开盒语义下评测。

### 6.3.1 合成数据集：alacarte-rectgen（可控输出密度）

合成数据采用 `alacarte-rectgen`：
- PyPI：https://pypi.org/project/alacarte-rectgen/
- 文档：https://dannhiroaki.github.io/Alacarte/
- 代码：https://github.com/DANNHIROAKI/Alacarte

该库生成两组 $d$ 维轴对齐半开盒集合 $R,S$，并使期望输出密度 $\alpha_{\mathrm{out}}$ 逼近目标：
$$
\alpha_{\mathrm{out}}=\frac{|J(R,S)|}{|R|+|S|},\quad
J(R,S)=\{(r,s)\in R\times S\mid r\cap s\neq\varnothing\}.
$$
生成过程以覆盖率参数控制尺度，并通过调参使 $\mathbb{E}[\alpha_{\mathrm{out}}]\approx \alpha^\star_{\mathrm{out}}$。

**生成配置（固定部分）**
除 $(N,d,\alpha_{\mathrm{out}})$ 外其余参数固定为：
- `universe=None`（默认 $[0,1)^d$）
- `volume_dist="normal"`, `volume_cv=0.25`
- `shape_sigma=0.5`
- `tune_tol_rel=0.01`
- `dtype=float32`

**规模定义**
- 两侧规模相同：$N:=|R|=|S|$，因此 $|R|+|S|=2N$。
- 合成数据主评测覆盖：$N\in[10^5,5\times 10^6]$、$d\in\{2,3,4,5,6\}$、$\alpha_{\mathrm{out}}\in\{0.1,0.3,1,3,10,30,100\}$。

**生成质量审计**
- 记录生成返回的 `info["alpha_expected_est"]` 与 `info["coverage"]`。
- 使用库提供的 pair-sampling 对 $R\times S$ 随机抽样估计 $\hat\alpha$，并记录相对误差
$$
\varepsilon_\alpha=\frac{|\hat\alpha-\alpha^\star_{\mathrm{out}}|}{\alpha^\star_{\mathrm{out}}}.
$$

---

### 6.3.2 建筑真实数据集：CMAB-Spatial-Join-0.08B

数据集：
- Hugging Face：https://huggingface.co/datasets/DannHiroaki/CMAB-Spatial-Join-0.08B
- Builder：https://github.com/DANNHIROAKI/CMAB-Spatial-Join-0.08B-Builder

该数据集以建筑物为对象，包含基础包围盒与扩张影响包围盒，并给出 4 个 workload level（影响范围强度不同）。数据以 Parquet 组织，按 level 与省份拆分为多个子集。

**字段（用于实验构造）**
- `xmin,ymin,xmax,ymax`：建筑基础 AABB
- `exmin,eymin,exmax,eymax`：扩张后的影响 AABB
- `d_m`：扩张尺度参数（作为影响强度记录）

**Join 场景输入（本章统一采用）**
- $R$：影响盒集合 $\{(exmin,eymin,exmax,eymax)\}$
- $S$：基础盒集合 $\{(xmin,ymin,xmax,ymax)\}$
- 维度：$d=2$
- 场景强度：level 越高影响盒越大，通常产生更大的结果集。

---

### 6.3.3 地理真实数据集：Geolife-Spatial-Join-0.15B

数据集：
- Hugging Face：https://huggingface.co/datasets/DannHiroaki/Geolife-Spatial-Join-0.15B
- Builder：https://github.com/DANNHIROAKI/Geolife-Spatial-Join-0.15B-Builder

该数据集基于 Geolife Trajectory 构造 encounter join：给定轨迹点 $p=(x,y,t)$，构造 encounter 盒
$$
b(p)=[x-\Delta d,\ x+\Delta d]\times [y-\Delta d,\ y+\Delta d]\times [t-\Delta t,\ t+\Delta t].
$$
提供 4 组阈值（group_1..4），对应不同 $(\Delta d,\Delta t)$：
- group_1：20m & 10min
- group_2：50m & 30min
- group_3：100m & 60min
- group_4：200m & 120min

数据坐标以厘米（cm）与毫秒（ms）存储为 int64。

**语义统一**
Geolife 数据的 encounter 盒以闭区间表达。为统一到本章半开语义，所有上界做一次整数刻度平移：
$$
R' = R + 1,
$$
从而将 $[L,R]$ 映射为 $[L,R')$，并与严格不等号的半开判定一致。

**Join 场景输入（本章统一采用）**
- 维度：$d=3$
- 坐标顺序：$(t,x,y)$（扫描维选 $t$）
- 两侧划分：按 `user_id` 奇偶划分
  - $R:=\{b(p)\mid user\_id\ \text{为偶}\}$
  - $S:=\{b(p)\mid user\_id\ \text{为奇}\}$
- 场景强度：group_1..group_4 四个阈值组分别作为 4 个测试场景。

---

### 6.3.4 CV 真实数据集：COCO-Spatial-Join-1.23B

数据集：
- Hugging Face：https://huggingface.co/datasets/DannHiroaki/COCO-Spatial-Join-1.23B
- Builder（参考实现与字段定义）：https://github.com/DANNHIROAKI/COCO-Spatial-Join-1B-Builder

该数据集基于 COCO（MS COCO 2017）目标检测标注与 RPN proposals 构建。对象以半开 3D 盒表示，其中 $z$ 维用于按图像切分（每张图像映射到 $[z,z+1)$ 的 slab），全局 join 等价于“每张图像内部 join 的总和”。

**字段（用于过滤与构造场景）**
- `x_min,x_max,y_min,y_max,z_min,z_max`
- `type,rank,score,category_id,coco_image_id`
  - `type=0`：annotation
  - `type=1`：proposal
  - proposal 在同一图像内按 `score` 形成 `rank`（从 0 开始）

**Join 场景输入**
- 维度：$d=3$
- 坐标顺序：$(z,x,y)$（扫描维选 $z$，使不同图像对象在扫描过程中天然分离）
- 任务 A（Detection Matching）：
  - $R$：annotations（`type=0`）
  - $S$：top-$k$ proposals（`type=1` 且 `rank<k`）
  - $k\in\{10,30,100,300\}$
- 任务 B（Proposal–Proposal Overlap）：
  - 取 top-$k$ proposals（`type=1` 且 `rank<k`）
  - 按 `rank` 奇偶分两侧：
    - $R=\{rank\ \text{even}\}$
    - $S=\{rank\ \text{odd}\}$
  - $k\in\{10,30,100,300\}$

---

### 6.3.5 数据统一化：坐标、ID 与扫描轴

为保证不同数据源在算法内部完全一致，本章统一：

1. **唯一 ID**：每个盒对象分配唯一 `id`（即便几何完全相同也视为不同元素）。
2. **扫描轴通过坐标重排实现**（算法固定扫描维为第 1 维）：
   - CMAB：$(x,y)$（扫描 $x$）
   - Geolife：$(t,x,y)$（扫描 $t$）
   - COCO：$(z,x,y)$（扫描 $z$）
3. **严格边界与重复坐标**：排序与边界比较统一采用键 $(x,\mathrm{id})$，并配合 END-before-START，避免贴边误计。

---

## 6.4 对照模型（7 个模型）与计时分解口径

### 6.4.1 模型集合

本章对比三类方法，总计 7 个模型：

1. **SJS（本文方法）**：事件块原语 `COUNT/REPORT/SAMPLE` 由“模式化分解 + 递归段树结构”实现，并组合三种采样框架：
   - **SJS-Enum**
   - **SJS-Adaptive**
   - **SJS-Sampling**

2. **RangeTree 基线**：以 $D=2(d-1)$ 维正交范围树实现事件块原语，并组合相同三种采样框架（仅替换原语实现）：
   - **RT-Enum**
   - **RT-Adaptive**
   - **RT-Sampling**

3. **KDS 基线（仅 Sampling）**：右侧集合做 $2d$ 维提升建立静态 kd-tree；先精确计算左侧度数 $w(r)$，再按度加权选择左端并在范围内均匀采样右端：
   - **KDS-Sampling**

### 6.4.2 三种采样框架（Enum / Adaptive / Sampling）

- **Enum**：单次扫描枚举并物化全部 $J$（通过 `REPORT(e)` 枚举每个事件块），之后在数组上做 i.i.d. 均匀抽样。时间与空间均为 $\Theta(|J|)$。
- **Adaptive**：两遍扫描。第一遍只 `COUNT` 得到每个事件块权重 $w_i$ 并完成位置指派；第二遍仅对被指派事件块调用 `SAMPLE(e_i,k_i)` 回填。
- **Sampling**：在 Adaptive 基础上引入预算 $B$，第一遍在预算内执行“小块全量缓存 + 预取采样”，尽量减少第二遍需要回填的事件块；当第一遍已覆盖全部回填需求时，第二遍直接跳过。

### 6.4.3 框架 Sampling 的预算参数（固定设定）

为保证可对照性，框架 Sampling 在所有数据集与所有配置下固定：

- 缓存预算：$B = 10^7$（以“可缓存的 partner 记录条数”计）
- 小块全量缓存阈值：$w_{\text{small}}=1024$
- 残差回填策略：
  - 对残差事件块 $i$，残差需求 $r_i$ 满足 $r_i \ge 0.25\,w_i$ 时，采用一次 `REPORT` 的流式有放回回填；
  - 否则调用 `SAMPLE(e_i,r_i)`。

---

## 6.5 实验一：鲁棒性测试 / 发烟测试 / Sanity Check

本实验以“强约束 + 差分对照 + 压力覆盖”的方式验证：实现严格遵守半开语义、事件序、事件块原语语义与全局 i.i.d. 均匀采样目标。

### 6.5.1 构建与最小发烟测试

对每次提交执行：

1. 编译三套实现（SJS / RangeTree / KDS），开启 `-O3` 与符号表。
2. 运行最小样例（$N=10^3,d=2,t=10^3$）检查：
   - 程序可执行且输出样本数量正确；
   - 所有输出对满足半开语义相交判定；
   - 统计信息（如 $W$、阶段耗时）可输出且无 NaN/溢出；
   - 在同一 `master_seed` 下重复运行两次得到逐字节一致的输出（用于锁定随机性与事件序确定性）。

### 6.5.2 半开语义与事件序一致性测试（贴边与同坐标）

构造手工用例覆盖边界情形：

- 贴边：$R_1(r)=L_1(s)$、$R_k(r)=L_k(s)$（任一维贴边都必须判定不相交）
- 极窄盒：$R_k(r)-L_k(r)$ 接近最小刻度（浮点与整数两类）
- 同坐标冲突：大量对象共享同一 $L_1$ 或 $R_1$

验证项：
- 相交判定严格使用 $<$；
- END-before-START 的处理使贴边不会进入同一事件块；
- 同坐标 tie-break（按 id）保证事件序完全确定。

### 6.5.3 事件块原语一致性测试（COUNT/REPORT/SAMPLE）

分别对 SJS 与 RangeTree 的事件块原语实现进行如下验证（随机选择事件与查询）：

1. **一致性**
   $$
   \texttt{COUNT}(e) = |\texttt{REPORT}(e)|
   $$
2. **无重复枚举**：`REPORT(e)` 输出元素不重复（每个 partner 恰好一次）。
3. **采样合法性**：`SAMPLE(e,k)` 输出的每个 partner 必属于 `REPORT(e)` 的集合。
4. **块内均匀性（小块精检）**
   - 选取 $w_e\le 200$ 的事件块；
   - 从该块采样 $10^6$ 次；
   - 对 partner 的出现频次做卡方检验（显著性水平 0.01），并记录统计量与 $p$ 值。

### 6.5.4 端到端分布真值验证（小规模）

在小规模数据上构造精确真值 $J$（暴力枚举或通过 Enum 物化），用于对“均匀与独立”的端到端验证。  
推荐设置：$N\le 5\times 10^4$，$d\in\{2,3,4\}$，$\alpha_{\mathrm{out}}\in\{1,10\}$，采样规模 $t=10^6$。

对每个模型输出序列 $Z_1,\dots,Z_t$ 执行：

#### (1) 几何正确性
验证每个输出对都属于真值集合 $J$（用哈希表或排序数组做成员查询）。

#### (2) 事件块边缘分布（聚合检验）
令 START 事件共有 $M$ 个，事件块权重为 $w_1,\dots,w_M$，$W=\sum_i w_i=|J|$。  
对输出序列记录事件索引序列 $I_1,\dots,I_t$（每个样本来自哪个事件块），理论概率为 $p_i=w_i/W$。

由于 $M$ 可能很大，直接对 $M$ 类做卡方检验会出现大量期望频数过小的问题。为此使用**分桶聚合**：  
- 以 $p_i$ 或 $w_i$ 为键，将事件划分为 $B$ 个桶（例如按对数区间或分位数分桶），并保证每个桶的理论期望次数不低于 50；  
- 对桶级命中次数做卡方检验，并报告桶划分方式、每桶期望次数最小值与 $p$ 值。

#### (3) 全局均匀性（哈希桶检验）
将真值集合 $J$ 按哈希函数 $h(r\_id,s\_id)\in\{0,\dots,B_u-1\}$ 划分为 $B_u$ 个桶。  
- 选择 $B_u$ 使每个桶的期望次数 $t\cdot |J_b|/|J|$ 不低于 50（必要时合并极小桶）；  
- 统计输出落入每个桶的次数，与理论桶概率 $|J_b|/|J|$ 做卡方检验。

该检验将“海量类别”的均匀性验证转化为“有限桶”的多项分布检验，满足检验前提且可解释。

#### (4) 独立性（相邻对检验 + 重复率检验）
**相邻对检验**：对输出序列的桶标号 $B_j=h(Z_j)$，构造相邻二元组 $(B_j,B_{j+1})$。在桶数较小（例如 64 或 128）时，统计二元组计数并与独立模型 $P(B_j=b,B_{j+1}=b')=P(B=b)P(B=b')$ 做卡方检验（确保期望频数充分）。  
**重复率检验**：记 $U$ 为 $t$ 个输出中的不同 pair 数量。在“从 $W=|J|$ 个元素均匀有放回抽样 $t$ 次”的模型下，
$$
\mathbb{E}[U]=W\left(1-\left(1-\frac{1}{W}\right)^t\right),
\qquad
\mathbb{E}[t-U]=t-\mathbb{E}[U].
$$
记录观测到的 $U$ 与 $t-U$，并与理论期望对比（必要时用泊松近似或自助法给出置信区间）。该统计量对非均匀或相关性引起的“过多重复/过少重复”具有敏感性。

### 6.5.5 差分测试（实现间对照）

在随机生成的小规模输入上进行差分对照：

- 对同一输入、同一 `master_seed`：
  - 对比 SJS 与 RangeTree 在每个事件的 `COUNT(e)` 是否一致；
  - 对比 `REPORT(e)` 的集合一致性（忽略输出顺序）；
  - 对比 KDS 在小规模情况下计算得到的 $W$ 是否与扫描法一致；
- 通过“逐事件对齐日志”定位任何偏差到具体事件、具体块与具体维度边界。

### 6.5.6 压力与退化场景（高强度稳定性）

覆盖三类极端输入：

1. **高重叠**：合成数据取 $\alpha_{\mathrm{out}}=100$，检验 `REPORT/SAMPLE` 吞吐与内存稳定性。
2. **高维**：合成数据取 $d=6$，检验结构构建、更新与查询在实现层面的边界稳定性。
3. **重复坐标/重复盒**：大量对象共享相同端点或完全相同几何参数，检验键 $(x,\mathrm{id})$ 与严格边界处理正确性。

---

## 6.6 实验二：合成数据集性能测试（7 模型对照）

本实验在可控合成数据上系统测量 7 个模型的端到端采样耗时，并拆解为细粒度阶段，展示算法随 $\alpha_{\mathrm{out}}$、$N$、$d$、$t$ 的可扩展性差异。

### 6.6.1 统一实验设置

- 数据生成：调用 `alacarte-rectgen` 的生成接口得到 $R,S$，记录 `info` 审计字段。
- 输入规模：两侧等大，$|R|=|S|=N$。
- 扫描维：合成数据维度对称，扫描维固定为第 1 维。
- Sampling 框架参数：固定为 $B=10^7$、$w_{\text{small}}=1024$（见 6.4.3）。
- 输出规模：除“变化 $t$”实验外固定 $t=100{,}000$。

### 6.6.2 参数网格（四组主实验）

四组主实验均对 7 模型完整运行并记录：

**(A) 变化 $\alpha_{\mathrm{out}}$：结果密度敏感性**
- 固定：$N=1{,}000{,}000$，$d=4$，$t=100{,}000$
- 取：$\alpha_{\mathrm{out}}\in\{0.1,0.3,1,3,10,30,100,300,1000\}$

**(B) 变化 $N$：输入规模可扩展性**
- 固定：$d=4$，$\alpha_{\mathrm{out}}=10$，$t=100{,}000$
- 取：
  $$
  N\in\{100{,}000,\ 200{,}000,\ 500{,}000,\ 1{,}000{,}000,\ 2{,}000{,}000,\ 5{,}000{,}000\}
  $$

**(C) 变化 $d$：维度效应**
- 固定：$N=1{,}000{,}000$，$\alpha_{\mathrm{out}}=10$，$t=100{,}000$
- 取：$d\in\{2,3,4,5,6\}$

**(D) 变化 $t$：输出规模可扩展性**

- 固定：$N=1{,}000{,}000$，$d=4$，$\alpha_{\mathrm{out}}=10$
- 取：
  $$
  t\in\{100{,}000,\ 300{,}000,\ 1{,}000{,}000,\ 3{,}000{,}000,\ 10{,}000{,}000,\ 30{,}000{,}000,\ 100{,}000{,}000,\ 300{,}000{,}000\}
  $$

### 6.6.3 采样耗时分解（统一口径）

#### 6.6.3.1 SJS / RangeTree（扫描线 + 事件块原语）的阶段定义

所有 SJS/RT 模型按以下阶段计时（差异体现在原语实现与框架调用路径）：

- **T0：事件准备**
  - 构造 START/END 事件数组
  - 按键 $(x,\text{type},\mathrm{id})$ 排序，且 END 严格先于 START
  - 预处理（坐标压缩 / 秩映射 / 索引初始化）

- **T1：扫描遍历**
  - **T1-upd**：活跃集更新（INSERT/DELETE）
  - **T1-qry**：事件块原语调用（COUNT / REPORT / SAMPLE）
  - **T1-meta**：缓存/堆/统计等元数据维护

- **T2：位置指派**
  - 计算 $W=\sum_i w_i$
  - 构建按 $w_i/W$ 的离散采样器
  - 对 $t$ 个位置独立抽取事件索引，生成每个事件的回填位置表 $L_i$

- **T3：回填与写回**
  - 对被选事件块调用 `SAMPLE` 或从缓存/预取样本取值，写回 $Z[1..t]$

- **T4：第二遍补齐（仅框架 Sampling 且存在残差时）**
  - 第二遍扫描，对残差事件块执行 `SAMPLE` 或流式 `REPORT` 并回填

#### 6.6.3.2 各框架的阶段组合

- **Enum**
  - T0 + T1（单遍扫描，`REPORT` 枚举物化全部 $J$）+ T2（在 `Pairs` 上均匀抽样）+ T3（输出写回）

- **Adaptive**
  - T0 + T1（Pass 1：全量 `COUNT`）+ T2（位置指派）+ T3（Pass 2：仅对 $k_i>0$ 的事件块 `SAMPLE(e_i,k_i)` 回填）

- **Sampling**
  - T0 + T1（Pass 1：`COUNT` + 小块 `REPORT` 缓存 + 预取 `SAMPLE`）+ T2 + T3（缓存/预取回填）+ T4（必要时第二遍补齐）

#### 6.6.3.3 KDS-Sampling（静态 kd-tree）的阶段定义

KDS-Sampling 按以下阶段计时：

- **K0：构建**
  - 构造 $2d$ 维点集 $P=\{p(s)\}$
  - 构建静态 kd-tree 索引

- **K1：度数计算**
  - 对每个 $r\in R_c$ 计算
    $$
    w(r)=|\{s\in R_{\bar c}\mid r\cap s\neq\varnothing\}|
    $$
  - 累加 $W=\sum_{r\in R_c} w(r)$

- **K2：左端加权采样器**
  - 构建按 $w(r)/W$ 的离散采样器（alias 或前缀和）

- **K3：生成 $t$ 个样本**
  - 按度选择左端 + 条件范围内均匀采样右端
  - 采用批量化：先抽取 $t$ 次左端并分组，对同一左端只做一次范围采样，再按位置写回（分布保持一致）

### 6.6.4 输出与图表规范

合成实验输出以下结果（每组配置至少一张图与一张表）：

- **总耗时与吞吐**：折线图（横轴为 $\alpha_{\mathrm{out}}$ / $N$ / $d$ / $t$）
- **分阶段耗时堆叠图**：展示 T0–T4 或 K0–K3 的时间组成
- **峰值内存表**：各模型 Peak RSS
- **可运行性汇总**：对 `TIMEOUT/OOM` 标注并保留阶段日志，用于解释 Enum 在大结果集下的不可运行性与各结构的内存边界

---

## 6.7 实验三：真实数据集 Join 场景与性能测试

本实验在三类真实数据上构造具有实际意义且结果规模较大的 join 场景，评测 7 模型的端到端采样耗时与分阶段组成。

### 6.7.1 统一设置

- 样本量两档：
  - $t=100{,}000$（在线分析级）
  - $t=1{,}000{,}000$（压力/吞吐级）
- Sampling 框架参数：$B=10^7$、$w_{\text{small}}=1024$
- 维度重排：按 6.3.5 固定（CMAB 扫描 $x$；Geolife 扫描 $t$；COCO 扫描 $z$）
- 对每个场景记录：$|R|,|S|,d,t,W$、阶段耗时、峰值内存、完成状态

---

### 6.7.2 CMAB（建筑）：影响范围邻接 Join

**场景含义**  
以“影响盒 $\bowtie$ 基础盒”刻画“某建筑影响范围内有哪些建筑”的邻接关系（覆盖分析、邻域查询、服务范围影响评估等）。

**场景集合（4 个 workload level）**
- CMAB-L1：level=1
- CMAB-L2：level=2
- CMAB-L3：level=3
- CMAB-L4：level=4（结果规模最大、重点压力场景）

**输入构造（每个 level 独立执行一次）**
- $R$：影响盒 $(exmin,eymin,exmax,eymax)$
- $S$：基础盒 $(xmin,ymin,xmax,ymax)$
- 维度：$d=2$

**评测输出**
- 7 模型在 CMAB-L1..L4 上的总耗时、分阶段耗时与峰值内存
- 同时输出各场景下 $W=|J|$（Adaptive/Sampling 由计数得到；Enum 由物化结果规模得到；KDS 由度数求和得到）

---

### 6.7.3 Geolife（地理）：时空 encounter Join

**场景含义**  
以 encounter 盒相交模拟“在给定空间与时间窗口内的相遇/同域活动”。

**场景集合（4 个阈值组）**
- Geo-G1：group_1（20m,10min）
- Geo-G2：group_2（50m,30min）
- Geo-G3：group_3（100m,60min）
- Geo-G4：group_4（200m,120min，结果规模最大、重点压力场景）

**输入构造（每个 group 独立执行一次）**
- 闭区间转半开：所有上界 $R\leftarrow R+1$（cm/ms）
- 两侧划分：
  - $R=\{user\_id\ \text{偶}\}$
  - $S=\{user\_id\ \text{奇}\}$
- 维度：$d=3/4$各跑一个版本，坐标顺序 $(t,x,y)/(t,x,y,z)$（扫描 $t$）

**评测输出**

- 7 模型在 Geo-G1..G4 的耗时与分阶段组成
- 重点分析 Geo-G4 上 Sampling 系列的端到端加速与阶段贡献（计数、预取、第二遍残差等）

---

### 6.7.4 COCO（CV）：检测匹配与 proposal 重叠 Join

**场景含义**
- Task A：annotation 与 proposal 的空间匹配（检测候选覆盖与匹配关系）
- Task B：proposal 与 proposal 的重叠关系（候选框冗余与聚集特征），通常产生更大的结果集，是重点压力场景

**场景集合**
- COCO-A10 / A30 / A100 / A300：Task A，$k\in\{10,30,100,300\}$
- COCO-B10 / B30 / B100 / B300：Task B，$k\in\{10,30,100,300\}$

**输入构造（每个 $k$ 独立执行一次）**
- 坐标顺序：$(z,x,y)$（扫描 $z$，按图像 slab 分解）
- Task A：
  - $R=\{type=0\}$（annotations）
  - $S=\{type=1\ \wedge\ rank<k\}$（top-$k$ proposals）
- Task B：
  - 先取 $T=\{type=1\ \wedge\ rank<k\}$（top-$k$ proposals）
  - 再按 rank 奇偶划分：
    - $R=\{p\in T\mid rank\ \text{even}\}$
    - $S=\{p\in T\mid rank\ \text{odd}\}$

**评测输出**
- 7 模型在 8 个 COCO 场景上的总耗时、分阶段耗时堆叠图与峰值内存
- 重点分析 Task B 在 $k=100,300$ 时的可扩展性与瓶颈阶段

---

## 6.8 实验记录与结果呈现格式

为保证实验过程与结果可复现、可追溯，每次运行输出一份结构化日志（JSON 或 CSV），至少包含：

- **实验标识**
  - `dataset`：`synthetic` / `cmab` / `geolife` / `coco`
  - `scenario`：如 `alpha=10`、`CMAB-L4`、`Geo-G4`、`COCO-B300`
  - `model`：7 模型之一
  - `run_id`：1..5（对应重复运行编号）

- **输入参数**
  - $|R|,|S|,d,t$
  - 框架参数：$B,w_{\text{small}}$（若适用）
  - 合成数据附加：$\alpha^\star_{\mathrm{out}}$、`alpha_expected_est`、$\hat\alpha$、$\varepsilon_\alpha$

- **结果标定**
  - $W=|J|$
  - `status`：`OK` / `TIMEOUT` / `OOM`

- **阶段耗时**
  - SJS/RT：`T0,T1_upd,T1_qry,T1_meta,T2,T3,T4`（不存在的阶段置 0）
  - KDS：`K0,K1,K2,K3`

- **资源占用**
  - `peak_rss_gb`
  - `peak_tmp_bytes`（如有）

- **随机性**
  - `master_seed`
  - `seed_derivation`（阶段与调用的派生规则版本号）

**结果呈现**
- 对合成实验：每组主实验（A/B/C/D）给出总耗时与吞吐折线图 + 分阶段堆叠图 + 内存表。
- 对真实数据：按场景给出总耗时与分阶段堆叠图，并同时报告 $W$ 以标定场景强度。
- 对 `TIMEOUT/OOM`：在图表中以缺失点或特殊标记呈现，并在附表中给出原因与阶段日志摘要（例如枚举导致内存爆炸、构建阶段超时等）。
