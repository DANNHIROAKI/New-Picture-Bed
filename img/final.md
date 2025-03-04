# 期末重点

**前言1：感谢20越杰朱学长，重点关注课件中的例子**

<img src="https://s2.loli.net/2023/12/30/wYVWbGc9DFzSrxI.png" alt="image-20231230195212866" style="zoom: 50%;" />  

**前言2：所有课件中的什么什么分类，老师的原话是这些都没意义，我觉得可以不用过分关注**

**前言3：太长的&特别重点的例题我都放到另一个PPT例题总结里了**

**前言4：奇怪的要求**

<img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/20b837e6170b7d73102f40f01336b78.jpg" alt="20b837e6170b7d73102f40f01336b78" style="zoom: 33%;" /> 

**前言5：感谢越杰徐学长**

<img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20240102195836480.png" alt="image-20240102195836480" style="zoom:50%;" /> 

**前言6：老师会捞**

<img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20240102201758670.png" alt="image-20240102201758670" style="zoom: 50%;" /> 

# 1. 绪论

> **1️⃣**人工智能是什么
>
> ==<font color='red'>**记下来**，曾经考过简答题，直接叫你阐述智能的三个观点，三言两语表达出核心就行</font>==
>
> 1. 思维理论：知识来自于思维，智能的核心是思考和推理能力
>
> 2. 知识阈值理论：智能取决于知识量及其一般化程度
> 3. 进化理论：有感知外界，适应环境的能力
>
> **2️⃣**智能的特征：具有感知能力，具有记忆和思维能力，能够学习和自适应，具有行为能力
>
> ==<font color='red'>**四点完全记下来**，曾经考过简答题，提问方式如：请阐述智能的具体特征</font>==
>
> **3️⃣**图灵测试：由能否在对话中欺骗人类，判断机器是否有人类智能；侧重知识阈值理论
>
> ==<font color='red'>了解</font>==
>
> **4️⃣**深蓝：IBM的超级计算机，首次击败国际象棋世界冠军；侧重思维理论
>
> ==<font color='red'>了解</font>==
>
> **5️⃣**图灵测试和深蓝哪个更难实现：
>
> 1. 图灵测试更难
> 2. 深蓝侧重思维理论，只需了解围棋的少量知识和规则，实质上没有智能
> 3. 图灵测试侧重的是知识阈值理论，需要有庞大的知识和理解能力
>
> ==<font color='red'>简答题重灾区，要搞明白如何解释这个问题</font>==
>
> **6️⃣**人工智能发展史
>
> 1. 孕育(1956前)：
>
> 2. 形成(1956-1969)：
>    - 1956：麦卡锡提出人工智能一词
>    - 1965：第一个专家系统诞生
>
> 3. 发展(1970后)：
>    - 1982：Hopfield网络提出
>    - 1986：BP网络提出
>    - 1995：SVM提出
>    - 2006：深度学习提出
>    - 2016：Alpha Go
>
> ==<font color='red'>了解一下，但老师又说不用记，看个人吧</font>==

# 2. 知识工程

> ## 2.1. 知识的概念
>
> > **1️⃣**知识的定义：有关信息关联在一起所形成的信息结构
> >
> > ==<font color='red'>了解</font>==
> >
> > **2️⃣**知识特性：相对正确，不确定性(随机/模糊/不完全/经验)，可表示可利用
> >
> > ==<font color='red'>**记下来**，非常中要，简答题重灾区</font>==
> >
> > **3️⃣**知识分类：常识/领域，事实(证据)/规则(知识)，确定/不确定
> >
> > ==<font color='red'>了解</font>==
>
> ## 2.2. 谓词逻辑知识表示
>
> > **1️⃣**一些重要的谓词等价式
> >
> > ```txt
> > P→Q       ⇔    ¬PvQ     //用于反证法
> > ```
> >
> > ==<font color='red'>不用刻意去背，用到了的(PPT例子中)再记下来(如上)</font>==
> >
> > **2️⃣**要重的永真蕴含式
> >
> > ==<font color='red'>不用刻意去背，用到了的(PPT例子中)再记下来(如上)</font>==
> >
> > **3️⃣**知识的谓词表示
> >
> >  比如人人都爱劳动：$\forall{x}(\text{Man{(x)}}\to{}\text{Love}(x,\text{labour}))$
> >
> > ==<font color='red'>了解一下，以前出过这样的题目</font>==
>
> ## 2.3. 产生式知识表示
>
> > **1️⃣**分为知识和规则的表示，规则表示为$P\to{}Q$或者`If P Then Q`
> >
> > ==<font color='red'>重点关注规则的表示法</font>==
> >
> > **2️⃣**产生式系统的构成(当然也是专家系统的原型)
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20240102133210542.png" alt="image-20240102133210542" style="zoom:40%;" /> 
> >
> > 1. 规则库：描述相应领域知识的产生式集合
> > 2. 事实库：存放已知事实和推导出的事实
> > 3. 控制系统：负责匹配规则的条件部分，匹配后执行规则，将执行结果加入数据库
> >
> > ==<font color='red'>**必考点**，简答题重灾区，出题方式为：请简述产生式系统的组成，以及各个模块的作用</font>==

# 3. 确定性推理

> ## 3.1. 知识匹配
>
> > **1️⃣**变量代换的定义：有限集$\theta{}=\{t_1/x_1,...,t_n/x_n\}$中
> >
> > 1. $t_i$为常量(也可能为变量or函数)，$x_i$为==互不相同==的变元，$/$表示左边换掉右边的
> > 2. <font color='red'>注意事项：$x_i$互不相同，$x_i$和$t_i$不能相同，$x_i$不能出现在另一个$t_j$中</font>
> > 3. 特例：$F\theta$表示用$\theta$替换去替换$F$中的变量
> >
> > ==<font color='red'>搞清楚这是怎么回事，主要是鲁滨逊归结原理要用到，但用到也只是简单的</font>==
> >
> > **2️⃣**代换的复合：对于$\theta{}=\{t_1/x_1,...,t_n/x_n\}$和$\lambda{}=\{u_1/y_1,...,u_m/y_m\}$
> >
> > 1. 令$\theta{}\lambda{}=\{t_1\lambda{}/x_1,...,t_n\lambda{}/x_n\}$，
> > 2. 合并为$(\theta{}\circ{}\lambda{})_{RAW}=\{\theta{}\lambda{},\lambda\}=\{t_1\lambda{}/x_1,...,t_n\lambda{}/x_n\,\,|\,\,u_1/y_1,...,u_m/y_m\}$
> > 3. 删除$\theta{}\lambda{}$中满足$t_i\lambda{}=x_i$的项，删除$\lambda$中满足$y_i∈\{x_1,...x_n\}$的
> > 4. 完成上述操作后$(\theta{}\circ{}\lambda{})_{RAW}\xrightarrow{形成复合}\theta{}\circ{}\lambda{}$
> >
> > ==<font color='red'>曾经出过计算题，掌握课件里的计算题，难度不会超过它</font>==
> >
> > **3️⃣**合一和最一般合一
> >
> > 1. 合一：使得公式集$F=\{F_1,...F_n\}$满足$F_1\theta{}=F_2\theta{}=...=F_n\theta{}$的==$\theta$==
> > 2. 最一般合一：公式集$F$的所有合一$\theta_i$都可以表示为$\theta{}_i=\sigma{}\circ{}\lambda{}_i$，$\lambda{}_i$为相应代换
> >
> > ==<font color='red'>以理解为重点，可能出简答题，要你解释什么是最一般合一</font>==
> >
> > ==<font color='red'>关注课件中求最一般合一的例子</font>==
>
> ## 3.2. 归结演绎推理
>
> > **1️⃣**什么是归结：利用反证法和归结原理来进行逻辑推理的方法
> >
> > 例如证明$P→Q(¬P\lor{}Q)$永真，只需证明$P\land{}¬Q$不可满足
> >
> > ==<font color='red'>可能作为名词解释题出现</font>==
> >
> > **2️⃣**子句和谓词公式：原子为此公式及其否定(aka文字)$\xrightarrow{互相\lor}$子句$\xrightarrow{一堆子句构成集合}$子句集
> >
> > ==<font color='red'>重点关注和理解，PPT中将谓词公式$\to$子句集的例子，曾经出过给出谓词公式，求出子句集的例子</font>==
> >
> > ==<font color='red'>PPT例子中，涉及到的永真蕴含和等价式都要记下来</font>==
> >
> > ==<font color='red'>谓词公式$\to$子句集是重点难点，例子我放后面了</font>==
> >
> > **3️⃣**子句集的意义：$P\to{}Q\iff{}\neg{}P\lor{Q}永真\iff{}P\land{}\neg{Q}永假\iff{}P\land{}\neg{Q}的子句集S不可满足$
> >
> > ==<font color='red'>了解即可</font>==
> >
> > <font color='red'>**4️⃣**Herbrand理论讲再多都不用看</font>
>
> ## 3.3. 鲁宾逊归结原理
>
> > **1️⃣**命题逻辑的归结
> >
> > 1. 互补文字：$P$和$\neg{P}$
> > 2. $C_1,C_2$是子句集中两个字句，$L_1$是$C_1$文字，$L_2$是$C_2$文字
> > 3. 如果$L_1L_2$互补，则分别在$C_1C_2$中删除二者
> > 4. 将现有$C_1C_2$余下部分析取，得到$C_{12}$
> > 5. $C_{12}$是$C_1C_2$的归结式，$C_1C_2$是$C_{12}$亲本子句，且有$C_1,C_2\to{}C_{12}$
> >
> > ==<font color='red'>了解即可</font>==
> >
> > **2️⃣**谓词逻辑的归结：==老师原话是说，PPT中这一个内容每个字都要看，每个例子都要看==
> >
> > 区别于命题逻辑，需要用最一般合一对变元代换后，再进行归结。示例：
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231121234844376.png" alt="image-20231121234844376" style="zoom:33%;" />  
> >
> > 字句本身内部若可以合一，则先内部合一
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20240102163617946.png" alt="image-20240102163617946" style="zoom:50%;" /> 
> >
> > ==<font color='red'>计算题重灾区</font>==
> >
> > **3️⃣**二元归结式：
> >
> > 1. 含义：$C_1,C_2$是无相同变元字句，$L_1,L_2$是其文字，$\sigma{}$是$L_1,\neg{}L_2$的最一般合一
> > 2. 则归结为：$C_{12} = (C_{1}\sigma - \{L_{1}{\sigma}\}) \cup (C_{2}\sigma - \{L_{2}{\sigma}\})$
> >
> > ==<font color='red'>真的需要看吗？不确定还是看一眼吧</font>== 
>
> ## 3.4. 归结反演
>
> > **1️⃣**是个啥东西：应用归结原理来证明定理$(P_1,P_2，...,P_n)\to{}Q$
> >
> > 1. $(P_1,P_2，...,P_n)\to{}Q$等价于$\neg{}(P_1\land{}P_2\land{},...,\land{}P_n)\lor{}Q$永真
> > 2. $\neg{}(P_1\land{}P_2\land{},...,\land{}P_n)\lor{}Q$永真等价于$(P_1\land{}P_2\land{},...,\land{}P_n)\land{}\neg{}Q$不可满足
> > 3. $(P_1\land{}P_2\land{},...,\land{}P_n)\land{}\neg{}Q$不可满足等价于==其子句集不可满足==
> > 4. 子句集的不可满足，可用归结来证明
> >
> > ==<font color='red'>了解，理解</font>==
> >
> > **2️⃣**归结反演的步骤：一般会把谓词公式的子句集告诉你？再让你去归结出空子句？
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231122001740141.png" alt="image-20231122001740141" style="zoom:50%;" /> 
> >
> > ==<font color='red'>至于归结反演的步骤，详见我后面给出的例子，这是计算题的重灾区</font>==
>
> ## 3.5. 归结策略
>
> > **1️⃣**归结的一般过程
> >
> > 1. $S=\{C_1,C_2,...,C_n\}$中任意两两归结，得到第一级归结式$S_1$
> > 2. 把$S,S_1$任意子句两两归结，得到二级归结式$S_2$
> > 3. 把$S,S_1,S_2$任意子句两两归结，得到三级归结式$S_3$
> > 4. 循环，直到出现空子句/不能继续归结，示例如下
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231122114814401.png" alt="image-20231122114814401" style="zoom:40%;" /> 
> >
> > ==<font color='red'>关注例子，原话是例子都要看，书上的例子也都要看</font>==
> >
> > **2️⃣**子句集归结策略：删除策略
> >
> > 1. 删除纯文字所在字句；假设$L$在子句集中，但$\neg{}L$不在子句集，那$L$就是纯文字
> >
> > 2. 删除重言式字句；一个字句中含有$P\lor{}\neg{P}$这样结构的，就是重言式字句
> >
> > 3. $C_1$包孕于$C_2$(aka$C_1\sigma{}\subseteq{C_2}$)时，删除$C_2$
> >
> > ==<font color='red'>应该要了解</font>==
> >
> > **3️⃣**支持集策略
> >
> > 1. 支持集：目标公式的否定得到的子句==(题目一般会告诉你)==+该字句归结出的后裔
> > 2. 支持集策略：每次归结，亲本子句中必须有一方属于支持集
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231122134944202.png" alt="image-20231122134944202" style="zoom: 33%;" /> 
> >
> > ==<font color='red'>要不要看老师说的含糊不清，我觉得还是关注一下例子吧</font>==
> >
> > **4️⃣**线性输入策略和祖先策略
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231122153225195.png" alt="image-20231122153225195" style="zoom:33%;" /> 
> >
> > 1. 线性输入策略：$C_1,C_2$归结，至少一个是初始子句集中的子句
> > 2. 祖先过滤策略：$C_1,C_2$归结，至少一个是初始子句集中的子句，或者$C_1$是$C_2$的祖先
> >
> > ==<font color='red'>这里可能会出**简答题**，比如：请简述祖先过滤策略。记下来</font>==
>
> ## 3.6. 用归结证明命题 
>
> > ==<font color='red'>重点关注ABC三人谁说谎了的那个例子，在PPT例子总结里</font>==

# 4. 不确定性推理

> ## 4.1. 主观Bayes方法
>
> > **1️⃣**知识(规则)的不确定性：产生式表示
> >
> > `IF E THEN (LS,LN) H(P(H))`
> >
> > ==<font color='red'>计算题重灾区，比如告诉你LS，LN等参数，然后要计算后验概率(E发生情况下，H发生的概率)</font>==
> >
> > ==<font color='red'>LS和LN的定义公式不用记，不会这么考</font>== 
> >
> > **2️⃣**证据(事实)的不确定性：对于证据$E$，由用户根据观测$S$，给出==动态强度==$P(E|S)$
> >
> > ==<font color='red'>这个重在理解</font>==
> >
> > **3️⃣**可信度$C(E|S)$​：以PROSPECTOR(矿勘专家系统)为例
> >
> > | $C(E \mid S) \in\{-5, \ldots, 5\}$ |                      含义                      |
> > | :--------------------------------: | :--------------------------------------------: |
> > |                 -5                 | 观测 $S$ 下证据 $E$ 不存在，即 $P(E \mid S)=0$ |
> > |                 5                  |   观测 $S$ 下证据 $E$ 存在， $P(E \mid S)=1$   |
> > |                 0                  |      $S, E$ 之间无关， $P(E \mid S)=P(E)$      |
> >
> > $$
> > \begin{flalign*}
> > & P(E|S) = \begin{cases}
> > \cfrac{C(E|S) + P(E) \times (5 - C(E|S))}{5}, & 0 \leq C(E|S) \leq 5 \\
> > \cfrac{P(E) \times (5 + C(E|S))}{5}, & -5 < C(E|S) < 0
> > \end{cases} &
> > \end{flalign*}
> > $$
> >
> > ==<font color='red'>了解即可</font>==
> >
> > **4️⃣**组合证据的不确定性
> >
> > 1. 二维组合
> >
> >    |                     |       最大最小法       |           概率法           |            有界法            |
> >    | :-----------------: | :--------------------: | :------------------------: | :--------------------------: |
> >    | $T(E 1 \land E 2 )$ | $min\{T(E_1),T(E_2)\}$ |      $T(E_1)*T(E_2)$       | $max{}\{0,T(E_1)+T(E_2)-1\}$ |
> >    |  $T(E1 \lor E 2 )$  | $max\{T(E_1),T(E_2)\}$ | $T(E1)+T(E2)－T(E1)*T(E2)$ |   $min\{1,T(E_1)+T(E_2)\}$   |
> >
> > 2. 多维组合
> >
> >    $P(E_1\land{}E_2\land{}...\land{}E_n|S)=min\{P(E_1|S),P(E_2|S),...,P(E_n|S)\}$
> >
> >    $P(E_1\lor{}E_2\lor{}...\lor{}E_n|S)=max\{P(E_1|S),P(E_2|S),...,P(E_n|S)\}$
> >
> >    $P(\neg{}E|S)=1-P(E|S)$
> >
> > ==<font color='red'>了解一下，其实也不难理解</font>==
> >
> > **5️⃣**几率函数和概率函数的转换：$\Theta{}(x)=\cfrac{P(x)}{1-P(x)}$，$P(x)=\cfrac{\Theta{}(x)}{1+\Theta{}(x)}$
> >
> > ==<font color='red'>一定牢记</font>==
> >
> > **6️⃣**不确定性的传递
> >
> > 1. 证据一定存在时：$P(E)=P(E|S)=1$，记住以下推导过程和结论
> >
> >    - 有Bayes可得：$P(H|E) = \cfrac{P(E|H) \times P(H)}{P(E)}$，$P(\neg H|E) = \cfrac{P(E|\neg H) \times P(\neg H)}{P(E)}$
> >    - 二者相除，结和$\text{LS}$几率函数定义，可得==$\Theta{}(H|E)=LS*\Theta{}(H)$==
> >
> >    ==<font color='red'>理解推导过程，记住结论，这是计算题的重灾区，比如告诉你$\text{LS}$和先验概率$P(H)$，求后验概率</font>==
> >
> > 1. 证据不存在时：$P(E)=P(E|S)=0$，记住以下推导过程和结论
> >
> >    - 有Bayes可得：$P(H|\neg{}E) = \cfrac{P(\neg{}E|H) \times P(H)}{P(\neg{}E)}$，$P(\neg H|\neg{}E) = \cfrac{P(\neg{}E|\neg H) \times P(\neg H)}{P(\neg{}E)}$
> >    - 二者相除，结和$\text{LN}$几率函数定义，可得==$\Theta{}(H|\neg{}E)=LN*\Theta{}(H)$==
> >
> >    ==<font color='red'>要求同上</font>==
> >
> > 1. 证据不确定时：$0<P(E|S)<1$
> >
> >    - $P(H|S) = P(H|E) * P(E|S) + P(H|\neg E) * P(\neg E|S)$
> >    - 主观Bayes的终极结论$P(H|S) =后面一堆$
> >
> >    ==<font color='red'>终极结论不用记，但是要明白它是通过分段差值得到的</font>==
> >
> > **7️⃣** LS和LN 的性质
> >
> > 1. $LS>1$说明证据$E$对结论$H$有利，$LN>1$说明证据$\neg{}E$对结论$H$有利
> >
> > 2. <font color='red'>记住：一般情况取$LS>1\,,LN<1$</font>
> >
> > ==<font color='red'>LS，LN的计算公式都不用记，理解LS和LN的含义，但是一定要记住红字部分</font>==
> >
> > **8️⃣**结论不确定性的合成：
> >
> > ```txt
> > if 前提E1 then (LS1,LN1) 结论H(p(H))---->对应观察S1
> > if 前提E2 then (LS2,LN2) 结论H(p(H))---->对应观察S2
> > ..........
> > if 前提En then (LSn,LNn) 结论H(p(H))---->对应观察Sn
> > ```
> >
> > $A_1,A_2,...,A_n\xrightarrow{合成算法}A$的不确定性为
> >
> > $\Theta(H|S_1,S_2,\ldots,S_n) = \cfrac{\Theta(H|S_1)}{\Theta(H)} \times \cfrac{\Theta(H|S_2)}{\Theta(H)} \times \cdots \times \cfrac{\Theta(H|S_n)}{\Theta(H)} \times \Theta(H)$
> >
> > ==<font color='red'>老师原话是：作为常识，了解一下，记下来(懵逼)</font>==
>
> ## 4.2. 带阈值的可信度模型
>
> > **1️⃣**表示不确定性：`IF E THEN H(CF(H,E),λ)`
> >
> > 1. $CF(E)∈{[0,1]}$：证据$E$的可信度，$CF(E)=1/0$表示证据绝对存在/不存在
> > 2. $CF(H,E)∈[0,1]$：证据$E$为真时$H$的可信度，$CF(H,E)=0,1$对应$P(H|E)=0,1$
> > 3. $\lambda{}∈[0,1]$：阈值，==只有证据$E$可信度$CF(E)\geq\lambda{}$时知识$E$才能被应用==
> > 4. ==$CF(H)=CF(H,E)\times{}CF(E)$==
> >
> > ==<font color='red'>计算题重灾区，比如以下的例题</font>==
> >
> > 给定CF(H,E)=0.8，CF(E)=0.7，λ=0.6问该规则可不可以启用。当然可以
> >
> > 启用后算出结论的可信度。证据0.7可信，证据成立后结论0.8可信，所以结论可信度为0.8*0.7=0.56
> >
> > **2️⃣**结论不确定性合成算法：`IF Ei THEN H(CF(H,Ei),λi)   i=1,2,...,n`
> >
> > 1. 前提：$\forall{i}∈\{1,2,...,n\},\,CF(E_i)\geq\lambda{}_i$
> >
> > 2. 结论$H$的可信度：
> >
> >    | 方法     | $CF(H)$等于                                                  |
> >    | -------- | ------------------------------------------------------------ |
> >    | 极大值法 | $\max\{CF_i(H)\}\,\,i=1,2,...,n$                             |
> >    | 加权和法 | $\cfrac{\sum\limits_{i=1}^{n} CF(H, E_i) \times CF(E_i)}{\sum\limits_{i=1}^{n} CF(H, E_i)}$ |
> >    | 有限和法 | $min\{\sum\limits_{i=1}^{n}CF_i(H),\,1\}$                    |
> >
> > ==<font color='red'>复习课老师对着这几个方法讲了贼久，最好还是看一下</font>==
>
> ## 4.3. 加权可信度模型
>
> > **1️⃣**知识的不确定性：
> >
> > 1. 表示：`IF E1(w1) AND E2(w2) AND ... AND  En(wn)  THEN  H (CF(H,E),λ)`
> > 2. 阈值$\lambda{}$：由专家给出
> > 3. 加权因子$\omega_i$：由专家给出，满足$\sum\limits_{i=1}^{n}\omega_i=1$
> >
> > **2️⃣**组合证据的不确定性：$CF(E)=\cfrac{\sum\limits_{i=1}^{n}[\omega_i*CF(E_i)]}{\sum\limits_{i=1}^{n}\omega_i}$
> >
> > **3️⃣**不确定性传递算法：当$CF(E)\geq\lambda{}$时$CF(H)=CF(H,E)\times{}CF(E)$
> >
> > ==<font color='red'>一定要了解，计算题重灾区，很容易出计算题</font>==
>
> ## 4.4. 前件带不确定性的可信度模型
>
> > **1️⃣**表示：
> >
> > 1. 表示1：$cf_i$是==专家==给出给出的子==条件==$E_i$的可信度，$cf_i'$是基于==观测==得到的子<mark>证据</mark>$E_i$的可信度
> >
> >    `IF E1(cf1) AND E2(cf2) AND ... AND  En(cfn)  THEN  H (CF(H,E),λ)`
> >
> > 2. 表示2：加上权值
> >
> >    `IF E1(cf1,w1) AND E2(cf2,w2) AND ... AND  En(cfn,wn)  THEN  H (CF(H,E),λ)`
> >
> > **2️⃣**不确定性匹配：
> >
> > 1. 无加权因子：$\sum\limits_{i=1}^{n} \max(0, c_{f_i} - c'_{f_i}) \leq \lambda$
> > 2. 有加权因子：$\sum\limits_{i=1}^{n} \omega_i \times \max(0, c_{f_i} - c'_{f_i}) \leq \lambda$
> >
> > **4️⃣**不确定性的传递算法
> >
> > 1. 无加权因子：$CF(H) = \left[ \prod\limits_{i=1}^{n} (1 - \max\{0, c_{f_i} - c'_{f_i}\}) \right] \times CF(H,E)$
> >
> > 2. ==有加权因子：$CF(H) = \left[ \prod\limits_{i=1}^{n} (1 - \omega_i\max\{0, c_{f_i} - c'_{f_i}\}) \right] \times CF(H,E)$==
> >
> > ==<font color='red'>一定记下来，很容易出计算题</font>==

# 5. 模糊推理

> ## 5.1. 模糊概念
>
> > **1️⃣**$\mu_A:U\to{[0,1]}$，将任意$u∈U$映射到$[0,1]$上的某个函数
> >
> > 1. 模糊集：$A=\{\mu_A(u),u∈U\}$称为$U$上的一个模糊集，
> > 2. $\mu_A$：定义在$U$上的模糊集$A$的<mark>隶属函数</mark>
> > 3. $\mu_A(u)$：$u$对模糊集$A$的隶属度
> >
> > **2️⃣**离散论域$U$的模糊集$A$，以下表示中应剔除$\mu_A(u_i)=0$的
> >
> > 1. 表示1：$A=\{\mu_A(u_1),\mu_A(u_2),...,\mu_A(u_n)\}$
> > 2. 表示2：$A=\mu_A(u_1)/u_1+\mu_A(u_2)/u_2+...+\mu_A(u_n)/u_n= \sum_{i=1}^{n} {\mu_A(u_i)}/{u_i}$
> > 3. 表示3：$A=\{\mu_A(u_1)/u_1,\mu_A(u_2)/u_2,...,\mu_A(u_n)/u_n\}=\bigcup_{i=1}^{n} {\mu_A(u_i)}/{u_i}$
> > 4. 表示4：$A=\{[\mu_A(u_1),u_1],[\mu_A(u_2),u_2],...,[\mu_A(u_n),u_n]\}$
> >
> > **3️⃣**$U$上所有模糊集表示为：$\mathcal{F}(U)$
> >
> > ==<font color='red'>了解</font>==
>
> ## 5.2. 模糊集的交并补运算
>
> > **1️⃣**$A \cup B : \mu_{A \cup B}(u) = \max\{\mu_{A}(u), \mu_{B}(u)\} = \mu_{A}(u) \vee \mu_{B}(u)$
> >
> > **2️⃣**$A \cap B : \mu_{A \cap B}(u) = \min\{\mu_{A}(u), \mu_{B}(u)\} = \mu_{A}(u) \wedge \mu_{B}(u)$
> >
> > **3️⃣**$\neg A : \mu_{\neg A}(u) = 1 - \mu_{A}(u)$
> >
> > ==例如==$\mu_{A}(1)=0.3/1\,,\mu_{B}(u)=0.4/1$则
> >
> > - $\mu_{A \cup B}(1)=\text{max}(0.3,0.4)/1=0.4/1$
> > - $\mu_{A \cap B}(1)=\text{min}(0.3,0.4)/1=0.3/1$
> > - $\mu_{\neg A}(1) = [1 - 0.3]/1=0.7/1$
> >
> > 对所有$u$求出$\mu$值，就得到了对应集合
> >
> > ==<font color='red'>会考计算题吧，以前有过给出AB求二者并集的题目</font>==
> >
> > ==<font color='red'>但是老师有说不要求，可能是想说不要求记下来，但要求会算吧</font>==
>
> ## 5.3. 模糊关系
>
> > **1️⃣**模糊关系的示例：$\mu_R(u_{1}, v_{1})$表示$u_1$信任$v_1$的程度
> > $$
> > \begin{flalign*}
> > R = \begin{bmatrix}
> > \mu_R(u_{1}, v_{1}) & \mu_R(u_{1}, v_{2}) & \cdots & \mu_R(u_{1}, v_{n}) \\
> > \mu_R(u_{2}, v_{1}) & \mu_R(u_{2}, v_{2}) & \cdots & \mu_R(u_{2}, v_{n}) \\
> > \vdots & \vdots & \ddots & \vdots \\
> > \mu_R(u_{m}, v_{1}) & \mu_R(u_{m}, v_{2}) & \cdots & \mu_R(u_{m}, v_{n})
> > \end{bmatrix}
> > \quad \to \quad
> > \left[ \begin{array}{ccc}
> > 1 & 0.3 & 0.8 \\
> > 0.9 & 1 & 0.6 \\
> > 0.7 & 0.5 & 1 \\
> > \end{array} \right]
> > \end{flalign*}
> > $$
> > ==<font color='red'>这个看看就行</font>==
> >
> > **2️⃣**模糊关系的合成：注意是==$\text{Max}\{min\{....\}\}$=
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231124004226145.png" alt="image-20231124004226145" style="zoom:44%;" /> 
> >
> > ==<font color='red'>非常重点，很容易就出题，一定记下这个例子</font>==
>
> ## 5.4. 模糊匹配
>
> > **1️⃣**含义：知识的$A$与证据的$A'$不一定完全相同
> >
> > 1. 规则：`IF (x_is_A) THEN (y_is_B)`
> > 2. 事实：`(x_is_A')`
> >
> > **2️⃣**如何衡量知识的$A$与证据的$A'$到底多匹配呢？
> >
> > 先算出语义距离
> >
> > 1. 海明距离：$d(A, B) = \cfrac{1}{n} \times \sum\limits_{i=1}^{n} |\mu_A(u_i) - \mu_B(u_i)|$
> > 2. 几何距离：$d(A, B) = \sqrt{\cfrac{1}{n} \times \sum\limits_{i=1}^{n} (\mu_A(u_i) - \mu_B(u_i))^2}$
> >
> > 得到匹配度$1-d$
> >
> > ==<font color='red'>一定要会算这两个距离，这是计算题重灾区，以前有考过告诉你两个模糊集(离散的)，要你算他们的匹配度</font>== 
>
> ## 5.5. 简单模糊推理
>
> > **1️⃣**简单模糊推理的基本规则
> >
> > 1. 对于规则`IF (x_is_A) THEN (y_is_B)`
> > 2. 需要构造一个$A,B$之间的模糊关系$R$，构造方法见下一条
> > 3. 然后通过R 与证据的合成求出结论
> >    - 如果已知证据`x_is_A'`，则==$B'=A'\circ{}R$==
> >    - 如果已知证据`x_is_B'`，则==$A'=R\circ{}B'$==
> >    - 可以想一种不存在的东西来记忆，比如 Ar(a)b Bar，另外证据必须在等号右边
> >
> > ==<font color='red'>极有可能考概念题，要能讲清楚模糊推理的基本过程</font>== 
> >
> > ==<font color='red'>也有可能考计算题，而且老师承诺一定会告诉你$R$</font>==
> >
> > **2️⃣**$R$的扎得构造方法：==只要求记住$R_m = (A \times B)\cup (\neg A \times V)$==，其他任何都不要了
> >
> > ==<font color='red'>这个一定出一道10分大题，关注例题(在PPT例题总结里)</font>==
> >
> > **3️⃣**$R$的曼德拉构造：$R_c = A \times B$
> >
> > ==<font color='red'>这个非常非常简单，了解一下就可以了，然后</font>==
> >
> > ==<font color='red'>みずもお方法一律不看</font>==
>
> ## 5.6. 各种模糊关系的性能分析
>
> > $$
> > \begin{cases}
> > A=\int\limits_U{}\mu_A(u)/u
> > \\
> > \text{very } A =\int\limits_U{}\mu_A^2(u)/u
> > \\
> > \text{more/less } A =\int\limits_U{}\mu_A^{0.5}(u)/u
> > \end{cases}
> > \\
> > \begin{cases}
> > \text{not } A =\int\limits_U{}1-\mu_A^{}(u)/u
> > \\
> > \text{not very } A =\int\limits_U{}1-\mu_A^{2}(u)/u
> > \\
> > \text{not more/less } A =\int\limits_U{}1-\mu_A^{0.5}(u)/u
> > \end{cases}
> > $$
> >
> > 示例
> >
> > | 状态             | 参数 | 处理方式    | 模糊集                                     |
> > | ---------------- | ---- | ----------- | ------------------------------------------ |
> > | 原始状态         | A    | \           | `1  0.80  0.60  0.40  0.20  0  0  0  0  0` |
> > | very             | A    | 平方        | `1  0.64  0.36  0.16  0.04  0  0  0  0  0` |
> > | more or less     | A    | 开根号      | `1  0.89  0.77  0.63  0.45  0  0  0  0  0` |
> > | not              | A    | 用1减       | `1  0.20  0.40  0.60  0.80  1  1  1  1  1` |
> > | not very         | A    | 用1减平方   | `1  0.36  0.64  0.84  0.96  1  1  1  1  1` |
> > | not more or less | A    | 用1减开根号 | `1  0.11  0.23  0.37  0.55  1  1  1  1  1` |
> >
> > ==<font color='red'>记住这个表格，会出计算题：给出A，然后计算下面的各种A</font>== 
>
> ## 5.7. 模糊三段论
>
> > **1️⃣**对于三个模糊关系，他们满足==$R(A,B)\circ{}R(B,C)=R(A,C)$==时，说明三段论成立
> >
> > ```txt
> > R1: IF x_is_A THEN y_is_B
> > R2: IF y_is_B THEN z_is_C
> > -------------------------
> > R3: IF x_is_A THEN z_is_C
> > ```
> >
> > **2️⃣**示例
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20240102194717568.png" alt="image-20240102194717568" style="zoom: 50%;" /> 
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20240102194800653.png" alt="image-20240102194800653" style="zoom:50%;" /> 
> >
> > 注意AB，BC，AC之间可构造不同的模糊关系
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20240102194945865.png" alt="image-20240102194945865" style="zoom:50%;" /> 
> >
> > ==<font color='red'>重点关注这个例子，一定亲笔做一下，这他妈是个原题</font>==
>
> ## 5.8. 多维模糊推理
>
> > ```txt
> > 知识：IF x1_is_A1 AND x2_is_A2 AND ... AND xn_is_An THEN y_is_B
> > 证据：   x1_is_A1'    x2_is_A2'    ...     xn_is_An'
> > ---------------------------------------------------------------
> > 结论：                                                   y_is_B'
> > ```
> >
> > **1️⃣**扎德方法：前提是$A_i$在相同的论域==，这个要记下来==
> >
> > 1. 所有$A_i$求交集得到$A$，所有$A_i’$求交集得到$A'$
> > 2. 求出$A,B$之间的模糊关系为$R(A,B)$
> > 3. 由$B'=A'\circ{}R(A,B)$，示例如下
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20240102200204265.png" alt="image-20240102200204265" style="zoom:50%;" /> 
> >
> > ==<font color='red'>一定要记下来，一定要了解基本思想，关注例子</font>==
> >
> > **2️⃣**つかもと方法：
> >
> > 1. 先求出每个$R(A_i,B),i=1,2,...,n$
> > 2. 求出$B_i'=A_i'\circ{}R(A_i,B),i=1,2,...,n$
> > 3. 最后$B' = B'_1 \cap B'_2 \cap \ldots \cap B'_n$
> >
> > ==<font color='red'>这个是要求的，也要记一记</font>==
>
> ## 5.9. 多重模糊推理
>
> > ```txt
> > 知识:IF x_is_A THEN y_is_B ELSE y_is_C
> > 证据:   x_is_A'
> > 结论：                          y_is_D
> > ```
> >
> > **1️⃣**通过A 、B 、C 构造模糊关系R
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20240102201419744.png" alt="image-20240102201419744" style="zoom: 67%;" /> 
> >
> > **2️⃣**求得结论D，==但是记住这个==
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20240102201442269.png" alt="image-20240102201442269" style="zoom: 67%;" /> 
> >
> > ==<font color='red'>了解一下这个规则，公式不用记，因为真要有这样的题目，R如何计算都会给出</font>== 

# 6. 搜索策略

> ==<font color='red'>和与或树有关的一律不看</font>==
>
> ## 6.1. 基本概念
>
> > **1️⃣**搜索：在知识库中找到可利用的知识，构造一条代价较小的推理路线，从而解决问题
> >
> > **2️⃣**状态空间：所有可能的状态+一切可用算符
> >
> > 1. 状态：描述问题求解过程中不同时刻的状态
> > 2. 算符：使得状态转换的操作
> >
> > ==<font color='red'>内化成常识</font>==
>
> ## 6.2. 搜索的一版过程即A*算法
>
> > **1️⃣**有关数据结构
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231125163352256.png" alt="image-20231125163352256" style="zoom:38%;" /> 
> >
> > 1. OPEN表：存放刚生成的结点，不同策略下结点的排列顺序不同
> >
> > 2. CLOSE表：存放将要扩展(生成子节点)的节点
> >
> > **2️⃣**算法流程：气泡中为$A^{*}$算法相较于一般搜索过程的具象，==这个一定会出大题==
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/搜索一般过程Axing.png" alt="搜索一般过程Axing" style="zoom: 10%;" />  
> >
> > 1. OPEN表中节点按估价函数$f(x)=g(x)+h(x)$小至大排序
> > 2. $g(x)$是初始$S_0\to{}x$路径的代价，$g^{*}(x)$是$S_0\to{}x$最小代价，$g(x)$可估计$g^{*}(x)$
> > 3. ==$h^{*}(x)$是$x\to{}目标节点$最小代价，$h(x)\leq{}h^{*}(x)$==
> >
> > **3️⃣**$A*$算法的性质：
> >
> > 1. 最优性：$h(x)$越大启发信息就越多，搜索效率越高
> >
> > 2. 单调性==(这个很重要)==：要求启发函数$h(x)$单调搜索代价会更低
> >
> >    - ==定义==：$h(S_g)=0$，$h(x_i)\leq{}h(x_j)+c(x_i,x_j)$
> >    - ==$g(x)=g^{*}(x)$是A*算法可操作的关键==
> >
> > ==<font color='red'>非常非常重要，A*算法一定会出大题，我认为应该补充一些A星算法的例题</font>== 
> >
> > ==<font color='red'>h(x)是什么，A*算法的单调性定义，背下来，可能会考</font>== 
>
> ## 6.3. 深度优先(左/首部)&广度优先(右/尾部)
>
> > ![image-20240102210959138](https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20240102210959138.png)
> >
> > ==<font color='red'>深度优先和广度优先一定要看懂</font>==
> >
> > ==<font color='red'>尤其关注PPT上那个重排九宫的例子，基本考试的时候就改几个数据，原封不动</font>==
> >
> > ==<font color='red'>有界深度优先不看</font>==
>
> ## 6.4. 代价树的广度优先搜索
>
> > **1️⃣**基本概念
> >
> > 1. 代价树：边上标有代价的树
> > 2. 代价：$g(x)$是初始结点$S_0$到结点$x$的代价，$c(x_1,x_2)=g(x_2)-g(x_1)$是两结点之间的代价
> >
> > **2️⃣**基本思想：==OPEN表中的节点按代价从小到大排序==，每次取出==代价最小==的结点
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/sx.png" alt="sx" style="zoom: 40%;" /> 
> >
> > ==<font color='red'>对应的深度优先不看</font>==
>
> ## 6.5. 启发式搜索，==一定要看==
>
> > **1️⃣**启发性信息：指导搜索，与问题有关的信息
> >
> > **2️⃣**估价函数：评估结点重要性，$f(x)=g(x)+h(x)$
> >
> > |      | $g(x)$                   | $h(x)$启发函数                 |
> > | ---- | ------------------------ | ------------------------------ |
> > | 含义 | 初始节点→x的代价         | x→目标结点，最优路径的预估代价 |
> > | 特点 | 有利于完备性，不利于效率 | 不利于完备性，有利于效率       |
> >
> > ➕重拍九宫中，$h(x)$是$x$的格局与目标节点格局不相同的牌数
> >
> > **3️⃣**全局择优搜索：一个结点扩展后，从OPEN表全体中，选一个估价最小的结点下一个考察
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/30750.jpg" alt="30750" style="zoom:33%;" />  
> >
> > ==<font color='red'>尤其关注PPT上那个重排九宫的例子，基本考试的时候就改几个数据，原封不动</font>==

# 7. 决策树

> ## 7.1. 决策树的构建
>
> > ```txt
> > 水果分类
> > 特征向量：{颜色，尺寸，形状，味道}
> > 四种属性：颜色(红绿蓝)，尺寸(大中小)，味道(酸甜)，形状(圆细)
> > 样本集合：一批水果，知道其特征向量+类别
> > 问题：一个新水果，观测得其特征向量，如何为其分类？
> > ```
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231126191336583.png" alt="image-20231126191336583" style="zoom: 50%;" /> 
> >
> > ==<font color='red'>对这个图要有映像</font>==
>
> ## 7.2. ID3算法：决策树的构建，核心在于分类精度
>
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/338.jpg" alt="338" style="zoom: 50%;" /> 
> >
> > ==<font color='red'>老师明确说了ID3算法的伪代码不用记，稍微关注一下吧</font>==
> >
> > ==<font color='red'>但需要注意，ID3算法的构建方向，就是最大化信息增益，提高了分类准确性</font>==
>
> ## 7.3. 熵和信息增益
>
> > **1️⃣**样本$D$的熵：$\text{Entropy}(D) = - \sum\limits_{i=1}^{c} P_i \log_2 (P_i)$
> >
> > 1. 样本来自$c$个不同类别
> > 2. $P_i$为$i$类样本所占比例
> > 3. 熵越小样本越纯净
> >
> > **2️⃣**信息增益：$\text{Gain}(D, A) = \text{Entropy}(D) - \sum\limits_{v ∈ \text{Values}(A)} [\cfrac{|D_v|}{|D|} \text{Entropy}(D_v)]$
> >
> > 1. $\text{Gain}(D, A)$：应用属性$A$于数据集$D$，导致的$D$期望熵的减少程度
> > 2. $\text{Values}(A)$：属性A 所有可能取值的集合
> > 3. $D_v=\{d|d∈D,d的属性A的值为v\}$
>
> ## 7.4. 决策树克服过学习/欠学习的方法
>
> > **1️⃣**阈值方法：
> >
> > 1. 决策树训练前设定阈值作为终止学习条件(比如信息增益小于某值)
> > 2. 学习过程中一个结点满足条件就终止学习
> >
> > **2️⃣**测试集方法：在学习过程的曲线上择优
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231126212005974.png" alt="image-20231126212005974" style="zoom: 50%;" /> 
> >
> > **3️⃣**规则后修剪
> >
> > 1. 先允许过度拟合
> > 2. 从根到叶的路径就是一条规则，把得到的决策树转化为规则集合
> > 3. 对每条规则，若删除这个规则的一个前件不会降低规则的分类精度，就删除此前件
> > 4. 按修建后规则的精度对所有规则排序，按照此顺序来对规则分类
> >
> > ==这几个一定要清楚==

# 8. 神经网络

> ==关注BP算法，反馈网络(真的吗)==，这个具体知识点，在知识总结里
>
> 关于以下几点
>
> ==过学习/欠学习问题，网络设计(输入/输出/隐藏层几个)==
>
> **1️⃣**隐藏层数目：过多容易过拟合，一般取样本数/10
>
> **2️⃣**过学习：对训练数据学习得太好，以至于学到了训练数据中的噪声特征，导致在新数据上表现不佳
>
> **3️⃣**欠学习：在训练数据上的表现就不佳，无法捕捉数据的复杂性
>
> 其他的看知识点总结吧

