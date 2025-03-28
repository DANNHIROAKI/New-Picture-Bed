# 知识点

# 1. 绪论

> ## 1.1. 什么是人工智能
>
> > ### 1.1.1. 什么是智能==(记下来)==
> >
> > > **1️⃣**思维理论：知识来自于思维，智能的核心是思考和推理能力
> > >
> > > **2️⃣**知识阈值理论：智能取决于知识量及其一般化程度
> > >
> > > **3️⃣** **进化理论**：有感知外界，适应环境的能力
> > >
> > > 🌱总结：**智能=智力+知识**
> > >
> > > 🔥补充：智能的具体特征——感知外界+记忆思维+**学习自适应**+行为能力
> >
> > ### 1.1.2. 智能的特征==(记下来)==
> >
> > > **1️⃣**具有感知能力
> > >
> > > **2️⃣**具有记忆和思维能力
> > >
> > > **3️⃣**能够学习和自适应
> > >
> > > **4️⃣**具有行为能力
> >
> > ### 1.1.2. 图灵&深蓝
> >
> > > **1️⃣**图灵测试：由能否在对话中欺骗人类，判断机器是否有人类智能；侧重知识阈值理论
> > >
> > > **2️⃣**图灵的梦想：创造能够模仿人类智能的机器或计算机系统
> > >
> > > **3️⃣**深蓝：IBM的超级计算机，首次击败国际象棋世界冠军；侧重思维理论
> > >
> > > **4️⃣**图灵的梦想更难实现
>
> ## 1.2. AI发展史及重要事件
>
> > **1️⃣**孕育(1956前)：
> >
> > **2️⃣**形成(1956-1969)：
> >
> > 1. 1956：麦卡锡提出人工智能一词
> > 2. 1965：第一个专家系统诞生
> >
> > **3️⃣**发展(1970后)：
> >
> > 1. 1982：Hopfield网络提出
> > 2. 1986：BP网络提出
> > 3. 1995：SVM提出
> > 4. 2006：深度学习提出
> > 5. 2016：Alpha Go
>
> ## 1.3. AI研究内容
>
> > **1️⃣**机器感知：模式识别
> >
> > **2️⃣**机器思维：知识表示/推理/搜索，神经网络
> >
> > **3️⃣**机器学习
> >
> > **4️⃣**机器行为：人工智能
>
> ## 1.4. 人工智能的应用
>
> > **1️⃣**专家系统，知识工程，机器学习，模式识别，NLP，机器定理证明
> >
> > **2️⃣**博弈论，智能决策，人工神经网络，机器人

# 2. 知识工程

> ## 2.1. 知识
>
> > **1️⃣**知识定义：有关信息关联在一起所形成的信息结构
> >
> > **2️⃣**==知识特性：相对正确，不确定性(随机/模糊/不完全/经验)，可表示可利用==
> >
> > **3️⃣**知识分类：常识/领域，事实/规则(因果)，确定/不确定
>
> ## 2.2. 知识的表示
>
> > ### 2.2.1. 概述
> >
> > > **1️⃣**什么是知识表示：用计算机可以接收的数据结构，来描述知识
> > >
> > > **2️⃣**知识的表示方法：==符号表示(逻辑性知识，有命题逻辑+谓词逻辑)==，神经网络表示
> >
> > ### 2.2.2. 符号逻辑表示：命题逻辑
> >
> > > **1️⃣**命题：有真假意义的语句，用大写字母表示
> > >
> > > **2️⃣**缺陷：无法描述客观事物的结构和逻辑关系；好比只能描述AB的对错，不能描述AB之间的关系
> >
> > ### 2.2.3. 符号逻辑表示：一阶谓词逻辑
> >
> > > 接近自然语言，容易在计算机上实现
> > >
> > > **1️⃣**谓词结构：谓词名+个体，如Teacher(Wang)，一般形式为$P(x_1,x_2...x_3)$
> > >
> > > 1. 个体：可以是变量，常量，函数(谓词)；变量个体必有定义域
> > > 2. 谓词名：含义由人指定；可以说Teacher(wang)是王老师，也可以说是王老湿
> > >
> > > **2️⃣**谓词有真假，又是需要确定变元才可判断真假
> > >
> > > **3️⃣**谓词连接词：非¬，析取∨，合取∧，蕴含→，等价$\Leftrightarrow$
> > >
> > > | P    | Q    | ¬P   | P∨Q  | P∧Q  | P→Q  | P$\Leftrightarrow$Q |
> > > | ---- | ---- | ---- | ---- | ---- | ---- | ------------------- |
> > > | T    | T    | F    | T    | T    | T    | T                   |
> > > | T    | F    | F    | T    | F    | F    | F                   |
> > > | F    | T    | T    | T    | F    | T    | F                   |
> > > | F    | F    | T    | F    | F    | T    | T                   |
> > >
> > > 有关蕴含关系
> > >
> > > 1. 蕴含实质上是一种因果，P→Q中P是前提Q是结论，P→Q为真表示因果关系成立
> > >
> > > 2. P→Q中P作为前提为F时，则认为P→Q=T，前提都不成立，姑且就认为你是对的
> > >
> > >    ==只要不是错误前提推出真确结论都认为蕴含关系成立==
> > >
> > > **4️⃣**谓词公式：按照以下规则得到的合式公式
> > >
> > > 1. 单个谓词(原子公式)
> > > 2. 单个合式公式的否定
> > > 3. 两个合式公式的∨，∧，→，$\Leftrightarrow$
> > > 4. 以上规则的有限步骤
> > >
> > > **5️⃣**谓词等价式&永真蕴含式(不用刻意记)
> > >
> > > ```)
> > > 等价式(等价号左右同真假)，取几个特殊的：
> > > (PvQ)vR   ⇔    Pv(QvR)
> > > Pv(Q^R)   ⇔    (PvQ)^(PvR)
> > > ¬(PvQ)    ⇔    ¬P^¬Q
> > > P→Q       ⇔    ¬PvQ
> > > P↔Q       ⇔    (P→Q)^(Q→P)
> > > 
> > > 永真蕴含式(=>表示左边为真右边也为真)，取几个特殊的：
> > > ¬P,PvQ    =>    Q(二者同时为真，则Q为真)
> > > P,P→Q     =>    Q
> > > ¬Q,P→Q    =>    ¬P
> > > P→Q,Q→R   =>    P→R
> > > (∃x)P(x)  =>    P(y)
> > > (∀x)P(x)  =>    P(a)
> > > ```
> > >
> > > **6️⃣**谓词推理规则
> > >
> > > 1. P规则：任何推理步可引入前提
> > >
> > > 2. T规则：前步骤的公式永真蕴含S，则S可引入为该步前提
> > >
> > > 3. CP规则：R+前提=>S，则前提=>R→S
> > >
> > > 4. 反证法：P=>Q，当且仅当P^¬Q不可满足
> > >
> > > **7️⃣**基于谓词的知识表示
> > >
> > > 1. 先定义谓词指出其特殊含义
> > > 2. 用连接词连接谓词：用∨∧连接表示事实，有→连接表示规则
> >
> > ### 2.2.4. 产生式表示法：模仿人的思考
> >
> > > **1️⃣**事实的表示
> > >
> > > 1. 三元组：(A, age, 40)A年龄为40
> > > 2. 四元组：用于表示<mark>不确定性</mark>，(gay, A, B, 0.8)A和B有80%可能是gay
> > >
> > > **2️⃣**规则表示：如果P则Q，P→Q(置信度为p)
> > >
> > > **3️⃣**产生式系统：综合数据库(事实库)$\xleftrightarrow{控制系统}$产生式规则库
> > >
> > > 1. 规则库：描述相应领域知识的产生式集合
> > > 2. 事实库：已知事实+推导出事实
> > > 3. 控制系统：负责规则匹配与执行管理
> > >
> > > **4️⃣**产生式系统求解步骤
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231120210147522.png" alt="image-20231120210147522" style="zoom:50%;" /> 
> > >
> > > **5️⃣**示例
> > >
> > > ```txt
> > > 数据库：{有毛,吃草,黑条纹}
> > > 规则库: R1：有毛→哺乳类
> > >        R3：哺乳类∧吃肉→食肉类
> > >        R4：哺乳类∧吃草→有蹄类
> > >        R5：食肉类∧黄褐色∧有斑点→猎狗
> > >        R6：有蹄类∧黑条纹→斑马                            
> > > 1.匹配到R1->更新事实库：{哺乳动物,有毛,吃草,黑条纹}
> > > 2.匹配到R4->更新事实库：{哺乳动物,有毛,吃草,黑条纹,有蹄类}
> > > 3.匹配到R8->更新事实库：{哺乳动物,有毛,吃草,黑条纹,有蹄类,斑马}
> > > ```

# 3. 确定性推理

> ## 3.1. 概述
>
> > ### 3.2.1. 推理机&推理种类
> >
> > > **1️⃣**推理机：程序实现的推理
> > >
> > > **2️⃣**推理种类：
> > >
> > > 1. 演绎推理(如几何证明)，归纳推理(个体到一般)，默认推理(假设某个缺省条件成立)
> > > 2. 确定性推理(基于确切知识)，不确定性推理(处理含糊信息)
> > > 3. 单调推理(推出结论单调增加)，非单调推理(比如缺省推理)
> >
> > ### 3.2.2. 推理方向
> >
> > > **1️⃣**正向推理(数据驱动)：
> > >
> > > 1. 用户已知事实+知识库(KB)中的适用知识$\to$可用知识集(KS)
> > > 2. 从KS选一个知识进行推理，推理结果存入数据库(DB)作为下次推理的已知事实
> > > 3. 循环往复
> > >
> > > **2️⃣**逆向推理：选一个假设目标→寻找支持假设所需证据，能全找到说明假设成立
> > >
> > > 其他还有：混合/双向推理
>
> ## 3.2. 知识匹配&变量代换
>
> > **0️⃣**知识匹配：判断两个知识模式是否一致，有确定(完全匹配)/不确定之分(有置信度)
> >
> > ### 3.2.1. 变量代换概念
> >
> > > **1️⃣**概念：
> > >
> > > 1. 简单的理解：用==变量表示规则==，用==常量换掉后就变成事实==
> > > 2. 定义：有限集$\theta{}=\{t_1/x_1,...,t_n/x_n\}$中$t_i$为常量(也可能为变量or函数)，$x_i$为==互不相同==的变元，$/$表示左边换掉右边的
> > > 3. 特例：用一组代换去换掉表达式中的变量，例如$F\theta$表示用$\theta$替换去替换$F$中的变量
> > >
> > > **2️⃣**代换规则：
> > >
> > > 1. $t_i\neq{x_i}$，两个相同还替换啥
> > > 2. $\forall{}t_j$中不含$\forall{}{x_i}$，例如$\{g(y)/\mathbf{x},f(\mathbf{x})/y\}$就不是代换
> >
> > ### 3.2.2. 代换的复合
> >
> > > **1️⃣**如何得到代换复合：对于$\theta{}=\{t_1/x_1,...,t_n/x_n\}$和$\lambda{}=\{u_1/y_1,...,u_m/y_m\}$
> > >
> > > 1. 令$\theta{}\lambda{}=\{t_1\lambda{}/x_1,...,t_n\lambda{}/x_n\}$，然后合并为$(\theta{}\circ{}\lambda{})_{RAW}=\{\theta{}\lambda{},\lambda\}$
> > > 2. 删除$\theta{}\lambda{}$中满足$t_i\lambda{}=x_i$的项，删除$\lambda$中满足$y_i\in{}\{x_1,...x_n\}$的
> > > 3. 完成上述操作后$(\theta{}\circ{}\lambda{})_{RAW}\xrightarrow{形成复合}\theta{}\circ{}\lambda{}$
> > >
> > > **2️⃣**复合的性质：$F(\theta{}\circ{}\lambda)=(F\theta{})\lambda\neq(F\lambda)\theta$
> > >
> > > **3️⃣**复合的示例：$\theta = \{f(y)/x, z/y\}$，$\lambda = \{a/x, b/y, y/z\}$
> > >
> > > 1. $\theta \circ \lambda =\{\{\theta{}\lambda{}\},\{\lambda\}\}= \{f(y)\lambda/x, z\lambda/y, a/x, b/y, y/z\}$
> > > 2. ==$f(y)\lambda$表示在$f(y)$中运用$\lambda$替换得到$f(b)$==，$z\lambda$表示在$z$中运用$\lambda$替换得到$y$
> > > 3. 化简为$\{f(b)/x, \cancel{y/y}, \cancel{a/x}, \cancel{b/y}, y/z\} = \{f(b)/x,y/z\}$
> >
> > ### 3.2.3. 公式的合一(特殊的代换)
> >
> > > **1️⃣**定义：公式集$F=\{F_1,...F_n\}$存在代换$\theta$使$F_1\theta{}=F_2\theta{}=...=F_n\theta{}$，则
> > >
> > > 1. $\theta$为公式$F$的一个合一，公式集合一==一般不唯一==
> > > 2. 公式集$F=\{F_1,...F_n\}$可合一
> > >
> > > **2️⃣**最一般合一(唯一)：$\sigma$为公式集$F$最一般合一$\iff$其它合一都可表示为$\theta{}_i=\sigma{}\circ{}\lambda{}_i$($\lambda{}_i$为相应代换)
> > >
> > > **3️⃣**最一般合一求法：发现差异，消除差异
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231121011535400.png" alt="image-20231121011535400" style="zoom:50%;" /> 
> > >
> > > **4️⃣**示例
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231121013536569.png" alt="image-20231121013536569" style="zoom:46%;" /> 
>
> ## 3.3. 归结演绎推理
>
> > :palm_tree:定理证明：证明P→Q(等价于¬PvQ)的永真，或者P^¬Q永假
> >
> > ### 3.3.1. 海伯论理论
> >
> > > **1️⃣**文字，子句，子句合取范式，子句集
> > >
> > > 1. 文字：==原子谓词公式==及其否定，如$¬P(x,f(x))，Q(x,g(x))$
> > > 2. 子句：任何==文字的析取==，如$C=¬P(x,f(x)) \lor{} Q(x,g(x))$，不包含任何文字的子句为空子句
> > > 3. 子句合取范式：$C_1\land{}C_2\land{}...\land{}C_n$
> > > 4. 子句集：$S=\{C_1,C_2,...,C_n\}$，仍和谓词公式都可化成子句集，途径是**等价**和**推理**
> > >
> > > **2️⃣**谓词公式化为子句集：==以$(\forall{x})[(\forall{y})P(x,y)→\neg{}(\forall{y})(Q(x,y)→R(x,y))]$为例==
> > >
> > > 1. 消去→和↔(用等价关系)
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231121145429852.png" alt="image-20231121145429852" style="zoom:33%;" /> 
> > >
> > > 2. 让¬紧靠谓词(用等价关系)
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231121150247762.png" alt="image-20231121150247762" style="zoom:34%;" /> 
> > >
> > > 3. 重命名变元，使不同量词约束的变元名字不同
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231121150931588.png" alt="image-20231121150931588" style="zoom:36.5%;" /> 
> > >
> > > 4. 消去存在量词
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231121152331150.png" alt="image-20231121152331150" style="zoom:43%;" /> 
> > >
> > > 5. 全称量词移到公式左边，示例中$\forall$已经在左边
> > >
> > > 6. 将公式化为Skolem标准形(用等价关系)
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231121153517775.png" alt="image-20231121153517775" style="zoom: 35%;" /> 
> > >
> > > 7. 消去全称量词，消去合取词，就得到了子句集
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231121154026611.png" alt="image-20231121154026611" style="zoom:42%;" /> 
> > >
> > > **3️⃣**子句集的意义
> > >
> > > 1. 不可满足性：$S=\{S_1,S_2,..,S_n\}$中，任何论域的任何解释，==都不可能让所有$S_i$都为真==
> > > 2. 定理：谓词公式$F$不可满足$\iff{}$其子句集$S$不可满足 
> > > 3. ==$P\to{}Q\iff{}\neg{}P\lor{Q}永真\iff{}P\land{}\neg{Q}永假\iff{}P\land{}\neg{Q}的子句集S不可满足$==
> > >
> > > **4️⃣**Herbrand理论
> > >
> > > 1. 判断子句集不可满足，需要对所有可能论域上的所有解释进行判定
> > > 2. Herbrand理论中，构造了Herbrand域，代替以上一句话中的所有域
> >
> > ### 3.3.2. 鲁滨逊归结原理
> >
> > > 原理概述：子句集$S$中有空子句必不可满足，==$S$中的子句进行归结后能得到空子句则也不可满足==
> > >
> > > #### 3.3.2.1. 命题逻辑中的归结原理
> > >
> > > > **1️⃣**互补文字：也就是$P$和$\neg{P}$，$P$可为谓词公式or命题
> > > >
> > > > **2️⃣**归结操作：归结式&亲本子句
> > > >
> > > > 1. $C_1,C_2$是子句集中的任意两个子句，$L_i$是$C_i$的文字
> > > > 2. 如果$L_1L_2$互补，则分别在$C_1C_2$中删除二者
> > > > 3. 将现有$C_1C_2$余下部分析取，得到$C_{12}$
> > > > 4. $C_{12}$是$C_1C_2$的归结式，$C_1C_2$是$C_{12}$亲本子句
> > > >
> > > > 例如$C_1=\neg{P\lor{Q}},C_2=\neg{Q}\lor{R}\iff{}C_{12}=\neg{P}\lor{R}$
> > > >
> > > > **3️⃣**定理：$C_1,C_2$是子句集$S$中的子句，则归结式$C_{12}$是$C_1,C_2$的逻辑结论，即$C_1\land{}C_2\Rightarrow{}C_{12}$
> > > >
> > > > 1. 推论1：$S\xrightarrow{C_{12}代替C_1,C_2}S_1$，则$S_1不可满足性\Rightarrow{}S不可满足性$
> > > > 2. 推论2：$S\xrightarrow{C_{12}加入S}S_2$，则$S_2不可满足性\iff{}S不可满足性$
> > > >
> > > > ➕证明子句集$S$不可满足：
> > > >
> > > > 1. 对句集进行归结子
> > > > 2. 把归结式进行加入/替换，得到$S_2$或$S_1$
> > > > 3. 证明二者之一不可满足，即可证明$S$不可满足(若是空集则$S$直接不可满足)
> > > >
> > > > **4️⃣**命题逻辑归结原理的完备性：子句集不可满足$\iff$必有一个从$S$到空子句的归结演绎
> > >
> > > #### 3.3.2.2. 谓词逻辑的归结原理
> > >
> > > > **1️⃣**规则：不能直接消去两个子句中互补的文字，需要先用最一般合一对变元代换，再归结
> > > >
> > > > **2️⃣**示例
> > > >
> > > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231121234844376.png" alt="image-20231121234844376" style="zoom:33%;" /> 
> >
> > ### 3.3.3. 归结反演：证明$Q$为$P_1,P_2，...,P_n$逻辑结论
> >
> > > **1️⃣**以下命题等价
> > >
> > > 1. $Q$是$P_1,P_2，...,P_n$的逻辑结论，即<mark>$(P_1,P_2，...,P_n)\to{}Q$</mark>
> > > 2. $\neg{}(P_1\land{}P_2\land{},...,\land{}P_n)\lor{}Q$永真
> > > 3. $(P_1\land{}P_2\land{},...,\land{}P_n)\land{}\neg{}Q$不可满足
> > > 4. <mark>$(P_1\land{}P_2\land{},...,\land{}P_n)\land{}\neg{}Q$</mark>的子句集不可满足
> > >
> > > **2️⃣**归结反演步骤：$F=\{P_1,P_2，...,P_n\}$为已知前提公式集，$Q$为目标公式(结论)
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231122001740141.png" alt="image-20231122001740141" style="zoom:50%;" /> 
> > >
> > > ➕示例
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231122003947821.png" alt="image-20231122003947821" style="zoom:50%;" /> 
> >
> > ### 3.3.4. 归结策略(看懂这个例子)
> >
> > > **1️⃣**归结的一般过程
> > >
> > > 1. $S=\{C_1,C_2,...,C_n\}$中任意两两归结，得到第一级归结式$S_1$
> > > 2. 把$S,S_1$任意子句两两归结，得到二级归结式$S_2$
> > > 3. 把$S,S_1,S_2$任意子句两两归结，得到三级归结式$S_3$
> > > 4. 循环，直到出现空子句/不能继续归结，示例如下
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231122114814401.png" alt="image-20231122114814401" style="zoom:40%;" /> 
> > >
> > > **2️⃣**删除策略：删除无用子句来缩小归结范围
> > >
> > > 1. 纯文字：
> > >
> > >    - 文字$L$在子句集中但$\neg{}L$不在，那么$L$就是纯文字
> > >    - 如$\{P,\neg{R},\neg{P\lor{Q}},Q\lor{W}\}$的$W$
> > >
> > >    - 纯文字所在子句删掉
> > >
> > > 2. 重言式：
> > >
> > >    - 子句包含互补文字时，该子句就为重言式子句
> > >    - 如$P\lor{}\neg{P}\lor{}Q$
> > >
> > >    - 重言式子句删掉
> > >
> > > 3. 包孕式：
> > >
> > >    - $\sigma$代换$C_1\sigma{}\subseteq{C_2}$，则$C_1$包孕于$C_2$($C_2$包孕了$C_1$)
> > >    - 如$C_1=P(x),C_2=P(y)\lor{}Q(z)$
> > >
> > >    - $C_1$包孕于$C_2$时删除$C_2$
> > >
> > > **3️⃣**支持集的限制策略：
> > >
> > > 1. 支持集：<mark>目标公式的否定得到的子句</mark>+其后裔(其与别的子句归结的结果)
> > > 2. 支持集策略：每次归结，亲本子句中必须有一方属于支持集
> > > 3. 支持集策略完备
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231122134944202.png" alt="image-20231122134944202" style="zoom: 33%;" /> 
> > >
> > > **4️⃣**线性输入的限制策略：归结的两个子句中，必须至少有一个是初始子句集中的子句
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231122153225195.png" alt="image-20231122153225195" style="zoom:33%;" /> 
> > >
> > > **5️⃣**线性输入策略$\xrightarrow{放宽}$祖先过滤策略：要求$C_1C_2$归结时
> > >
> > > 1. 二者至少之一是初始子句集中的子句
> > > 2. 二者中一个是另外一个的祖先
> >
> > ### 3.3.5. 用归结原理求取问题(这里的例子可以不用看)
> >
> > > **1️⃣**步骤概览
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231122161305735.png" alt="image-20231122161305735" style="zoom:58%;" />
> > > **2️⃣**示例：A,B,C 三人，$T(x)$表示$x$说真，证明$\neg{}T(A)$
> > >
> > > 1. 有人全说谎，有人全说真
> > > 2. A说：BC说谎
> > > 3. B说：AC说谎
> > > 4. C说：AB至少一人说谎
> > >
> > > :dagger:前提化为子句集
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231122164643111.png" alt="image-20231122164643111" style="zoom:25%;" /> 
> > >
> > > 🎅整理，简化子句集
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231122170334006.png" alt="image-20231122170334006" style="zoom: 43%;" /> 
> > >
> > > 😍开始归结
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231122170829291.png" alt="image-20231122170829291" style="zoom:50%;" /> 
> > >
> > > 🥰推广：题目改为，判断谁是说真话，引入Answer谓词
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231122172621503.png" alt="image-20231122172621503" style="zoom: 40%;" /> 

# 4. 不确定性推理

> ## 4.1. 概念
>
> > **1️⃣**不确定性推理：不确定<mark>事实(证据)</mark>$\xrightarrow{不确定的规则(知识)}$不确定结论
> >
> > 1. <mark>知识静态强度</mark>：知识不确定性，用数值表示，由专家系统给出(所以静态)
> > 2. <mark>知识动态强度</mark>：证据不确定性，用数值表示
> >
> > **2️⃣**不确定性匹配算法
> >
> > 1. 匹配度：证据和知识的前提(已知事实)的相似程度
> > 2. 匹配算法：求这个匹配度
> >
> > **3️⃣**不确定性推理方法的分类
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231122213332471.png" alt="image-20231122213332471" style="zoom:70%;" /> 
>
> ## 4.2. 基本概率方法
>
> > **1️⃣**经典概率方法：`if 前提E then 结论H`
> >
> > 1. 大量统计可得$E$发生后$H$发生的概率，记为<mark>$P(H|E)$</mark>(后验概率)
> >
> > 2. <mark>$P(H|E)$</mark>就是证据$E$下结论$H$的可信度
> >
> > **2️⃣**逆概率方法：前提$E$对应有很多结论$H_1,H_2,...,H_n$
> >
> > ```txt
> > if 前提E then 结论H1,H2,...,Hn
> > =>
> > if 前提E then 结论H1
> > if 前提E then 结论H2
> > .......
> > if 前提E then 结论Hi
> > .......
> > if 前提E then 结论Hn
> > ```
> >
> > $P(H_i | E) = \cfrac{P(H_i) \times P(E | H_i)}{\sum\limits_{j=1}^{n} P(H_j) \times P(E | H_j)}, \quad i = 1,2,...,n$(这个公式不用记)
>
> ## 4.3. 主观Bayes方法
>
> > ### 4.3.1. 不确定性的表示
> >
> > > **1️⃣**知识的不确定性：`IF E THEN (LS,LN) H(P(H))`
> > >
> > > 1. $P(H)$：结论$H$的先验概率，由专家根据经验给出
> > >
> > > 2. ==静态强度==$LS,LN$：由专家给出，==这两个表达式不用记==
> > >
> > >    - $LS=\cfrac{P(E|H)}{P(E|\neg{}H)}\in{[0,\infty{})}$：充分性度量，指出$E$对$H$的支持程度
> > >
> > >    - $LN=\cfrac{P(\neg{}E|H)}{P(\neg{}E|\neg{}H)}=\cfrac{1-P(E|H)}{1-P(E|\neg{}H)}\in{[0,\infty{})}$：必要性度量，指出$\neg{}E$对$H$的支持程度
> > >
> > > **2️⃣**证据的不确定性
> > >
> > > 1. ==动态强度==$P(E|S)$：对于证据$E$，由用户根据观测$S$，给出$P(E|S)$
> > >
> > > 2. 可信度$C(E|S)$：用来代替$P(E|S)$，其在PROSPECTOR(矿勘专家系统)中
> > >    $$
> > >    \begin{flalign*}
> > >    & P(E|S) = \begin{cases}
> > >    \cfrac{C(E|S) + P(E) \times (5 - C(E|S))}{5}, & 0 \leq C(E|S) \leq 5 \\
> > >    \cfrac{P(E) \times (5 + C(E|S))}{5}, & -5 < C(E|S) < 0
> > >    \end{cases} &
> > >    \end{flalign*}
> > >    $$
> > >
> > >    | $C(E \mid S) \in\{-5, \ldots, 5\}$ | 含义                                           |
> > >    | :--------------------------------: | :--------------------------------------------- |
> > >    |                 -5                 | 观测 $S$ 下证据 $E$ 不存在，即 $P(E \mid S)=0$ |
> > >    |                 5                  | 观测 $S$ 下证据 $E$ 存在， $P(E \mid S)=1$     |
> > >    |                 0                  | $S, E$ 之间无关， $P(E \mid S)=P(E)$           |
> > >
> > > **3️⃣**结论的不确定性：二者关系间==4.3.2.2.**3️⃣**==
> > >
> > > 1. 后验概率$P(H|E)$：给出证据$E$情况下结论$H$为真的概率
> > > 2. 结论可信度$P(H|S)$：考虑了观察$S$后结论$H$为真的概率
> > >
> > > **3️⃣**组合证据的不确定性
> > >
> > > 1. 二维组合
> > >
> > >    |                     |       最大最小法       |           概率法           |            有界法            |
> > >    | :-----------------: | :--------------------: | :------------------------: | :--------------------------: |
> > >    | $T(E 1 \land E 2 )$ | $min\{T(E_1),T(E_2)\}$ |      $T(E_1)*T(E_2)$       | $max{}\{0,T(E_1)+T(E_2)-1\}$ |
> > >    |  $T(E1 \lor E 2 )$  | $max\{T(E_1),T(E_2)\}$ | $T(E1)+T(E2)－T(E1)*T(E2)$ |   $min\{1,T(E_1)+T(E_2)\}$   |
> > >
> > > 2. 多维组合
> > >
> > >    $P(E_1\land{}E_2\land{}...\land{}E_n|S)=min\{P(E_1|S),P(E_2|S),...,P(E_n|S)\}$
> > >
> > >    $P(E_1\lor{}E_2\lor{}...\lor{}E_n|S)=max\{P(E_1|S),P(E_2|S),...,P(E_n|S)\}$
> > >
> > >    $P(\neg{}E|S)=1-P(E|S)$
> >
> > ### 4.3.2. 不确定的传递
> >
> > > #### 4.3.2.1. 前置概念
> > >
> > > > **1️⃣**什么是不确定的传递：证据不确定性+知识不确定性$\xrightarrow[用遗传算法算出算出]{遗传给}$结论的不确定性
> > > >
> > > > **2️⃣**==几率函数$\Theta{}(x)=\cfrac{P(x)}{1-P(x)}$==：$[0,1]\to[0,\infty]$，$\Theta{}(x)P(x)$单调性相同
> > >
> > > #### 4.3.2.2. 三种情况(推导详见课本P92)==计算题重灾区==
> > >
> > > > **1️⃣**证据一定存在$P(E)=P(E|S)=1$(不论观测与否，证据都成立)
> > > >
> > > > - $P(H|E) = \cfrac{LS \times P(H)}{(LS - 1) \times P(H) + 1}$
> > > > - 或者==$\Theta{}(H|E)=LS*\Theta{}(H)$==
> > > >
> > > > **2️⃣**证据一不存在$P(E)=P(E|S)=0$(不论观测与否，证据都不成立)
> > > >
> > > > - $P(H|\neg E) = \cfrac{LN \times P(H)}{(LN - 1) \times P(H) + 1}$
> > > > - 或者==$\Theta{}(H|\neg{}E)=LN*\Theta{}(H)$==
> > > >
> > > > **3️⃣**证据可能存在$0<P(E|S)<1$时：$P(H|S) = P(H|E) * P(E|S) + P(H|\neg E) * P(\neg E|S)$
> > > >
> > > > 1. $P(E|S)=1$时，$P(H|S)=P(H|E)$
> > > >
> > > >    观察后一定能获得证据，所以基于观察的可信度和基于证据的可信度就都一样了
> > > >
> > > > 2. $P(E|S)=0$时，$P(H|S)=P(H|\neg E)$
> > > >
> > > >    观察后一定不能获得证据，所以基于观察的可信度等于基于无证据的可信度
> > > >
> > > > 3. $P(E|S)=P(E)$时，$P(H|S)=P(H)$
> > > >
> > > >    不论观察与否，证据的可信度都不变，所以不论观察与否结论的可信度也不会变
> > > >
> > > > 4. 其余情况：
> > > >    $$
> > > >    \begin{flalign*}
> > > >    P(H|S) = 
> > > >    \begin{cases} 
> > > >    P(H|\neg E) + \cfrac{P(H) - P(H|E)}{P(E)} \times P(E|S), & 0 \leq P(E|S) < P(E) \\
> > > >    P(H) + \cfrac{P(H|E) - P(H)}{1 - P(E)} \times [P(E|S) - P(E)], & P(E) \leq P(E|S) \leq 1
> > > >    \end{cases}
> > > >    \end{flalign*}
> > > >    $$
> > > >
> > > >    - 这是$EH$公式，$C(E/S)$换掉$P(E/S)$便是$CP$公式
> > > >    - 这个公式不用记，只需要知道这是分段差值法求来的
> > > >
> > >
> > > #### 4.3.2.3. LS和LN 的性质
> > >
> > > > **1️⃣**$LS>1$说明证据$E$对结论$H$有利，$LN>1$说明证据$\neg{}E$对结论$H$有利
> > > >
> > > > **2️⃣**==记住==：一般情况取$LS>1\,,LN<1$
> >
> > ### 4.3.3. 结论不确定性的合成
> >
> > > **1️⃣**概念：$A_1,A_2,...,A_n\xrightarrow{合成算法}A$的不确定性到底是多少
> > >
> > > ```txt
> > > if 前提E1 then (LS1,LN1) 结论H(p(H))---->对应观察S1
> > > if 前提E2 then (LS2,LN2) 结论H(p(H))---->对应观察S2
> > > ..........
> > > if 前提En then (LSn,LNn) 结论H(p(H))---->对应观察Sn
> > > ```
> > >
> > > **2️⃣**==$\Theta(H|S_1,S_2,\ldots,S_n) = \cfrac{\Theta(H|S_1)}{\Theta(H)} \times \cfrac{\Theta(H|S_2)}{\Theta(H)} \times \cdots \times \cfrac{\Theta(H|S_n)}{\Theta(H)} \times \Theta(H)$==
> ## 4.4. 可信度模型
>
> > ### 4.4.1. 阈值可信度模型==(计算题重灾区)==
> >
> > >  **1️⃣**表示不确定性：`IF E THEN H(CF(H,E),λ)`
> > >
> > >  1. $CF(E)\in{[0,1]}$：证据$E$的可信度，$CF(E)=1/0$表示证据绝对存在/不存在
> > >  2. $CF(H,E)\in{}[0,1]$：证据$E$为真时$H$的可信度，$CF(H,E)=0,1$对应$P(H|E)=0,1$
> > >  3. $\lambda{}\in{}[0,1]$：阈值，==只有证据$E$可信度$CF(E)\geq\lambda{}$时知识$E$才能被应用==
> > >
> > >  **2️⃣**组合证据的不确定性
> > >
> > >  $CF(E_1\land{}E_2\land{}...\land{}E_n)=min\{CF(E_1),CF(E_2),...,CF(E_n)\}$
> > >
> > >  $CF(E_1\lor{}E_2\lor{}...\lor{}E_n)=max\{CF(E_1),CF(E_2),...,CF(E_n)\}$
> > >
> > >  **3️⃣**==不确定性的传递算法：$CF(E)\geq\lambda{}$时$CF(H)=CF(H,E) × CF(E)$==
> > >
> > >  **4️⃣**结论不确定性合成算法：`IF Ei THEN H(CF(H,Ei),λi)   i=1,2,...,n`
> > >
> > >  1. 前提：$\forall{i}\in\{1,2,...,n\},\,CF(E_i)\geq\lambda{}_i$
> > >
> > >  2. 结论$H$的可信度：
> > >
> > >     | 方法     | $CF(H)$等于                                                  |
> > >     | -------- | ------------------------------------------------------------ |
> > >     | 极大值法 | $\max\{CF_i(H)\}\,\,i=1,2,...,n$                             |
> > >     | 加权和法 | $ \cfrac{\sum\limits_{i=1}^{n} CF(H, E_i) \times CF(E_i)}{\sum\limits_{i=1}^{n} CF(H, E_i)}$ |
> > >     | 有限和法 | $min\{\sum\limits_{i=1}^{n}CF_i(H),\,1\}$                    |
> > >     | 递推法   | $C_0=0$,$C_k = C_{k-1} + (1 - C_{k-1}) \times CF(H, E_k) \times CF(E_k)$ |
> >
> > ### 4.4.2. 加权可信度模型==(计算题重灾区)==
> >
> > >  **1️⃣**知识的不确定性：
> > >
> > >  1. 表示：`IF E1(w1) AND E2(w2) AND ... AND  En(wn)  THEN  H (CF(H,E),λ)`
> > >  2. 阈值$\lambda{}$：由专家给出
> > >  3. 加权因子$\omega_i$：由专家给出，满足$\sum\limits_{i=1}^{n}\omega_i=1$
> > >
> > >  **2️⃣**组合证据的不确定性：$CF(E)=\cfrac{\sum\limits_{i=1}^{n}[\omega_i*CF(E_i)]}{\sum\limits_{i=1}^{n}\omega_i}$
> > >
> > >  **3️⃣**不确定性传递算法：当$CF(E)\geq\lambda{}$时$CF(H)=CF(H,E)\times{}CF(E)$
> >
> > ### 4.4.3. 前件带不确定性的可信度模型==(计算题重灾区)==
> >
> > > **1️⃣**表示：
> > >
> > > 1. 表示1：$cf_i$是子==条件==$E_i$的可信度由专家给出，子<mark>证据</mark>$E_i$的可信度$cf_i'$
> > >
> > >    `IF E1(cf1) AND E2(cf2) AND ... AND  En(cfn)  THEN  H (CF(H,E),λ)`
> > >
> > > 2. 表示2：加上权值
> > >
> > >    `IF E1(cf1,w1) AND E2(cf2,w2) AND ... AND  En(cfn,wn)  THEN  H (CF(H,E),λ)`
> > >
> > > **3️⃣**不确定性匹配
> > >
> > > 1. 无加权因子：$\sum\limits_{i=1}^{n} \max(0, c_{f_i} - c'_{f_i}) \leq \lambda$
> > > 2. 有加权因子：$\sum\limits_{i=1}^{n} \omega_i \times \max(0, c_{f_i} - c'_{f_i}) \leq \lambda$
> > >
> > > **4️⃣**不确定性的传递算法
> > >
> > > 1. 无加权因子：$CF(H) = \left[ \prod\limits_{i=1}^{n} (1 - \max\{0, c_{f_i} - c'_{f_i}\}) \right] \times CF(H,E)$
> > >
> > > 2. ==有加权因子：$CF(H) = \left[ \prod\limits_{i=1}^{n} (1 - \omega_i\max\{0, c_{f_i} - c'_{f_i}\}) \right] \times CF(H,E)$==
>
> ## 4.5. 模糊推理
>
> > ### 4.5.1. 模糊理论模糊集
> >
> > > #### 4.5.1.1. 模糊集
> > >
> > > > **1️⃣**$\mu_A:U\to{[0,1]}$，将任意$u\in{}U$映射到$[0,1]$上的某个函数
> > > >
> > > > 1. 模糊集：$A=\{\mu_A(u),u\in{}U\}$称为$U$上的一个模糊集，
> > > > 2. $\mu_A$：定义在$U$上的模糊集$A$的<mark>隶属函数</mark>
> > > > 3. $\mu_A(u)$：$u$对模糊集$A$的隶属度
> > > >
> > > > **2️⃣**离散论域$U$的模糊集$A$，以下表示中应剔除$\mu_A(u_i)=0$的
> > > >
> > > > 1. 表示1：$A=\{\mu_A(u_1),\mu_A(u_2),...,\mu_A(u_n)\}$
> > > > 2. 表示2：$A=\mu_A(u_1)/u_1+\mu_A(u_2)/u_2+...+\mu_A(u_n)/u_n= \sum\limits_{i=1}^{n} {\mu_A(u_i)}/{u_i}$
> > > > 3. 表示3：$A=\{\mu_A(u_1)/u_1,\mu_A(u_2)/u_2,...,\mu_A(u_n)/u_n\}=\bigcup_{i=1}^{n} {\mu_A(u_i)}/{u_i}$
> > > > 4. 表示4：$A=\{[\mu_A(u_1),u_1],[\mu_A(u_2),u_2],...,[\mu_A(u_n),u_n]\}$
> > > >
> > > > **3️⃣**连续论域$U$的模糊集$A=\int\limits_{u\in{}U} {\mu_A(u)}/{u}$
> > > >
> > > > **4️⃣**$U$上所有模糊集表示为：$\mathcal{F}(U) = \{\mathcal{A} | \mu_{A} : U \to [0,1]\}$或$F(U) = \{\mu_{A} | \mu_{A} : U \to [0,1]\}$
> > >
> > > #### 4.5.1.2. 模糊集的运算
> > >
> > > > **1️⃣**$B$包含于$A$：$A,B\in{}\mathcal{F}(U)\,,\forall{}u\in{}U\,,\mu_B(u)\leq{}\mu_A(u)\to{}B\subseteq{}A$
> > > >
> > > > **2️⃣**==$A$，$B$的交并补：==
> > > >
> > > > 1. $A \cup B : \mu_{A \cup B}(u) = \max\{\mu_{A}(u), \mu_{B}(u)\} = \mu_{A}(u) \vee \mu_{B}(u)$
> > > >
> > > >    例如$\mu_{A}(1)=0.3/1\,,\mu_{B}(u)=0.4/1$则$\mu_{A \cup B}(1)=\text{max}(0.3,0.4)/1=0.4/1$
> > > >
> > > > 2. $A \cap B : \mu_{A \cap B}(u) = \min\{\mu_{A}(u), \mu_{B}(u)\} = \mu_{A}(u) \wedge \mu_{B}(u)$
> > > >
> > > > 3. $\neg A : \mu_{\neg A}(u) = 1 - \mu_{A}(u)$
> > > >
> > >
> > > #### 4.5.1.3. 模糊关系
> > >
> > > > **1️⃣**笛卡尔乘积：了解即可
> > > >
> > > > 1. $A_i$是$U_i$上的模糊集
> > > > 2. $A_1A_2...,A_n$的笛卡尔乘积：
> > > >
> > > > $$
> > > > A_1 \times A_2 \times \cdots \times A_n = \int\limits_{U_1 \times U_2 \times \cdots \times U_n} (\mu_{A_1}(u_1) \wedge \mu_{A_2}(u_2) \wedge \cdots \wedge \mu_{A_n}(u_n)) \, d(u_1, u_2, \ldots, u_n)
> > > > $$
> > > >
> > > > 3. 笛卡尔乘积是$U_1\times{}U_2\times{}...\times{}U_n$上的一个模糊集
> > > >
> > > > **2️⃣**$n$元模糊关系：了解即可
> > > >
> > > > 1. 基于$U_1\times{}U_2\times{}...\times{}U_n$论域
> > > >
> > > > 2. $R = \int\limits_{U_1 \times U_2 \times \cdots \times U_n} \mu_{R}(u_1,u_2,...,u_n)/(u_1, u_2, \ldots, u_n)$
> > > >
> > > > 3. 二元模糊关系：基于$U\times{}V$，当二者都为有限论域时，模糊关系可表示为举证
> > > >
> > > >    例如$U=V=\{u_1,u_2,u_3\}$表示信任关系则有：
> > > >
> > > > $$
> > > > \begin{flalign*}
> > > > R = \begin{bmatrix}
> > > > \mu_R(u_{1}, v_{1}) & \mu_R(u_{1}, v_{2}) & \cdots & \mu_R(u_{1}, v_{n}) \\
> > > > \mu_R(u_{2}, v_{1}) & \mu_R(u_{2}, v_{2}) & \cdots & \mu_R(u_{2}, v_{n}) \\
> > > > \vdots & \vdots & \ddots & \vdots \\
> > > > \mu_R(u_{m}, v_{1}) & \mu_R(u_{m}, v_{2}) & \cdots & \mu_R(u_{m}, v_{n})
> > > > \end{bmatrix}
> > > > \quad \to \quad
> > > > \left[ \begin{array}{ccc}
> > > > 1 & 0.3 & 0.8 \\
> > > > 0.9 & 1 & 0.6 \\
> > > > 0.7 & 0.5 & 1 \\
> > > > \end{array} \right]
> > > > \end{flalign*}
> > > > $$
> > > >
> > > > ==**3️⃣**模糊关系的合成(重点)==
> > > >
> > > > 1. $R_1,R_2$分别是$U\times{}V,V\times{}W$的模糊关系，其合成即为$R_1\circ{}R_2$
> > > >
> > > > 2. $\mu_{R_1 \circ R_2}(u, w) = \bigvee_\limits{v \in V} \{ \mu_{R_1}(u, v) \wedge \mu_{R_2}(v, w) \}$
> > > >
> > > > 3. 示例：
> > > >
> > > >    <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231124004226145.png" alt="image-20231124004226145" style="zoom:44%;" /> 
> > >
> > > #### 4.5.1.4. 模糊逻辑
> > >
> > > > **1️⃣**含义：含有模糊概念、模糊数据的语句
> > > >
> > > > **2️⃣**形式：`x_is_A`，$A$是模糊概念(模糊集)；比如`张三 is 如存在的`
> > >
> > > #### 4.5.1.5. 模糊匹配
> > >
> > > > **1️⃣**形式：`IF (x_is_A) THEN (y_is_B)  `，证据和结论都用模糊命题表示，$AB$是模糊概念
> > > >
> > > > **2️⃣**核心问题：条件的$A$与证据的$A'$不一定完全相同，例如
> > > >
> > > > ```txt
> > > > IF x_is_小(知识) THEN y_is_大(结论)
> > > > x_is_微(证据)
> > > > ```
> > > >
> > > > **3️⃣**==匹配度：计算两个模糊概念(集)之间的相似程度：计算题重灾区==
> > > >
> > > > 1. 海明距离：
> > > >    $$
> > > >    \begin{cases}
> > > >    d(A, B) = \cfrac{1}{n} \times \sum\limits_{i=1}^{n} |\mu_A(u_i) - \mu_B(u_i)|\\d(A, B) = \cfrac{1}{b-a} \int_{a}^{b} |\mu_A(u) - \mu_B(u)| du
> > > >    \end{cases}
> > > >    $$
> > > >    
> > > > 2. 欧几里得距离：
> > > >    $$
> > > >    d(A, B) = \sqrt{\frac{1}{n} \times \sum\limits_{i=1}^{n} (\mu_A(u_i) - \mu_B(u_i))^2}
> > > >    $$
> > > >
> > > > **4️⃣**复合条件的模糊匹配
> > > >
> > > > 1. 条件：`E= x1_is_A1 AND x2_is_A2 AND...AND xn_is_An`
> > > >
> > > > 2. 证据：`x1_is_A1',x2_is_A2',...,xn_is_An'`，匹配度$\delta_{match}(A_i，A_i')\,i=1,2,...,n$
> > > >
> > > > 3. 整个条件与证据的匹配度
> > > >
> > > >    $\delta_{match}(E,E')=min\{\delta_{match}(A_i,A_i')\,i=1,2,...,n\}$
> > > >
> > > >    $\delta_{match}(E,E')=\prod\limits_{i=1}^{n}\delta{}_{match}(A_i,A_i')$
> > >
> > > #### 4.5.1.6. 模糊推理的基本模式
> > >
> > > > **1️⃣**<mark>模糊假言推理</mark>
> > > >
> > > > ```txt
> > > > 知识：IF  x_is_A  THEN  y_is_B
> > > > 证据：    x_is_A'
> > > > 结论：                  y_is_B'
> > > > ```
> > > >
> > > > **2️⃣**<mark>模糊拒取式推理</mark>
> > > >
> > > > ```txt
> > > > 知识：IF  x_is_A  THEN  y_is_B
> > > > 证据：                  y_is_B'
> > > > 结论：    x_is_A'            
> > > > ```
> > > >
> > > > **3️⃣**模糊推理方法(扎德)：由`IF (x_is_A) THEN (y_is_B)  `求出$AB$之间的模糊关系$R$，通过$R$与相应证据合成求出模糊结论
> >
> > ### ==4.5.2. 简单模糊推理==
> >
> > > #### 4.5.2.1. 概述`IF (x_is_A) THEN (y_is_B)  `
> > >
> > > > **1️⃣**什么是模糊推理：知识中只含简单条件，不带可信度因子
> > > >
> > > > **2️⃣**推理方法：
> > > >
> > > > 1. 构造$AB$间的模糊关系$R$
> > > > 2. 如果已知证据`x_is_A'`且$AA'$模糊匹配，则==$B'=A'\circ{}R$==
> > > > 3. 如果已知证据`x_is_B'`且$BB'$模糊匹配，则==$A'=R\circ{}B'$==记住这两个公式，Bar和Arab
> > >
> > > #### 4.5.2.2. 构造模糊关系$R$的方法==：必出一道计算题==
> > >
> > > > **1️⃣**扎德方法：
> > > >
> > > > 1. 概述：$A\in\mathcal{F}(U),B\in\mathcal{F}(V)$表示为$A = \int_{U} {\mu_A(u)}/{u}, \quad B = \int_{V} {\mu_B(u)}/{u}$
> > > >
> > > > $$
> > > > \begin{cases}
> > > > 极大/小规则：R_m = (A \times B) \cup (\neg A \times V) = \int\limits_{U \times V} [\mu_A(u) \land \mu_B(v)] \vee [1 - \mu_A(u)] / (u, v)\\
> > > > \\
> > > > 算术规则：R_a = (\neg A \times V) \oplus (U \times B) = \int\limits_{U \times V} 1 \land [1 - \mu_A(u) + \mu_B(v)]/(u, v)
> > > > \end{cases}
> > > > $$
> > > >
> > > > ==只要求记住$R_m = (A \times B)\cup (\neg A \times V)$==，其他仍和都不要了
> > > >
> > > > 2. 推理方法
> > > >
> > > >    - ==若已知证据`x_is_A'`则$B'_m=A'\circ{R_m}$或$B'_a=A'\circ{R_a}$==
> > > >    - ==若已知证据`y_is_B'`则$A'_m=R_m\circ{}B'$或$A'_a=R_a\circ{}B'$==
> > > >
> > > > 3. 示例
> > > >
> > > >    题目与核心公式
> > > >
> > > >    <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231124145847188.png" alt="image-20231124145847188" style="zoom: 33%;" /> 
> > > >
> > > >    模糊集补集的求法
> > > >
> > > >    <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231124151457307.png" alt="image-20231124151457307" style="zoom:33%;" /> 
> > > >
> > > >    求笛卡尔乘积
> > > >
> > > >    <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231124153246939.png" alt="image-20231124153246939" style="zoom: 28%;" /> 
> > > >
> > > >    求笛卡尔并——(笛卡尔交就是取二者之一的较小者)
> > > >
> > > >    <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/7HNDE%257O4%7DF5FJW2%5B63F2%5DK.png" alt="7HNDE%7O4}F5FJW2[63F2]K" style="zoom:29%;" /> 
> > > >
> > > >    求笛卡尔有界和
> > > >
> > > >    <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/LD%25ZSAD%40%7D%7BTVI0Q6%7BRGQ%7E%7DL.png" alt="LD%ZSAD@}{TVI0Q6{RGQ~}L" style="zoom: 29%;" /> 
> > > >
> > > >    合成模糊关系
> > > >
> > > >    <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231124161240979.png" alt="image-20231124161240979" style="zoom:29%;" /> 
> > > >
> > > > **2️⃣**Mamdani方法：
> > > >
> > > > $R_c = A \times B = \int\limits_{U \times V} \mu_A(u) \land \mu_B(v) /(u,v)$则$B_c' = A' \circ R_c$或者$A_c' = R_c \circ B'$
> > > >
> > > > **3️⃣**Mizumoto方式<font color='red'>**：一律不用看**</font>
> > > >
> > > > 1. $R_s,R_g$
> > > >    $$
> > > >    \begin{cases}
> > > >    R_s = [A \times V \underset{s}{\Rightarrow}
> > > >     U \times B] = \int_{U \times V} \left[ \mu_A(u) \overset{s}{\rightarrow} \mu_B(v) \right] (u, v)
> > > >     \\
> > > >    \mu_A(u) \overset{s}{\rightarrow} \mu_B(v) =
> > > >    \begin{cases}
> > > >    1, & \mu_A(u) \leq \mu_B(v) \\
> > > >    0, & \mu_A(u) > \mu_B(v)
> > > >    \end{cases}
> > > >    \end{cases}
> > > >    $$
> > > >
> > > >    $$
> > > >    \begin{cases}
> > > >    R_g = [A \times V \underset{g}{\Rightarrow} U \times B] = \int_{U \times V} \left[ \mu_A(u) \overset{g}{\rightarrow} \mu_B(v) \right] (u, v)
> > > >    \\
> > > >    \mu_A(u) \overset{g}{\rightarrow} \mu_B(v) =
> > > >    \begin{cases}
> > > >    1, & \mu_A(u) \leq \mu_B(v) \\
> > > >    \mu_B(v), & \mu_A(u) > \mu_B(v)
> > > >    \end{cases}
> > > >                         
> > > >    \end{cases}
> > > >    $$
> > > >
> > > >    - 若已知证据`x_is_A'`则$B'_s=A'\circ{R_s}$或$B'_g=A'\circ{R_g}$
> > > >
> > > >      若已知证据`y_is_B'`则$A'_s=R_s\circ{}B'$或$A'_g=R_g\circ{}B'$
> > > >
> > > >    - 当$U,V$离散时：
> > > >      $$
> > > >      R_{s \text{ or } g} = \left[ \begin{array}{ccc}
> > > >      \mu_A(u_1) \overset{s \text{ or } g}{\rightarrow} \mu_B(v_1) & \cdots & \mu_A(u_1) \overset{s \text{ or } g}{\rightarrow} \mu_B(v_{\text{max}}) \\
> > > >      \vdots & \ddots & \vdots \\
> > > >      \mu_A(u_{\text{max}}) \overset{s \text{ or } g}{\rightarrow} \mu_B(v_1) & \cdots & \mu_A(u_{\text{max}}) \overset{s \text{ or } g}{\rightarrow} \mu_B(v_{\text{max}})
> > > >      \end{array} \right]
> > > >      $$
> > > >
> > > > 2. 其他关系的表示
> > > >    $$
> > > >    \begin{cases}
> > > >     R_{sg} = [(A \times B) \underset{s}{\Rightarrow} (U \times B) ]\cap [(\neg A \times B) \underset{g}{\Rightarrow} (U \times \neg B)]
> > > >    = \int\limits_{U \times V} \left[ \mu_A(u) \overset{s}{\rightarrow} \mu_B(v) \right] \land \left[ (1 - \mu_A(u)) \overset{g}{\rightarrow} (1 - \mu_B(v)) \right] \, /(u, v) \\
> > > >                         
> > > >     R_{gg} = [(A \times B) \underset{s}{\Rightarrow} (U \times B) ]\cap [(\neg A \times B) \underset{g}{\Rightarrow} (U \times \neg B)]
> > > >    = \int\limits_{U \times V} \left[ \mu_A(u) \overset{g}{\rightarrow} \mu_B(v) \right] \land \left[ (1 - \mu_A(u)) \overset{g}{\rightarrow} (1 - \mu_B(v)) \right] \, /(u, v) \\
> > > >                         
> > > >     R_{gs} = [(A \times B) \underset{g}{\Rightarrow} (U \times B) ]\cap [(\neg A \times B) \underset{g}{\Rightarrow} (U \times \neg B)]
> > > >    = \int\limits_{U \times V} \left[ \mu_A(u) \overset{g}{\rightarrow} \mu_B(v) \right] \land \left[ (1 - \mu_A(u)) \overset{s}{\rightarrow} (1 - \mu_B(v)) \right] \, /(u, v) \\
> > > >                         
> > > >     R_{ss} = [(A \times B) \underset{g}{\Rightarrow} (U \times B) ]\cap [(\neg A \times B) \underset{g}{\Rightarrow} (U \times \neg B)]
> > > >    = \int\limits_{U \times V} \left[ \mu_A(u) \overset{s}{\rightarrow} \mu_B(v) \right] \land \left[ (1 - \mu_A(u)) \overset{s}{\rightarrow} (1 - \mu_B(v)) \right] \, /(u, v)
> > > >    \end{cases}
> > > >    $$
> > > >
> > > >
> > > > 3. 示例
> > > >
> > > >    <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231124164310801.png" alt="image-20231124164310801" style="zoom:33%;" /> 
> > > >
> > > >    求解：
> > > >
> > > >    <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231124171212732.png" alt="image-20231124171212732" style="zoom:45%;" /> 
> > >
> > > #### 4.5.2.3. 各种模糊关系的性能分析
> > >
> > > >==:🚗:Pre：这里这几个公式要记下来，会出计算题==，最好能有示例(见我的考点&例题总结)
> > > >$$
> > > >\begin{cases}
> > > >A=\int\limits_U{}\mu_A(u)/u
> > > >\\
> > > >\text{very } A =\int\limits_U{}\mu_A^2(u)/u
> > > >\\
> > > >\text{more/less } A =\int\limits_U{}\mu_A^{0.5}(u)/u
> > > >\end{cases}
> > > >\\
> > > >\begin{cases}
> > > >\text{not } A =\int\limits_U{}1-\mu_A^{}(u)/u
> > > >\\
> > > >\text{not very } A =\int\limits_U{}1-\mu_A^{2}(u)/u
> > > >\\
> > > >\text{not more/less } A =\int\limits_U{}1-\mu_A^{0.5}(u)/u
> > > >\end{cases}
> > > >$$
> > > >
> > > ><font color='red'>**后面的都不用看了**</font>
> > > >
> > > >**1️⃣**模糊假言推理
> > > >
> > > >```txt
> > > >原则1
> > > >知识: IF  x_is_A        THEN  y_is_B  
> > > >证据:     x_is_A
> > > >结论:                         y_is_B 
> > > >
> > > >原则2
> > > >知识: IF  x_is_A        THEN  y_is_B  
> > > >证据:     x_is_very_A
> > > >结论:                         y_is_B or very_B
> > > >
> > > >原则3
> > > >知识: IF  x_is_A        THEN  y_is_B  
> > > >证据:     x_is_more/less_A
> > > >结论:                         y_is_B or more/less_B
> > > >
> > > >原则4
> > > >知识: IF  x_is_A        THEN  y_is_B  
> > > >证据:     x_is_not_A
> > > >结论:                         y_is_not_B or unknown
> > > >```
> > > >
> > > >**2️⃣**模糊拒取式推理
> > > >
> > > >```txt
> > > >原则5
> > > >知识: IF  x_is_A        THEN  y_is_B  
> > > >证据:                         y_is_not_B
> > > >结论:     x_is_not_A                           
> > > >
> > > >原则6
> > > >知识: IF  x_is_A        THEN  y_is_B  
> > > >证据:                         y_is_not_very_B
> > > >结论:     x_is_not_very_A
> > > >
> > > >原则7
> > > >知识: IF  x_is_A        THEN  y_is_B  
> > > >证据:                         y_is_not_more/less_B
> > > >结论      x_is_not_more/less_B
> > > >
> > > >原则4
> > > >知识: IF  x_is_A        THEN  y_is_B  
> > > >证据:                         y_is_B
> > > >结论      x_is_A or x_is_unknown
> > > >```
> >
> > ### 4.5.3. 模糊三段论推理：
> >
> > > ==这里PPT对应一个例子看一下，因为是原题==
> > >
> > >  ```txt
> > > R1: IF x_is_A THEN y_is_B
> > > R2: IF y_is_B THEN z_is_C
> > > -------------------------
> > > R3: IF x_is_A THEN z_is_C
> > > ```
> > >
> > > **1️⃣**$ABC$的论域是$UVW$，$R(A,B)\,\,R(B,C)\,\,R(A,C)$是由模糊知识得到的模糊关系
> > >
> > > **2️⃣**当满足<mark>$R(A,B)\circ{}R(B,C)=R(A,C)$</mark>时，模糊三段论成立，模糊关系$R$满足三段论
> > >
> > > **3️⃣**满足关系：
> > > $$
> > > \begin{cases}
> > > Y'=X_1'\circ{}R(X_1,X_2)\circ{}R(X_2,X_3)\circ{}...\circ{}R(R_n,Y)
> > > \\
> > > 在此处为：
> > > \begin{cases}
> > > B'=A'\circ{}R(A,B)
> > > \\
> > > C'=B'\circ{}R(B,C)=A'\circ{}R(A,B)\circ{}R(B,C)=A'\circ{}R(A,C)
> > > \end{cases}
> > > \end{cases}
> > > $$
> > >
> >
> > ### 4.5.4. 多维模糊理论
> >
> > > ```txt
> > > 知识：IF x1_is_A1 AND x2_is_A2 AND ... AND xn_is_An THEN y_is_B
> > > 证据：   x1_is_A1'    x2_is_A2'    ...     xn_is_An'
> > > ---------------------------------------------------------------
> > > 结论：                                                   y_is_B'
> > > ```
> > >
> > > $A_i,A_i'\in\mathcal{F}(U_i)\,\,\,B,B'\in\mathcal{F}(V)\,\,\,i=1,2,...,n$
> > >
> > > **1️⃣**扎德方法：前提是$A_i$在相同的论域==，这个要记下来==
> > >
> > > 1. 求出所有$A_i,i=1,2,...,n$和$A_i’,i=1,2,...,n$的交集(最小值)为$A$和$A'$
> > > 2. 求出$A,B$之间的模糊关系为$R(A,B)$
> > > 3. 由$B'=A'\circ{}R(A,B)$，示例如下
> > >
> > > **2️⃣**つかもと方法，==这个也要记下来==
> > >
> > > 1. 先求出每个$R(A_i,B),i=1,2,...,n$
> > > 2. 求出$B_i'=A_i'\circ{}R(A_i,B),i=1,2,...,n$
> > > 3. 最后$B' = B'_1 \cap B'_2 \cap \ldots \cap B'_n$
> >
> > ### 4.5.5. 其他模糊推理
> >
> > > **1️⃣**多重模糊推理：==这个需要了解==
> > >
> > > ```txt
> > > IF x_is_A1 THEN y_is_B1 ELSE
> > > IF x_is_A2 THEN y_is_B2 ELSE
> > > ...........
> > > IF x_is_An THEN y_is_Bn
> > > ```
> > >
> > > 简单情况
> > >
> > > ```\
> > > 知识:IF x_is_A THEN y_is_B ELSE y_is_C
> > > 证据:   x_is_A'
> > > 结论：                          y_is_D
> > > ```
> > >
> > > $A_i\in\mathcal{F}(U_i)\,\,\,B,C,D\in\mathcal{F}(V)$
> > >
> > > <img src="https://s2.loli.net/2023/12/31/pU4zoDTHLfw6Y13.png" alt="image-20231231031301211" style="zoom:50%;" /> 
> > >
> > > 以下公式(多重模糊推理中的模糊关系)不用记，考到了会直接给你，PPT中往下复杂的公式都同理
> > >
> > > <img src="https://s2.loli.net/2023/12/31/yxc8GToA3KfDwPq.png" alt="image-20231231031557716" style="zoom:50%;" /> 
> > >
> > > **2️⃣**带可信度因子的模糊推理
> > >
> > > ```txt
> > > 知识:IF x_is_A THEN y_is_B  CF1
> > > 证据:   x_is_A'             CF2
> > > 结论:               y_is_B' CF
> > > ```
> > >
> > > 可信度计算：
> > > $$
> > > CF = 
> > > \begin{cases}
> > > \delta_{\text{match}}(A, A') \times CF_1 \times CF_2, & \text{for the first case}, \\
> > > \delta_{\text{match}}(A, A') \times \min\{CF_1, CF_2\}, & \text{for the second case}, \\
> > > \delta_{\text{match}}(A, A') \times \max\{0, CF_1 + CF_2 - 1\}, & \text{for the third case}, \\
> > > \min\{\delta_{\text{match}}(A, A'), CF_1, CF_2\}, & \text{for the fourth case}.
> > > \end{cases}
> > > $$
> > >

# 5. 搜索策略：

> ## 5.1. 基本概念
>
> > ### 5.1.1. 什么是搜索
> >
> > > **1️⃣**定义：在知识库中找到可利用的知识，构造一条代价较小的推理路线，从而解决问题
> > >
> > > **2️⃣**分类
> > >
> > > 1. 盲目搜索：按预设控制策略搜索，
> > > 2. 启发式搜索：搜索过程加入与问题有关的启发信息
> > >
> > > :put_litter_in_its_place:问题求解过程可看作是搜索过程
> >
> > ### 5.1.2. 问题及其求解过程的表示：状态空间法
> >
> > > **1️⃣**状态：描述问题求解过程中不同时刻的状态，是一个向量$S_k=(S_{k_0},S_{k_1},...)$
> > >
> > > **2️⃣**算符：使问题从一个状态转为另一个状态的操作，比如产生是系统的产生式规则
> > >
> > > **3️⃣**状态空间：所有可能的状态+一切可用算符
> > >
> > > **4️⃣**求解过程：
> > >
> > > 1. 过程：表示出所有状态，定义一组算符→不断把算符作用给状态→得到目标状态
> > > 2. 解：从初始状态→目标状态，所采用的算符序列。算符最少的话就是最优解
> >
> > ### 5.1.3. AND/OR树表示法：==这个不用看==
> >
> > > **1️⃣**复杂问题的处理
> > >
> > > 1. 分解：把复杂问题分解为若干子问题，子问题再继续分解，形成与或数
> > > 2. 等价变换：把原问题化为若干容易求解的新问题吗，形成与或数
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231125160201461.png" alt="image-20231125160201461" style="zoom: 59%;" /> 
> > >
> > > **2️⃣**基本概念
> > >
> > > 1. 本源问题：不能再分解的，可以直接解的子问题
> > > 2. 端节点：没有子节点的结点
> > > 3. 终止结点：本源问题对应的结点
> > > 4. 可解结点：是<mark>终止结点</mark>or是<mark>或结点且子节点中至少一个可解</mark>or是<mark>与结点且子节点全部可解</mark>
> > > 5. 解树：可解节点构成的子树，是表示问题解决路径的树形结构
>
> ## 5.2. 状态空间的搜索策略
>
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/1700900938496.jpg" alt="1700900938496" style="zoom:45%;" /> 
> >
> > ### 5.2.1. 状态空间搜索概述
> >
> > > **1️⃣**有关数据结构：
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231125163352256.png" alt="image-20231125163352256" style="zoom:38%;" /> 
> > >
> > > 1. OPEN表：存放刚生成的结点，不同策略下结点的排列顺序不同
> > >
> > > 2. CLOSE表：结点扩展前要将其放入CLOSE表
> > >
> > >    PS.结点扩展：用算符作用于节点，生成一系列结点
> > >
> > > **2️⃣**搜索过程：给定初始/最终结点，找到最佳路径，==详见$A*$算法==
> >
> > ### 5.2.2. 广度优先搜索
> >
> > > **1️⃣**概述
> > >
> > > 1. 基本思想：初始节点开始，逐层扩展并考察(是否为目标节点)结点，一层结点全扩展考察完，再扩展下一层的结点
> > > 2. OPEN表中的结点按照进入顺序排列
> > > 3. 特点是：盲目性大，但一定能得到最短解<font color='cornflowerblue'>(反之深度优先不一定得到的是最短解)</font>
> > >
> > > **2️⃣**广度优先搜索过程：注意以下将其子节点放入OPEN表尾部，<mark>尾部改为首部就是**<font color='red'>深度优先</font>**了</mark>
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231125193758728.png" alt="image-20231125193758728" style="zoom:50%;" />  
> > >
> > > **3️⃣**示例：从$S_0\to{}S_g$的重排九宫，==这个例子一定看明白==
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231125193843463.png" alt="image-20231125193843463" style="zoom:59%;" /> 
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231125193937300.png" alt="image-20231125193937300" style="zoom:33%;" /> 
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/1.png" alt="1" style="zoom: 51%;" /> 
> >
> > ### 5.2.3. 有解深度优先算法==(了解深度优先就行了，这个不用看)==
> >
> > > **1️⃣**基本思想：深度优先搜索的时候设置界限$d_m$，搜索达到深度界限后，不管找没找到都换个分支
> > >
> > > **2️⃣**流程：
> > >
> > >  <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/acd.png" alt="acd" style="zoom:33%;" /> 
> > >
> > > **3️⃣**特点：由于深度限制不一定能找到解，所以深度设置也困难，找到的也不一定是最优的
> > >
> > > **4️⃣**改进1：先设置一较小$d_m$，当到达了$d_m$深度还未找到，就加深深度
> > >
> > > **5️⃣**改进2：增加一个表$R$
> > >
> > > 1. 每找到一个目标节点就，把他放进$R$并且设置$d_m$为其路径长
> > > 2. 然后以$d_m$有限度的深度优先搜索
> > > 3. 不断找到目标节点不断更新$d_m$，最后一定得到最优解
> >
> > ### 5.2.4. 代价树的广度/深度优先搜索
> >
> > > ==**(广度优先看，深度优先不看)**==
> > >
> > > **1️⃣**基本概念
> > >
> > > 1. 代价树：边上标有代价的树
> > > 2. 代价：$g(x)$是初始结点$S_0$到结点$x$的代价，$c(x_1,x_2)=g(x_2)-g(x_1)$是两结点之间的代价
> > >
> > > **2️⃣**基本思想：OPEN表中的节点按代价从小到大排序，每次取出代价最小的结点
> > >
> > > **3️⃣**特点：代价树的广度优先一定可以求得最优解
> > >
> > > **4️⃣**过程
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/sx.png" alt="sx" style="zoom: 33%;" /> 
> >
> > ### 5.2.5. 启发式搜索，==一定要看==
> >
> > > **1️⃣**启发性信息：指导搜索，与问题有关的信息
> > >
> > > **2️⃣**估价函数：评估结点重要性，$f(x)=g(x)+h(x)$
> > >
> > > |      | $g(x)$                   | $h(x)$启发函数                 |
> > > | ---- | ------------------------ | ------------------------------ |
> > > | 含义 | 初始节点→x的代价         | x→目标结点，最优路径的预估代价 |
> > > | 特点 | 有利于完备性，不利于效率 | 不利于完备性，有利于效率       |
> > >
> > > ➕重拍九宫中，$h(x)$是$x$的格局与目标节点格局不相同的牌数
> > >
> > > **3️⃣**局部择优搜索==(不看)==：一个结点扩展后，按$f(x)$对子节点计算估计价值，选择价值最小的下一个考察
> > >
> > > 深度优先+代价树深度优先，均为局部择优搜索的特例
> > >
> > > **4️⃣**全局择优搜索：一个结点扩展后，从OPEN表全体中，选一个估价最小的结点下一个考察
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/30750.jpg" alt="30750" style="zoom:33%;" /> 
> >
> > ### 5.2.6. A*算法：这一部分全部PPT都要看
> >
> > > **1️⃣**算法思想：基于一般搜索过程，加上以下限制($S_0$是初始节点)
> > >
> > > 1. OPEN表中节点按估价函数$f(x)=g(x)+h(x)$小至大排序
> > > 2. $g(x)$是$S_0\to{}x$路径的代价，$g^{*}(x)$是$S_0\to{}x$最小代价，$g(x)$可估计$g^{*}(x)$
> > > 3. $h^{*}(x)$是$x\to{}目标节点$最小代价，$h(x)$是$h^{*}(x)$下界
> > >
> > > **2️⃣**算法流程：气泡中为$A^{*}$算法相较于一般搜索过程的具象，==这个一定会出大题==
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/搜索一般过程Axing.png" alt="搜索一般过程Axing" style="zoom:10%;" /> 
> > >
> > > **3️⃣**$A*$算法的性质：
> > >
> > > 1. 可纳性==(这个不用看)==：$A*$算法可纳
> > >
> > >    - 可解状态空间图：初始→目标节点有路径存在
> > >    - 可纳性：搜索算法可在有限步终止，并找到最优秀解
> > >
> > > 2. 最优性：$h(x)$越大启发信息就越多，搜索效率越高
> > >
> > > 3. 单调性==(这个很重要)==：要求启发函数$h(x)$单调搜索代价会更低
> > >
> > >    单调的条件：$h(S_g)=0$，$h(x_i)\leq{}h(x_j)+c(x_i,x_j)$
> > >

# 6. 机器学习==(重点看6.2，6.1扫一遍吧)==

> ## 6.1. 概述
>
> > ### 6.1.1. 什么是ML
> >
> > > **1️⃣**机器学习定义：机器模拟人类学习
> > >
> > > **2️⃣**每次释义：KDD知识发现，DM数据挖掘，ML机器学习
> > >
> > > **3️⃣**机器学习系统基本结构
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231126123729231.png" alt="image-20231126123729231" style="zoom:50%;" /> 
> > >
> > > **4️⃣**学习环节的一般过程：确定模型→收集数据→清洗数据→提取特征→训练→获得知识
> >
> > ### 6.1.2. ML分类
> >
> > > **1️⃣**有无指导：有监督，无监督，强化学习
> > >
> > > **2️⃣**学习方法：类比，返利，解释学习
> > >
> > > **3️⃣**推理策略：演绎，归纳，类比学习
> >
> > ### 6.1.3. ML基本问题
> >
> > > ML：在现有观察基础上求得学习函数$L:S(数据空间)\to{}Z(目标空间)$
> > >
> > > | 问题类型 | 分类                                                         | 预测                                                         | 聚类            | 联想                                                         | 优化                                                         |
> > > | :------- | :----------------------------------------------------------- | :----------------------------------------------------------- | :-------------- | :----------------------------------------------------------- | :----------------------------------------------------------- |
> > > | $L$      | 分类器                                                       | 回归曲线(面)                                                 | 聚类函数        | 自身映射                                                     | $S\to{}Max\{d[F(S)]\}$                                       |
> > > | $Z$      | 已知离散值                                                   | 连续值                                                       | 未知离散值      | $Z=S$                                                        | 数据空间的函数$F(S)$，其度量为$d[F(S)]$                      |
> > > | 训练数据 | $<D,C>,D\subset{}S$                                          | $<D,R>,D\subset{}S$                                          | $D,D\subset{}S$ | $D,D\subset{}S$                                              | $D,D\subset{}S$                                              |
> > > | 图示     | <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231126124830491.png" alt="image-20231126124830491" style="zoom: 33%;" /> | <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231126125819266.png" alt="image-20231126125819266" style="zoom: 33%;" /> | 同分类的        | <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231126130507957.png" alt="image-20231126130507957" style="zoom: 33%;" /> | <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231126130624342.png" alt="image-20231126130624342" style="zoom:33%;" /> |
> > > | 特点     | 监督学习(目标已知)                                           | 已知模型类别求参数                                           | 无监督学习      | \                                                            | \                                                            |
> > > | 目的     | 分为已知类                                                   | 预测                                                         | 分为未知类      | 发现不同数据间的依赖关系                                     | 寻找使$d[F(S)]$到达极值的方法                                |
> > > | 方法     | SVM，Bayes，决策树                                           | 神经网络，回归                                               | 各种聚类方法    | 反馈神经网络                                                 | 遗传算法                                                     |
> >
> > ### 6.1.4. 学习成果的评估
> >
> > > **1️⃣**评估原则：模型泛化能力，算法复杂度
> > >
> > > 1. 鲁棒性：健壮性，系统处理非正常数据的能力
> > > 2. 适应性：对不同数据，模型本身需要做多少人工调整
> > > 3. 简洁性，可解释性：有限选择更简单的假设
> > >
> > > **2️⃣**测试集：
> > >
> > > 1. 保留法：从数据集中，抽取部分为测试集，剩下为训练集
> > > 2. 交叉验证：数据集划分为k个不相交的子集，选一个作测试，其他做训练
> > > 3. 随机法：随机抽一部分数据做训练集，好处是可以无限重复
> > >
> > > **3️⃣**学习有效性指标 
> > >
> > > |         指标         |                             公式                             |                             解释                             |
> > > | :------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
> > > |         误差         | $\operatorname{Error}(T)=\sum\limits_{i=1}^{\|T\|} P_i\left\|E_i-L_i\right\|$ 或者 <br> $\operatorname{Error}(T)=\sum\limits_{i=1}^{\|T\|} P_i\left(\sum\limits_{j=1}^d\left(E_{i j}-L_{i j}\right)^2\right)$ |        理想结果与学习结果的插值的加权平均or方 <br> 差        |
> > > |        正确率        | $\operatorname{Accuracy}(T)=\frac{\left\|T_{\text {Enror }<\varepsilon}\right\|}{\|T\|}$ ，错误率=1-正确率 |                    被正确处理的数据的比例                    |
> > > |    复合指 <br> 标    | $\operatorname{precision}(T)=\frac{\|a\|}{\|a+b\|}, \operatorname{recall}(T)=\frac{\|a\|}{\|a+d\|}$ <br> $\operatorname{Accuracy}(T)=\frac{\|a+c\|}{\|a+b+c+d\|}$ |                            见下图                            |
> > > | $F_\beta$ 度 <br> 量 | $F_\beta(T)=\frac{\left(\beta^2+1\right) \operatorname{precision}(T) \times \operatorname{recall}(T)}{\beta^2 \operatorname{precision}(T)+\operatorname{recall}(T)}$ | 精度和召回率的调和平均数， $\beta$ 为精度较召 <br> 回率权重  |
> > > |    宏平均 <br> 法    | $\operatorname{precision}_{\text {macro }}(T)=\frac{1}{k} \sum\limits_{i=1}^k \operatorname{precision}_i(T)$ <br> $\operatorname{recall}_{\text {macro }}(T)=\frac{1}{k} \sum\limits_{i=1}^k \operatorname{recall}_i(T)$ | 将测试集分为 $\mathrm{k}$ 个目标类别，计算每类精度 <br> 和召回率，再求平均 |
> > > |    微平均 <br> 法    | $\operatorname{precision}_{\text {micro }}(T)=\frac{\sum\limits_{i=1}^k\left\|a_i\right\|}{\sum\limits_{i=1}^k\left\|a_i\right\|+\sum\limits_{i=1}^k\left\|b_i\right\|}$ <br> $\operatorname{recall}_{\text {micro }}(T)=\frac{\sum\limits_{i=1}^k\left\|a_i\right\|}{\sum\limits_{i=1}^k\left\|a_i\right\|+\sum\limits_{i=1}^k\left\|d_i\right\|}$ | 把整个测试集看作单分类问题，一次性计算 <br> 所有个体样本指标的平均值 |
> > >
> > > ➕复合指标的参数含义
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231126163644444.png" alt="image-20231126163644444" style="zoom: 50%;" /> 
> >
>
> ## 6.2. 决策树学习
>
> > ### 6.2.1. 决策树的表示
> >
> > > ```txt
> > > 水果分类
> > > 特征向量：{颜色，尺寸，形状，味道}
> > > 四种属性：颜色(红绿蓝)，尺寸(大中小)，味道(酸甜)，形状(圆细)
> > > 样本集合：一批水果，直到其特征向量+类别
> > > 问题：一个新水果，观测得其特征向量，如何为其分类？
> > > ```
> > >
> > > ==关注这个图==
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231126191336583.png" alt="image-20231126191336583" style="zoom: 50%;" /> 
> > >
> > > **1️⃣**结构：把实例从根排列到叶子来分类
> > >
> > > 1. 叶结点：实例所属分类
> > > 2. 结点：实例对应的一个属性
> > >
> > > **2️⃣**实例分类流程：测试根节点属性→按实例属性向下移动→在新结点重复以上过程
> >
> > ### 6.2.2. ID3算法
> >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/338.jpg" alt="338" style="zoom: 50%;" /> 
> > >
> > > #### 6.2.2.1. 算法思想：自顶向下贪婪遍历可能的决策树空间
> > >
> > > > **1️⃣**核心：<font color='red'>如何选取每个结点测试的属性</font>
> > > >
> > > > **2️⃣**主要步骤
> > > >
> > > > 1. 根节点放什么属性：根据每个实例属性单独分类训练样例的能力，选<font color='red'>分类能力最好</font>的属性为根
> > > > 2. 根属性的每个可能值产生一个分支，顺着分支把样例分配过去
> > > > 3. 重复以上过程
> > >
> > > #### 6.2.2.2. 信息增益：==重点关注==
> > >
> > > > <mark>度量属性区分训练样例的能力</mark>
> > > >
> > > > **1️⃣**样本$D$的熵：$\text{Entropy}(D) = - \sum\limits_{i=1}^{c} P_i \log_2 (P_i)$
> > > >
> > > > 1. 样本来自$c$个不同类别
> > > > 2. $P_i$为$i$类样本所占比例
> > > > 3. 熵越小样本越纯净
> > > >
> > > > **2️⃣**信息增益：$\text{Gain}(D, A) = \text{Entropy}(D) - \sum\limits_{v \in \text{Values}(A)} [\frac{|D_v|}{|D|} \text{Entropy}(D_v)]$
> > > >
> > > > 1. $\text{Gain}(D, A)$：应用属性$A$于数据集$D$，导致的$D$期望熵的减少程度
> > > > 2. $\text{Values}(A)$：属性A 所有可能取值的集合
> > > > 3. $D_v=\{d|d\in{}D,d的属性A的值为v\}$
> > >
> > > #### 6.2.2.3. 决策树结点划分的三种情况
> > >
> > > > **1️⃣**结点的样本同属一类：绝对纯净，不可再分
> > > >
> > > > **2️⃣**结点样本不属同一类
> > > >
> > > > 1. 不可再分情形：划分后，子节点不可能具有更低的平均熵
> > > > 2. 可再分情形：划分后，子节点具有更低的平均熵
> > >
> > > #### 6.2.2.4. 确定叶节点&过学习欠学习：==重点关注，尤其是阈值方法==
> > >
> > > > **1️⃣**是否要极致追求叶节点纯净度？(把所有可划分的结点全部划分到底，成为不可划分子节点)
> > > >
> > > > 1. 划分到底：会过学习(过拟合)
> > > > 2. 不分到底：可能会欠学习
> > > >
> > > > **2️⃣**平衡过学习/欠学习：测试集方法
> > > >
> > > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231126212005974.png" alt="image-20231126212005974" style="zoom: 50%;" /> 
> > > >
> > > > **3️⃣**阈值方法：
> > > >
> > > > 1. 决策树训练前，设定阈值作为终止学习条件(信息增益<阈值，结点样本量/全体样本比例<阈值)
> > > > 2. 学习过程中一个结点满足条件，就终止学习
> > > >
> > > > **4️⃣**叶节点分类：叶结点的样本中，$i$类样本最多，就称之为$i$类叶节点
> > >
> > > #### 6.2.2.5. 决策树的规则后修剪：重点关注四个步骤
> > >
> > > > **1️⃣**先允许过度拟合
> > > >
> > > > **2️⃣**从根到叶的路径就是一条规则，把得到的决策树转化为规则集合
> > > >
> > > > **3️⃣**对每条规则，若删除这个规则的一个前件不会降低规则的分类精度，就删除此前件
> > > >
> > > > **4️⃣**按修建后规则的精度对所有规则排序，按照此顺序来对规则分类
> > >
> > > #### 6.2.2.6. 过拟合问题
> > >
> > > > **1️⃣**假设空间$H$训练集$D$(太长是因为训练集选太小了)
> > > >
> > > > **2️⃣**对假设$h\in{}H$，$\exists{}h'\in{H}$​使得
> > > >
> > > > 1. 训练集$D$上$h$错误率更小
> > > > 2. 全体可能数据上$h'$错误率更小
> > > >
> > > > **3️⃣**那么假设$h$就过拟合了训练集D
> > > >
> > > > ➕信息增益度量问题：信息增益度量会偏向于<mark>有较多可能值</mark>的属性

# 7. 神经网络：

> ==和决策树一样，关注一下过学习/欠学习的知识==
>
> 关注网络的设计，几个输入单元几个输出端元，隐藏层多少个神经元等等==(简答题)==
>
> ## 7.1. 神经网络单元形态
>
> > ### 7.1.1. MP模型
> >
> > > **1️⃣**结构
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231229234940457.png" style="zoom: 50%;" /> 
> > >
> > > | 符号       | 名称                            | 值                                                 |
> > > | ---------- | ------------------------------- | -------------------------------------------------- |
> > > | $w_{ij}$   | 连接权($i,j$两神经元的连接强度) | 可调整                                             |
> > > | $u_i$      | $i$神经元活跃值(状态)           | $U_i = \sum\limits_{k=1}^{n} w_{ki}x_k - \theta_i$ |
> > > | $y_i$      | $i$神经元的输出                 | $y_i=f(u_i)$                                       |
> > > | $\theta_i$ | $i$神经元的阈值                 | 可调整                                             |
> > >
> > > **2️⃣**神经元的输出函数$f(u_i)$
> > >
> > > 1. 采用阶跃函数：$f(u_i)=1\,(u_i>0)$
> > > 2. 采用非线性函数：特点是真负无穷都可饱和
> > >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231230000522795.png" alt="image-20231230000522795" style="zoom:50%;" /> 
> >
> > ### 7.1.2. PDP(并行分布处理)模型
> >
> > > **1️⃣**PDP单元的结构
> > >
> > > 1. 输入：把所有输入加权就和就是其输入
> > > 2. 激活状态：当前状态+输入→激活值
> > > 3. 每个单元都有输出
> > >
> > > **2️⃣**单元之间：PDP有很多单元，按照一定模式连接，系统通过经验修改连接强度学习规则
> > >
> > > **3️⃣**原理：通过连接强度的变化来产生学习
>
> ## 7.2. 神经网络拓扑形态
>
> > **1️⃣**前馈网络：每层直接收上一层输入，如BP
> >
> > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231230001544907.png" alt="image-20231230001544907" style="zoom:50%;" /> 
> >
> > **2️⃣**反馈网络：从输出层到输入层有反馈，如Hopfield
> >
> > <img src="https://s2.loli.net/2023/12/30/jR7DXYwfF6Cnqyh.png" alt="image-20231230001734746" style="zoom: 67%;" /> 
> >
> > **3️⃣**竞争网络：同一层的神经元之间也互相连接，如ART/SOM
> >
> > <img src="https://s2.loli.net/2023/12/30/iplSGrW7dvwUOHD.png" alt="image-20231230002028737" style="zoom: 60%;" /> 
> >
> > **4️⃣**全连接网络：疯子
>
> ## 7.3. 神经网络的学习
>
> > **1️⃣**学习概念：改变神经网络内部表示，使得输入输出向好的方向发展
> >
> > **2️⃣**学习分类：
> >
> > 1. 按如何改变网络结构：改变权值，改变拓扑，二者兼具
> > 2. 按确定性：确定性学习/不确定性学习
> >
> > **3️⃣**学习规则
> >
> > 1. Hebb：两个神经元同时处于兴奋状态，则二者连接权值加强
> > 2. $\delta$：根据期望值和实际值差值，来调整权值，核心在于调整量$\propto$差值
> > 3. W-H：根据期望值和实际值均方差，来调整权值
> > 4. 竞争学习：一层中输出最大的为胜者神经元，只调整胜者神经元有关权值
>
> ## 7.4. BP算法
>
> > ### 7.4.1. 结构
> >
> > > <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231230001544907.png" alt="image-20231230001544907" style="zoom:50%;" /> 
> > >
> > > **1️⃣**输入：所有输入信号和，$\text{net}_j = \sum\limits_{i=1}^{n} w_{ij}x_i = \sum\limits_{i=1}^{n} w_{ji}y_i$
> > >
> > > **2️⃣**输出：可以是$y_j = f(\text{net}_j) = \cfrac{1}{1 + e^{-\text{net}_j}}$，也可以是$\text{tanh}$，必须是处处可微分的
> > >
> > > **2️⃣**拓扑结构：由输入维度，隐藏层数量，输出维度决定
> >
> > ### 7.4.2. 算法过程
> >
> > > **1️⃣**初始化：设置权值为小随机数
> > >
> > > **2️⃣**向前传播：输入一个样本$X_p$，输出$Y_p$
> > >
> > > **3️⃣**向后传播：
> > >
> > > 1. 算出与预期的误差$E_p = \cfrac{1}{2} \sum\limits_{i=1}^{m} (y_{pi} - o_{pi})^2$
> > >
> > > 2. 按极小化误差的方式调整权矩阵
> > >
> > >    注意这样会导致后出现的样本更影响结果，solution是样本输入完后才累计调整一次
> > >
> > > 3. 累计误差$E=E+E_p$最后$E$就是整个样本集误差
> > >
> > > **4️⃣**重复2-3步，直到误差足够小
> >
> > ### 7.4.3. 权值调整：使权值沿误差函数的负梯度方向改变
> >
> > > **1️⃣**上一层第$i$个神经元$\xleftrightarrow{连接}$输出层第$j$个神经元
> > >
> > > $\Delta w_{ij} = \eta \cdot \delta_j \cdot o_i$其中$\delta_j = (t_j - o_j) \cdot f'(net_j)$
> > >
> > > | 参数            | 含义                                   |
> > > | --------------- | -------------------------------------- |
> > > | $\Delta w_{ij}$ | 第$i$个神经元到第$j$个神经元的权值变化 |
> > > | $\eta$          | 学习率，决定权值调整的步长             |
> > > | $\delta_j$      | 输出层第$j$个神经元的误差项            |
> > > | $t_j$           | 输出层第$j$个神经元的目标输出          |
> > > | $o_j$           | 输出层第$j$个神经元的实际输出          |
> > > | $f'(net_j)$     | 输出层第$j$个神经元激活函数的导数      |
> > > | $o_i$           | 上一层（隐藏层）第$i$个神经元的输出    |
> > >
> > > **2️⃣**$k-1$层第$i$个神经元$\xleftrightarrow{连接}$$k$层第$j$个神经元：
> > >
> > > $\Delta w_{ij}^{(k)} = \eta \cdot \delta_j^{(k)} \cdot o_i^{(k-1)}$其中$\delta_j^{(k)} = \left( \sum \delta_m^{(k+1)} w_{jm}^{(k+1)} \right) \cdot f'(net_j^{(k)})$
> > >
> > > | 参数                                   | 含义                                                         |
> > > | -------------------------------------- | ------------------------------------------------------------ |
> > > | $\Delta w_{ij}^{(k)}$                  | $k-1$层第$i$个神经元到$k$层第$j$个神经元的权值变化           |
> > > | $\delta_j^{(k)}$                       | $k$层第$j$个神经元的误差项                                   |
> > > | $\sum \delta_m^{(k+1)} w_{jm}^{(k+1)}$ | $k+1$层的误差项与权值的加权和，表示$k$层第$j$个神经元对后续层误差的贡献 |
> > > | $f'(net_j^{(k)})$                      | $k$层第$j$个神经元激活函数的导数                             |
> > > | $o_i^{(k-1)}$                          | $k-1$层第$i$个神经元的输出                                   |
> >
> > ### 7.4.4. BP算法注意事项
> >
> > > **1️⃣**激活函数：满足非线性/有界/连续/可微 
> > >
> > > **2️⃣**目标向量的选择：向量的每个分量，不超过激活函数的值域
> > >
> > > **3️⃣**权系数的初始化：随机化均匀赋值
> > >
> > > **4️⃣**学习速率$\eta$：一般取0.1，过大会导致不收敛。过小会收敛慢
> > >
> > > **5️⃣**误差函数局部极小值：当训练过程被困在局部极小值时，需要对权值重新初始化，然后重新训练
> > >
> > > **6️⃣**学习曲线：
> > >
> > > <img src="https://s2.loli.net/2023/12/30/dSQmT1sPfiKx7ck.png" alt="image-20231230131913242" style="zoom:50%;" /> 
> > >
> > > 1. 样本划分：训练集(调整权值)，确认集(确定最佳权值)，测试集(评估实际性能)
> > > 2. Epoch：所有样本都参与了一次对权值的更新，就算一次Epoch
> > >
> > > **7️⃣**过拟合：对训练集分类能力强，对测试集却很差，solution是在确认集平均误差函数最小时终止
> > >
> > > **8️⃣**隐藏层数目：过多容易过拟合，一般取样本数/10

# PS. 补充内容

> ## PS.1. 模式识别
>
> > ### PS.1.1. 什么是模式识别(eg.计算机认字)
> >
> > > **1️⃣**含义：对观察到的物理对象识别分类
> > >
> > > **2️⃣**模式：存储于计算机内的有关物理对象的观测信息
> > >
> > > **3️⃣**如何让机器模式识别：采集数据+模式识别算法
> > >
> > > **4️⃣**模式识别的意义：对外界感知和识别是智能的基础
> >
> > ### PS.1.2. 模式识别系统设计：鱼类分类为例
> >
> > > **1️⃣**获取观测量：比如获得+分割图像，得到单个鱼的图像
> > >
> > > **2️⃣**特征提取：
> > >
> > > 1. 如鱼的长短，暗亮等得到<mark>特征向量</mark>
> > > 2. 特征向量所有可能取值得到特征空间
> > > 3. 样本$<x,y>$中$x$表示特征向量$y$表示样本类别，如$y=-1/1$表示A类/B类鱼
> > >
> > > **3️⃣**分类器训练：根据样本集，构造一个判定函数$d(x)$来实现分类
> > >
> > > 1. 设计函数是要求A/B类鱼样本有，$d(x)<0/d(x)>0$，$d_{归一化}(x)=-1/d_{归一化}(x)=1$
> > >
> > > 2. $\{x|d(x)=0\}$为分类面
> > >
> > > 3. 本例中采用线性模型，则如下图，这里分类器复杂度极低(出现欠学习，泛化性低)
> > >
> > >    <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231127014413374.png" alt="image-20231127014413374" style="zoom:50%;" /> 
> > >
> > > **4️⃣**分类器测试：用新样本测试分类器性能，看能不能正确分类得到误差，随后进行优化
> > >
> > > **5️⃣**分类器设计的考虑因素：
> > >
> > > 1. 两类鱼特征向量的真实分布
> > >
> > > 2. 训练集测试集越多越好
> > >
> > > 3. 分类器不能太复杂，不然会过拟合导致泛化程度降低(如下图)，复杂度应该适中且遵从以下原则
> > >
> > >    - Occam razor原则：误差相近选择复杂度低的分类器
> > >    - 用学习曲线方法选择
> > >
> > >    <img src="https://raw.githubusercontent.com/DANNHIROAKI/New-Picture-Bed/main/img/image-20231127015556965.png" alt="image-20231127015556965" style="zoom:50%;" /> 



