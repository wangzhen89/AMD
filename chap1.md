# (PART) 热身 {-}

# 最简单的情况：具有同质误差的完全随机设计结构中的单向处理结构 {#chap1}

> "Randomization is too important to be left to chance." - J. D. Petruccelli

[^treatment]: 本文将 "treatment" 统一翻译为 “处理” ，不论语境。

假设一名实验者想比较几种不同处理[^treatment]的效应，例如不同药物对人们心率的影响或几种不同小麦品种的产量。通常，分析此类实验数据的第一步是使用一种称为单向方差分析模型的统计方法来刻画数据。单向方差分析所基于的模型是统计学领域中最有用的模型之一。许多实验情况只是该模型的特殊情况。其他看起来要复杂得多的模型通常可以被认为是单向模型。本章分为几个部分。在前两部分中，定义了单向模型，并讨论了其参数的估计。在 \@ref(sec1-3) 节和\@ref(sec1-5)节中，提供了处理效应特定线性组合的推断程序。在 \@ref(sec1-7) 节和 \@ref(sec1-9) 节中，我们介绍了开发检验统计量的两种基本方法。这两种方法在本书的其余部分被广泛使用。最后，在 \@ref(sec1-11) 中，我们讨论了使用上述技术可用的计算机分析。使用一个示例来演示每一节中描述的概念和计算。

## 模型定义和假设 {#sec1-1}

[^population]: 本文将 "population" 统一翻译为 “总体”，不论语境。

假设含 $N$ 个实验单元的样本是从可能的实验单元总体[^population]中完全随机地选取的。实验单元被定义为施加处理并独立观察的基本单元。在第 \@ref(chap4) 章和第 \@ref(chap5) 章中可以找到对实验单元更完整的描述。

为了比较 $t$ 种不同处理的效应，将 $N$ 个实验单元的样本随机分为 $t$ 组，使第 $i$ 组中有 $n_i$ 个实验单元，其中 $i=1,2,\dots,t$ 且 $N=\sum_{i=1}^tn_i$. 将实验单元随机分为 $t$ 组应消除任何系统偏差。也就是说，随机化应该确保在施加处理之前，$t$ 组实验单元在性质上是相似的。最后，$t$ 个处理中的每个处理应随机分配到每组实验单元。等价地，可以使用一些随机化装置将实验单元随机分配到 $t$ 个处理组，例如将 $n_1$ 个标签置于表示处理 1 的碗中，将 $n_2$ 个标签置于表示处理 2 的碗中，......，将 $n_t$ 个标签置于表示处理 t 的碗中。混合标签，然后从碗中随机选择标签，以确定分配给每个实验单元的处理。这种在碗中使用标签的过程显然可以使用随机数生成的软件来执行。

设 $y_{ij}$ 表示分配给第 $i$ 种处理的第 $j$ 个实验单元的反应。值 $y_{11},y_{12},\dots,y_{1n_1}$ 可以被认为是来自均值为 $m_1$ 且方差为 $s^2_1$ 的总体的大小为 $n_1$ 的随机样本，值 $y_{21},y_{22},\dots,y_{2n_1}$ 可以被认为是来自均值为 $m_2$ 且方差为 $s^2_2$ 的总体的大小为 $n_2$ 的随机样本，对于 $i = 3,4,\dots,t$ 也是如此。参数 $\mu_i$ 和 $\sigma^2_i$ 表示对实验单元的整个总体施加处理 $i$ 时的总体均值和总体方差。

本章考虑了最简单的情况，即假设各处理之间的方差是同质的或相等的，即 $\sigma^2_1= \sigma^2_2=\cdots = \sigma^2_t$. 也就是说，假设对实验单元施加第 $i$ 个处理可能会影响响应的均值，但不会影响响应的方差。第 \@ref(chap2) 章讨论了等方差假设以及不等方差的方差分析。

一个好的统计分析的基本目标是估计模型的参数并对它们进行推断，推断的方法通常包括检验假设和构造置信区间。

有几种方法可以为上述情况下的数据编写模型。第一个要使用的模型称为 $\mu_i$ 模型或**均值模型** (**means model**)。均值模型是：

$$y_{ij}=\mu_{i}+\varepsilon_{ij}\quad i=1,2,\dots,t,\,j=1,2,\dots,n_{i}$$

其中假设了

```{=tex}
\begin{equation} 
  \varepsilon_{ij}\sim \text{i.i.d.}\,N(0,\sigma^{2})\quad i=1,2,\ldots,t,\,j=1,2,\ldots,n_{i}
  (\#eq:1-1)
\end{equation}
```
记号 $\varepsilon_{ij}\sim \text{i.i.d.}\,N(0,\sigma^{2})$ 贯穿本书。它表示 $\varepsilon_{ij}\,(i=1,2,\dots,t;\,j=1,2,\dots,n_i)$ 独立同分布，每个 $\varepsilon_{ij}$ 的抽样分布是均值为 $0$，方差为 $s^2$ 的正态分布。

## 参数估计 {#sec1-2}

统计分析最重要的方面是得到每个实验单元的误差方差 $\sigma^2$ 的良好估计. 误差方差衡量实验的准确性------- $s^2$ 越小，实验越准确。在任何实验或研究中，如果不了解实验误差方差，就不能做出任何在统计上有效的推论。

在上述情况下，当 $n_i > 1$ 时，第 $i(=1,2,\dots,t)$ 个样本提供 $\sigma^2$ 的估计。从第 $i$ 种处理的数据中获得的 $\sigma^2$ 的估计值为

$$\hat{\sigma}_i^2=\sum_{j=1}^{n_j}\frac{\left(y_{ij}-\overline{y}_{i\cdot}\right)^2}{n_i-1}$$

这是 $\sigma^2$ 的一个无偏估计。其中

$$\overline{y}_{i\cdot}=\frac{\sum_{j=1}^{n_i}y_{ij}}{n_i}$$ 第 $i$ 个处理的 $\sigma^2$ 的估计值是 $\sigma^2_i$ ，基于 $n_i - 1$ 个自由度，且 $(n_i - 1)\hat\sigma_i^2/\sigma^2$ 的抽样分布是具有 $n_i - 1$ 个自由度的卡方分布。

$\sigma^2$ 的这些 $t$ 个独立估计值的加权平均提供了如下情况下 $\sigma^2$ 可能的最佳估计值：其中方差的每个估计值由其相应的自由度加权。$\sigma^2$ 的最佳估计为

$$\hat{\sigma}^2=\sum_{i=1}^t(n_i-1)\hat{\sigma}_i^2\left/\sum_{i=1}^t(n_i-1)\right.$$

出于计算目的，每个方差乘以其权重可以表示为

$$(n_{i}-1)\hat{\sigma}_{i}^{2}=\sum_{i=1}^{t}\left(y_{ij}-\overline{y}_{i\cdot}\right)^{2}=\sum_{i=1}^{t}y_{ij}^{2}-n_{i}\overline{y}_{i\cdot}^{2}=\sum_{i=1}^{t}y_{ij}^{2}-(y_{i.})^{2}/n_{i}=SS_{i}$$ 其中 $y_{i\cdot}=\sum_{j=1}^{n_{i}}y_{ij}$. 则方差的合并估计 (pooled estimate) 为

$$\hat{\sigma}^2=\frac{SS_1+SS_2+\cdots+SS_t}{(n_1-1)+(n_2-1)+\cdots+(n_t-1)}=\frac{\sum_{i=1}^tSS_i}{N-t}$$ 方差 $\sigma^2$ 的合并估计值基于 $N-t$ 个自由度，并且 $(N-t)\hat\sigma^2/\sigma^2$ 的抽样分布是具有 $N-t$ 个自由度的卡方分布，即，$(N-t)\hat\sigma^2/\sigma^2\sim \chi^2_{N-t}$.

每个 $\mu_i$ 的最佳估计是 $\hat{\mu}_{i}=\bar{y}_{i\cdot},\,i=1,2,\dots,t$.

在式 \@ref(eq:1-1) 给出的假设下，$\hat\mu_i$ 的抽样分布是正态分布，具有均值 $\mu_i$ 和方差 $\sigma^2/n_i$. 即，

```{=tex}
\begin{equation} 
\hat{\mu}_i\sim N{\left(\mu_i,\frac{\sigma^2}{n_i}\right)}\quad i=1,2,\ldots,t
(\#eq:1-2)
\end{equation}
```
使用 $\hat\mu_i,\hat\sigma_i^2$ 的抽样分布，那么

```{=tex}
\begin{equation}
t_{i}=\frac{\hat{\mu}_{i}-\mu_{i}}{\sqrt{\hat{\sigma}^{2}/n_{i}}}\sim t_{\Tiny{N-t}}\quad i=1,2,\ldots,t
(\#eq:1-3)
\end{equation}
```
也就是说，$t_i$ 的抽样分布是具有 $N-t$ 个自由度的 $t$ 分布。此外，$\mu_1,\mu_2,\dots,\mu_t$ 和 $\sigma_i^2$ 在统计上是独立的。

## 线性组合的推断：检验与置信区间 {#sec1-3}

本节提供均值模型中参数线性函数的假设检验和置信区间。上一节中的结果可用于检验有关单个 $\mu_i$ 的假设。这些结果也可用于检验有关 $\mu_i$ 线性组合的假设或构造 $\mu_i$ 线性组合的置信区间。

对于涉及几种处理的实验，研究者选择研究中的处理是因为存在需要研究的有趣假设。这些有趣的假设形成了研究的目的。涉及处理的假设很可能涉及方法的特定线性组合。这些线性组合使研究者得以比较不同处理的效应，或者等价地，比较不同处理或总体的平均值。如下所示，其中 $c_1,c_2,\dots,c_t$ 和 $a$ 为某组已知常数：

$$H_{01}{:}\,\sum_{i=1}^{t}c_{i}\mu_{i}=a\,\mathrm{~vs.~}\, H_{a1}{:}\mathrm{~(not~}H_{01}{:})\\
H_{02}{:}\,\mu_{1}=\mu_{2}=\cdots=\mu_{t}\,\mathrm{~vs.~}\,H_{a2}{:}\mathrm{~(not~}H_{02}{:})\\
H_{03}{:}\,\mu_i=\mu_{i^\prime}\text{ for some }i\neq i^{\prime}\mathrm{~vs.~}H_{a3}{:}\mathrm{~(not~}H_{03}{:})$$

对于 $H_{01}$ 给出的线性组合，可以证明：

```{=tex}
\begin{equation}
\frac{\sum_{i=1}^tc_i\hat{\mu}_i-\sum_{i=1}^tc_i\mu_i}{\sqrt{\hat{\sigma}^2\sum_{i=1}^tc_i^2/n_i}}\thicksim t_{\Tiny{N-t}}
(\#eq:1-4)
\end{equation}
```
这个结果可以用来推断形式为 $\sum_{i=1}^tc_i\mu_i$ 的线性组合。由于 $H_{03}$ 中的假设可以写为 $H_{03}:\,\mu_i-\mu_{i^\prime}=0$，它是 $H_{01}$ 的特殊情况，其中 $c_i=1,c_{i^\prime}=-1$，且 $c_k=0$ 如果 $k\ne i$ 或 $i^\prime$. $H_{02}$ 的一个检验在 \@ref(sec1-5) 节给出。$\sum_{i=1}^tc_i\hat{\mu}_i$ 的标准误估计为

```{=tex}
\begin{equation}
\widehat{s.e.}\left(\sum c_i\hat{\mu}_i\right)=\sqrt{\hat{\sigma}^2\sum\frac{c_i^2}{n_i}}
(\#eq:1-5)
\end{equation}
```
为了检验 $H_{01}{:}\,\sum_{i=1}^{t}c_{i}\mu_{i}=a\,\mathrm{~vs.~}\, H_{a1}{:}\mathrm{~(not~}H_{01}{:})$，计算 $t$ 统计量：

```{=tex}
\begin{equation}
t_c=\frac{\sum c_i\hat{{\mu}}_i-a}{\widehat{s.e.}\left(\sum c_i\hat{\mu}_i\right)}
(\#eq:1-6)
\end{equation}
```

[^alpha=100]: 原文为 $\alpha=100\%$，应为笔误。
[^criticalpoint]: 原文为critical point，后同。

如果 $|t_{c}| >t_{\alpha/2,v}$，其中 $v=N-t$，那么 $H_{01}$ 在 $\alpha100\%$[^alpha=100]显著性水平下被拒绝，其中 $t_{\alpha/2,v}$ 是自由度为 $v$ 的 $t$ 分布的上 $\alpha/2$ 分位数[^criticalpoint]。$\sum_{i=1}^{t}c_{i}\mu_{i}$ 的 一个 $(1-\alpha)100\%$ 置信区间为


```{=tex}
\begin{equation}
\sum c_i\hat{\mu}_i\pm t_{\alpha/2,v}\widehat{s.e.}\left(\sum c_i\hat{\mu}_i\right)
(\#eq:1-7)
\end{equation}
```

## 示例：任务和脉搏率 {#sec1-4}

表 \@ref(tab:table1-1) 中的数据来自一项实验，该实验旨在确定六种不同的工作任务如何影响工人的脉搏率。在这项实验中，$78$ 名男性工人被随机分配到六个不同的组，每组 $13$ 名工人。每组工人都接受了执行指定任务的培训。在培训后的选定一天，在工人完成指定任务 $1$ 小时后测量他们的脉搏率。不幸的是，在训练过程中，一些个体退出实验，使得一些组包含少于 $13$ 个个体。记录的数据表示 $20s$ 内的心脏搏动次数，其中有 $N = 68$ 个观测值，和为 $y = 2197$.

<table class="table table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:table1-1)六项任务的脉搏数据和摘要信息</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%201.1.png">
</td>
  </tr>
</tbody>
</table>

对于任务数据，$\sigma^2$ 的最佳估计为

$$\hat{\sigma}^{2}=\sum_{i=1}^{6}SS_{i}\left(N-t\right)=1,916.0761/62=30.9045$$

基于 $62$ 个自由度。$\mu_i$ 的最佳估计分别为 $\hat\mu_1 = 31.923,\hat\mu_2 = 31.083,\hat\mu_3 = 35.800,\hat\mu_4 = 38.000,\hat\mu_5 = 29.500$ 和 $\hat\mu_6 = 28.818$.

为了说明，假设研究者有兴趣回答以下关于任务均值线性组合的问题：

a) 检验 $H_0:\,\mu_3=30\,\mathrm{~vs.~}\,H_a:\,\mu_3\ne30$.
b) 为 $\mu_1$ 构造 $95\%$ 置信区间。
c) 检验 $H_0:\,\mu_4=\mu_5\,\mathrm{~vs.~}\,H_a:\,\mu_4\ne\mu_5$.
d) 检验 $H_0:\,\mu_1=(\mu_2+\mu_3+\mu_4)/3\,\mathrm{~vs.~}\,H_a:\,\mu_1\ne(\mu_2+\mu_3+\mu_4)/3$.
e) 为 $4\mu_1-\mu_3-\mu_4-\mu_5-\mu_6$ 构造 $90\%$ 置信区间。

这些问题可以通过应用本节的结果来回答。

**a的结果**。通过代入式 \@ref(eq:1-6) 获得用于检验 $H_0:\,\mu_3=30$ 的 $t$ 统计量

$$t_{c}=\frac{\hat{\mu}_{3}-30}{\widehat{s.e.}(\hat{\mu}_{3})}=\frac{\hat{\mu}_{3}-30}{\sqrt{\hat{\sigma}^{2}/n_{3}}}=\frac{35.8-30.0}{\sqrt{30.9045/10}}=3.30$$

该 $t$ 值的显著性概率是 $\hat \alpha = \operatorname{Pr}\{|t_c| > 3.30\} = 0.0016$，其中 $\operatorname{Pr}\{|t_c| > 3.30\}$ 是具有 $62$ 个自由度的 $t$ 分布中 $3.30$ 右侧的面积加上 $-3.30$ 左侧的面积。上述 $\hat\alpha$ 值是由电脑计算出来的，但也可由一些特别的手持计算器计算出来。没有电脑或计算器的读者，应比较 $t_c = 3.30$ 与 $t_{\alpha/2,62}$ 来选择 $\alpha$.

**b的结果**。$\mu_1$ 的 $95\%$ 置信区间由下式给出

\begin{aligned}
\hat{\mu}_{1}\pm t_{0.025,62}\widehat{s.e.}(\hat{\mu}_{1})& =31.923\pm2.00\sqrt{30.9045/13}  \\
&=31.923\pm2.00\times1.542
\end{aligned}

因此 $\mu_1$ 的 $95\%$ 置信区间为 $28.839<\mu_1<35.007$，我们有 $95\%$ 的把握认为这个区间包含了未知的真实值 $\mu_1$.

**c的结果**。为检验 $H_0:\,\mu_4=\mu_5$，令 $l_1=\mu_4-\mu_5$，那么 $\hat l_1=\hat \mu_4-\hat\mu_5=38.0-29.5=8.5$。由于 $c_1=c_2=c_3=c_6=0,c_4=1,c_5=-1$，那么

$$\widehat{s.e.}(\hat{l}_{1})=\sqrt{\hat{\sigma}^{2}\sum_{i=1}^{6}c_{i}^{2}/n_{i}}=\sqrt{30.9045{\left(\frac{1}{10}+\frac{1}{12}\right)}}=2.380$$

用于检验 $H_0:\,\mu_4=\mu_5$ 的 $t$ 统计量为

$$t_c=\frac{8.5}{2.380}=3.57$$

该检验的显著性概率为 $\hat\alpha=0.0007$.

**d的结果**。检验 $H_0:\,\mu_1=(\mu_2+\mu_3+\mu_4)/3$ 等价于检验 $H_0:\,\mu_1-\frac{1}{3}\mu_2-\frac{1}{3}\mu_3-\frac{1}{3}\mu_4=0$ 或检验 $H_0:\,3\mu_1-\mu_2-\mu_3-\mu_4=0$. 选择最后一个版本，计算稍微容易一些，并且 $t_c$ 检验统计量的值相对于常数乘数是不变的。

令 $l_2=3\mu_1-\mu_2-\mu_3-\mu_4$，那么

$$\hat{l}_2=3\hat{\mu}_1-\hat{\mu}_2-\hat{\mu}_3-\hat{\mu}_4=3(31.923)-31.083-35.8-38.0=-9.114$$

$\hat l_2$ 的标准误估计为

$$\widehat{s.e.}(\hat{l}_2)=\sqrt{30.9045{\left(\frac9{13}+\frac1{12}+\frac1{10}+\frac1{10}\right)}}=5.491$$

用于检验 $H_0:\,3\mu_1-\mu_2-\mu_3-\mu_4=0$ 的 $t$ 统计量为

$$t_c=\frac{-9.114}{5.491}=-1.66$$

该检验的显著性概率为 $\hat\alpha=0.1020$.

**e的结果**。令 $l_{3}=4\mu_{1}-\mu_{3}-\mu_{4}-\mu_{5}-\mu_{6}$. 那么 $\hat{l}_{3}=-4.426$ 且 $\widehat{s.e.}(\hat{l}_{3})=7.0429$. $l_3$ 的 $90\%$ 置信区间为

$$\hat{l}_{3}\pm t_{0.05,62}\widehat{s.e.}\left(\hat{l}_{3}\right)=-4.426\pm1.671\times7.043=-4.426\pm11.769$$

因此 $90\%$ 置信区间为 $-16.195<4\mu_1-\mu_3-\mu_4-\mu_5-\mu_6<7.343$.

## 几个线性组合的同时检验 {#sec1-5}

在许多情况下，研究人员希望检验关于处理效应或均值的几个**线性组合的同时假设** (**simultaneous hypothesis**)。例如，涉及处理均值的 $k$ 个**线性独立线性组合** (**linearly independent linear combinations**) 的一般假设可以表示为：

\begin{equation}
H_0{:}\quad\begin{aligned}c_{11}\mu_1+c_{12}\mu_2+&\cdots+c_{1t}\mu_t=a_1\\c_{21}\mu_1+c_{22}\mu_2+&\cdots+c_{2t}\mu_t=a_2\\&\vdots\\c_{k1}\mu_1+c_{k2}\mu_2+&\cdots+c_{kt}\mu_t=a_k\end{aligned}\quad\mathrm{~vs.~}\quad H_a{:}(\operatorname{not}H_0)
(\#eq:1-8)
\end{equation}

本节中给出的结果是用向量和矩阵来说明的。然而，对于有矩阵处理软件的读者来说，向量和矩阵的知识并不是真正必要的，因为大多数计算机甚至允许新手用户轻松地进行矩阵计算。

式 \@ref(eq:1-8) 中的假设可以用矩阵表示法写为：

\begin{equation}
H_0{:}\,\boldsymbol{C}\boldsymbol\mu=\boldsymbol a\quad\mathrm{~vs.~}\quad H_a{:}\,\boldsymbol{C}\boldsymbol\mu\neq\boldsymbol a
(\#eq:1-9)
\end{equation}

其中

\begin{equation}
\boldsymbol{C}=\begin{bmatrix}c_{11}&c_{12}&\cdots&c_{1t}\\c_{21}&c_{22}&\cdots&c_{2t}\\\vdots&\vdots&\ddots&\vdots\\c_{k1}&c_{k2}&\cdots&c_{kt}\end{bmatrix},\quad\boldsymbol{\mu}=\begin{bmatrix}\mu_{1}\\\mu_{2}\\\vdots\\\mu_{t}\end{bmatrix},\quad\mathrm{and}\quad\boldsymbol{a}=\begin{bmatrix}a_{1}\\a_{2}\\\vdots\\a_{k}\end{bmatrix}
(\#eq:1-10)
\end{equation}

假设 $\boldsymbol C$ 中的 $k$ 行是线性独立的，这意味着 $\boldsymbol C$ 中没有一行可以表示为其余行的线性组合。如果 $\boldsymbol C$ 中的 $k$ 行不是线性独立的，那么总是可以选择 $\boldsymbol C$ 的子集：包含所有线性独立的行，它们包含关于所需假设的所有必要信息。

例如，假设你有三种处理，你希望检验

$$H_0{:}\,\mu_1-\mu_2=0,\,\mu_1-\mu_3=0\quad\mathrm{and}\quad\mu_2-\mu_3=0$$

相应的 $\boldsymbol C$ 矩阵为

$$\boldsymbol C=\begin{bmatrix}1&-1&0\\1&0&-1\\0&1&-1\end{bmatrix}$$

但是 $\boldsymbol C$ 的第三行是第二行和第一行之间的差，因此这三行不是线性独立的。此时，一个等价的假设可以表述为 $H_{0}\colon\,\mu_{\mathrm{l}}-\mu_{2}=0\mathrm{~and~}\mu_{\mathrm{l}}-\mu_{\mathrm{3}}=0$，因为如果 $\mu_1-\mu_2=0$ 且 $\mu_1-\mu_3=0$，那么 $\mu_2-\mu_3$ 必为 $0$.  下面的讨论使用 $\boldsymbol C$ 的行是线性独立的假设。

将样本均值向量表示为 $\hat{\boldsymbol \mu}$，则 $\hat{\boldsymbol \mu}$ 抽样分布的矩阵表示为

$$\left.\hat{\boldsymbol{\mu}}\sim {N_t}(\boldsymbol{\mu},{\sigma}^2\boldsymbol{D})\quad\mathrm{~其中~}\quad \boldsymbol D=\left[\begin{array}{ccccc}1/n_1&0&\cdots&0\\0&1/n_2&\cdots&0\\\vdots&\vdots&\ddots&\vdots\\0&0&\cdots&1/n_t\end{array}\right.\right]$$

此式读作：$t×1$ 向量 $\hat{\boldsymbol \mu}$ 的元素的联合抽样分布为 $t$ 元正态分布，其均值为 $\boldsymbol \mu$，方差和协方差由矩阵 $\boldsymbol \sigma^2\boldsymbol D$ 中的元素给出。$\boldsymbol \sigma^2\boldsymbol D$ 第 $i$ 个对角元为 $\hat\mu_i$ 的方差，非对角元 $(i,\,j)$ 为 $\hat\mu_i,\hat\mu_j$ 的协方差。

$\boldsymbol C \hat{\boldsymbol \mu}$ 的抽样分布为

$$\boldsymbol C\hat{\boldsymbol{\mu}}\sim N_k(\boldsymbol C\boldsymbol{\mu},\sigma^2 \boldsymbol {CDC}^{\prime})$$

偏离 $H_0$ 的平方和或用于检验 $H_0{:}\,\boldsymbol C \boldsymbol \mu =\boldsymbol a$ 的平方和由下式给出

\begin{equation}
SS_{H0}=(\boldsymbol C\boldsymbol{\hat{\mu}}-\boldsymbol{a})^{\prime}(\boldsymbol{C}\boldsymbol D\boldsymbol{C}^{\prime})^{-1}(\boldsymbol{C}\hat{\boldsymbol{\mu}}-\boldsymbol{a})
(\#eq:1-11)
\end{equation}

基于 $k$ 个自由度，即 $\boldsymbol C$ 线性无关行的数量。利用正态性假设，$SS_{H0}/\sigma^2$ 的抽样分布为自由度为 $k$ 的非中心卡方分布，如果 $H_0$ 为真，则 $SS_{H_0}/\sigma^2\sim\chi^2_k$. 用于检验 $H_0$ 的统计量为

$$F_c=\frac{SS_{H0}/k}{\hat{\sigma}^2}$$

若 $F_c>F_{\alpha,k,N-t}$，其中 $F_{\alpha,k,N-t}$ 是分子分母自由度分别为 $k$ 和 $N-t$ 的 $F$ 分布的上 $\alpha$ 分位数，假设 $H_0{:}\,\boldsymbol C\boldsymbol \mu=\boldsymbol \alpha$. 这里给出的结果是 Graybill (1976) 中定理 6.3.1 的特殊情况。

当 $H_0$ 为真，$\frac{SS_{H0}/k}{\hat{\sigma}^2}$ 为 $\sigma^2$ 的无偏估计，然后将其与 $\hat\sigma^2$ 进行比较，而无论 $H_0$ 是否为真，$\hat\sigma^2$ 又是 $\sigma^2$ 的无偏估计。因此若 $H_0$ 为真，上面给出的 $F$ 统计量应该接近 $1$。若 $H_0$ 为假， $SS_{H0}/k$ 是下式的无偏估计

$$\sigma^2+\frac1k(\boldsymbol C\boldsymbol\mu-\boldsymbol a)^{\prime}(\boldsymbol C\boldsymbol D\boldsymbol C^{\prime})^{-1}(\boldsymbol C\boldsymbol\mu-\boldsymbol a)$$

因此，若 $H_0$ 为假，则 $F$ 统计量的值应大于 $1$ 。若计算的 $F$ 统计量显著大于 $1$ ，则拒绝假设 $H_0$.

## 示例：任务和脉搏率（续） {#sec1-6}

以下是 \@ref(sec1-4) 节示例中的信息总结，包含六项任务中每项任务的样本量和平均值。

| 任务 $i$ | $1$ | $2$ | $3$ | $4$ | $5$ | $6$ |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| $n_i$ | $13$ | $12$ | $10$ | $10$ | $12$ | $11$|
| $\bar y_{i\cdot}$ | $31.9231$ | $31.0833$ | $35.8000$ | $38.0000$ | $29.5000$ |  $28.8182$ |

方差的合并估计为 $\hat\sigma^2 = 30.9045$，基于 $62$ 个自由度。与均值向量估计抽样分布相关的 $\boldsymbol D$ 矩阵为

$$\boldsymbol D=\begin{bmatrix}\frac1{13}&0&0&0&0&0\\0&\frac1{12}&0&0&0&0\\0&0&\frac1{10}&0&0&0\\0&0&0&\frac1{10}&0&0\\0&0&0&0&\frac1{12}&0\\0&0&0&0&0&\frac1{11}\end{bmatrix}$$

假设研究者有兴趣同时检验以下假设 (simultaneously testing the following hypothesis)，该假设涉及任务均值的两个线性组合：

$$H_0{:}\,\mu_4-\mu_5=4\,\mathrm{~and~}\,3\mu_1-\mu_2-\mu_3-\mu_4=0\,\mathrm{~vs.~}\,H_a{:}\,(\operatorname{not}H_0)$$

提供 $SS_{H0}$ 值所需的初步计算如下：

$$\begin{aligned}
\boldsymbol{C\hat{\mu}}-\boldsymbol{a}& =\begin{bmatrix}8.5-4\\-9.114-0\end{bmatrix}=\begin{bmatrix}4.500\\-9.114\end{bmatrix}  \\
\boldsymbol C\boldsymbol D \boldsymbol C^\prime& =\begin{bmatrix}\dfrac{1}{10}+\dfrac{1}{12}&-\dfrac{1}{10}\\-\dfrac{1}{10}&\dfrac{9}{13}+\dfrac{1}{12}+\dfrac{1}{10}+\dfrac{1}{10}\end{bmatrix}  \\
&=\begin{bmatrix}0.1833&-0.1000\\-0.1000&0.9756\end{bmatrix} \\
(\boldsymbol C \boldsymbol D \boldsymbol C^{\prime})^{-1}& =\begin{bmatrix}5.7776&0.5922\\0.5922&1.0856\end{bmatrix} 
\end{aligned}$$

以及

$$SS_{H0}=(\boldsymbol C\hat{\boldsymbol \mu}-\boldsymbol a)^{\prime}(\boldsymbol C\boldsymbol D\boldsymbol C^{\prime})^{-1}(\boldsymbol C\hat{\boldsymbol \mu}-\boldsymbol a)=158.602$$

自由度为 $2$. 检验统计量为

$$F_c=\frac{158.602/2}{30.9045}=2.566$$

该 $F$ 统计量的显著性概率为 $\hat{\alpha}=\mathrm{Pr}\{F>2.566\}=0.0850$.

## 检验所有均值相等 {#sec1-7}

大多数研究人员感兴趣的第一个假设通常是检验平均值是否同时相等。假设为 $H_{0}\colon\,\mu_{1}=\mu_{2}=\cdots=\mu_{t}\,\mathrm{~vs.~}H_{a}{:}\,(\operatorname{not}H_{0})$。检验等均值假设可使用两个基本程序。对于本章讨论的特定情况，这两个程序给出相同的统计检验。然而，对于大多数混乱数据的情况（对于除单向外的处理结构），这两个程序可以给予不同的检验。第一个程序将在本节介绍，而第二个程序将在 \@ref(sec1-9) 节介绍。

等均值假设 $H_{0}\colon\,\mu_{1}=\mu_{2}=\cdots=\mu_{t}$ 与下列形式的假设等价：$H_{0}\,\colon\mu_{1}-\mu_{2}=0,\mu_{1}-\mu_{3}=0,\ldots,\mu_{1}-\mu_{t}=0$，或与涉及 $\mu_i$ 的 $t - 1$ 个线性独立线性组合的任何其它假设等价。对应于 $t - 1$ 个成对差异集合的 $\boldsymbol C$ 矩阵和 $\boldsymbol a$ 向量为：

$$\boldsymbol C=\begin{bmatrix}1&-1&0&0&\cdots&0\\1&0&-1&0&\cdots&0\\1&0&0&-1&\cdots&0\\1&0&0&0&\cdots&0\\\vdots&\vdots&\vdots&\vdots&\ddots&\vdots\\1&0&0&0&\cdots&-1\end{bmatrix}\quad\mathrm{and}\quad\boldsymbol{a}=\begin{bmatrix}0\\0\\0\\0\\\vdots\\0\end{bmatrix}$$

对应于 $\mu_i$ 的以下 $t - 1$ 个线性无关线性组合的集合 $H_{0}{:}\,\mu_{1}-\mu_{2}=0,\mu_{1}+\mu_{2}-2\mu_{3}=0,\mu_{1}+\mu_{2}+\mu_{3}-3\mu_{4}=0,\ldots,\mu_{1}+\mu_{2}+\cdots-\left(t-1\right)\mu_{t}=0$ 的 $\boldsymbol C$ 和 $\boldsymbol a$ 向量为：

$$\boldsymbol C=\begin{bmatrix}1&-1&0&0&\cdots&0\\1&1&-2&0&\cdots&0\\1&1&1&-3&\cdots&0\\\vdots&\vdots&\vdots&\vdots&\ddots&0\\1&1&1&1&\cdots&t-1\end{bmatrix}\quad\mathrm{and}\quad\boldsymbol{a}=\begin{bmatrix}0\\0\\0\\\vdots\\0\end{bmatrix}$$

存在许多其他矩阵使得 $\boldsymbol C \boldsymbol \mu = \boldsymbol 0$ 当且仅当 $\mu_1 = \mu_2 = \cdots = \mu_t$；然而所有这些矩阵产生相同的偏离 $H_0$ 的平方和以及相同的自由度 $t - 1$，因此产生相同的 $F$ 统计量。对于这种特殊情况，式 \@ref(eq:1-11) 总可以简化为：

\begin{equation}
SS_{H0:\,\mu_{1}=\mu_{2}=\cdots=\mu_{t}}=\sum_{i=1}^{t}n_{i}(\overline{y}_{i.}-\overline{y}_{..})^{2}=\sum_{i=1}^{t}\left(\frac{y_{i.}^{2}}{n_{i}}\right)-\frac{y_{..}^{2}}{N}
(#eq:1-12)
\end{equation}

## 示例：任务和脉搏率（续） {#sec1-8}

对于 \@ref(sec1-4) 节中的任务和脉搏数据，使用式 \@ref(eq:1-11) 和 \@ref(eq:1-12) 可计算 $SS_{H0:\,\mu_{1}=\mu_{2}=\cdots=\mu_{t}}$.

使用式 \@ref(eq:1-12) 可得到

$$\begin{aligned}
SS_{H0}& =\frac{415^2}{13}+\frac{373^2}{12}+\frac{358^2}{10}+\frac{380^2}{10}+\frac{354^2}{12}+\frac{317^2}{11}-\frac{2197^2}{68}  \\
&=694.4386
\end{aligned}$$

自由度为 $t - 1 = 5$. $F_c$ 统计量的值为

$$F_c=\frac{694.4386/5}{30.9045}=4.49$$

并且显著性概率为 $\hat \alpha = 0.0015$.

然后，使用式 \@ref(eq:1-11)，矩阵 $\boldsymbol C$，向量 $\boldsymbol a$ 和矩阵 $\boldsymbol D$ 为

$$\boldsymbol C=\begin{bmatrix}1&-1&0&0&0&0\\1&0&-1&0&0&0\\1&0&0&-1&0&0\\1&0&0&0&-1&0\\1&0&0&0&0&-1\end{bmatrix},\quad \boldsymbol a=\begin{bmatrix}0\\0\\0\\0\\0\end{bmatrix}$$

$$\boldsymbol D=\begin{bmatrix}\frac{1}{13}&0&0&0&0&0\\0&\frac{1}{12}&0&0&0&0\\0&0&\frac{1}{10}&0&0&0\\0&0&0&\frac{1}{10}&0&0\\0&0&0&0&\frac{1}{12}&0\\0&0&0&0&0&\frac{1}{11}\end{bmatrix}$$

接下来，计算式 \@ref(eq:1-11) 中的各个量

$$\boldsymbol{C\hat{\boldsymbol{\mu}}}-\boldsymbol{a}=\begin{bmatrix}0.844\\-3.877\\-6.077\\2.423\\3.105\end{bmatrix}\quad\mathrm{and}\quad\boldsymbol{C}\boldsymbol D\boldsymbol{C}^{\prime}=\begin{bmatrix}\frac{25}{156}&\frac1{13}&\frac1{13}&\frac1{13}&\frac1{13}\\\frac1{13}&\frac{23}{130}&\frac1{13}&\frac1{13}&\frac1{13}\\\frac1{13}&\frac1{13}&\frac{23}{130}&\frac1{13}&\frac1{13}\\\frac1{13}&\frac1{13}&\frac1{13}&\frac{25}{156}&\frac1{13}\\\frac1{13}&\frac1{13}&\frac1{13}&\frac1{13}&\frac{24}{143}\end{bmatrix}$$

$\boldsymbol{CDC'}$ 的逆为

$$(\boldsymbol{CDC'})^{-1}=\begin{bmatrix}9.882&-1.765&-1.765&-2.118&-1.941\\-1.765&8.529&-1.471&-1.765&-1.618\\-1.765&-1.471&8.529&-1.765&-1.618\\-2.118&-1.765&-1.765&9.882&-1.941\\-1.941&-1.618&-1.618&-1.941&9.221\end{bmatrix}$$

最后，平方和的值为


$$SS_{H0}=(\boldsymbol C\hat{\boldsymbol{\mu}}-\boldsymbol{a})^{\prime}(\boldsymbol C\boldsymbol D\boldsymbol{C}^{\prime})^{-1}(\boldsymbol{C}\hat{\boldsymbol{\mu}}-\boldsymbol{a})=694.4386$$

这与使用式 \@ref(eq:1-12) 计算的平方和相同。

显然，如果必须手工计算，这个公式并不容易使用。然而，在许多混乱数据的情况下，为了获得检验有意义假设的统计量，像这样的公式是必要的。幸运的是，通过利用计算机，可以为特定假设构造 $\boldsymbol C$矩阵，然后可以让计算机进行繁琐的计算。

## 比较两种模型的一般方法：条件误差原理 {#sec1-9}

计算检验统计量的第二个过程是比较两个模型的拟合度。在本节中，比较的两个模型是 $y_{ij} = \mu_i + \varepsilon_{ij}$，这叫做一般 (general) 或**未缩减模型** (**unreduced model**)，$y_{ij} = \mu_i + \varepsilon_{ij}$，这是如果 $H_0:\,\mu_1 = \mu_2 = \cdots = \mu_t = \mu$（假设）成立时得到的模型。第一个模型称为**完全模型** (**full model**) 或**无约束模型** (**unrestricted model**)，第二个模型称为**缩减模型** (**reduced model**) 或**约束模型** (**restricted model**).

**条件误差原理** (**principle of conditional error**) 用于比较两个模型，其中一个模型是通过对另一个模型的参数施加限制而获得的。该原理非常简单，要求获得完全模型和缩减模型的残差或误差平方和。令 $ESS_F$ 表示拟合完全模型得到的误差平方和，$ESS_R$ 表示拟合缩减模型得到的误差平方和。则由于假设的限制而导致的，或由于偏离零假设的平方和为 $SS_{H0}=ESS_R-ESS_F$. $ESS_R$ 和 $ESS_F$ 的自由度由数据集中的观测总数与待估计的（基本）参数数量之差给出 [**基本参数** (**essential parameters**) 将在第 \@ref(chap6) 章中讨论]。分别用 $df_R$ 和 $df_F$ 表示与 $ESS_R$ 和 $ESS_F$ 对应的自由度。对应于 $SS_{H0}$ 的自由度为 $df_{H0} = df_R - df_F$. 用于检验 $H_0$ 的 $F$ 统计量由下式给出：

$$F_c=\frac{SS_{H0}/df_{H0}}{ESS_F/df_F}$$

若 $F_{c}>F_{\alpha,df_{{H0}},df_{F}}$，则在显著性水平下拒绝 $H_0$.

对于上面讨论的情况，$y_{ij} = \mu_i + \varepsilon_{ij}$ 是完全模型而 $y_{ij} = \mu + \varepsilon_{ij}$ 是缩减模型。完全模型的误差平方和为

$$ESS_{F}=\sum_{i=1}^{t}\sum_{j=1}^{n_{i}}(y_{ij}-\overline{y}_{i\cdot})^{2}=(N-t)\hat{\sigma}^{2}$$

其自由度 $df_F=N-t$. 缩减模型的误差平方和为

$$ESS_{_R}=\sum_{i=1}^{t}\sum_{j=1}^{n_{i}}(y_{_{ij}}-\overline{y}_{\cdot\cdot})^{2}$$

其自由度 $df_R=N-1$. 因此，由于偏离 $H_0$ 而产生的平方和为：

$$SS_{H0{:}\,\mu_1=\mu_2=\cdots=\mu_t}=ESS_{R}-ESS_{F}=\sum_{i=1}^{t}n_{i}(\overline{y}_{i\cdot}-\overline{y}_{\cdot\cdot})^{2}$$

其自由度为 $t-1$. 这与式 \@ref(eq:1-12) 得到的平方和相同。

在检验情况下感兴趣的平方和通常被放在一个叫做**方差分析** (**analysis of variance**, **ANOVA**)表的表格中，这样的表格通常有类似于表 \@ref(tab:table1-2)  的形式，在 "source of variation" 栏下的条目被分组为多个集合。在给定的情况下，每组只使用一个标签，选择完全由实验者决定。

<table class="table table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;border-bottom: 0;">
<caption>(\#tab:table1-2)单向模型检验均值相等的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%201.2.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> df=degrees of freedom 表示自由度, SS=sum of square 表示平方和, MS=mean square 表示均方. 这些标准缩写贯穿本书</td></tr></tfoot>
</table>

## 示例：任务和脉搏率（续） {#sec1-10}

条件误差原则也被称为模型比较程序，其过程相当灵活。例如，如果对检验任务和脉搏率数据的假设感兴趣，例如 $H_{0}\colon\,\mu_{1}=\mu_{2}=\mu_{3}\,\mathrm{~vs.~}\,H_{a}\colon\,(\mathrm{not}\,H_{0})$，那么在 $H_0$ 的条件下模型有如下形式

$$\begin{array}{ll}y_{ij}=\mu_0+\varepsilon_{ij}&\mathrm{~for~}i=1,2,3\\y_{ij}=\mu_i+\varepsilon_{ij}&\mathrm{~for~}i=4,5,6\end{array}$$

这样的模型可以使用大多数软件包来拟合，其中包含一个定性或分类变量，当任务为 $1,2,3$ ，其值为 $0$，当任务为 $4,5,6$ ，其值等于任务值.

## 计算机分析 {#sec1-11}

将条件误差原理应用于 \@ref(sec1-4) 节的任务和脉搏率数据，以提供等均值假设的检验，$H_{0}\colon\,\mu_{1}=\mu_{2}=\mu_{3}=\mu_{4}=\mu_{5}=\mu_{6}\,\mathrm{~vs.~}\,H_{a}\colon(\mathrm{not}\,H_{0})$. 完全模型的误差平方和为 $ESS_F=1916.076$, $df_F=62$. 缩减模型的误差平方和为 $ESS_R=73593-(2197)^2/68=2610.545$, $df_R=67$. 因此，$SS_{H0}=2610.545-1916.076=694.439$, $df_{H0}=67-62=5$. 汇总这些计算的方差分析表如表 \@ref(tab:table1-3) 所示。

<table class="table table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:table1-3)任务和脉搏率数据检验均值相等的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%201.3.png">
</td>
  </tr>
</tbody>
</table>

本章最后介绍了一些关于使用计算机和统计计算软件包（如SAS^®^、BMDP^®^、SYSTAT^®^、JMP^®^和SPSS^®^）的说明。前面章节中提供的所有方法和公式都可以在大多数计算机上轻松使用。如果计算机使用 MATLAB、SAS-IML 或 APL 等编程语言，则所需的矩阵计算很容易按照前面章节中给出的矩阵公式进行。SAS、JMP、BMDP、SYSTAT 和 SPSS 每个都包含程序，使用户能够生成自己的处理的线性组合，以检验假设。此外，这些包都提供了方差分析表、处理均值及其标准误的分析。表 \@ref(tab:table1-4) 包含 SAS-GLM 代码，其中包含检验任务和脉搏率数据所述假设所需的 estimate 语句和 contrast 的语句。estimate 语句用于评估均值的一个线性组合，提供的结果是对比的估计、其标准误估计并生成 $t$ 统计量及其相应的显著性水平。contrast 语句用于评估均值的一个或多个线性组合，提供的结果是平方和、自由度和 $F$ 统计量。对于SAS-GLM中的 estimate 和 contrast 语句，假设中的 $\boldsymbol a$ 的唯一值为零，也就是说，只能检验均值的线性组合等于零。

<table class="table table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:table1-4)Proc GLM 代码，用于拟合任务和脉搏率数据，并提供文本中描述的分析所需的 estimate 和 contrast 语句</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%201.4.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table1-5) 包含 SAS-IML 代码，为 \@ref(sec1-6) 节中检验的假设提供计算。通过用矩阵语言构造代码，可以得到 $\boldsymbol C \boldsymbol \mu = \boldsymbol a$ 形式的任何假设的检验。

<table class="table table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:table1-5)Proc IML代码执行 1.6 节中任务和脉搏率数据的计算</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%201.5.png">
</td>
  </tr>
</tbody>
</table>

## 结束语 {#sec1-12}

本章对单向方差分析模型进行了分析。提供了对不同处理的效应进行统计推断的一般程序，并对同质误差的情况进行了说明。介绍了获得实验设计模型的统计分析的两个基本程序。这些程序在本书的其余部分被广泛用于描述更复杂的模型设计的实验和更混乱的数据情况。还给出了同时比较所有处理效应均值的检验。这种检验可以被认为是统计分析的初始步骤。完成数据集分析所应使用的程序取决于是否拒绝处理均值相等的假设。

## 练习 {#sec1-13}

## R 代码 {#sec1-14}


```r
# Chap 1 ----
library(gmodels)

data <- data.frame(
  task = factor(rep(1:6, times=c(13,12,10,10,12,11)), levels = 1:6), 
  y = c(27,31,26,32,39,37,38,39,30,28,27,27,34,
        29,28,37,24,35,40,40,31,30,25,29,25,
        34,36,34,41,30,44,44,32,32,31,
        34,34,43,44,40,47,34,31,45,28,
        28,28,26,35,31,30,34,34,26,20,41,21,
        28,26,29,25,35,34,37,28,21,28,26)
)

m <- aov(y ~ -1 + task, data) # "-1" 删去截距项，得到均值模型，IMPORTANT!
n <- aov(y ~ task, data)      # 保留截距项，得到效应模型，第一个水平被置零（置零限制，R 默认将第一个水平置零）

## a) test: mu3=30 ----
# 使用均值模型 m
c1 <- c(0,0,1,0,0,0)
names(c1) <- paste0("task",1:6) # 名字要与向量中的元素对应，如果向量的元素按顺序对应于因子的水平，则可以不赋名字
c11 <- c(1,0,0,0,0,0)
names(c11) <- paste0("task",c(3,1,2,4,5,6)) # 或者这样是等价的
estimable(m,c1,beta0 = 30)
```

```
##               beta0 Estimate Std. Error  t value DF    Pr(>|t|)
## (0 0 1 0 0 0)    30     35.8   1.757966 3.299267 62 0.001609292
```

```r
estimable(m,c11,beta0 = 30)
```

```
##               beta0 Estimate Std. Error  t value DF    Pr(>|t|)
## (0 0 1 0 0 0)    30     35.8   1.757966 3.299267 62 0.001609292
```

```r
# 可以为线性组合添加标签
c1 <- rbind("c1"=c(0,0,1,0,0,0))
estimable(m,c1,beta0 = 30)
```

```
##    beta0 Estimate Std. Error  t value DF    Pr(>|t|)
## c1    30     35.8   1.757966 3.299267 62 0.001609292
```

```r
# 效应模型 n 只给一个例子
d1 <- c(1,0,1,0,0,0) # 第一个 1 表示水平 1 的效应，随后是各水平相较水平 1 的效应
estimable(n,d1,beta0 = 30)
```

```
##               beta0 Estimate Std. Error  t value DF    Pr(>|t|)
## (1 0 1 0 0 0)    30     35.8   1.757966 3.299267 62 0.001609292
```

```r
## b) 95% CI for mu1 ----
c2 <- rbind("c2"=c(1,0,0,0,0,0))
estimable(m,c2,conf.int = 0.95)
```

```
##    Estimate Std. Error  t value DF Pr(>|t|) Lower.CI Upper.CI
## c2 31.92308   1.541838 20.70455 62        0 28.84099 35.00517
```

```r
## c) test: mu3=mu4 ----
c3 <- c(0,0,1,-1,0,0)
c31 <- c(0,0,-1,1,0,0) # 或等价地
c3 <- rbind("c3"=c3,
            "c31"=c31)
estimable(m,c3,beta0 = 0)
```

```
##     beta0 Estimate Std. Error   t value DF  Pr(>|t|)
## c3      0     -2.2    2.48614 -0.884906 62 0.3796264
## c31     0      2.2    2.48614  0.884906 62 0.3796264
```

```r
## d) test: mu1=(mu2+mu3+mu4)/3 ----
c4 <- c(1,-1/3,-1/3,-1/3,0,0)
c41 <- -c4
c42 <- c(3,-1,-1,-1,0,0)
c4 <- rbind("c4"=c4, # 三个都等价
            "c41"=c41,
            "c42"=c42)
estimable(m,c4,beta0 = 0)
```

```
##     beta0  Estimate Std. Error  t value DF  Pr(>|t|)
## c4      0 -3.038034   1.830351 -1.65981 62 0.1020029
## c41     0  3.038034   1.830351  1.65981 62 0.1020029
## c42     0 -9.114103   5.491052 -1.65981 62 0.1020029
```

```r
# 事实上这是一个对比
library(emmeans)
contrastc4 <- emmeans(m, ~ task)
c4 <- list("c4"=c(3,-1,-1,-1,0,0))
test(contrast(contrastc4,c4))
```

```
##  contrast estimate   SE df t.ratio p.value
##  c4          -9.11 5.49 62  -1.660  0.1020
```

```r
## e) 90% CI for 4mu1-mu3-mu4-mu5-mu6 ----
c5 <- rbind("c5"=c(4,0,-1,-1,-1,-1))
estimable(m,c5,conf.int = .90)
```

```
##     Estimate Std. Error    t value DF  Pr(>|t|)  Lower.CI Upper.CI
## c5 -4.425874   7.042869 -0.6284192 62 0.5320376 -16.18609 7.334338
```

```r
# 事实上这也是一个对比
contrastc5 <- emmeans(m, ~ task)
c5 <- list("c5"=c(4,0,-1,-1,-1,-1))
confint(contrast(contrastc5,c5),level = .90)
```

```
##  contrast estimate   SE df lower.CL upper.CL
##  c5          -4.43 7.04 62    -16.2     7.33
## 
## Confidence level used: 0.9
```

```r
## sec 1.6, https://bookdown.org/wangzhen/AMD/chap1.html#sec1-6 ----

c <- rbind("a"=c(0,0,0,1,-1,0),
           "b"=c(3,-1,-1,-1,0,0))
estimable(m,c,beta0 = c(4,0),joint.test=TRUE) # 需指定联合检验（同时检验）
```

```
##    X2.stat DF Pr(>|X^2|)
## 1 5.132057  2 0.07684013
```

```r
# 书上用的是 F 检验，这里用的是 Wald 检验
# 改用 car::linearHypothesis
library(car)
linearHypothesis(m,
                 c("task4-task5=4",
                   "3*task1-task2-task3-task4=0"), 
                 test = 'F') # 指定 F 检验
```

```
## Linear hypothesis test
## 
## Hypothesis:
## task4 - task5 = 4
## 3 task1 - task2 - task3 - task4 = 0
## 
## Model 1: restricted model
## Model 2: y ~ -1 + task
## 
##   Res.Df    RSS Df Sum of Sq     F  Pr(>F)  
## 1     64 2074.7                             
## 2     62 1916.1  2     158.6 2.566 0.08498 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# 事实上 linearHypothesis 也可以做上面所有的检验

## sec 1.7 - 1.11, https://bookdown.org/wangzhen/AMD/chap1.html#sec1-7 ----

linearHypothesis(m,
                 c("task2-task1=0",
                   "task3-task2=0",
                   "task4-task3=0",
                   "task5-task4=0",
                   "task6-task5=0"), 
                 test = 'F')
```

```
## Linear hypothesis test
## 
## Hypothesis:
## - task1  + task2 = 0
## - task2  + task3 = 0
## - task3  + task4 = 0
## - task4  + task5 = 0
## - task5  + task6 = 0
## 
## Model 1: restricted model
## Model 2: y ~ -1 + task
## 
##   Res.Df    RSS Df Sum of Sq      F   Pr(>F)   
## 1     67 2610.5                                
## 2     62 1916.1  5    694.44 4.4941 0.001471 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
