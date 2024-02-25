
# 使用均值模型分析具有不等子类数的均衡双向处理结构 {#chap9}

第 \@ref(chap7) 章和第 \@ref(chap8) 章考虑了相同样本量的情况，其中每个处理组合的观测次数相同。第 \@ref(chap13)-\@ref(chap15) 章考虑了一些处理组合缺失的情况，但本章以及第 \@ref(chap10)-\@ref(chap12) 章假设观测到了每种处理组合，并且每种组合观测了多次。

## 模型定义和假设 {#sec9-1}

如第 \@ref(sec7-1-1) 节所示，令 $\mu_{ij}$ 为当处理 T 的可能性 $i$ 和处理 B 的可能性 $j$ 均应用于同一实验单元时的预期响应。本章假设观测到的响应 $Y_{ijk}$ 可以通过以下方式建模

\begin{equation}
Y_{ijk}=\mu_{ij}+\varepsilon_{ijk},\quad i=1,2,\ldots,t,\quad j=1,2,\ldots,b,k=1,2,\ldots,n_{ij}
(#eq:9-1)
\end{equation}

满足理想条件

$$\varepsilon_{ijk}\sim i.i.d.N(0,\sigma^2),\quad i=1,2,\ldots,t,\quad j=1,2,\ldots,b,\quad k=1,2,\ldots,n_{ij}$$

以及对于每个 $i$ 和 $j$，$n_{ij}>0$.

## 参数估计 {#sec9-2}

若将 $bt$ 个处理组合视为 $bt$ 个不同的处理，那么第 \@ref(chap1)-\@ref(chap3) 章中讨论的单向处理结构的所有内容也适用于双向处理结构。对于不均衡的数据问题，这通常是分析数据最佳和最简单的方法。均值模型中参数的最佳估计是

\begin{equation}
\hat{\mu}_{ij}=\frac1{n_{ij}}\sum_{k=1}^{n_{ij}}y_{ijk}=\bar{y}_{ij\cdot},i=1,2,\ldots,t,j=1,2,\ldots,b
(#eq:9-2)
\end{equation}

以及

\begin{equation}
\hat{\sigma}^2=\frac1{N-bt}\sum_{ijk}(y_{ijk}-\overline{y}_{ij\cdot})^2
(#eq:9-3)
\end{equation}

其中 $N=n_{\cdot}$. $\hat \mu_{ij},\hat \sigma^2$ 的抽样分布为

$$\hat{\mu}_{ij}\thicksim N(\mu_{ij},\sigma^2/n_{ij})\quad\mathrm{for~}i=1,2,\ldots,t,j=1,2,\ldots,b$$

以及

$$(N-bt)\hat{\sigma}^2/\sigma^2\thicksim\chi^2(N-bt)$$

此外，所有的 $\hat \mu_{ij},\hat \sigma^2$ 独立分布，像之前章节一样。

当数据不均衡时，实验者通常想回答与数据均衡时相同的问题。回想一下，这些问题是：

1. 两组处理发生交互吗？
2. T 处理如何影响响应？
3. B 处理如何影响响应？

根据均值模型的参数，这些问题可以表述为假设，相应的假设是：

$$\begin{array}{l}H_{T\times B}\colon\mu_{ij}-\mu_{i'j}-\mu_{ij'}+\mu_{i'j'}=0&\text{ for all }i\neq i'\text{ and }j\neq j'\\H_T\colon\bar{\mu}_{1\cdot}=\bar{\mu}_{2\cdot}=\cdots=\bar{\mu}_{t\cdot}\\H_B\colon\bar{\mu}_{\cdot1}=\bar{\mu}_{\cdot2}=\cdots=\bar{\mu}_{\cdot b}\end{array}$$

**检验上述假设应被视为分析任何双向实验的第一步**。通常会有明确的对比，直接解决研究人员感兴趣的其他问题。假设 $H_{T\times B},H_T,H_B$ 经常被检验，以帮助选择适当的多重比较程序来解决可能感兴趣的其他问题。

如式 \@ref(sec7-2) 所示，始终存在参数 $\mu,\tau_i,\beta_j,\gamma_{ij},i=1,2,\cdots,t,j=1,2,\cdots,b$ 使得 $\mu_{ij}$ 可以在效应模型中表示为

$$\mu_{ij}=\mu+\tau_i+\beta_j+\gamma_{ij},\quad i=1,2,\ldots,t,\quad j=1,2,\ldots,b$$

许多实验者更喜欢使用效应模型来查看处理组合方式的表示。这可能是因为许多统计学家鼓励实验者考虑这样的模型。因此，许多现有的计算机软件引导数据分析师使用效应模型的表示。本书考虑了这两类模型。本章考虑均值模型，第 \@ref(chap10) 章考虑效应模型。第 \@ref(sec1-5) 节和第 \@ref(sec1-6) 节介绍了计算检验统计量的不同程序。对于单向情况，所有这些方法都会产生相同的检验统计量。事实上，对于匀称的 (well-balanced) 双向实验也总是如此。然而，不均衡数据集的情况并不总是如此。匀称意味着每种处理组合都有相同数量的观测结果。本章使用矩阵程序来获取检验统计量，第 \@ref(chap10) 章使用模型拟合程序。

## 检验所有均值是否相等 {#sec9-3}

考虑表 \@ref(tab:table9-1) 中的数据。数据来自于在完全随机设计结构中进行的小型双向处理结构实验。首先，计算表 \@ref(tab:table9-1) 中数据的双向单元格均值和边际均值。表 \@ref(tab:table9-2) 给出了这些均值，其中**行边际均值** (row marginal mean) 定义为给定行单元格均值的均值，**列边际均值** (column marginal mean) 定义为给定列单元格均值的均值。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table9-1)一个不均衡的双向实验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%209.1.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table9-2)表 \@ref(tab:table9-2) 中数据的单元平均值和边际平均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%209.2.png">
</td>
  </tr>
</tbody>
</table>

应注意的是，在不均衡双向实验中，有两种不同的方法可以计算 T 均值和 B 均值。表 \@ref(tab:table9-1) 给出了通过行（列）综总和并除以行（列）观测数量计算得出的平均值。表 \@ref(tab:table9-2) 给出了分两步计算的平均值。首先计算每个 T × B 组合的单元格均值，然后计算每行和列的单元格均值。除非实验是匀称的，否则这两种方法通常会给出不同的答案。例如，表 \@ref(tab:table9-1) 中的 $T_1$ 平均值为 22.75，而表 \@ref(tab:table9-2) 中的 $T_1$ 平均值为 23. 因此，必须解决的一个问题是在计算处理主效应均值时应使用这两种方法中的哪一种。应使用哪种方法将在 \@ref(sec9-5) 节中讨论。

表 \@ref(tab:table9-1) 中数据的实验误差方差为

$$\hat{\sigma}^2=\frac1{N-bt}\sum_{ijk}(y_{ijk}-\overline{y}_{ij.})^2=20/10=2$$

自由度为 $N-bt=16-6=10$.

接下来，将该实验视为具有六种处理的单向处理结构，并检验所有六种处理组合均值是否彼此相等。也就是说，考虑检验

$$\begin{aligned}H_0\colon\mu_{11}=\mu_{12}=\mu_{13}=\mu_{21}=\mu_{22}=\mu_{23}\end{aligned}$$

使用式 \@ref(eq:8-1) 可得

$$SSH_0=\frac{60^2}3+\frac{50^2}2+\frac{72^2}3+\frac{52^2}2+\frac{69^2}3+\frac{96^2}3-\frac{399^2}{16}=238.9375$$

基于 5 个自由度。因此用于检验 $H_0$ 的 $F$ 统计量为

$$F_c=\frac{238.9375/5}2=23.89$$

在 $\hat\alpha = 0.00003$ 水平上显著。因此，$H_0$ 将被拒绝，因此六种不同处理组合均值之间存在显著差异。

## 交互作用和主效应假设 {#sec9-4}

在上一节中，确定了六种处理组合均值之间存在显著差异。现在有必要看看差异出现在哪里。第一步，考虑表 \@ref(tab:table9-1) 中给出的数据是否存在显着的 T × B 交互作用并检验

$$\begin{aligned}H_{T\times B}\colon\mu_{ij}-\mu_{i'j}-\mu_{ij'}+\mu_{i'j'}=0&&\text{ for all }i\neq i'&\mathrm{~and~}j\neq j'\end{aligned}$$

这可以通过利用第 \@ref(sec1-4) 节中讨论的矩阵程序来开发检验统计量来实现。假设 $H_{T\times B}$ 为真当且仅当

$$\mu_{11}-\mu_{12}-\mu_{21}+\mu_{22}=0$$

以及

$$\mu_{11}-\mu_{13}-\mu_{21}+\mu_{23}=0$$

反过来，这些声明为真当且仅当 $\boldsymbol C\boldsymbol \mu=0$，其中

$$\boldsymbol C=\begin{bmatrix}1&-1&0&-1&1&0\\1&0&-1&-1&0&1\end{bmatrix}\quad\mathrm{and}\quad \boldsymbol \mu^{\prime}=[\mu_{11}~\mu_{12}~\mu_{13}~\mu_{21}~\mu_{22}~\mu_{23}]$$

根据式 \@ref(eq:1-11)，由于 $\boldsymbol D=\mathrm{Diag}[1/3,1/2,1/3,1/2,1/3,1/3]$

$$SSH_{T \times B}=[\boldsymbol C\hat{\boldsymbol{\mu}}]^{\prime}[\boldsymbol C\boldsymbol D\boldsymbol C^{\prime}]^{-1}[\boldsymbol C\hat{\boldsymbol{\mu}}]=[-8\quad2]\begin{bmatrix}\frac{10}6&\frac56\\\frac56&\frac96\end{bmatrix}^{-1}\begin{bmatrix}-8\\2\end{bmatrix}$$

因此 $SSH_{T\times B}=776(6/65)=71.631$ 且基于 2 个自由度。相应的 $F$ 统计量为 $F_c=(71.631/2)/2=17.91$ 且基于自由度 2 和 10，在 $\hat \alpha=0.0005$ 水平下显著。其他 $\boldsymbol C$ 矩阵也可以用来检验无交互作用的假设，但所有这些矩阵都会产生相同的检验统计量。读者可以验证这是真的。

接下来考虑检验预期行边际均值的相等性。这里出于说明目的，尽管由于显著的 T × B 交互作用，这种检验在这里可能不合适。合适的假设为 $H_T{:{\bar{\mu}}_{1\cdot}}=\bar{{\mu}}_{2\cdot}$. 注意到 $H_T$ 为真当且仅当 $\boldsymbol C\boldsymbol \mu=\boldsymbol 0$，其中 $\boldsymbol C=[1\,\,1\,\,1\,\,-1\,\,-1\,\,-1]$. 使用式 \@ref(eq:1-11) 可得

$$SSH_T=[\boldsymbol C\hat{\boldsymbol{\mu}}]^{\prime}[\boldsymbol C\boldsymbol D\boldsymbol C^{\prime}]^{-1}[\boldsymbol C\hat{\boldsymbol{\mu}}]=[-12][14/6]^{-1}[-12]=61.714$$

基于 1 个自由度，且相应的 $F$ 统计量为 $F_{c}=(61.714/1)/2=30.857$，在 $\hat \alpha=0.00024$ 水平下显著。也可以取 $\boldsymbol C=[\frac{1}{3}\,\,\frac{1}{3}\,\,\frac{1}{3}\,\,-\frac{1}{3}\,\,-\frac{1}{3}\,\,-\frac{1}{3}]$。读者应该验证 $\boldsymbol C$ 的第二个选择是否会产生与上面使用的 $\boldsymbol C$ 相同的 $F$ 检验统计量。

最后，考虑通过检验 $H_B\colon\bar{\mu}_{\cdot1}=\bar{\mu}_{\cdot2}=\bar{\mu}_{\cdot3}$ 来测试预期列边际均值的相等性。注意，$H_B$ 为真当且仅当 $\boldsymbol C\boldsymbol \mu=\boldsymbol 0$ 时，其中

$$\boldsymbol C=\begin{bmatrix}1&-1&0&1&-1&0\\1&0&-1&1&0&-1\end{bmatrix}$$

使用式 \@ref(eq:1-11) 可得

$$\begin{aligned}SSH_B&=[C\hat{\boldsymbol{\mu}}]^{\prime}[\boldsymbol C\boldsymbol D\boldsymbol C^{\prime}]^{-1}[\boldsymbol C\hat{\boldsymbol{\mu}}]=[-2\,\,-10]\begin{bmatrix}\frac{10}6&\frac56\\\frac56&\frac96\end{bmatrix}^{-1}\begin{bmatrix}-2\\-10\end{bmatrix}\\\\&=[-2\,\,-10]\frac6{65}\begin{bmatrix}9&-5\\-5&10\end{bmatrix}\begin{bmatrix}-2\\-10\end{bmatrix}=77.169\end{aligned}$$

基于 2 个自由度。因此，相应的 $F$ 统计量为 $F_c=19.29$，在 $\hat \alpha=0.0037$ 水平下显著。

表 \@ref(tab:table9-3) 中给出的方差分析表总结了上述检验。请注意，此分析结果与第 \@ref(chap7) 章讨论的均衡情形下的结果之间存在一些差异。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table9-3)表 \@ref(tab:table9-1) 中数据的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%209.3.png">
</td>
  </tr>
</tbody>
</table>

1. 对于均衡的数据，总有
$$SS_T+SS_B+SS_{T\times B}=SS_{\mu_{11}=\mu_{12}=\cdot\cdot\cdot=\mu_{23}}$$
对于不均衡的数据通常不然。
2. 对于均衡的数据，$SS_T,SS_B,SS_{T\times B}$ 在统计学上是独立的；对于不均衡的数据，这通常是不正确的。

人们不必过于担心 T, B 和 T × B 的平方和之和不等于检验所有均值相等的平方和，或者 T, B 和 T × B 的平方和在统计上不独立。之所以在此提出，是因为它们确实成立，而并非因为它们是必须解决的问题。

还有其他一些平方和经常与双向处理结构的分析有关，第 \@ref(chap10) 章将讨论其中两种常用的平方和。

## 总体边际均值 {#sec9-5}

实验者常常对比较每个主效应的可能性感兴趣。在均衡的情况下，人们可以比较 $\bar{\mu}_{1\cdot},\bar{\mu}_{2\cdot},\ldots,\bar{\mu}_{t\cdot}$。如前所述，这些均值在均衡和不均衡的情况下都被称为总体边际均值。

$\bar{\mu}_{i\cdot}$ 的最佳估计是

\begin{equation}
\hat{\bar{\mu}}_{i\cdot}=\frac1b\sum_{j=1}^b\hat{\mu}_{ij}=\bar{\hat{\mu}}_{i\cdot}\quad i=1,2,\ldots,t
(#eq:9-4)
\end{equation}

$\hat{\bar{\mu}}_{i\cdot}$ 的标准误估计为

\begin{equation}
\widehat{s.e.}(\hat{\bar{\mu}}_{i\cdot})=\frac{\hat{\sigma}}b\sqrt{\sum_{j=1}^b\frac1{n_{ij}}},\quad i=1,2,\ldots,t
(#eq:9-5)
\end{equation}

$\bar{\mu}_{\cdot j}$ 的最佳估计是

\begin{equation}
\hat{\bar{\mu}}_{\cdot j}=\frac1t\sum_{i=1}^t\hat{\mu}_{ij}=\bar{\hat{\mu}}_{\cdot j},\quad j=1,2,\ldots,b
(#eq:9-6)
\end{equation}

\begin{equation}
\widehat{s.e.}(\hat{\bar{\mu}}_{\cdot j})=\frac{\hat{\sigma}}t\sqrt{\sum_{i=1}^t\frac1{n_{ij}}},\quad j=1,2,\ldots,b
(#eq:9-7)
\end{equation}

需注意，在不均衡数据问题中，通常 $\hat{\bar{\mu}}_{i\cdot}$ 不同于 $\bar y_{i\cdot}$ 且 $\hat{\bar{\mu}}_{\cdot j}$ 不同于 $\bar y_{\cdot j}$. 在这个例子中，$\bar{y}_{1\cdot\cdot}=22.75,\hat {\bar \mu}_{1\cdot}=23$.

$\hat{\bar{\mu}}_{i\cdot},\hat{\bar{\mu}}_{\cdot j}$ 分别为 $\bar{\mu}_{i\cdot},\bar{\mu}_{\cdot j}$ 的无偏估计。注意到 $\bar{y}_{i\cdot\cdot},\bar{y}_{\cdot j\cdot}$ 分别为

$$\tilde{\mu}_{i\cdot}=\left(\sum_{j=1}^bn_{ij}\mu_{ij}\right)\Bigg/n_{i\cdot},\quad\tilde{\mu}_{\cdot j}=\left(\sum_{i=1}^tn_{ij}\mu_{ij}\right)\Bigg/n_{\cdot j}$$

的无偏估计。即，$\bar{y}_{i\cdot\cdot}$ 使用第 $i$ 行单元格样本量作为权重，提供该行单元格均值参数加权平均的无偏估计。类似的，$\bar{y}_{\cdot j \cdot}$ 使用第 $j$ 列单元格样本量作为权重，提供该列单元格均值参数加权平均的无偏估计。因此主效应均值的估计有两种：$\hat{\bar{\mu}}_{i\cdot},\hat{\bar{\mu}}_{\cdot j}$ 和 $\bar{y}_{i\cdot\cdot},\bar{y}_{\cdot j \cdot}$。当使用软件包来分析数据时，确定它计算的是哪种是极其重要的。

对于表 \@ref(tab:table9-1) 中的数据，如表 \@ref(tab:table9-2) 所示，总体边际均值的估计为 $\hat{\bar{\mu}}_{1\cdot}=23,\hat{\bar{\mu}}_{2\cdot}=27,\hat{\bar{\mu}}_{\cdot1}=23,\hat{\bar{\mu}}_{\cdot2}=27,\hat{\bar{\mu}}_{\cdot3}=28$. 这些估计的标准误估计为

$$\begin{aligned}
&\widehat{s.e.}(\hat{\bar{\mu}}_{1\cdot}) =\frac{\sqrt{2}}3\sqrt{\frac13+\frac12+\frac13}=0.51  \\
&\widehat{s.e.}(\hat{\bar{\mu}}_{2\cdot}) =\frac{\sqrt2}3\sqrt{\frac12+\frac13+\frac13}=0.51  \\
&\widehat{s.e.}(\hat{\bar{{\mu}}}_{\cdot1}) =\frac{\sqrt{2}}2\sqrt{\frac13+\frac12}=0.65  \\
&\widehat{s.e.}(\hat{\bar{\mu}}_{\cdot2}) =\frac{\sqrt2}2\sqrt{\frac12+\frac13}=0.65 \\
&\widehat{s.e.}(\hat{\bar{\mu}}_{\cdot3}) =\frac{\sqrt2}2\sqrt{\frac13+\frac13}=0.58
\end{aligned}$$

关于总体边际均值线性组合的推断，如 $\sum_ic_i\bar{\mu}_{i\cdot}$ 或 $\sum_jd_j\bar{\mu}_{\cdot j}$，可以证明

\begin{equation}
\frac{\sum_ic_i\hat{\bar{\mu}}_{i\cdot}-\sum_ic_i\bar{\mu}_{i\cdot}}{\frac{\hat{\sigma}}b\sqrt{\sum_ic_i^2\left(\sum_j\frac1{n_{ij}}\right)}}\thicksim t(v)
(#eq:9-8)
\end{equation}

以及

\begin{equation}
\frac{\sum_jc_j\hat{\bar{\mu}}_{\cdot j}-\sum_jc_j\bar{\mu}_{\cdot j}}{\frac{\hat{\sigma}}t\sqrt{\sum_jc_j^2\left(\sum_i\frac1{n_{ij}}\right)}}\thicksim t(v)
(#eq:9-9)
\end{equation}

方程 \@ref(eq:9-8) 和 \@ref(eq:9-9) 中的公式可以作为方程 \@ref(eq:1-4) 的特殊情况获得。例如，检验 $H_T{:{\bar{\mu}}_{1\cdot}}=\bar{{\mu}}_{2\cdot}$ 的 $t$ 统计量为

$$t_c=\frac{\hat{\bar{\mu}}_{1\cdot}-\hat{\bar{\mu}}_{2\cdot}}{\frac\partial b\sqrt{\sum_j\frac1{n_{1j}}+\sum_j\frac1{n_{2j}}}}=\frac{23-27}{\frac{\sqrt{2}}3\sqrt{\left(\frac13+\frac12+\frac13\right)+\left(\frac12+\frac13+\frac13\right)}}=\frac{-4}{\frac{1.414}3\sqrt{\frac73}}=\frac{-4}{0.72}=-5.55$$

在 $\hat \alpha=0.00024$ 水平下显著。$\hat{\bar{\mu}}_{1\cdot}-\hat{\bar{\mu}}_{2\cdot}$ 的一个 95% 置信区间为

$$\begin{aligned}
\hat{\bar{\mu}}_{\cdot1}-\hat{\bar{\mu}}_{\cdot2}\mp t_{\frac\alpha2,v}\frac{\hat{\sigma}}t\sqrt{\sum_i\frac1{n_{i1}}+\sum_j\frac1{n_{i2}}}& =23-24\mp t_{.025,v}\cdot\frac{1.414}2\sqrt{\left(\frac13+\frac12\right)+\left(\frac12+\frac13\right)}  \\
&=-1\mp(2.228)(0.91)=-1\mp2.03
\end{aligned}$$

## 同时推断与多重比较 {#sec6-5}

在每个处理组合观测次数不等的双向实验中，很少有好的程序可以用于进行多重比较。如果想比较所有双向单元格均值对子，那么第 \@ref(chap3) 章中讨论的任何技术都可以简单地使用，只需将双向实验视为单向处理结构实验。在这种情况下，读者应参阅第 \@ref(sec3-2) 节中给出的建议。如果希望对总体边际均值进行多重比较，建议使用基于方程 \@ref(eq:9-8) 和 \@ref(eq:9-9) 的 $t$ 检验。如果用于比较相应边际均值的相应 $F$ 检验是显著的，则使用给定的显著性水平。如果 $F$ 检验不显著，那么仍然建议使用这些 $t$ 检验。然而，在这种情况下，应该使用 Bonferroni's method，并声称只有当计算出的显著水平小于 $\alpha/p$ 时，两个总体边际均值才显著不同，其中 $\alpha$ 是选定的实验误差率，$p$ 是在收集数据之前考虑的比较次数。如果确定数据中存在交互作用，那么人们可能会想要比较在其他处理的每种可能性下，某一种处理的效果。也就是说，人们可能希望将每一行内的单元格均值相互比较，并将每一列内的单元格均值相互比较。有 $bt(t + 1)/2 + bt(b + 1)/2$ 个 这样的两两比较。如果比较所有均值的 $F$ 检验显著，则可以对所有此类成对比较使用实际显著性水平。这种用法相当于 Fisher's LSD 程序。如果比较所有均值的 $F$ 检验不显着，则可以采用 Bonferroni's method. 对于数据窥探和计划外比较，应该使用 Scheffé's procedure. \@ref(sec9-4) 节中描述的检验统计量可以通过许多统计计算包自动获得。由于这些包采用了效应模型，感兴趣的读者应该参见第 \@ref(sec10-7) 节。

## 结束语 {#sec9-7}

本章是考虑对子类数量不等的双向处理结构进行分析的七章中的第一章。本章所介绍的分析方法是基于均值模型得出的。其中一个重要的假设是，所有处理组合至少被观测过一次。检验主效应和相互作用假设的程序是作为第一章中介绍的一般技术的特例而获得的。本章定义了总体边际均值，并给出了对总体边际均值进行推断的程序。第 \@ref(chap10) 章利用效应模型回答了类似的问题；然而，作者希望所有读者都能对均值模型和效应模型感到同样舒适。

## 练习 {#sec9-8}
