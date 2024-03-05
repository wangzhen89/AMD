
# 方差分量的估计方法 {#chap19}

> "By a small sample, we may judge of the whole piece." - Miguel de Cervantes from Don Quixote

对于一般随机效应模型，估计方差分量有几种方法。当设计均衡（每个单元格样本量相等且没有缺失的单元格）时，其中的一些程序会产生相同的估计，而当设计不均衡时，则会产生不同的估计。本章讨论的四种技术分别是**矩法** (method of moments)、**最大似然** (maximum likelihood, ML)、**受限**或**残差最大似然** (restricted, or residual maximum likelihood, REML) 和 **MIVQUE**. 矩法产生无偏估计，最大似然和 REML 估计是一致的，并且具有最大似然估计通常的大样本量特性，而 MIVQUE 法产生的估计在二次无偏估计类中具有最小方差。当设计均衡且方差分量的解均为正时，矩法、REML 和 MIVQUE 估计是相同的。当设计不均衡时，矩法估计最容易计算，而其他三种方法则需要迭代算法。另一方面，最大似然、REML 和 MIVQUE 法提供的估计比矩法具有更好的性质。REML 通常是估计方差分量的首选方法。

## 矩法 {#sec19-1}

自从 Eisenhart (1947) 将随机效应模型命名为 MODEL II 以来，矩法一直被用于获取得方差分量的估计。许多研究人员在接下来的 20 年里致力于矩法的研究，推导了估计，并开发了检验假设的方法以及构造方差分量置信区间的方法（参见  Searle, 1987; Graybill, 1976; Henderson 1984; Searle et al. 1992; Burdick and Graybill, 1992
等提供了很好的参考文献列表）。在本节中，讨论矩法估计过程的广义版本。

一般随机效应模型可以写为

$$\boldsymbol y=\boldsymbol j_N\boldsymbol\mu+\boldsymbol Z_1\boldsymbol u_1+\boldsymbol Z_2\boldsymbol u_2+\cdots+\boldsymbol Z_r\boldsymbol u_r+\boldsymbol \varepsilon $$

其中

\begin{align}
E(\boldsymbol{u}_i)&=\boldsymbol{0},\quad i=1,2,\ldots,r\\\operatorname{Var}(\boldsymbol{u}_i)&={\sigma}_i^2\boldsymbol{I}_{t_i},\quad i=1,2,\ldots,r\\E(\boldsymbol{\varepsilon})&=\boldsymbol{0},\operatorname{Var}(\boldsymbol{\varepsilon})={\sigma}_{\boldsymbol{\varepsilon}}^2\boldsymbol{I}_N
(#eq:19-1)
\end{align}

$\boldsymbol u_1,\boldsymbol u_2,\cdots,\boldsymbol u_r,\boldsymbol \varepsilon$ 为独立的随机向量。

用于估计方程 \@ref(eq:19-1) 的一般随机效应模型的方差分量的矩法技术包括以下步骤：

1. 计算模型中方差分量的平方和及其相应的均方。
2. 根据方差分量评估每个均方的期望；这些期望必须不涉及 $\boldsymbol \mu$（或任何其他固定效应参数），并且每个方差分量必须包括在至少一个均方的期望中。
3. 令均方的期望等于均方的观测值，从而生成方差分量的线性方程组（用方程组中的方差分量解替换方差分量参数）。
4. 求解所得方程组以获得每个方差分量的估计。

矩法的一个问题是方差分量的一些估计可能具有负值。当方差分量的解为负时，将方差分量的估计设置为零（将估计保留在参数空间中）。

当模型 \@ref(eq:19-1) 的随机效应 $\boldsymbol u_1,\boldsymbol u_2,\cdots,\boldsymbol u_r,\boldsymbol \varepsilon$ 联合独立且服从正态分布时，并且当平方和彼此独立分布时，所得方差分量的估计为**最小方差无偏的** (minimum variance unbiased)。如果 $\boldsymbol u_1,\boldsymbol u_2,\cdots,\boldsymbol u_r,\boldsymbol \varepsilon$ 具有与正态分布相同的前四阶矩，则估计是**最小方差二次无偏的** (minimum variance quadratic unbiased) (Graybill, 1976, p. 632). 矩法技术不需要为了获得估计而假设正态性。在不假设正态性，或不假设随机向量的分布具有与正态分布相同前四阶矩的情况下，这些估计所拥有的唯一已知属性是：它们是无偏的。然而，当解为负时将估计设置为零的过程意味着估计不再是无偏的。

矩法的关键是确定如何计算平方和，然后评估所得均方期望。这些主题在第 \@ref(chap18) 章中进行了讨论。

如果模型具有 r+1 个方差分量 $\sigma_\varepsilon^2,\sigma_1^2,\sigma_2^2,\ldots,\sigma_r^2$，则需要 r+1 个平方和（或均方）。设 $Q_0=\boldsymbol y^\prime \boldsymbol A_0 \boldsymbol y,Q_1=\boldsymbol y^\prime \boldsymbol A_1 \boldsymbol y,\cdots,Q_r=\boldsymbol y^\prime \boldsymbol A_r \boldsymbol y$，表示具有各自期望值的平方和，这些期望值为

$$\begin{aligned}E(Q_i)=b_{i0}\sigma_\varepsilon^2+b_{i1}\sigma_1^2+b_{i2}\sigma_\gamma^2+\cdots+b_{ir}\sigma_r^2,\quad i=0,1,2,\ldots,r\end{aligned}$$

令每个平方和等于其期望值，在方差上加上波浪号（~），将解表示为

$$Q_i=b_{i0}\tilde{\sigma}_\varepsilon^2+b_{i1}\tilde{\sigma}_1^2+b_{i2}\tilde{\sigma}_2^2+\cdots+b_{ir}\tilde{\sigma}_r^2,\quad i=0,1,2,\ldots,r$$

或以矩阵表示法

$$\begin{bmatrix}Q_0\\Q_1\\Q_2\\\vdots\\Q_r\end{bmatrix}=\begin{bmatrix}b_{00}&b_{01}&b_{02}&\cdots&b_{0r}\\b_{10}&b_{11}&b_{12}&\cdots&b_{1r}\\b_{20}&b_{21}&b_{22}&\cdots&b_{2r}\\\vdots&\vdots&\vdots&\ddots&\vdots\\b_{r0}&b_{r1}&b_{r2}&\cdots&b_{rr}\end{bmatrix}\begin{bmatrix}\tilde{\sigma}_\varepsilon^2\\\tilde{\sigma}_1^2\\\tilde{\sigma}_2^2\\\vdots\\\tilde{\sigma}_r^2\end{bmatrix}$$

或 $\boldsymbol Q=\boldsymbol B \tilde \sigma^2$.

如果矩阵 $\boldsymbol B$ 的秩为 r+1，那么所有的方差分量都是可估的。如果 $\boldsymbol B$ 的秩小于 r+1，则不是所有的方差分量都是可估的，并且只有方差分量的一些线性组合是可估的。假设 $\boldsymbol B$ 的秩为 r+1，则方程组的解为

$$\tilde{\sigma}^2=\boldsymbol B^{-1}\boldsymbol Q{=}\boldsymbol C\boldsymbol Q\mathrm{~(say)}$$

在不将值限制在参数空间的情况下获得该解；也就是说，有些解可能是负的。$\sigma_i^2$ 的解由 $\tilde \sigma_i^2$ 表示，估计由 $\hat \sigma_i^2$ 表示，其中

$$\hat{\sigma}_i^2=\begin{cases}\tilde{\sigma}_i^2&\mathrm{~if~}\tilde{\sigma}_i^2>0\\0&\mathrm{~if~}\tilde{\sigma}_i^2\leq0&\end{cases}\quad i=0,1,2,\ldots,r$$

在许多计算平方和的模型和方法中，矩阵 $\boldsymbol B$ 是三角阵，因此可以在不求 $\boldsymbol B$ 的逆的情况下获得解。

每个解都是观测到的平方和 $Q_0,Q_1,Q_2,\cdots,Q_r$ 的线性组合

$$\tilde{\sigma}_i^2=c_{i0}Q_0+c_{i1}Q_1+c_{i2}Q_2+\cdots+c_{ir}Q_r\quad i=0,1,2,\ldots,r$$

其中 $c_i'=[c_{i0}\quad c_{i1}\quad c_{i2}\quad \cdots\quad  c_{ir}]$ 是 $\boldsymbol C=\boldsymbol B^{-1}$ 的第 $i$ 行。$\tilde \sigma_i^2$ 的方差为

$$\mathrm{Var}(\tilde{\sigma}_i^2)=\mathrm{Var}(c_{i0}Q_0+c_{i1}Q_1+c_{i2}Q_2+\cdots+c_{ir}Q_r)$$

当 $Q_i,i=0,1,\cdots,r$ 不相关时，$\tilde \sigma_i^2$ 的方差为

$$\operatorname{Var}(\tilde{\sigma}_i^2)=c_{i0}^2\operatorname{Var}(Q_0)+c_{i1}^2\operatorname{Var}(Q_1)+c_{i2}^2\operatorname{Var}(Q_2)+\cdots+c_{ir}^2\operatorname{Var}(Q_r)$$

Searle (1971, Chapter 11) 总结了几种模型的矩法估计及其方差。对于大多数均衡模型（假设随机变量分布的矩对应于正态分布的前四阶矩），矩法估计是方差分量的**一致最小方差无偏估计** (uniformly minimum variance unbiased estimators) (Graybill, 1976). 因此，对于几乎均衡的模型，矩法估计应该具有相当好的性质。

### 应用。示例 19.1：不均衡单向模型 {#sec19-1-1}

示例 [18.2](#sec18-2-1) 的不均衡单向随机效应模型为

$$y_{ij}=\mu+u_i+\varepsilon_{ij}\quad i=1,2,\ldots,t\mathrm{~and~}j=1,2,\ldots,n_i$$

其中 $u_i$ 不相关且均值为 0 方差为 $\sigma^2_1$，$\varepsilon_{ij}$ 不相关且均值为 0 方差为 $\sigma^2_{\varepsilon}$，并且 $u_i$ 与 $\varepsilon_{ij}$ 不相关。可以使用的两个平方和是组内平方和 (sum of squares within) $Q_0$ 或 SSW ，以及组间平方和 (sums of squares between) $Q_1$ 或 SSB，其中

$$\begin{aligned}Q_0&=\sum_{i=1}^t\sum_{j=1}^{n_i}(y_{ij}-\bar{y}_{i\cdot})^2=\sum_{i=1}^t\sum_{j=1}^{n_i}y_{ij}^2-\sum_{i=1}^tn_i\bar{y}_{i\cdot}^2=SSW\\\\\\Q_1&=\sum_{i=1}^tn_i(\bar{y}_{i\cdot}-\bar{y}_{\cdot\cdot})^2=\sum_{i=1}^tn_i\bar{y}_{i\cdot}^2-\left(\sum_{i=1}^tn_i\right)\bar{y}_{\cdot}^2=SSB\end{aligned}$$

$Q_0,Q_1$ 的期望值在第 \@ref(chap18) 章中评估为

$$\begin{aligned}E(Q_0)&=(N-t){\sigma}_\varepsilon^2\\\\E(Q_1)&=(t-1){\sigma}_\varepsilon^2+\left(N-\frac{{\sum}_{i=1}^t{n}_i^2}N\right){\sigma}_1^2\end{aligned}$$

其中 $N=\sum_{i=1}^tn_i$.

通过令平方和等于其期望值而获得的方程为

$$\begin{aligned}&Q_0=(N-t)\tilde{{\sigma}}_\varepsilon^2\\\\&Q_1=(t-1)\tilde{{\sigma}}_\varepsilon^2+\left(N-\frac{\sum_{i=1}^tn_i^2}N\right)\tilde{{\sigma}}_1^2\end{aligned}$$

或以矩阵表示法

$$\begin{bmatrix}Q_0\\Q_1\end{bmatrix}=\begin{bmatrix}N-t&0\\(t-1)&N-\frac{\sum_{i=1}^tn_i^2}N\end{bmatrix}\begin{bmatrix}\tilde{\sigma}_\varepsilon^2\\\tilde{\sigma}_1^2\end{bmatrix}$$

方程组也可以通过令观测到的均方等于它们的期望来生成，因为均方是通过将每个方程除以其相应的自由度来获得的。解是相同的。所得到的涉及均方的待解方程组为

$$\begin{aligned}\frac{Q_0}{N-t}&=\tilde{\sigma}_\varepsilon^2\\\frac{Q_1}{t-1}&=\tilde{\sigma}_\varepsilon^2+\frac{\left(N-\frac{\sum_{i=1}^tn_i^2}N\right)}{t-1}\tilde{\sigma}_1^2\end{aligned}$$

此方程组的解为

$$\begin{aligned}\widetilde{\sigma}_\varepsilon^2&=\frac{Q_0}{N-t}\\\tilde{\sigma}_1^2&=\frac{Q_1-(t-1)\tilde{\sigma}_\varepsilon^2}{N-\frac{\sum_{i=1}^tn_i^2}N}\end{aligned}$$

矩法估计为

$$\hat{\sigma}_\varepsilon^2=\tilde{\sigma}_\varepsilon^2$$

且

$$\hat{{\sigma}}_1^2=\begin{cases}\tilde{{\sigma}}_1^2&\mathrm{~if~}\tilde{{\sigma}}_1^2>0\\0&\mathrm{~if~}\tilde{{\sigma}}_1^2\leq0&\end{cases}$$

### 示例 19.2：单向随机效应模型中的小麦品种 {#sec19-1-2}

实验者从一批小麦品种中随机选取了四个品种，进行了一项实验，以评估小麦植株在抽穗前受到昆虫破坏的程度。实验设计结构为完全随机设计，每个品种有四个重复或区（区为实验单元）。由于环境条件的影响，一些区 (plots) 被毁坏（因雨水过多而淹没）。在小麦植株开始抽穗的前一天，实验者从每个区中随机选取了 20 株植物，并使用 0 到 10 的评分量表评估每株植物受到的昆虫破坏程度，其中 0 表示无破坏，10 表示严重破坏。因此，每个区上测得的响应是 20 株植物评分的平均值。数据如表 \@ref(tab:table19-1) 所示。计算平方和 $Q_0$ 和 $Q_1$ 所需的计算，以及由此产生的平方和、均方、期望均方、方程组、求解结果和方差分量估计均列在表 \@ref(tab:table19-1) 中。从方差分量的估计获得的信息是，一个品种内区之间的方差约为 0.056，而品种总体的方差约为 0.067. 将随机选择的品种种植到随机选择的区的方差是两个方差分量的总和，即

$$\hat{\sigma}_{\mathrm{Damage}}^2=\hat{\sigma}_{\varepsilon}^2+\hat{\sigma}_{\mathrm{Var}}^2=0.056+0.067=0.123$$

<table>
<caption>(\#tab:table19-1)示例 19.2 小麦品种虫害的数据和计算</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2019.1.png">
</td>
  </tr>
</tbody>
</table>

组内相关性 (intraclass correlation) 的估计为

$$\hat{\rho}=\frac{\hat{\sigma}_\mathrm{Var}^2}{\hat{\sigma}_\varepsilon^2+\hat{\sigma}_\mathrm{Var}^2}=\frac{0.067}{0.123}=0.545$$

当用非单向或非完全随机设计结构的处理结构进行实验时，没有普遍接受的技术来获得平方和，从而得出方差分量的估计。第 \@ref(chap18) 章中提出的计算平方和的方法用于估计双向随机效应模型的方差分量。

### 示例 19.3：表 18.2 中的双向设计数据 {#sec19-1-3}

表 \@ref(tab:table19-2) 中的数据是图 \@ref(fig:figure18-1) 的观测，其中通过综合法评估了几种类型的平方和的期望值。值 $Q_0,Q_1,Q_2$ 和 $Q_3$ 对应于 Henderson's method I 的平方和为

$$\begin{aligned}
&Q_0 =SSError=30.6666  \\
&Q_1 =SSA=\sum_{i=1}^2\frac{y_{i\cdot\cdot}^2}{n_{i\cdot}}-\frac{y_{\cdot\cdot\cdot}^2}{n_{\cdot\cdot}}=0.6428  \\
&Q_2 =SSB=\sum_{j=1}^3\frac{y_{\cdot j\cdot}^2}{n_{\cdot j}}-\frac{y_{\cdot\cdot\cdot}^2}{n_{\cdot\cdot}}=15.2143  \\
&Q_3 =SSAB=\sum_{i=1}^2\sum_{j=1}^3\frac{y_{ij\cdot}^2}{n_{ij}}-\sum_{i=1}^2\frac{y_{i\cdot\cdot}^2}{n_{i\cdot}}-\sum_{j=1}^3\frac{y_{\cdot j\cdot}^2}{n_{\cdot j}}+\frac{y_{\cdot\cdot\cdot}^2}{n_{\cdot\cdot}}=108.6905 
\end{aligned}$$

<table>
<caption>(\#tab:table19-2)示例 19.3 的双向随机效应处理结构数据</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2019.2.png">
</td>
  </tr>
</tbody>
</table>

这些平方和的期望值是：

$$\begin{aligned}
&E(Q_0)=8{\sigma}_\varepsilon^2 \\
&E(Q_{1}) =\sigma_\varepsilon^2+0.1429\sigma_b^2+2.4286\sigma_c^2+7.0\sigma_a^2  \\
&E(Q_2) =2\sigma_\varepsilon^2+9.286\sigma_b^2+4.77\sigma_c^2+0.20\sigma_a^2  \\
&E(Q_3) =2\sigma_\varepsilon^2+6.37\sigma_c^2 
\end{aligned}$$

令平方和的值等于其各自的预期值，从而得到以下方程组

$$\begin{aligned}
\text{30.6666}& =8\tilde{\sigma}_\varepsilon^2  \\
\text{0.6428}& =\tilde{{\sigma}}_\varepsilon^2+0.1429\tilde{{\sigma}}_b^2+2.4286\tilde{{\sigma}}_c^2+7.0\tilde{{\sigma}}_a^2  \\
\text{15.2143}& =2\tilde{{\sigma}}_\varepsilon^2+9.286\tilde{{\sigma}}_b^2+4.77\tilde{{\sigma}}_c^2+0.20\tilde{{\sigma}}_a^2  \\
\text{108.6905}& =2\tilde{{\sigma}}_\varepsilon^2+6.37\tilde{{\sigma}}_c^2 
\end{aligned}$$

解为

$$\begin{aligned}
&\tilde{{\sigma}}_\varepsilon^2 =3.83325  \\
&\tilde{{\sigma}}_c^2 =15.8593  \\
&\tilde{{\sigma}}_b^2 =-10.8817  \\
&\tilde{{\sigma}}_a^2 =-8.2523 
\end{aligned}$$

提供方差分量的无偏估计。因为上面的一些值是负值，所以方差分量的矩法估计的最终被认为是

$$\begin{aligned}
&\hat{{\sigma}}_\varepsilon^2 =3.83325  \\
&\hat{{\sigma}_c^2} \text{=15.8593}  \\
&\hat{{\sigma}_b^2} =0.00  \\
&\hat{{\sigma}_a^2} =0.00 
\end{aligned}$$

示例 [19.3]({#sec19-1-3}) 的结果指出了使用矩法技术时经常遇到的一个问题：它可以产生方差分量的负解，这对于方差分量是不可接受的。以下列出了通过拟合常数法或 Henderson's method III 平方和（SAS^®^ I 型）和SAS III 型平方和（这些平方和的期望值是通过综合法评估的，并在表 \@ref(tab:table19-4) 中给出）生成的方程组求解得出的示例 [19.3]({#sec19-1-3}) 的矩法解和相应的估计值。

I 型解为

$$\begin{aligned}
&\tilde{{\sigma}}_\varepsilon^2 =3.8333  \\
&\tilde{\sigma}_c^2 =\frac{109.1451-2(3.8333)}{4.5178}=22.4620  \\
&\tilde{{\sigma}}_b^2 =\frac{14.5797-2(3.8333)-4.6252(22.4620)}{9.1429}=-10.5877  \\
&\tilde{\sigma}_{a}^{2} =\frac{0.6429-3.833-0.1428(-10.71322)-2.4284(22.4620)}{7.00}=-8.0329 
\end{aligned}$$

以及 I 型估计为

$$\begin{aligned}
&\hat{{\sigma}}_\varepsilon^2 =3.8333  \\
&\hat{\sigma}_c^2 =22.4620  \\
&\hat{{\sigma}_b^2} =0  \\
&\hat{{\sigma}_a^2} =0 
\end{aligned}$$

III 型解为

$$\begin{aligned}
&\tilde{{\sigma}}_\varepsilon^2 =3.8333  \\
&\tilde{{\sigma}}_c^2 =\frac{109.1451-2(3.8333)}{4.5178}=22.4620  \\
&\tilde{{\sigma}}_b^2 =\frac{8.9098-2(3.8333)-4.5178(22.4620)}{9.0353}=-11.1080  \\
&\tilde{{\sigma}}_a^2 =\frac{0.677-3.833-2.2500(22.4620)}{6.75}=-8.0305 
\end{aligned}$$

[^type3solution]: 原文为 “the type III solutions are”，应为笔误。

III 型估计[^type3solution]为

$$\begin{aligned}
&\hat{{\sigma}}_\varepsilon^2 =3.8333  \\
&\hat{{\sigma}_c^2} =22.4620  \\
&\hat{{\sigma}_b^2} =0  \\
&\hat{\sigma}_a^2 =0 
\end{aligned}$$

当方差分量的解为负时，标准过程是将相应的估计值设置为零。为了证明这一过程的结果，单向随机效应模型的组内均方 (mean square within) 和组间均方 (mean square between) 的期望均方可以表示为

$$\begin{aligned}E(MSWithin)&=\sigma_\varepsilon^2\\E(MSBetween)&=\sigma_\varepsilon^2+c\sigma_u^2\end{aligned}$$

$\sigma^2_u$ 的解为

$$\tilde{\sigma}_u^2=\frac{MSBetween-MSWithin}c$$

在正态性假设下，MSWithin 和 MSBetween 是独立的随机变量。如果 $\sigma^2_u=0$，那么两个期望均方都等于 $\sigma^2_\varepsilon=0$，因此，$\sigma^2_u$ 的解为负的概率约为 0.50，这取决于自由度。如果分子和分母的自由度相等，那么

$$P(MSBetween<MSWithin\mid\sigma_u^2=0)=0.50$$

随着 $\sigma^2_u$ 变大，得到负解的概率降低。因此，当解为负时，将 $\hat\sigma^2_u$ 设置为零是合理的（更详细的讨论请参见 Searle et al. 1992）。第 \@ref(chap20) 章讨论了由矩法估计构造的置信区间和假设检验。

## 最大似然 {#sec19-2}

在统计学中，估计分布参数最常用的技术是最大似然法 (method of maximum likelihood). 该过程使用假设的观测分布，并构造似然函数，该函数是数据和未知模型参数的函数。最大似然估计是在参数空间中使似然函数值最大化的参数值。在实践中最大化似然函数的 log~e~. 等价地，可以确定参数空间中最小化 -2log~e~似然函数的参数值。一般随机效应模型 \@ref(eq:19-1) 的参数空间为

$$\{-\infty<\mu<+\infty,\quad0<\sigma_i^2<+\infty,\quad i=1,2,\ldots,k;\quad0<\sigma_\varepsilon^2<\infty\}$$

对于一般随机效应模型 \@ref(eq:19-1)，观测向量的分布为

$$
\boldsymbol y\thicksim N(\boldsymbol j_n\boldsymbol\mu,\sigma_\varepsilon^2\boldsymbol I_n+\sigma_1^2\boldsymbol Z_1\boldsymbol Z_1^{\prime}+\sigma_2^2\boldsymbol Z_2\boldsymbol Z_2^{\prime}+\cdots+\sigma_k^2\boldsymbol Z_k\boldsymbol Z_k^{\prime})\quad\mathrm{or}\quad y\thicksim N (\boldsymbol j_n\boldsymbol \mu,\boldsymbol \Sigma)
$$

观测值的似然函数为

$$L(\boldsymbol\mu,\sigma_\varepsilon^2,\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2|y)=(2\pi)^{-n/2}|\boldsymbol{\Sigma}|^{-1/2}\exp[-\frac12(\boldsymbol{y}-\boldsymbol{j}_n\boldsymbol\mu)^{\prime}\boldsymbol{\Sigma}^{-1}(\boldsymbol{y}-\boldsymbol{j}_n\boldsymbol\mu)]$$

-2log~e~似然函数为

$$\begin{aligned}\ell(\boldsymbol\mu,\sigma_e^2,\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2|y)&=-2\log_e[L(\boldsymbol\mu,\sigma_e^2,\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2|y)]\\&=n\log_e(2\pi)+\log_e(|\boldsymbol{\Sigma}|)+(\boldsymbol{y}-\boldsymbol{j}_n\boldsymbol{\mu})^{\prime}\boldsymbol{\Sigma}^{-1}(\boldsymbol{y}-\boldsymbol{j}_n\boldsymbol{\mu})\end{aligned}$$

在参数空间上最小化 $\ell(\boldsymbol\mu,\sigma_e^2,\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2|y)$ 的过程通常需要利用似然方程的迭代过程，其中涉及到关于模型 $\ell(\boldsymbol\mu,\sigma_e^2,\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2|y)$ 每个参数的一阶导数或一阶和二阶导数。

当数据来自均衡设计时（相等的 $n$ 以及无缺失单元格），通过将似然函数的一阶导数等于零而生成的似然方程组通常可以显式求解。由于方差分量解中的一些值可能为负值，因此所得到的解不受参数空间的限制。对于某些均衡模型，可以证明，当似然方程中 $\sigma^2_i$ 的解为负时，$\sigma^2_i$ 的最大似然估计为 $\hat\sigma^2_i=0$ (Searle, 1971). 对于不均衡设计，需要迭代技术，其中估计过程应将方差分量的估计值限制在参数空间中。

### 示例 19.4：均衡单向模型的最大似然解 {#sec19-2-1}

均衡单向随机效应模型可表示为

$$y_{ij}=\mu+u_i+\varepsilon_{ij},\quad i=1,2,\ldots,t\quad\mathrm{and~}j=1,2,\ldots,n$$

其中 $\mu_i\sim i.i.d.\,N(0,\sigma^2_u),\varepsilon_{ij}\sim i.i.d.\,N(0,\sigma^2_\varepsilon)$，且 $\mu_i,\varepsilon_{ij}$ 独立，或 $\boldsymbol y\sim N(\boldsymbol j_n\otimes \boldsymbol j_t)\mu,\boldsymbol \Sigma]$ 其中 $\boldsymbol \Sigma=\sigma^2_u \boldsymbol j_n\otimes \boldsymbol j_t   +\sigma^2_\varepsilon \boldsymbol I_n\otimes \boldsymbol I_t$.

[^kronecker]: 设 $\boldsymbol A, \boldsymbol B$ 分别为 m × n 和 p × q 阶矩阵，则 $\boldsymbol A \otimes \boldsymbol B=\begin{bmatrix}a_{11}\boldsymbol B&\cdots&a_{1n}\boldsymbol B\\\vdots&\ddots&\vdots\\a_{m1}\boldsymbol B&\cdots&a_{mn}\boldsymbol B\end{bmatrix}$. Knronecker 积 "$\otimes$" 可作用于两任意大小的矩阵，那么也包括向量。

符号 $A \otimes B$ 表示 $A,B$ 两矩的 Kronecker 积[^kronecker]或直积 (direct product)  (Graybill, 1976). 协方差阵 $\boldsymbol \Sigma$ 可以表示为

$$\boldsymbol{\Sigma}={\sigma}_\varepsilon^2\left[\left(\boldsymbol{I}_n-\frac1n\boldsymbol{J}_n\right)\otimes\boldsymbol{I}_t\right]+({\sigma}_\varepsilon^2+n{\sigma}_u^2)\left[\left(\frac1n\boldsymbol{J}_n\right)\otimes\boldsymbol{I}_t\right]$$

其中 ${\sigma}_\varepsilon^2$ 和 ${\sigma}_\varepsilon^2+n{\sigma}_u^2$ 分别为 $\boldsymbol \Sigma$ 和 $\left[\left(\boldsymbol{I}_n-\frac1n\boldsymbol{J}_n\right)\otimes\boldsymbol{I}_t\right]$ 的特征根，以及 $\left[\left(\frac1n\boldsymbol{J}_n\right)\otimes\boldsymbol{I}_t\right]$ 为正交幂等矩阵。

协方差阵的逆为

$$\boldsymbol{\Sigma}^{-1}=\frac1{{\sigma}_\varepsilon^2}{\left[\left(\boldsymbol{I}_n-\frac1n\boldsymbol{J}_n\right)\otimes\boldsymbol{I}_t\right]}+\frac1{{\sigma}_\varepsilon^2+n{\sigma}_u^2}{\left[\left(\frac1n\boldsymbol{J}_n\right)\otimes\boldsymbol{I}_t\right]}$$

使用协方差的逆矩阵这种表示，可以获得以下结果

$$|\boldsymbol{\Sigma}|=(\sigma_\varepsilon^2)^{t(n-1)}(\sigma_\varepsilon^2+n\sigma_u^2)^t$$

或

$$\log_\mathrm{e}|\boldsymbol{\Sigma}|=t(n-1)\log_\mathrm{e}(\sigma_\varepsilon^2)+t\log_\mathrm{e}(\sigma_\varepsilon^2+n\sigma_u^2)$$

以及

$$(\boldsymbol y-\boldsymbol j_{nt}\mu)^{\prime}\boldsymbol{\Sigma}^{-1}(\boldsymbol y-\boldsymbol j_{nt}\mu)=\frac{nt(\bar{y}_{\cdot\cdot}-\mu)^2}{\sigma_\varepsilon^2+n\sigma_u^2}+\frac{SSE}{\sigma_\varepsilon^2}+\frac{SSU}{\sigma_\varepsilon^2+n\sigma_u^2}$$

其中

$$SSE=\sum_{i=1}^t\sum_{j=1}^n(y_{ij}-\bar{y}_{i\cdot})^2\quad\mathrm{~and~}\quad SSU=n\sum_{i=1}^t(\bar{y}_{i\cdot}-\bar{y}_{\cdot\cdot})^2$$

使用这些表达式，$-2\log_{\mathrm{e}}(L(\mu,\sigma_{u}^2,\sigma_{\varepsilon}^2|\boldsymbol y)$ 可以写为

$$\begin{aligned}
\ell(\mu,{\sigma}_u^2,{\sigma}_\varepsilon^2|\boldsymbol{y}) =tn\log_\mathrm{e}(2\pi)+t(n-1)\log_\mathrm{e}(\sigma_\varepsilon^2)+t\log_\mathrm{e}(\sigma_\varepsilon^2+n\sigma_u^2) \\
+\frac{nt(\bar{y}_{\cdot\cdot}-\mu)^2}{\sigma_\varepsilon^2+n\sigma_u^2}+\frac{SSE}{\sigma_\varepsilon^2}+\frac{SSU}{\sigma_\varepsilon^2+n\sigma_u^2} 
\end{aligned}$$

通过关于三个参数 $\mu,\sigma_u^2,\sigma_\varepsilon^2$ 对 $\ell(\mu,\sigma_u^2,\sigma_\varepsilon^2|\boldsymbol y)$ 求导，然后将导数设置为零来获得似然方程。将导数设置为零，在求解参数时评估的导数为

$$\begin{aligned}
&\frac{\partial\ell(\mu,\sigma_u^2,\sigma_\varepsilon^2|\boldsymbol y)}{\partial\mu} =\frac{-2nt(\bar{y}_{\cdot\cdot}-\tilde{\mu})}{\sigma_\varepsilon^2+n\sigma_u^2}=0  \\
&\frac{\partial\ell(\mu,\sigma_u^2,\sigma_\varepsilon^2|\boldsymbol y)}{\partial\sigma_\varepsilon^2} =\frac{t(n-1)}{\tilde{\sigma}_\varepsilon^2}+\frac t{\tilde{\sigma}_\varepsilon^2+n\tilde{\sigma}_u^2}-\frac{nt(\bar{y}_{\cdot\cdot}-\tilde{\mu})^2}{(\tilde{\sigma}_\varepsilon^2+n\tilde{\sigma}_u^2)^2}-\frac{SSE}{(\tilde{\sigma}_\varepsilon^2)^2}-\frac{SSU}{(\tilde{\sigma}_\varepsilon^2+n\tilde{\sigma}_u^2)^2}=0  \\
&\frac{\partial\ell(\mu,\sigma_u^2,\sigma_\varepsilon^2|\boldsymbol y)}{\partial\sigma_u^2} =\frac{nt}{\tilde{\sigma}_\varepsilon^2+n\tilde{\sigma}_u^2}-\frac{n^2t(\bar{y}_{\cdot\cdot}-\tilde{\mu})^2}{(\tilde{\sigma}_\varepsilon^2+n\tilde{\sigma}_u^2)^2}-\frac{nSSU}{(\tilde{\sigma}_\varepsilon^2+n\tilde{\sigma}_u^2)^2}=0 
\end{aligned}$$

最大似然方程的解为

$$\tilde{\mu}=\bar{y}_{\cdot\cdot    },\quad\tilde{\sigma}_\varepsilon^2=\frac{SSE}{t(n-1)}=MSError,\quad\mathrm{and}\quad\tilde{\sigma}_u^2=\frac1n{\left[\frac{SSU}t-MSError\right]}$$

因此，最大似然估计为

$$\hat{\mu}=\tilde{\mu}=\bar{y}_{\cdot\cdot    },\quad\hat{\sigma}_\varepsilon^2=\tilde{\sigma}_\varepsilon^2=\frac{SSE}{t(n-1)}=MSError$$

以及

$$\hat{\sigma}_u^2=\begin{cases}\tilde{\sigma}_u^2&\mathrm{~if~}&\tilde{\sigma}_u^2\geq0\\0&\mathrm{~if~}&\tilde{\sigma}_u^2<0&\end{cases}$$

当 ${\sigma}_u^2$ 的估计为零时，则通过合并 SSE 和 SSU 及其自由度来重新计算 ${\sigma}_\varepsilon^2$ 的估计，从而获得

$$\hat{\sigma}_\varepsilon^2=\frac{SSE+SSU}{tn-1}$$

如果恰好存在负的组内相关性，那么通过合并获得的 ${\sigma}_\varepsilon^2$ 的估计将低估方差。必须对假设及其适当性进行仔细研究，因为如果 $u_i$ 水平内的实验单元之间存在一定程度的竞争，则负相关将是合适的。在这种情况下，协方差阵可以表示为

$$\boldsymbol{\Sigma}={\sigma}_A^2{\rho}[c{J}_n\otimes\boldsymbol{I}_t]+{\sigma}_A^2(1-{\rho})[\boldsymbol{I}_n\otimes\boldsymbol{I}_t]$$

其中 $\sigma_{\varepsilon}^{2}=\sigma_{A}^{2}\left(1-\rho\right)$ 以及 $\sigma_{u}^{2}=\sigma_{A}^{2}\rho$.

已经开发了几种计算算法来最大化似然函数，从而提供模型参数的最大似然估计（参数空间中的值）(Hemmerle and Hartley, 1973; Corbeil and Searle, 1976).  最大似然估计的大样本量方差可以通过求二阶导数矩阵的逆来获得，其中二阶导数在最大似然估计处。Searle (1971) 以及 Searle et al. (1992) 报道了几个设计的实验的最大似然估计及其方差。

当使用计算机软件来构造这些模型时，实验者应彻底研究所使用的算法，并确定其性质，即它是否总是产生有意义的估计，以及它是否在参数空间上使似然函数最大化。

通过使用 SAS-Mixed，我们获得了示例 [19.2](#sec19-1-2) 和示例 [19.3](#sec19-1-3) 中方差分量的最大似然估计。示例 [19.2](#sec19-1-2) 中 $\sigma^2_\varepsilon$ 和 $\sigma^2_u$ 的最大似然估计值分别为 $\hat\sigma^2_\varepsilon= 0.05749$ 和 $\sigma^2_u = 0.04855$，如表 \@ref(tab:table19-3) 所示。示例 [19.3](#sec19-1-3) 中模型参数的最大似然估计为 $\hat\sigma^2_\varepsilon= 3.5989,\hat\sigma^2_c= 8.5604,\hat\sigma^2_a = 0$ 以及 $\hat\sigma^2_b = 0$，如表 \@ref(tab:table19-4) 所示。SAS-Mixed 中的 ML 算法在参数空间上进行最大化，如表 \@ref(tab:table19-4) 中 $\hat\sigma^2_a,\hat\sigma^2_b$ 设置为零。

<table>
<caption>(\#tab:table19-3)Proc Mixed 代码来计算示例 19.2 中数据的方差分量和均值的最大似然估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2019.3.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table19-4)Proc Mixed 代码来获得示例 19.3 中双向随机效应模型方差分量的最大似然估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2019.4.png">
</td>
  </tr>
</tbody>
</table>

## 受限或残差最大似然估计 {#sec19-3}

通过最大化似然函数中不包括任何混合效应的部分，或者通过在从模型中去除固定效应之后最大化残差的似然函数，来获得受限或残差最大似然估计 (Restricted or residual maximum likelihood estimates, REML). 对于本章中的模型，只有一个固定效应参数 $\mu$. 这也相当于在给定的总体样本均值下查看一组平方和的条件分布。该过程通过将似然函数分解为两部分来完成，其中一部分涉及固定效应参数，另一部分仅涉及方差分量。通过关于方差分量对 -2log~e~残差似然函数进行微分，并将其设置为零，从而得到 REML 方程。

对于一般随机效应模型 \@ref(eq:19-1)，观测值的似然函数为

$$L(\mu,\sigma_\varepsilon^2,\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2|\boldsymbol y)=(2\pi)^{-n/2}\left|\boldsymbol{\Sigma}\right|^{-1/2}e^{\left[-(1/2)(\boldsymbol{y}-\boldsymbol{j}_n{\mu})^{\prime}\boldsymbol{\Sigma}^{-1}(\boldsymbol{y}-\boldsymbol{j}_n{\mu})\right]}$$

并且 -2log~e~似然函数为

$$\begin{aligned}
\ell(\mu,\sigma_\varepsilon^2,\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2|\boldsymbol y) =&-2\log_e[L(\mu,\sigma_\varepsilon^2,\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2|\boldsymbol y)]  \\
=&\,n\log_\mathrm{e}(2\pi)+\log_\mathrm{e}(|\boldsymbol{\Sigma}|)+(\boldsymbol{y}-\boldsymbol{j}_n{\mu})^{\prime}\boldsymbol{\Sigma}^{-1}(\boldsymbol{y}-\boldsymbol{j}_n{\mu}) \\
=&\,\ell(\mu,\sigma_\varepsilon^2,\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2|\bar{y}) \\
&+\ell(\sigma_\varepsilon^2,\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2| SSE,SSU_1,SSU_2,\ldots,SSU_k)
\end{aligned}$$

其中 $SSE,SSU_1,SSU_2,\ldots,SSU_k$ 表示不依赖于 $\mu$ 的独立平方和的集合，$\ell(\sigma_\varepsilon^2,\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2| SSE,SSU_1,SSU_2,\ldots,SSU_k)$ 是残差似然函数。残差似然方程的解提供了方差分量的 REML 估计。

### 示例 19.5：均衡单向模型的 REML 解 {#sec19-3-1}

使用示例 [19.4](#sec19-2-1) 中描述的均衡单向随机效应模型，-2log~e~似然函数可以表示为充分统计量 $\bar y_{\cdot\cdot},SSE,SSU$ 的函数

$$\begin{aligned}
\ell(\mu,\sigma_u^2,\sigma_\varepsilon^2|\boldsymbol y) =&\,tn\log_{\mathrm{e}}(2\pi)+t(n-1)\log_{\mathrm{e}}(\sigma_{\varepsilon}^2)+t\log_{\mathrm{e}}(\sigma_{\varepsilon}^2+n\sigma_{u}^2)  \\
&+\frac{nt(\bar{y}_{\cdot\cdot}-\mu)^2}{\sigma_\varepsilon^2+n\sigma_u^2}+\frac{SSE}{\sigma_\varepsilon^2}+\frac{SSU}{\sigma_\varepsilon^2+n\sigma_u^2} \\
=&\left[\frac{nt(\bar{y}_{\cdot\cdot}-\mu)^2}{\sigma_\varepsilon^2+n\sigma_u^2}+\log_{\mathrm{e}}(2\pi)+\log_{\mathrm{e}}(\sigma_\varepsilon^2+n\sigma_u^2)\right] \\
&+\left[(tn-1)\log_{\mathrm{e}}(2\pi)+t(n-1)\log_e(\sigma_e^2)+(t-1)\log_{\mathrm{e}}(\sigma_\varepsilon^2+n\sigma_u^2)+\frac{SSE}{\sigma_\varepsilon^2}+\frac{SSU}{\sigma_\varepsilon^2+n\sigma_u^2}\right] \\
=&\,\ell(\mu,\sigma_u^2,\sigma_\varepsilon^2\mid\bar{y}_{\cdot\cdot})+\ell(\sigma_u^2,\sigma_\varepsilon^2\mid SSE,SSU)
\end{aligned}$$

残差似然函数为 $\ell(\sigma_u^2,\sigma_\varepsilon^2|SSE,SSU)$，其中

$$\begin{aligned}
\ell(\sigma_u^2,\sigma_\varepsilon^2|SSE,SSU) =&\left[(tn-1)\log_{\mathrm{e}}(2\pi)+t(n-1)\log_{\mathrm{e}}(\sigma_{\varepsilon}^2)\right. \\
&\left.+(t-1)\log_{\mathrm{e}}(\sigma_\varepsilon^2+n\sigma_u^2)+\frac{SSE}{\sigma_\varepsilon^2}+\frac{SSU}{\sigma_\varepsilon^2+n\sigma_u^2}\right]
\end{aligned}$$

通过关于两个参数 $\sigma^2_u,\sigma^2_\varepsilon$ 对 $\ell(\sigma_u^2,\sigma_\varepsilon^2|SSE,SSU)$ 进行微分，然后将导数设置为零，可以获得受限最大似然方程。将导数设置为零，在求解参数时评估的导数为

$$\begin{aligned}\frac{\partial\ell({\sigma}_u^2,{\sigma}_\varepsilon^2|SSE,SS{U})}{\partial{\sigma}_\varepsilon^2}&=\frac{t(n-1)}{{\tilde{\sigma}}_\varepsilon^2}+\frac{t-1}{{\tilde{\sigma}}_\varepsilon^2+n{\tilde{\sigma}}_u^2}-\frac{SS{E}}{({\tilde{\sigma}}_\varepsilon^2)^2}-\frac{SS{U}}{({\tilde{\sigma}}_\varepsilon^2+n{\tilde{\sigma}}_u^2)^2}=0\\\frac{\partial\ell({\sigma}_u^2,{\sigma}_\varepsilon^2|{SSE},{SSU})}{\partial{\sigma}_u^2}&=\frac{n(t-1)}{{\tilde{\sigma}}_\varepsilon^2+n{\tilde{\sigma}}_u^2}-\frac{nSS{U}}{\left({\tilde{\sigma}}_\varepsilon^2+n{\tilde{\sigma}}_u^2\right)^2}=0\end{aligned}$$

残差最大似然方程的解是

$$\tilde{\sigma}_\varepsilon^2=\frac{SSE}{t(n-1)}=MSError,\quad\mathrm{and}\quad\tilde{\sigma}_u^2=\frac1n{\left[\frac{SSU}{t-1}-MSError\right]}=\frac1n[MSU-MSError]$$

残差最大似然估计为

$$\hat{\sigma}_\varepsilon^2=\tilde{\sigma}_\varepsilon^2=\frac{SSE}{t(n-1)}=MSError$$

以及

$$\hat{\sigma}_u^2=\begin{cases}\tilde{\sigma}_u^2&\mathrm{~if~}\tilde{\sigma}_u^2\geq0\\0&\mathrm{~if~}\tilde{\sigma}_u^2<0&\end{cases}$$

当 ${\sigma}_u^2$ 的估计为零时，则通过合并 SSE 和 SSU 及其自由度来重新计算 ${\sigma}_\varepsilon^2$ 的估计，从而获得

$$\hat{\sigma}_\varepsilon^2=\frac{SSE+SSU}{tn-1}$$

表 \@ref(tab:table19-5) 和 \@ref(tab:table19-6) 分别包含用于提取示例 [19.2](#sec19-1-2) 和示例 [19.3](#sec19-1-3) 中数据集方差分量的 REML 估计的 SAS-Mixed 代码和结果。


<table>
<caption>(\#tab:table19-5)Proc Mixed 代码来计算示例 19.2 中数据的方差分量和均值的受限最大似然估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2019.5.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table19-6)Proc Mixed 代码来获得示例 19.3 中双向随机效应模型方差分量的受限最大似然估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2019.6.png">
</td>
  </tr>
</tbody>
</table>

## MIVQUE 法 {#sec19-4}

Rao (1971) 描述了获得方差分量的**最小方差二次无偏估计** (minimum variance quadratic unbiased estimators, MIVQUE) 的一般过程。对于第 \@ref(sec18-2) 节的一般随机效应模型，方差的线性组合

$$\theta=C_0\sigma_\varepsilon^2+C_1\sigma_1^2+C_2\sigma_2^2+\cdots+C_k\sigma_k^2$$

的 MIVQUE 是观测值的二次函数，其关于 $\theta$ 是无偏的，并且在 $\theta$ 的二次无偏估计类中具有最小方差。因此，方差分量的 MIVQUE 估计具有最小方差性质，而矩法估计通常不具有。每个单独的方差分量都可以被选作可能要估计的参数。选择 $C_0=1,C_1=C_2=\cdots=C_k=0$ 提供 $\theta=\sigma^2_\varepsilon$. $C_i$ 值的其它选择将提供 $\theta=\sigma^2_i$ 以及方差分量的其它线性组合。

### 方法说明 {#sec19-4-1}

$\theta$ 的估计必须是 $\boldsymbol y$ 的二次函数，因此，对于某个矩阵 $\boldsymbol A$，$\theta$ 的估计的形式为 $\boldsymbol y^\prime\boldsymbol A\boldsymbol y$. $\boldsymbol y^\prime\boldsymbol A\boldsymbol y$ 的期望是

$$E(\boldsymbol y^{\prime}\boldsymbol A\boldsymbol y)=\mathrm{tr}(\boldsymbol{\Sigma}\boldsymbol A)+\mu^2\boldsymbol{j}_n^{\prime}\boldsymbol A\boldsymbol{j}_n$$

根据假设，$E(\boldsymbol y^{\prime}\boldsymbol A\boldsymbol y)=\theta$。由于期望不取决于 $\mu$，因此必须选择 $\boldsymbol A$ 以满足 $\mu^2\boldsymbol{j}_n^{\prime}\boldsymbol A\boldsymbol{j}_n=0$. 在正态性条件下，当 $\mu^2\boldsymbol{j}_n^{\prime}\boldsymbol A\boldsymbol{j}_n=0$ 时，$\boldsymbol y^{\prime}\boldsymbol A\boldsymbol y$ 的方差为

$$\mathrm{Var}(\boldsymbol{y}^{\prime}A\boldsymbol{y})=2\mathrm{~tr}[\boldsymbol{\Sigma}A]^2$$

因此，$\theta$ 的 MIVQUE 是 $\boldsymbol y^{\prime}\boldsymbol A\boldsymbol y$，其中 $ \boldsymbol A$ 被选择为使得 $\operatorname{tr}[\boldsymbol{\Sigma}\boldsymbol A]=\theta $ 并且 $\operatorname{tr}[\boldsymbol{\Sigma}\boldsymbol A]^2$ 在如下参数空间中达到最小化

$$\{0<\sigma_\varepsilon^2<\infty,\quad0<\sigma_1^2<\infty,\quad0<\sigma_2^2<\infty,\ldots,0<\sigma_k^2<\infty\}$$

Rao (1971) 证明了$\boldsymbol \sigma^2=[\sigma_\varepsilon^2,\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2]'$ 的 MIVQUE 是 $\hat{\boldsymbol \sigma}^2=\boldsymbol S^{-1}\boldsymbol f$，其中 $\boldsymbol S$ 是具有如下元素的 (k+1) × (k+1) 的矩阵

$$s_{ii^{\prime}}=\mathrm{tr}[\boldsymbol X_i\boldsymbol X_i^{\prime}\boldsymbol R\boldsymbol X_{i^{\prime}}\boldsymbol X_{i^{\prime}}]\quad i,i^{\prime}=0,1,2,\ldots,k$$


$\boldsymbol f$ 是具有如下元素的 (k+1) × 1 的向量

$$\boldsymbol f_i=\boldsymbol y^{\prime}\boldsymbol R\boldsymbol X_i\boldsymbol X_i^{\prime}\boldsymbol R\boldsymbol y,\quad i=0,1,2,\ldots,k$$

其中

$$\boldsymbol R=\boldsymbol\Sigma^{-1}[\boldsymbol I_n-\boldsymbol j_n(\boldsymbol j_n^{\prime}\boldsymbol\Sigma^{-1}\boldsymbol j_n)^{-1}\boldsymbol j_n^{\prime}]\boldsymbol \Sigma^{-1}$$

  $\hat\sigma^2$ 的解取决于 $\boldsymbol \Sigma$ 的元素，这些元素是未知方差分量的函数。为了计算 $\sigma^2$ 的 MIVQUE，必须将一些常数代入 $\boldsymbol \Sigma$ 中的 $\sigma_\varepsilon^2,\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2$ 。对于这组常数，$\sigma^2$ 的估计是 MIVQUE（并且是 $\boldsymbol y$ 的二次函数）。为了使 $\hat\sigma^2$ 是 $\sigma^2$ 的 MIVQUE，$\boldsymbol R$ 的元素必须不依赖于数据向量。一些软件使用 1 作为残差方差的值，使用 0 作为其他方差和协方差的值。

通常，最好将接近真实值的 $\sigma_\varepsilon^2,\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2$ 的值代入 $\boldsymbol \Sigma$. 一种可能的程序是从其他实验中获得值。使用固定值（不依赖于当前数据）作为非迭代过程中 $\sigma^2$ 的起始值，或使用零迭代 (zero iterations) 的过程，提供了一个称为 MIVQUE0 的解，其中 0 表示未执行任何迭代。Swallow and Monahan (1984) 使用矩法估计作为方差的初始值，并使用 MIVQUE-A 来描述方差分量的最终估计。另一种方法是使用迭代程序 (Brown, 1976)，通过使用方差分量的一些初始值，例如 $\sigma_{\varepsilon0}^2,\sigma_{10}^2,\sigma_{20}^2,\ldots,\sigma_{k0})^2$ 来开始该过程。使用这些初始值来评估 $\boldsymbol \Sigma$，并获得 $\hat\sigma^2_{(0)}$。这里，$\hat\sigma^2_{(0)}$ 取决于为 $\sigma^2_{(0)}$ 选择的值。然后使用 $\hat\sigma^2_{(0)}$ 来评估 $\boldsymbol \Sigma$，以获得第二次迭代估计 $\hat\sigma^2_{(1)}$. 继续迭代过程，直到一次迭代与下一次迭代之间几乎没有变化。因为 $\boldsymbol \Sigma$ 的元素是 $\boldsymbol y$ 的函数，所得到的迭代 MIVQUE 值不再是 $\boldsymbol y$ 的二次函数。可以将 $\hat\sigma^2_{i}$ 的最终估计，例如在在第 m+1 步时，称为给定先前值 $\hat\sigma^2_{(m)}$ 的 MIVQUE. 对于均衡模型，Swallow and Searle (1978) 已经表明方程可以简化从而可以获得显式解。此解与矩法提供的解相同。当存在不等样本量和/或空单元格时，迭代过程可能是合适的方法。Swallow and Searle (1978) 的一项模拟研究表明，REML, ML 和矩法比 MIVQUE0 提供更好的方差分量估计。

MIVQUE 的值（无论是以 $\sigma_i^2$ 的常数值评估还是在给定的前一步中评估）都是 $\boldsymbol y$ 的二次型的线性组合。因此，可以评估方差，因为二次型     $\boldsymbol y^{\prime}\boldsymbol A\boldsymbol y$ 的方差为 $2\, \text{tr}(\boldsymbol B\boldsymbol \Sigma)^2$. Swallow and Searle (1978) 展示了如何使用这些表达式来获得不均衡单向模型的估计的方差。他们计算并比较了具有不同总体数量、样本量和方差的非均衡单向模型的 MIVQUE 和矩法估计的方差。MIVQUE 估计的方差是根据 $\sigma_u^2$ 和 $\sigma^2_\varepsilon$ 的真实值在估计过程中进行评估的。从矩法获得的 $\sigma^2_\varepsilon$ 的估计与从 MIVQUE 法获得的估计非常相似，其中 MIVQUE 的方差不大于矩法估计的方差的 4％.

对于相当均衡的模型（$n_i$ 差异不大），$\sigma^2_u$ 的 MIVQUE 的方差不大于矩法估计的方差的 10%. 对于许多不均衡样本量，$\sigma^2_u$ 的 MIVQUE 方差比相应的矩法估计小 60%. 如果在估计过程中使用除 $\sigma^2_u$ 和 $\sigma^2_\varepsilon$ 真实值之外的值，则两种方法的估计的方差将更为相似。SAS-Mixed 具有一个 `MIVQUE0` 选项用于 MIVQUE 法，是一种非迭代方法。

本节最后以一个使用 MIVQUE 估计的不均衡单向设计的例子结束。

### 应用。示例 19.6：MIVQUE 用于不均衡单向设计 {#sec19-4-2}

Swallow and Searle (1978) 提出的方程用于一般单向模型，然后将其应用于示例 [19.2]({#sec19-1-2}) 中的数据。模型为

$$y_{ij}=\mu+u_i+\varepsilon_{ij}\quad i=1,2,\ldots,t\mathrm{~and~}j=1,2,\ldots,n_i$$

其中 $\boldsymbol{u}\thicksim{N}(0,\sigma_u^2\boldsymbol{I}_t),\boldsymbol{\varepsilon}\thicksim{N}(0,\sigma_\varepsilon^2\boldsymbol{I}_N)$ 且 $\boldsymbol{u},\boldsymbol{\varepsilon}$ 是独立随机变量。$N=\sum_{i=1}^tn_i$.

定义

$$k_i=\frac{n_i}{\sigma_{\varepsilon0}^2+n_i\sigma_{u0}^2}\quad\mathrm{and}\quad K=\frac1{\sum_{i=1}^tk_i}$$

矩阵 $\boldsymbol \Sigma$ 的元素是

$$\begin{aligned}
S_{11}& =\sum_{i=1}^tk_i^2-2K\sum_{i=1}^tk_i^3+K^2\left(\sum_{i=1}^tk_i^2\right)^2  \\
S_{12}& =\sum_{i=1}^t\frac{k_i^2}{n_i}-2K\sum_{i=1}^t\frac{k_i^3}{n_i}+K^2\biggl(\sum_{i=1}^tk_i^2\biggr)\biggl(\sum_{i=1}^t\frac{k_i^2}{n_i}\biggr)  \\
S_{22}& =\frac{N-t}{\sigma_{\varepsilon0}^4}+\sum_{i=1}^t\frac{k_i^2}{n_i^2}-2K\sum_{i=1}^t\frac{k_i^3}{n_i^2}+K^2\left(\sum_{i=1}^t\frac{k_i^2}{n_i}\right)^2 
\end{aligned}$$

向量 $\boldsymbol f$ 的元素是

$$\begin{aligned}f_1&=\sum_{i=1}^tk_i^2\bigg(\bar{y}_i-K\sum_{i=1}^tk_i\bar{y}_{i\cdot}\bigg)^2\\f_2&=\frac{\sum_{i=1}^t\sum_{j=1}^{n_i}y_{ij}^2-\sum_{i=1}^tn_i\bar{y}_{i\cdot}^2}{\sigma_{\epsilon0}^4}+\sum_{i=1}^t\frac{k_i^2\left(\bar{y}_{i\cdot}-K\sum_{i=1}^tk_i\bar{y}_{i\cdot}\right)^2}{n_i}\end{aligned}$$

对于给定的 $\sigma^2_{u0}$ 和 $\sigma^2_{\varepsilon0}$ 的值，$\sigma^2_{u}$ 和 $\sigma^2_{\varepsilon}$ 的 MIVQUE 估计为 $\hat{\sigma}^2=\boldsymbol S^{-1}\boldsymbol f$，或

$$\begin{aligned}\hat{\sigma}_\varepsilon^2&=\frac{s_{11}f_2-s_{12}f_1}c\\\hat{\sigma}_u^2&=\frac{s_{22}f_1-s_{12}f_2}c\end{aligned}$$

其中 $c=s_{11}s_{22}-s_{12}^2$. 表 \@ref(tab:table19-1) 中数据的两个方差分量的 MIVQUE 估计在表 \@ref(tab:table19-7) 中，表 \@ref(tab:table19-2) 中数据的四个方差分量的 MIVQUE 估计在表 \@ref(tab:table19-8) 中。这些估计是使用 SAS Mixed 的非迭代解获得的。估计的方差为

$$\mathrm{Var}(\hat{\sigma}_\varepsilon^2)=\frac{2s_{11}}c,\quad\mathrm{Var}(\hat{\sigma}_u^2)=\frac{2s_{2}}c,\quad\mathrm{and}\quad\mathrm{Cov}(\hat{\sigma}_\varepsilon^2,\hat{\sigma}_u^2)=\frac{-2s_{12}}c$$

<table>
<caption>(\#tab:table19-7)Proc Mixed 代码来计算示例 19.2 中数据的方差分量和均值的 MIVQUE0 估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2019.7.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table19-8)Proc Mixed 代码来获得示例 19.3 中双向随机效应模型的方差分量的 MIVQUE0 估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2019.8.png">
</td>
  </tr>
</tbody>
</table>

估计对 $\sigma^2_{u0}$ 和 $\sigma^2_{\varepsilon0}$ 的选择不是很敏感，但方差确实取决于 $\sigma^2_{u0}$ 和 $\sigma^2_{\varepsilon0}$ 初始值的选择。表 \@ref(tab:table19-9) 包含表 \@ref(tab:table19-1) 中数据的 MIVQUE0 估计及其方差，其中使用了所列起始值。当 $\sigma^2_{u0}$ 和 $\sigma^2_{\varepsilon0}$ 的起始值与 $\hat\sigma^2_{u}$ 和 $\hat\sigma^2_{\varepsilon}$ 的估计相距甚远时，估计的方差变化非常大。如果使用迭代程序，则解收敛于 $\hat{\sigma}_\varepsilon^2=0.057003$ 和 $\hat{\sigma}_u^2=0.073155$，其中 $\mathrm{Var}(\hat{\sigma}_\varepsilon^2)=0.000721,\mathrm{Var}(\hat{\sigma}_u^2)=0.005694$ 以及 $\mathrm{cov}(\hat{\sigma}_\varepsilon^2,\hat{\sigma}_u^2)=-0.000235$. 迭代过程从几个值开始（$\sigma^2_{u0}=2,\sigma^2_{\varepsilon0}=1$ 以及 $\sigma^2_{u0}=50,\sigma^2_{\varepsilon0}=1000$ 等）。所有起始值的选择在四次迭代中收敛到上述值。经过两次迭代，估计相当稳定，但方差仍在变化。

<table>
<caption>(\#tab:table19-9)示例 19.2 所选起始值的 MIVQUE0 方差分量估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2019.9.png">
</td>
  </tr>
</tbody>
</table>

## 使用 JMP 估计方差分量 {#sec19-5}

方差分量的估计可以使用 JMP 软件 (SAS Institute, Inc., 2005) 的 fit model 选项获得。图 \@ref(fig:figure19-1) 给出了示例 [19.2](#sec19-1-2) 的数据集，显示在 JMP 数据表中（从 SAS 数据集导入）。在 Analyze 菜单上，选择 fit model，如图 \@ref(fig:figure19-2) 所示。在 fit model 界面上，选择 Damage 作为 Y 变量，选择 variety 作为模型效应。使用 attributes 菜单来指定 variety 是一个随机效应。默认的估计方法是 REML，但可以选择 EMS 来提供使用 III 型平方和的矩法估计。单击 run model 按钮以获取图 \@ref(fig:figure19-3) 中的结果。方差分量的估计及其标准误估计与表 \@ref(tab:table19-5) 中 SAS 的结果相似。主要区别是使用了 Wald 法而不是 Satterthwaite 近似来计算 variety 方差分量的置信区间（详见第 \@ref(chap20) 章）。

<div class="figure" style="text-align: center">
<img src="figure/figure 19.1.png" alt="示例 19.2 数据集的 JMP 表" width="744" />
<p class="caption">(\#fig:figure19-1)示例 19.2 数据集的 JMP 表</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 19.2.png" alt="示例 19.2 的 JMP fit model 表" width="925" />
<p class="caption">(\#fig:figure19-2)示例 19.2 的 JMP fit model 表</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 19.3.png" alt="示例 19.2 的 JMP REML 结果" width="925" />
<p class="caption">(\#fig:figure19-3)示例 19.2 的 JMP REML 结果</p>
</div>

图 \@ref(fig:figure19-4) 是示例 [19.3](#sec19-1-3) 数据的 JMP 数据表。fit model 界面如图 \@ref(fig:figure19-5) 所示，其中 row, col 和 row × col 被选择为随机效应，并选择 REML 方法进行估计。方差分量的 REML 估计如图 \@ref(fig:figure19-6) 所示，其中的结果与表 \@ref(tab:table19-6) 中 SAS 的结果相似。fit model 界面有另一个选项，可以在其中选中 "unbounded variance components" 框。此选项不会将解限制在参数空间中（类似于在 SAS Mixed 中使用 `unbounded` 选项）。图 \@ref(fig:figure19-7) 选中了 "unbounded variance components" 框，解如图 \@ref(fig:figure19-8) 所示。row 和 col 方差分量的解是负的，Wald 法用于计算置信区间（残差除外）。JMP fit model 过程为具有随机效应的模型提供了适当的分析，SAS Mixed也是如此。

<div class="figure" style="text-align: center">
<img src="figure/figure 19.4.png" alt="示例 19.3 数据集的 JMP 表" width="738" />
<p class="caption">(\#fig:figure19-4)示例 19.3 数据集的 JMP 表</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 19.5.png" alt="用于示例 19.3 REML 估计的 JMP fit model 界面" width="782" />
<p class="caption">(\#fig:figure19-5)用于示例 19.3 REML 估计的 JMP fit model 界面</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 19.6.png" alt="示例 19.3 的 JMP REML 结果" width="897" />
<p class="caption">(\#fig:figure19-6)示例 19.3 的 JMP REML 结果</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 19.7.png" alt="示例 19.3 框选 &quot;unbounded variance components&quot; 的JMP fit model 界面" width="817" />
<p class="caption">(\#fig:figure19-7)示例 19.3 框选 "unbounded variance components" 的JMP fit model 界面</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 19.8.png" alt="示例 19.3 的 JMP REML 的无界结果" width="885" />
<p class="caption">(\#fig:figure19-8)示例 19.3 的 JMP REML 的无界结果</p>
</div>

## 结束语 {#sec19-6}

本章介绍了四种方法，矩法、最大似然法、残差最大似然法和 MIVQUE 法，用于获得随机效应模型方差分量的估计。使用两个示例，一个单向随机效应模型和一个双向随机效应模型来演示每种估计方法。当数据均衡且提供正解时，REML, MIVQUE0 和 矩法的估计是相同的。当数据集不均衡时，每种方法都会产生不同的估计。使用 SAS-MIXED 和 JMP 进行了计算。

## 练习 {#sec19-7}
