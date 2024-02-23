
# 具有异质误差的完全随机设计结构中的单向处理结构 {#chap2}

在本章中，考虑分配给实验单元的处理可能影响响应的方差以及均值的情况。从单向均值模型 $y_{ij} = \mu_i + \varepsilon_{ij}$ 开始，其中 $i = 1,2,\cdots,t;\,j =1,2,\cdots,n_i$. 在第 \@ref(chap1) 章中，假设实验误差都具有相同的方差；也就是说，预期处理可能会改变抽样总体的均值，但不会改变方差。当处理影响方差和均值时，本章描述了一些分析数据的方法。在这种情况下，实验者想要回答的关于平均值的问题类型与第 \@ref(chap1) 章中的类似。也就是说，

1) 所有均值都相等吗？
2) 是否可以在均值之间进行成对比较？
3) 是否可以检验形如 $\sum_{i=1}^tc_i\mu_i=a$ 的假设并为 $\sum_{i=1}^tc_i\mu_i$ 构造置信区间？

此外，关于方差还存在可能感兴趣的问题，例如

1) 所有方差都相等吗？
2) 是否存在组内方差相等而组间方差不相等的处理分组？

在回答有关模型均值的问题之前，必须获得处理方差的适当描述。

[^non-normality_of_the_error_variances]: 原文为: heterogeneous error variances pose a much more serious problem when ignored than non-normality of the error variances，应为笔误，此处译者做了修改。

**方差齐性检验** (**tests of homogeneity of variances**) 用于回答有关各个处理的数据方差的问题。如果有两种处理，当方差不等时比较均值的问题通常被称为 Behrens–Fisher 问题。此外，异质误差方差被忽略比数据的非正态性带来了更严重的问题[^non-normality_of_the_error_variances]。第一章中的程序在非正态性方面是稳健的，但在异质误差方差方面不是那么稳健。在之前考虑的分析中，假设总体方差都相等，这在许多情况下是合理的假设。当方差不相等时，分析数据的一种方法是简单地忽略它们不相等的事实，并计算在方差相等的情况下计算的相同的 $F$ 统计量或 $t$ 检验。也许令人惊讶的是，模拟研究表明，这些通常的检验都非常好，特别是当样本量都相等或几乎相等时。此外，如果较大的样本量对应于具有较大方差的处理或总体，那么用等方差假设计算的检验也相当好。事实上，通常的检验非常好，以至于许多统计学家甚至不建议对等方差进行检验。其他人试图找到一种能稳定处理方差的变换，即通过变换数据使处理方差相等。当方差不相等时，有一些技术可以在不等方差模型的框架内对均值进行比较。

以下小节描述了单向模型的处理方差相等性的检验程序，以及当方差不相等时分析处理均值的程序。当怀疑通常的技术时，应使用这些程序。接下来描述不等方差模型。

## 模型定义和假设 {#sec2-1}

不等方差模型为

\begin{equation}
y_{ij}=\mu_{i}+\varepsilon_{ij}\quad\mathrm{for}\quad i=1,2,\ldots,t,\,j=1,2,\ldots,n_{i}\quad\mathrm{and}\quad\varepsilon_{ij}\sim\text{independent}\,N(0,\sigma_{i}^{2})
(#eq:2-1)
\end{equation}

其中 $\varepsilon_{ij}\sim\text{independent}\,N(0,\sigma_{i}^{2})$ 表示误差 $\varepsilon_{ij}$ 都是独立地服从正态分布，并且每个正态分布的方差取决于 $i$，对于每个总体或处理可能不同。

## 参数估计 {#sec2-2}

模型中参数的最佳估计值：

$$\hat{\mu}_i=\sum_{j=1}^{n_i}y_{ij}/n_i=\bar{y}_{i\cdot},\quad i=1,2,\ldots,t$$

和

$$\hat\sigma_i^2=\frac{\sum_{j=1}^{n_i}(y_{ij}-\overline{y}_i)^2}{n_i-1},\quad i=1,2,\ldots,t$$

与参数估计相关的抽样分布为

$$\hat{\mu}_i\thicksim\text{independent}\,N(\mu_i,\sigma_i^2/n_i),\quad i=1,2,\ldots,t$$

和

$$\frac{\left(n_i-1\right)\hat\sigma_i^2}{\sigma_i^2}\thicksim\text{independent }\,\chi_{n_i-1}^2,\quad i=1,2,\ldots,t$$

这些抽样分布被用作建立方差相等性检验的基础，并在方差不相等时提供均值分析。

## 方差齐性检验 {#sec2-3}

本节描述了检验等方差假设的五个程序

$$H_{0}\colon\,\sigma_{1}^{2}=\sigma_{2}^{2}=\cdots=\sigma_{t}^{2}\,\mathrm{~vs.~}\,H_{a}\colon(\mathrm{not~}H_{0}:)$$

在尝试分析均值之前，应研究等方差假设。如果没有足够的证据得出方差不相等的结论，则可以使用第 \@ref(chap1) 章中的等方差模型来研究均值。如果有足够的证据认为方差不相等，则应使用 \@ref(sec2-5) 节中描述的程序来提供不等方差框架中的均值分析。建议当不等方差假设在 $\alpha\le0.01$ 时被拒绝时使用不等方差模型。

### Hartley's *F*-Max Test {#sec2-3-1}

所描述的第一个检验被称为 Hartley's *F*-max检验 (Hartley, 1950). 该检验要求所有样本具有相同的大小，即 $n_1 = n_2 = \cdots= n_t$. 该检验基于统计量

$$F_{\max}=\frac{\max_i\{\hat\sigma_i^2\}}{\min_i\{\hat\sigma_i^2\}}$$

[^choice_of_alpha]: 原文为: so you are going to reject the null hypothesis more often than specified by the choice of $\alpha$.

[^nottoounnequl]: 原文为: not too unequal.

附表 A.1 提供了 $\alpha = 0.05$ 和 $0.01$ 时 $F_{\max}$ 的百分位点。若 $F_{\max} > F_{\max,\alpha,v,k}$，其中 $v = n - 1$ 为与 $k$ 个处理方差中的每一个相关的自由度，则拒绝零假设 $H_0$. 如果 $n_i$ 不都相等，则可以通过取 $n=\max_i\{n_i\}-1$ 来获得 $H_0\,\mathrm{~vs.~} H_a$ 的“自由” (liberal)检验。这个检验是自由的，因为它假设所有处理都有相同的（最大）样本量，所以你会比选定的显著性水平 $\alpha$ 更频繁地拒绝零假设[^choice_of_alpha]。当样本量相差不大[^nottoounnequl]时，此程序提供了合理的检验。它还可以防止人们进行通常的方差分析，即使方差分析不合适的可能性极小。{#sec2-4} 节中有一个说明该检验使用的示例。

### Bartlett's Test {#sec2-3-2}

方差齐性检验的第二个检验是 (Bartlett, 1937) 提出的检验，其优点是不要求 $n_i$ 相等。Bartlett 检验统计量为

\begin{equation}
U=\frac1C{\left[v\log_\mathrm{e}(\hat\sigma^2)-\sum_{i=1}^tv_i\log_\mathrm{e}(\hat\sigma_i^2)\right]}
(#eq:2-2)
\end{equation}

其中

$$v=n_i-1,\quad v=\sum_{i=1}^tv_i,\quad\hat\sigma^2=\sum_{i=1}^tv_i\hat\sigma_i^2/v $$

且

$$C=1+\frac1{3(t-1)}{\left[\sum_{i=1}^t\frac1{v_i}-\frac1v\right]}$$

如果 $U>\chi_{\alpha,t-1}^2$，则拒绝等方差假设。方差齐性的前两个检验的缺点之一是它们对偏离正态性以及偏离等方差假设非常敏感。

### Levene's Test {#sec2-3-3}

Levene (1960) 建议对来自单向均值或效应模型的残差绝对值进行单向方差分析。残差绝对值由 $z_{ij}=|y_{ij}-\overline{y}_{i\cdot}|,i=1,2,\ldots,t;j=1,2,\ldots,n_{i}$ 给出。方差分析的 $F$ 检验提供残差绝对值的处理均值相等性检验。如果均值不等，则有证据表明，一种处理的残差平均而言大于另一种处理的残差。残差可以指导哪些方差不相等，多重比较检验（见第 \@ref(chap3) 章）可以用来对这些均值进行成对比较。Levene's Test 的一个改进是在方差分析中使用残差平方。

### Brown and Forsythe's Test {#sec2-4-4}

Brown and Forsythe (1974)  使用 Levene 的过程，并通过对每个处理的观察值与中位数的偏差的绝对值进行单向方差分析来对其进行修改。中位数偏差的绝对值由 $u_{ij}=|y_{ij}-y_{i\mathrm{med}}|,i=1,2,\ldots,t;j=1,2,\ldots,n_{i}$ 给出。方差分析的 $F$ 检验提供了偏差绝对值的处理均值相等性检验。如果均值不等，则有证据表明，一种处理的偏差平均而言大于另一种处理的偏差。偏离中位数的绝对值的均值可以提供关于哪些方差不等的指导，因为多重比较检验可以用于在这些均值之间进行成对比较。当数据的分布不对称时，这种利用中位数偏差的方法提供了比 Levene's method 更大的功效。

### O'Brien's Test {#sec2-3-5}

O'Brien (1979) 计算了一种得分

\begin{equation}
r_{ij}=[(w+n_i-2)n_i(y_{ij}-\overline{y}_{i\cdot})^2-w\hat\sigma_i^2(n_i-1)]/[(n_i-1)(n_i-2)]
(#eq:2-3)
\end{equation}

其中 $w$ 为权重参数。该程序是对计算出的得分进行方差分析。当 $w=0.5$ 时，得分的均值是样本方差 $\hat \sigma^2_i$，因此得分均值的比较是数据方差的比较。

还有其他几种程序可用于使用参数和非参数方法来检验方差相等或尺度参数相等 [如 Conover et al. (1981) 和 Olejnik & Algina (1987)]。McGaughey (2003) 提出了一种使用数据深度 (data depth) 概念的检验，并将该程序应用于单变量和多变量总体。数据深度超出了本书的范围。

### 一些建议 {#sec2-3-6}

Conover (1981) 以及 Olejnik & Algina (1987) 对方差检验的同质性进行了模拟研究，其中包括上述检验以及许多其他检验。研究表明，没有一种检验对所有情况都是稳健和具有最大功效的。根据 Olejnik 和 Algina，Levene's Test 是 Conover 等人研究的更好的检验之一。O'Brien's Test 似乎提供了一个合适的大小检验 (size test)，而不会损失太多功效。当分布具有重尾 (heavy tails) 时，Brown–Forsyth's Test 似乎更好。根据他们的结果，我们提出以下建议：

1) 如果分布具有较重的尾部，请使用 Brown–Forsyth's test.
2) 如果分布有些偏斜，请使用 O'Brien's test.
3) 如果数据几乎是正态分布的，那么任何检验都是合适的，包括 Bartlett's 和 Hartley's tests.

Levene's 和 O'Brien's tests 可以很容易地进行定制，用于涉及多个因素的设计实验，包括协方差分析 (Milliken and Johnson, 
2002). Levene's, O'Brien's 和 Brown–Forsyth's tests 在正态分布数据方面几乎与 Bartlett's 和 Hartley's tests 一样好，在非正态分布的数据方面优于它们。Conover et al. 以及 Olejnik 和 Algina 讨论了一些非参数检验，但它们更难计算，并且上述推荐的检验几乎同样有效。下面是一个示例，其中演示了用于方差相等的每个检验。

## 示例：药物和错误 {#sec2-4}

[^subject]: 本文将 "subject" 统一翻译为“个体”，不论语境。

表 \@ref(tab:table2-1) 中的数据来自在两种可能的药物影响下对个体[^subject]进行的配对联想学习任务 (paired-association learning task) 实验。第1组为对照组（无药物），第2组给药1，第3组给药2，第4组同时给两种药。每组数据的样本量、总和、中位数、平均值和方差见表 \@ref(tab:table2-1)。

<table class="table table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:table2-1)配对联想学习任务实验数据</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%202.1.png">
</td>
  </tr>
</tbody>
</table>

*F*-max 统计量为 $F_\max=16.286/1.867=8.723$. 取 $k=t=4,v=7$ 的自由 $5\%$ 临界点从表 A.1 得到。临界点为 $8.44$，由于 $8.723>8.44$，以 $0.05$ 的显著性水平拒绝 $H_{0}\colon\,\sigma_{1}^{2}=\sigma_{2}^{2}=\cdots=\sigma_{t}^{2}\,\mathrm{~vs.~}\,H_{a}\colon(\mathrm{not~}H_{0}:)$，但不能在 $0.01$ 的水平拒绝。

Bartlett's test 的计算为：

$$C=1+\frac1{3\times3}{\left(\frac16+\frac15+\frac17+\frac17-\frac1{25}\right)}$$

且

$$\hat\sigma^2=\frac{6(16.2857)+5(1.8667)+7(9.6964)+7(2.7860)}{25}=7.7769$$

因此

$$\begin{aligned}
U& =\frac1C{\left(v\log_e\hat\sigma^2-\sum_{i=1}^4v_i\log_e\hat\sigma_i^2\right)}  \\
&=\frac1{1.068}[25\log_{\mathrm{e}}(7.7769)-6\log_{\mathrm{e}}(16.2857)-5\log_{\mathrm{e}}(1.8667) \\
&-7\log_{\mathrm{e}}(9.6964)-7\log_{\mathrm{e}}(2.7860)] \\
&=7.8111
\end{aligned}$$

与 $U$ 相关的渐近抽样分布是基于自由度为 $3$ 的卡方分布，检验的显著性水平为 $0.0501$，可以再次得出结论，方差在大约 $5\%$ 的显著性水平上是不相等的。

Levene's test 的计算从残差或观测值与处理均值的偏差的计算开始。接下来，如表 \@ref(tab:table2-2) 所示计算残差的绝对值。最后，对残差的这些绝对值进行单向方差分析。所得到的 $F$ 统计量的值为 $6.97$，其基于自由度 $3,25$。观察到的 Levene's test 的显著性水平为 $0.0015$. 在进行方差分析之前，可以通过对表 \@ref(tab:table2-2) 中的项目进行平方来获得 Levene's test 的平方偏差或平方残差版本。在这种情况下，$F$ 统计量的值为 $7.36$，观察到的显著水平为 $0.0011$（也基于自由度 $3,25$）。

<table class="table table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:table2-2)计算 Levene's test 的 $z_{ij} = |y_{ij}-y_{i\cdot}|$ 值，其中 $y_{ij}$ 来自表 \@ref(tab:table2-1)</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%202.2.png">
</td>
  </tr>
</tbody>
</table>

Brown–Forsyth's test 统计量是通过计算观察值与处理中位数偏差的绝对值（中位数见表 \@ref(tab:table2-1)）获得的。表 \@ref(tab:table2-3) 包含了与中位数偏差的绝对值。接下来，单向方差分析提供的 $F$ 统计量为 $5.49$，并且观察到的显著性水平是 $0.0049$（也基于自由度 $3,25$）。

<table class="table table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:table2-3)观测结果相对于处理中位数偏差的绝对值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%202.3.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table2-4) 包含使用式 \@ref(eq:2-3) 计算的 $r_{ij}$ 值，其中 $w = 0.5$. 通过进行方差分析获得 O'Brien's test 统计量。$F$ 统计量的值为 $6.30$，观察到的显著性水平为 $0.0025$. 使用 $w = 0.7$（计算未显示）的 $F$ 统计量值为 $5.90$，观察到的显著性水平为 $0.0035$. 
每个 O'Brien's test 都基于自由度 $3,25$.

<table class="table table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:table2-4)使用 $w=0.5$ 进行 O'Brien's Test 的得分</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%202.4.png">
</td>
  </tr>
</tbody>
</table>

每个检验统计量表明，有足够的证据得出方差不相等的结论。残差绝对值的组均值见表 \@ref(tab:table2-5)。这些处理的残差绝对值均值之间的成对比较见表 \@ref(tab:table2-6)。无药物和药物2的残差绝对值的平均值没有差异，药物1和同时给药没有差异，但是这两组之间存在差异。可以使用具有两个方差的简单模型来继续分析治疗均值。使用简单方差模型将提高一些关于均值检验的功效。两个方差模型和相应的均值比较将遵循使用四个方差的分析的讨论。

<table class="table table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:table2-5)残差绝对值均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%202.5.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:table2-6)残差绝对值组均值的成对比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%202.6.png">
</td>
  </tr>
</tbody>
</table>

## 关于线性组合的推断 {#sec2-5}

当方差 $\sigma_i^2$ 相差很大而无法应用第 \@ref(chap1) 章中讨论的检验和置信区间时，本节讨论了处理均值的任意线性组合 $\sum_{i=1}^tc_i\mu_i$ 的假设检验和置信区间的构造问题。如果方差相等假设在 $0.01$ 或 $1\%$ 水平被拒绝，建议使用本节和下一节中的程序。如果没有足够的证据证明方差不等，则可以使用第 \@ref(chap1) 章中的结果对处理均值进行推断。

$\sum_{i=1}^tc_i\mu_i$ 的最佳估计为 $\sum_{i=1}^tc_i\hat \mu_i$，其抽样分布为

$$\sum_{i=1}^tc_i\hat{\mu}_i\sim N{\left(\sum_{i=1}^tc_i\mu_i,\sum_{i=1}^tc_i^2\sigma_i^2/n_i\right)}$$

因此

$$z=\frac{\sum_{i=1}^{t}c_{i}\hat{{\mu}}_{i}-\sum_{i=1}^{t}c_{i}\mu_{i}}{\sqrt{\sum_{i=1}^{t}c_{i}^{2}\sigma_{i}^{2}/n_{i}}}\sim N(0,1)$$

当方差未知且不等时，一个明显可用于推断 $\sum_{i=1}^tc_i\mu_i$ 的统计量为

$$z=\frac{\sum_{i=1}^tc_i\hat{\mu}_i-\sum_{i=1}^tc_i{\mu}_i}{\sqrt{\sum_{i=1}^tc_i^2\hat\sigma_i^2/n_i}}$$

如果非零 $c_i$ 对应的 $n_i$ 都非常大，则可以合理地假设 $Z$ 具有近似 $N(0,1)$ 分布，因此 $Z$ 可以用于推断 $\sum_{i=1}^tc_i\mu_i$. 此时 $\sum_{i=1}^tc_i\mu_i$ 的近似 $(1 - \alpha)100\%$ 置信区间由下式提供：

$$\sum_{i=1}^tc_i\hat{\mu}_i\pm z_{\alpha/2}\sqrt{\sum_{i=1}^tc_i^2\hat\sigma_i^2/n_i}$$

如果$|z| > z_{\alpha/2}$，则在显著性水平 $\alpha$ 处拒绝 $H_0$.

在其他情况下，注意 $z$ 可以写为

$$z=\frac{\left(\sum_{i=1}^tc_i\hat{{\mu}}_i-\sum_{i=1}^tc_i\mu_i\right)\Big/\sqrt{\sum_{i=1}^tc_i^2\sigma_i^2/n_i}}{\sqrt{\sum_{i=1}^tc_i^2\hat{\sigma}_i^2/n_i}\Big/\sqrt{\sum_{i=1}^tc_i^2\sigma_i^2/n_i}}$$

$z$ 的分子具有标准正态分布，$z$ 的分子和分母独立。$z$ 的分布可以近似为 $t(v)$ 分布，如果 $v$ 的值可以使得

$$V=v\times\frac{\sum_{i=1}^tc_i^2\hat\sigma_i^2\big/n_i}{\sum_{i=1}^tc_i^2\sigma_i^2\big/n_i}$$

近似服从 $\chi^2(v)$. 为了在方差不等时得到 $V$ 的分布的良好卡方近似，选择具有与 $V$ 相同的前两阶矩的卡方分布。也就是说，为了在方差不等的情况下找到 $v$，求出 $v$ 使得 $V$ 的矩等于 $\chi^2(v)$ 分布的前两阶矩（这就是所谓的 Satterthwaite's method ）。这导致自由度为：

$${v}=\frac{\left(\sum_{i=1}^tc_i^2\sigma_i^2/n_i\right)^2}{\sum_{i=1}^t\left[c_i^4\sigma_i^4/n_i^2(n_i-1)\right]}$$

不幸的是，由于 $v$ 取决于 $\sigma_1^2,\sigma_2^2,\ldots,\sigma_t^2$，因此无法精确确定。通常的程序是通过下式来估计 $v$：

\begin{equation}
\hat{{v}}=\frac{\left(\sum_{i=1}^tc_i^2\hat{\boldsymbol\sigma}_i^2/n_i\right)^2}{\sum_{i=1}^t\left[c_i^4\hat{\sigma}_i^4/n_i^2(n_i-1)\right]}
(#eq:2-4)
\end{equation}

综上所述，拒绝 $H_0\colon\sum_{i=1}^tc_i\mu_i=a\,\mathrm{~vs.~}H_a{:}\sum_{i=1}^tc_i\mu_i\neq a$，如果

$$\left|t_c\right|=\frac{\left|\sum_{i=1}^tc_i\hat{\mu}_i-a\right|}{\sqrt{\sum_{i=1}^tc_i^2\hat\sigma_i^2/n_i}}>t_{\alpha/2,\hat{v}}$$

其中 $\hat v$ 由式 \@ref(eq:2-4) 确定。$\sum_{i=1}^tc_i\mu_i$ 近似 $(1 - \alpha)100\%$ 置信区间由下式给出

$$\sum_{i=1}^tc_i\hat{\mu}_i\pm t_{\alpha/2,\hat{v}}\sqrt{\sum_{i=1}^tc_i^2\hat\sigma_i^2/n_i}$$

不幸的是，每当人们想要检验一个新的假设或构造另一个置信区间时，自由度 $\hat v$ 必须重新估计。可以证明 $n_*-1\leq\hat{v}\leq t(n^*-1)$ 其中 $n_*=\min\{n_1,n_2,\ldots,n_t\}$ 且 $n^*=\max\{n_1,n_2,\ldots,n_t\}$. 因此，如果 $|t_c|>t_{\alpha/2,n_*-1}$，则可以保证 $|t_c|>t_{\alpha/2,\hat{v}}$；如果 $|t_c|<t_{\alpha/2,t(n^*-1)}$，则可以保证 $|t_c|<t_{\alpha/2,\hat{v}}$. 此时我们可以避免计算 $\hat v$. 当 $t_{\alpha/2,t(n^*-1)}<|\left.t_c\right|<t_{\alpha/2,n^*-1}$，为了确定是否应该拒绝或不拒绝被检验的零假设，必须计算 $\hat v$ 的值。为了获取置信区间，通常需要计算 $\hat v$. 接下来，用药物误差的示例来证明前面的结果。

## 示例：药物和错误（续） {#sec2-6}

考虑表 \@ref(tab:table2-1) 中的数据，假设实验者有兴趣回答以下问题：

1) 平均而言，药物对学习有任何影响吗？
2) 当给予两种药物时，个体犯的错误是否比只给予一种药物时更多？
3) 两种药物对错误数量的影响是否不同？

为了回答第一个问题，可以检验三个药物组的均值等于对照组均值的假设。

$$H_{01}{:}\,l_1=\mu_1-\frac{(\mu_2+\mu_3+\mu_4)}3=0\mathrm{~vs~}H_{a1}{:}\,l_1\neq0$$

该线性组合的估计为

$$\hat{l}_1=\hat{\mu}_1-\frac{\hat{\mu}_2+\hat{\mu}_3+\hat{\mu}_4}3=4.571-\frac13(34.042)=-6.776$$

 $\hat l_1$ 相应的标准误估计为

$$\begin{aligned}
s.e.(\hat{l}_1)& =\sqrt{\sum_{i=1}^4\left(\frac{c_i^2\hat\sigma_i^2}{n_i}\right)}  \\
&=\sqrt{\frac{\hat\sigma_1^2}7+\frac19{\left(\frac{\hat\sigma_2^2}6\right)}+\frac19{\left(\frac{\hat\sigma_3^2}8\right)}+\frac19{\left(\frac{\hat\sigma_4^2}8\right)}}\\
&=\sqrt{2.535}\\
&=1.592
\end{aligned}$$

与此标准误估计相关的近似自由度通过使用

$$\sum_{i=1}^4\frac{c_i^4\hat\sigma_i^4}{n_i^2(n_i-1)}=0.9052$$

使得

$$\hat{v}=\frac{(2.535)^2}{0.9052}=7.10$$

检验统计量的值为 $t_c=-6.776/1.992=-4.256$，观察到的显著性水平为 $\hat \alpha= 0.0038$.

$l_1$ 的 $95\%$ 置信区间为

$$\hat{l}_1\pm t_{\alpha/2,\hat{v}}\times\widehat{s.e.}(\hat{l}_1)=-6.776\pm(2.365)(1.592)$$

简化为

$$-10.54<\mu_1-\frac{\mu_2+\mu_3+\mu_4}3<-3.01$$

下一个检验是看给予两种药物组的均值是否等于两种药物组的均值的均值。也就是检验

$$H_{02}\colon l_2=\mu_4-\frac{\mu_2+\mu_3}2=0\,\mathrm{~vs.~}H_{a2}{:}\,l_2\neq0$$

该线性组合的估计为

$$\hat l_2=\hat{\mu}_4-\frac{\hat{\mu}_2+\hat{\mu}_3}2=3.6042$$

其估计标准误为

$$\begin{aligned}
\widehat{s.e.}(\hat{l}_{2})& =\sqrt{\sum_{i=1}^4\left(\frac{c_i^2\hat\sigma_i^2}{n_i}\right)}  \\
&=\sqrt{\frac14\!\left(\frac{\hat\sigma_2^2}6\right)+\frac14\!\left(\frac{\hat\sigma_3^2}8\right)+\left(\frac{\hat\sigma_4^2}8\right)}\\
&=\sqrt{0.7290}\\
&=0.8538
\end{aligned}$$

检验统计量的值为 $t_c = 3.6042/0.8538 = 4.221$，因为 $|t_c|>t_{0.005,5}$ 它在 $\alpha=0.01$ 水平下显著。此时 $\hat v$ 的值不需要用 $n_ *-1$ 作为近似自由度来计算，计算出的 $\hat v$ 的值为 $16.8$ ，这是构造 $l_2$ 置信区间所需要的

最后，为了检验假设以查看两种药物的均值是否不同，检验 $H_0{:}\,l_3 = \mu_2 - \mu_3 = 0 \,\mathrm{~vs.~} H_a{:}\,l_3 = \mu_2 - \mu_3 \ne 0$. 该线性组合的估计值为 $\hat l_3 = \hat\mu_2 - \hat\mu_3=3.042$，其标准误估计为

$$\widehat{s.e.}(\hat{l}_3)=\sqrt{\sum_{i=1}^4\left(\frac{c_i^2\hat\sigma_i^2}{n_i}\right)}=\sqrt{\left(\frac{\hat\sigma_2^2}{6}\right)+\left(\frac{\hat\sigma_3^2}8\right)}=\sqrt{1.523}=1.234$$

近似自由度的计算通过使用：

$$\sum_{i=1}^4\frac{c_i^4\hat\sigma_i^4}{n_i^2(n_i-1)}=0.229$$

使得

$$\hat{v}=\frac{(1.523)^2}{0.229}=10.1$$

因此，$t_c = 3.042/1.234 = 2.465$，观察到的显著性水平为 $\hat\alpha = 0.0334$.

## 自由度的一般 Satterthwaite 近似 {#sec2-7}

与估计标准误差相关的自由度的 Satterthwaite 近似由下式获得

$$v=\frac{2*(E\{[\widehat{s.e.}(\hat{l})]^2\})^2}{\mathrm{Var}\{[\widehat{s.e.}(\hat{l})]^2\}}$$

[^partial]: 原文分子没有 $\partial$，应为笔误，此处已修改。

其中 $[\widehat{s.e.}(\widehat{l})]^2$ 用于估计 $E[s.e.(\hat{l})]^2$，且 $\mathrm{Var}[\widehat{s.e.}(\hat{l})]^2$ 的值由 $\sum_{i=1}^t\left[c_i^4\sigma_i^4/n_i^2(n_i-1)\right]$ 估计。对于更复杂的模型，$\mathrm{Var}[\widehat{s.e.}(\widehat{l})]^2$ 可以通过使用一阶泰勒级数 (Kendall and Stuart, 1952) 近似为 $\boldsymbol q'\boldsymbol M \boldsymbol q$，其中 $\boldsymbol M$ 是方差估计的渐近协方差矩阵估计，向量 $\boldsymbol q$ 的元素是 $E[s.e.(\hat{l})]^2$ 关于单个方差的一阶导数，即[^partial]

$$q_{i}=\frac{\partial E[(s.e.(\hat{l})]^{2}}{\partial\sigma_{i}^{2}},\quad i=1,2,\ldots,t$$

根据各处理方差的估计值评估 $q_i$ (Montgomery and Runger, 1993, 1994). 当每个样本的数据呈正态分布时，则

$$\frac{(n_i-1)\hat\sigma_i^2}{\sigma_i^2}$$

为中心卡方随机变量。因此 $E(\hat\sigma_i^2)=\sigma_i^2$ 且 $\operatorname{Var}(\hat\sigma_i^2)=2\sigma_i^4/(n_i-1)$. 令感兴趣的线性组合为 $l=\sum_{i=1}^{t}c_{i}\mu_{i}$，其方差为 $\sigma_{l}^{2}=\sum_{i=1}^{t}c_{i}^{2}\sigma_{i}^{2}/n_{i}$，则 $\sigma_{l}^{2}$ 关于 $\hat \sigma_{i}^{2}$ 的偏导数为

$$\frac{\partial\sigma_l^2}{\partial\sigma_i^2}=\frac{c_i^2}{n_i}$$

使用泰勒级数一阶近似获得的 $\sigma_{l}^{2}$ 的近似方差为

$$\mathrm{Var}(\sigma_l^2)=\sum_{i=1}^t\left[\left[\frac{c_i^2}{n_i}\right]^2\left[\frac{2\sigma_i^4}{n_i-1}\right]\right]^2$$

下一步是用提供近似自由度的相应样本估计值替换总体方差

$$\hat{v}=\frac{2*\left(E\{[\widehat{s.e.}(\hat{l})]^2\}\right)^2}{\mathrm{Var}\{[\widehat{s.e.}(\hat{l})]^2\}}=\frac{\left(\sum_{i=1}^tc_i^2\hat\sigma_i^2\right)^2}{\sum_{i=1}^tc_i^4\hat\sigma_i^4/[n_i^2(n_i-1)]}$$

与上述 Satterthwaite 近似相同。

## 比较所有均值 {#sec2-8}

如前所述，当方差不等时，通常的 $F$ 检验非常稳健，前提是样本量几乎相等，或者假设较大的样本量对应于方差差异较大的总体样本。本节提供了两个等均值假设的附加检验。等均值假设的第一个检验 $H_{0}\colon\,\mu_{1}=\mu_{2}=\cdots=\mu_{t}\,\mathrm{~vs.~}H_{a}\colon\,(\operatorname{not}H_{0}{:})$ ，由 Welch (1951) 给出，称为 Welch's test。定义权重 $W_{i}=n_{i}/\hat\sigma_{i}^{2}$，令 $\bar{y}^{*}=\sum_{i=1}^{t}W_{i}\bar{y}_{i\cdot}\big/\sum_{i=1}^{t}W_{i}$ 为样本均值的加权平均值，令

$$\Lambda=\sum_{i=1}^{t}\frac{(1-W_{i}/W_{\cdot})^{2}}{n_{i}-1}$$

其中 $W_\cdot=\sum_{i=1}^tW_{i\cdot}$. 那么 Welch's test 统计量为

\begin{equation}
F_c=\frac{\sum_{i=1}^tW_i\frac{(\bar{y}_i.-\bar{y}^*)}{(t-1)}}{1+2(t-1)\Lambda/(t^2-1)}
(#eq:2-5)
\end{equation}

它具有近似的 $F$ 分布，分子分母自由度分别为 $v_{1}=t-1,\,v_{2}=(t^{2}-1)/3\Lambda$. 因此当 $F_{c}>F_{\alpha,v_{1},v_{2}}$，拒绝零假设 $H_{0}\colon\,\mu_{1}=\mu_{2}=\cdots=\mu_{t}$. 式 \@ref(eq:2-5) 的分子也可以计算为：$[\sum_{i=1}^{t}W_{i}\bar{y}_{i\cdot}^{2}-\overline{W}\bar{y}^{*2}]/(t-1)$. 该程序使用 \@ref(sec2-7) 节中的数据进行了演示，初步计算见表 \@ref(tab:table2-7).

<table class="table table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:table2-7)计算 Welch's test 的量</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%202.7.png">
</td>
  </tr>
</tbody>
</table>


根据以上信息计算 $W_\cdot=7.341,\bar{y}^{*}=11.724$

$$\Lambda=\frac{(1-0.430/7.341)^2}6+\frac{(1-3.214/7.341)^2}5+\frac{(1-0.825/7.341)^2}7+\frac{(1-2.872/7.341)^2}7=0.376$$

且 $\sum_{i=1}^{t}W_{i}\left.\bar{y}_{i\cdot}^{2}-\overline{W}\right.\bar{y}^{*2}=1050.8069-1009.0954=43.7114$.

Welch's test 统计量的值为

$$F_c=\frac{41.7114/3}{1+2\times2\times0.376/15}=\frac{13.9038}{1.1003}=12.6355$$

自由度为 $v_1=3,v_2=15/(3×0.376) = 13.283$. 观察到的 $F_c$ 相应的显著性概率为 $\hat \alpha=0.00035$. 出于比较目的，通常的 $F$ 统计量 $F_c=14.91$，自由度为 $3,25$. Welch's test 可以使用 SAS^®^GLM 通过在 MEANS 语句中指定 `Welch` 选项来获得。表 \@ref(tab:table2-8) 包含 GLM 代码，用于提供方差相等的 BF 检验和均值相等的 Welch 检验。部分重要结果见表 \@ref(tab:table2-8) 的第二部分。其他方差相等的检验可以通过指定 `O'Brien`、`Levene` 或 `Bartlett` 来获得。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table2-8) SAS-GLM 代码提供 Brown-Forsythe 方差齐性检验和 Welch 方差齐性检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%202.8-1.png">
</td>
  </tr>
  <tr>
   <td style="text-align:center;">  <img src="table/table%202.8-2.png">
</td>
  </tr>
  <tr>
   <td style="text-align:center;">  <img src="table/table%202.8-3.png">
</td>
  </tr>
</tbody>
</table>

检验处理均值相等的第二个程序是将检验一组处理均值 $\mu_i$ 线性组合假设的过程加以推广。
设一个假设是由 $r$ 个独立的 $\mu_i$ 的线性组合组成的，例如 $H_{0}\colon\,\sum_{i=1}^{t}c_{1i}\mu_{i}=0,\sum_{i=1}^{t}c_{2i}\mu_{i}=0,\ldots,\sum_{i=1}^{t}c_{ri}\mu_{i}=0\,\mathrm{~vs.~}H_{a}:(\mathrm{not}\,H_{0})$. 设 $\boldsymbol C$ 是一个矩阵，其第 $k$ 行是第 $k$ 个线性组合的系数，若假定每个总体或处理的数据均服从正态分布，则处理均值向量的联合抽样分布为 $\hat{\boldsymbol{\mu}}\sim N[\boldsymbol{\mu},\boldsymbol V]$，其中 $\boldsymbol V$ 是一个对角阵，其第 $i$ 个对角元素为 $\sigma^2_i /n_i$. 线性组合集 $\boldsymbol C \boldsymbol \mu$ 的联合抽样分布为 $\boldsymbol{C\hat{\boldsymbol{\mu}}}\sim N[\boldsymbol{C\mu},\boldsymbol{CV}\boldsymbol{C}^{\prime}]$. 与零假设的偏差导致的平方和为 $SSH_{0} = [\boldsymbol C\hat{\boldsymbol{\mu}}]^{\prime}[\boldsymbol C\hat{\boldsymbol V}\boldsymbol C^{\prime}]^{-1}[\boldsymbol C\hat{\boldsymbol{\boldsymbol \mu}}]$，其渐近分布为具有 $r$ 个自由度的卡方分布。近似小样本量统计量为 $F_c = SSH_ 0/r$，近似分布为 $F$，自由度为 $r,v$，其中 $v$ 需要近似 (Fai and Cornelius, 1996; SAS Institute, Inc., 1999, p. 2118). 近似自由度的计算开始于对 $\boldsymbol C  \hat {\boldsymbol V} \boldsymbol C' = \boldsymbol Q \boldsymbol D \boldsymbol Q'$ 执行谱分解，其中 $\boldsymbol D$ 是具有 $\boldsymbol C  \hat {\boldsymbol V} \boldsymbol C'$ 的特征根作为对角元素的 $r×r$ 对角阵，并且 $\boldsymbol Q$ 是 $\boldsymbol C  \hat {\boldsymbol V} \boldsymbol C'$ 的对应特征向量的 $r×r$ 正交矩阵。令 $\boldsymbol z_k'$ 为 $\boldsymbol Q \boldsymbol C$ 的第 $k$ 行，并令

$$\nu_k=\frac{2(d_k)^2}{\boldsymbol b_k'\boldsymbol M\boldsymbol b_k}$$

其中 $d_k$ 是 $\boldsymbol D$ 的第 $k$ 个对角元素，$\boldsymbol b_k$ 包含 $\boldsymbol z_k'\boldsymbol V \boldsymbol z_k$ 关于 $\boldsymbol V$ 中每个方差参数的偏导数在方差估计处的值，$\boldsymbol M$ 是方差向量的渐近协方差。令

$$S=\sum_{k=1}^r\frac{v_k}{v_k-2}I[v_k>2]$$

其中 $I[v_k>2]$ 是一个指示函数，当 $v_k > 2$ 时值为 $1$，否则为 $0$. $F_c$ 分布的近似分母自由度为

$$v=\begin{cases}\frac{2S}{S-r}&\mathrm{if~}S>r\\0&\mathrm{if~}S\leq r&\end{cases}$$

上述过程可用于通过选择一组 $t - 1$ 个线性独立的 $\mu_i$ 的对比来提供等均值假设的检验。

SAS-Mixed 程序实现了一种对分母自由度的近似，该分母自由度与分子中具有多个自由度的近似 $F$ 统计量相关。SAS-Mixed 可用于拟合每个处理组方差不等的模型，或使用 REPEATED 语句并指定 `GROUP =` 选项。表 \@ref(tab:table2-9) 中的 Mixed 代码将不等方差模型拟合至表 \@ref(tab:table2-1) 中的数据。REPEATED 语句用于指定要为每个处理估计不同的方差（组的每个值）。三个 Estimate 语句用于提供与 \@ref(sec2-6) 节中的三个问题相对应的计算。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table2-9)将不等方差模型拟合到表 \@ref(tab:table2-1) 中的数据的 SAS-Mixed 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> proc mixed cl covtest data=task;\
class group;\
model errors=group/ddfm=kr;\
repeated/group=group;\
estimate “part(1)” group -1 -1 -1 3/divisor=3 cl alpha=0.05;\
estimate “part(2)” group 2 -1 -1 0/divisor=2 cl alpha=0.05;\
estimate “part(3)” group 0 1 -1 0/cl alpha=0.05;\
lsmeans group/diff cl; </td>
  </tr>
</tbody>
</table>

Mixed 程序的结果见表 \@ref(tab:table2-10)，其中协方差参数估计值是四个处理方差的估计值，拟合统计量中的 AIC 是 Akaike Information Criteria (Akaike, 1974)，模型似然比检验提供了等方差假设的检验，固定效应的 III 型检验使用第二个统计量和分母的相应近似自由度提供了等均值假设的检验，估计值包含与 \@ref(sec2-6) 节中的三个问题对应的结果，其中提供了 $t$ 统计量、近似分母自由度和 $95\%$ 置信区间。表 \@ref(tab:table2-11) 包含了处理均值估计及其相应的标准误估计。名义自由度是与它们各自的方差相对应的自由度。表 \@ref(tab:table2-11) 的第二部分包含处理均值的成对比较，包括每次比较的近似分母自由度。该模型可以通过对药物1和两种药物使用一个方差，对药物2和无药物使用一种方差来简化。这可以通过定义一个变量来实现，比如T，对于药物1和两种药物都是1，对于其他两种处理是0。然后在模型指定中将T放在class语句中，并使用 `Repeated/Group=T;`。两个方差的估计值分别为 $2.4028$ 和 $12.7376$，$AIC$ 为 $126.4$，这比四方差模型的 $AIC$ 值小，表明两方差模型足以描述数据。使用在模型指定中具有较少方差的模型为相应的标准误提供了更多的自由度，从而提供了关于模型中固定效应的假设的更有效的检验。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table2-10)将不等方差模型拟合到表 \@ref(tab:table2-1) 中的数据的结果</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%202.10-1.png">
</td>
  </tr>
  <tr>
   <td style="text-align:center;">  <img src="table/table%202.10-2.png">
</td>
  </tr>
  <tr>
   <td style="text-align:center;">  <img src="table/table%202.10-3.png">
</td>
  </tr>
  <tr>
   <td style="text-align:center;">  <img src="table/table%202.10-4.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table2-11)利用不等方差模型估计药物组均值和配对比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%202.11-1.png">
</td>
  </tr>
  <tr>
   <td style="text-align:center;">  <img src="table/table%202.11-2.png">
</td>
  </tr>
</tbody>
</table>

## 结束语 {#sec2-9}

总之，为了比较所有均值，建议采用以下方法：

1) 如果方差齐性检验在 $1\%$ 的水平上不显著，则进行通常的方差分析检验。
2) 如果方差齐性检验在 $1\%$ 的水平上显著，则使用 Welch's test 或混合模型检验以及相应的近似分母自由度。
3) 如果方差齐性检验在 $1\%$ 的水平上显著，则使用 $AIC$ 来确定是否可以使用更简单或更少数量的方差来充分描述数据，以增加与均值有关检验的功效。

许多教科书和文章都写过关于对数据进行转换以实现相等的处理方差，从而可以使用通常的方差分析来比较处理。由于能够建立一个不等方差模型来提供均值的标准误估计和均值的比较，许多情况下不需要使用变换。不必使用变换来实现相等方差的一个主要好处是，平均值的单位是测量单位，从而简化了解释。

本章讨论了具有异质误差的单向方差分析模型的统计分析。讨论包括确定误差方差齐性的几个统计检验，以及何时使用每个检验的建议。讨论了在发现误差方差不齐时，对不同处理的效应进行统计推断的适当程序，还回顾了说明这些程序使用的示例。

## 练习 {#sec2-10}
