
# 分析子类数不等的大型均衡双向实验 {#chap11}

> "The science of statistics is the chief instrumentality through which the progress of civilization is now measured, and by which its development hereafter will be largely controlled." - S. N. D. North

在本章中，我们提出了一种方法，用于获得均衡处理结构实验的近似分析，这些实验在每个单元格中具有不等数量的观测结果。

## 可行性问题 {#sec11-1}

我们通常建议使用通用计算包，如第 \@ref(sec10-7) 节中讨论的计算包，来分析不均衡的数据集，其中每个处理组合至少观察一次。然而，在实践中出现了许多这样做可能不可行的情况，特别是在发展中国家，人们可能有几个处理因素，每个因素都有几个不同的水平。在这些情况下，精确的程序需要进行多次矩阵求逆，而要求逆的矩阵的大小可能超过研究人员和可用计算机的能力。

例如，考虑一个有四个因素的实验，每个因素都发生在五个水平上。为了获得精确的分析，必须能够对几个大矩阵求逆。在这种情况下，为了进行精确分析，需要求逆的最大矩阵的大小是不均衡情况下的 624×624 矩阵。显然，这样的实验并不罕见。此外，一些大型实验可能需要良好的初始起始值，以便迭代算法能够收敛到合理的解。本章中的方法可以让研究人员在难以使程序收敛的情况下找到良好的初始起始值。

分析这些类型的混乱数据情况的另一种方法称为未加权均值法 (methods of unweighted means)（Bancroft 1968）。第 \@ref(sec11-2) 节对双向处理结构的方法进行了描述。虽然人们很少遇到无法通过一般计算程序进行分析的双向实验，但对于双向实验，未加权均值法很容易被讨论和理解，并且这种讨论很容易推广到更大的实验情况。在继续之前，我们想指出的不是测量到的观测结果的数量使一般程序不可行，而是正在研究的处理组合的数量。

## 未加权均值法 {#sec11-2}

基本上，未加权均值法通过使用各种处理组合的观测均值来近似对应于每个效应的平方和。这里给出的公式是用于等样本量的公式；然而，如果样本量只是略有不同，它们仍然可以给出相当准确的结果。正确的平方和可以通过使用 III 型分析获得，本章中给出的公式近似于 III 型平方和。

考虑第 \@ref(chap9) 章和第 \@ref(chap10) 章中描述的双向模型，其中因素 T 有 t 种可能性，因素 B 有 b 种可能性。令 $\mu_{ij}$ 代表人们期望在随机选择的实验单元上观测的响应，该实验单元接收 $T_i$ 与 $B_j$ 交叉的组合。令

$$\hat{\mu}_{ij}=\bar{y}_{ij\cdot}\quad i=1,2,\ldots,t,j=1,2,\ldots,b$$

并令

$$\hat{\bar{\mu}_{i\cdot}}=\frac1b\sum_{j=1}^b\hat{\mu}_{ij}\quad\mathrm{~and~}\quad\hat{\bar{\mu}}_{\cdot j}=\frac1t\sum_{i=1}^t\hat{\mu}_{ij}$$

对于检验 $H_{0T}\colon\bar{\mu}_{1\cdot}=\bar{\mu}_{2\cdot}=\cdots=\bar{\mu}_{t\cdot}$，我们计算

$$\begin{aligned}SST=b\sum_{i=1}^t(\hat{\bar{\mu}}_{i\cdot}-\hat{\bar{\mu}}_{\cdot\cdot})^2=\sum_{i=1}^tb\hat{\bar{\mu}}_{i\cdot}^2-bt\hat{\bar{\mu}}_{\cdot\cdot}^2\end{aligned}$$

基于 t-1 个自由度。对于检验 $H_{0B}{\colon}\bar{{\mu}}_{\cdot1}=\bar{{\mu}}_{\cdot2}=\cdots=\bar{{\mu}}_{\cdot b}$，我们计算

$$SSB=t\sum_{j=1}^b(\hat{\bar{\mu}}_{\cdot j}-\hat{\bar{\mu}}_{\cdot\cdot})^2=t\sum_{j=1}^b\hat{\bar{\mu}}_{\cdot j}^2-bt\hat{\bar{\mu}}_{\cdot\cdot}^2$$

基于 b-1 个自由度。对于所有 $i\ne i^\prime$ 和 $j\ne j^\prime$ 的检验 $H_{0T\times B}{:{\mu}_{ij}-{\mu}_{i'j}-{\mu}_{ij'}+{\mu}_{i'j'}=0}$，我们计算

$$SST\times B=\sum_{i=1}^t\sum_{j=1}^b(\hat{\mu}_{ij}-\hat{\bar{\mu}}_{i\cdot}-\hat{\bar{\mu}}_{\cdot j}+\hat{\bar{\mu}}_{\cdot\cdot})^2=\sum_{i=1}^t\sum_{j=1}^b\hat{\mu}_{ij}^2-b\sum_{i=1}^t\hat{\bar{\mu}}_{i\cdot}^2-t\sum_{j=1}^b\hat{\bar{\mu}}_{\cdot j}^2+bt\hat{\bar{\mu}}_{\cdot\cdot}^2$$

基于 (t-1)(b-1) 个自由度。由于上述平方和是根据平均值计算的，因此通常的误差平方和

$$SSError=\sum_{i=1}^t\sum_{j=1}^{b}\sum_{k=1}^{n_{ij}}(y_{ijk}-\bar{y}_{ij\cdot})^2$$

需要调整。因为 $\hat\mu_{ij}$ 的方差由 $\sigma^2/n_{ij}$ 给出而不是由 $y_{ijk}$ 的方差 $\sigma^2$ 给出。该调整是通过将误差平方和除以

$$\tilde{n}=\left[\frac1{bt}{\left(\sum_{i=1}^{t}\sum_{j=1}^{b}\frac1{n_{ij}}\right)}\right]^{-1}$$

单元格样本量的调和平均值以给出。

$\tilde{n}$ 是样本量的调和平均值，它是 $n_{ij}$ 的一种可能平均值。误差自由度仍然是 N - bt. 表 \@ref(tab:table11-1) 给出了未加权均值分析的方差分析表。当单元格样本量相差不大时，此分析会产生 $F$ 分布的合理近似值。通常的建议是，如果样本量的变化不超过 2 倍，则该分析是可以接受的。

<table>
<caption>(\#tab:table11-1)未加权均值的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2011.1.png">
</td>
  </tr>
</tbody>
</table>

## 同时推断与多重比较 {#sec11-3}

T 边际均值和 B 边际均值的定义与第 \@ref(chap9) 章和第 \@ref(chap10) 章相同。因此，T 边际均值由 $\bar{{\mu}}_{i\cdot},i=1,2,\ldots,t$ 给出，和B的边际平均值由 $\bar{{\mu}}_{\cdot j},j=1,2,\ldots,b$ 给出。这些边际均值的最佳估计分别为

$$\hat{\bar{{\mu}}}_{i\cdot}=\frac1b\sum_{j=1}^b\hat{{\mu}}_{ij},\quad i=1,2,\ldots,t$$

以及

$$\hat{\bar{{\mu}}}_{\cdot j}=\frac1b\sum_{i=1}^t\hat{{\mu}}_{ij},\quad i=1,2,\ldots,t$$

$\hat{\bar{\mu}}_{i\cdot}$ 标准误的精确估计 (exact estimated standard error) 为

$$\widehat{s.e.}(\hat{\bar{\mu}}_{i\cdot})=\frac{\hat{\sigma}}b\sqrt{\sum_{j=1}^b\frac1{n_{ij}}}$$

$\hat{\bar{\mu}}_{\cdot j}$ 标准误的精确估计 (exact estimated standard error) 为

$$\widehat{s.e.}(\hat{\bar{\mu}}_{\cdot j})=\frac{\hat{\sigma}}t\sqrt{\sum_{i=1}^t\frac1{n_{ij}}}$$

其中

$$\hat{\sigma}=\sqrt{\frac{\sum_{i=1}^{t}\sum_{j=1}^{b}\sum_{k=1}^{n_{ij}}(y_{ijk}-\bar y_{ij\cdot})^2}{N-tb}}=\sqrt{ErrorMS}$$

如果单元格样本量相差不大，则这些标准误可以分别用 $\hat \sigma/\sqrt{}(b\tilde n)$ 和 $\hat \sigma/\sqrt{}(t\tilde n)$ 近似。$\hat{\bar{\mu}}_{i\cdot}-\hat{\bar{\mu}}_{i'\cdot}$ 的标准误估计为

$$\widehat{s.e.}(\hat{\bar{\mu}}_{i\cdot}-\hat{\bar{\mu}}_{i'\cdot})=\frac{\hat{\sigma}}{b}\sqrt{\sum_{j=1}^b\frac1{n_{ij}}+\sum_{j=1}^b\frac1{n_{i'j}}}$$

可用 $\hat \sigma\sqrt{}2/\sqrt{}(b\tilde n)$ 近似。类似地，

$$\widehat{s.e.}(\hat{\bar{\mu}}_{\cdot j}-\hat{\bar{\mu}}_{\cdot j'})=\frac{\hat{\sigma}}{t}\sqrt{\sum_{i=1}^t\frac1{n_{ij}}+\sum_{i=1}^t\frac1{n_{ij'}}}$$

可用 $\hat \sigma\sqrt{}2/\sqrt{}(t\tilde n)$ 近似。

接下来考虑对于 $i\ne i^\prime$ 和 $j\ne j^\prime$ 的 2×2 交互对比 $\mu_{ij}-\mu_{i'j}-\mu_{ij'}+\mu_{i'j'}$. 这种交互对比的最佳估计为 $\hat\mu_{ij}-\hat\mu_{i'j}-\hat\mu_{ij'}+\hat\mu_{i'j'}$，其标准误估计为

$$\widehat{s.e.}(\hat{\mu}_{ij}-\hat{\mu}_{i^{\prime}j}-\hat{\mu}_{ij^{\prime}}+\hat{\mu}_{i^{\prime}j^{\prime}})=\hat{\sigma}\sqrt{\frac1{n_{ij}}+\frac1{n_{i^{\prime}j}}+\frac1{n_{ij^{\prime}}}+\frac1{n_{i^{\prime}j}}}$$

如果样本量相差不大，可用 $\hat \sigma/\sqrt{}(4/\tilde n)$ 近似。

对于单元格均值参数的上述函数中的每一个，检验统计量由下式给出

$$t_c=\frac{\mathrm{estimate}}{\text{estimated standard error}}$$

其中如果 $|t_c|>t_{\alpha/2,N-tb}$，则拒绝相应的假设。此外，一个 $(1-\alpha)100\%$ 置信区间为

$$\mathrm{estimate}\pm t_{\alpha/2,N-tb}(\text{estimated standard error})$$

可以使用第 \@ref(chap3) 章中描述的同时推断程序，并鼓励研究人员这样做。这些想法可以推广到单元格均值和/或边际均值的更一般的对比中。如果简单地用 $\tilde n$ 代替 $n$，则可以近似第 \@ref(chap8) 章中的所有公式。

## 未加权均值的示例 {#sec11-4}

<table>
<caption>(\#tab:table11-2)表 \@ref(tab:table9-1) 中的单元格均值和边际均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2011.2.png">
</td>
  </tr>
</tbody>
</table>

作为示例，请再次考虑表 \@ref(tab:table9-1) 中的数据。回想一下，误差平方和为 20，自由度为 10，并且均值表在表 \@ref(tab:table11-2) 中重复。$\tilde n$ 的值为

$$\tilde{n}=\left[\frac1{bt}\left(\sum_{i=1}^t\sum_{j=1}^b\frac1{n_{ij}}\right)\right]^{-1}=\left[\frac1{(3)(2)}\left(\frac13+\frac12+\frac13+\frac12+\frac13+\frac13+\frac13\right)\right]^{-1}=\left(\frac7{18}\right)^{-1}=2.5714$$

接着注意到

$$\begin{aligned}
&\sum_{i,j}{\hat{\mu}}_{ij}^2 =20^2+25^2+\cdots+28^2=3830  \\
&\sum_{i,j}\hat{\bar{\mu}}_{i\cdot}^2 =23^2+27^2=1258  \\
&\sum_j\hat{\bar{\mu}}_{\cdot j}^2 =23^2+24^2+28^2=1889 \\
& \hat{\bar{\mu}}_{\cdot \cdot}^2=(25)^2=625
\end{aligned}$$

因此

$$\begin{aligned}SST&=3(1258)-6(625)=24\\SSB&=2(1889)-6(625)=28\\SST \times B &= 3830 - 3(1258) - 2(1889) + 6(625) = 28\end{aligned}$$

方差分析表见表 \@ref(tab:table11-3)。请注意，表 \@ref(tab:table11-3) 中的 $F$ 统计量与表 \@ref(tab:table9-3) 中给出的精确 $F$ 统计量非常相似。

<table>
<caption>(\#tab:table11-3)表 \@ref(tab:table9-1) 中数据未加权均值方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2011.3.png">
</td>
  </tr>
</tbody>
</table>

## 计算机分析 {#sec11-5}

尽管我们不推荐这样做，但未加权均值分析所需的统计量可以在 Excel^®^ 和其他电子表格程序中获得。通过将实验视为单向实验并利用均值模型，可以最有效地获得未调整误差平方和。要获得未加权均值分析所需的 T, B 和 T × B 平方和，可以按照以下步骤操作。

1. 获取每种处理组合的单元格均值。
2. 获取 T 边际均值、B 边际均值和总体均值。
3. 计算单元格均值、T 边际均值和 B 边际均值的方差。分别记为 $S_{TB}^2,S_{B}^2,S_{T}^2$.
4. 那么

\begin{align}
\text{error } SS& =\sum_{i,j,k}y_{ijk}^2-(bt-1)S_{TB}^2+bt\hat{\bar{\mu}}_{\cdot\cdot}^2  \\
SST& =b(t-1)S_T^2  \\
SSB& =t(b-1)S_B^2 \\
SST\times B&=(bt-1)S_{TB}^2-b(t-1)S_T^2-t(b-1)S_B^2
(#eq:11-1)
\end{align}

如果有可以执行单向方差分析的程序，以及一个可以针对均衡数据进行双向方差分析的程序，那么她可以通过将双向处理结构视为使用均值模型的单向处理结构，从单向方差分析程序中获取误差平方和。通过用双向方差分析程序分析单元格均值，可以获得 T, B和 T×B 的平方和。表 \@ref(tab:table11-4) 给出了 SAS^®^ 代码，该代码提供执行未加权均值分析所需的统计量。如果读者希望运行此代码，那么将由她自己决定。

<table>
<caption>(\#tab:table11-4)提供未加权均值分析所需统计量的 SAS 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2011.4.png">
</td>
  </tr>
</tbody>
</table>

## 结束语 {#sec11-6}

在本章中，我们介绍了一种方法，用于对涉及大量不同处理组合的实验进行满意的统计分析，其中每个处理组合至少观察一次。这些技术对于可能无法使用高级统计软件程序的人员来说是有用的。
