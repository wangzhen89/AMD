
# 使用效应模型分析子类数不等的均衡双向处理结构 {#chap10}

> "It is a capital mistake to theorize before one has data." - Sir Arthur Conan Doyle

第 \@ref(chap9) 章讨论了使用均值模型分析具有不等子类数的双向处理结构。本章考虑在相同情况下使用效应模型。所有可以使用效应模型回答的问题也可以使用均值模型回答，反之亦然。之所以讨论效应模型，是因为当使用统计计算软件包分析双向处理结构时，效应模型常常是一个重要的工具。因为统计软件通常会自动产生主效应和交互效应的检验统计量，以及边际均值、双向均值及其标准误估计。

## 模型定义 {#sec10-1}

与均值模型 \@ref(eq:9-1) 对应的效应模型定义如下

\begin{equation}
y_{ijk}=\mu+\tau_i+\beta_j+\gamma_{ij}+\varepsilon_{ijk},\quad i=1,2,\ldots,t;~j=1,2,\ldots,b;~k=1,2,\ldots,n_{ij}
(#eq:10-1)
\end{equation}

其中 ${\varepsilon}_{ijk}\thicksim i.i.d.N(0,{\sigma}^2)$.

## 参数估计和 I 型分析 {#sec10-2}

除了第 \@ref(chap9) 章介绍的那些平方和，在分析具有不等子类数量的双向处理结构时，还会涉及其他一些平方和。本章将探讨其中两种。首先是一系列平方和，它们涉及通过推广第 \@ref(sec3-2) 节所述的模型比较方法，以序贯 (sequential) 的方式将双向效应模型逐步拟合到观测数据上。常用的一种序列步骤如下

1. 第一步. 拟合 $y_{ijk}=\mu+\varepsilon_{ijk}$，将其残差平方和表示为 $RSS_1$.
2. 第二步. 拟合 $y_{ijk}=\mu+\tau_i+\varepsilon_{ijk}$，将其残差平方和表示为 $RSS_2$.
3. 第三步. 拟合 $y_{ijk}=\mu+\tau_i+\beta_j+\varepsilon_{ijk}$，将其残差平方和表示为 $RSS_3$.
4. 第四步. 拟合 $y_{ijk}=\mu+\tau_i+\beta_j+\gamma_{ij}+\varepsilon_{ijk}$，将其残差平方和表示为 $RSS_4$.

$RSS_i$ 为第 $i(i=1,2,3,4)$ 步拟合模型之后的残差平方和。$RSS_1$ 与 $RSS_2$ 之间的差记作 $R(\tau|\mu)$，称为调整 $\mu$ 后 $\tau$ 引起的缩减。即 $R(\tau|\mu)=RSS_1-RSS_2$. 这个缩减给出了通过考虑将 $\tau_i$ 也包含在内的模型，可以减小第一步中模型的残差平方和。$R(\tau|\mu)$ 越大，模型中包含 $\tau_i$ 就越重要。因此，$R(\tau|\mu)$ 处理 T 不同可能性的效应的度量。$R(\beta|\mu,\tau)=RSS_2-RSS_3$ 称为同时调整 $\mu$ 和 $\tau$ 后 $\beta$ 引起的缩减。它表示通过在模型中进一步包含 $\beta_j$，能够额外减少第二步模型残差平方和的量。$R(\beta|\mu,\tau)$ 是对处理 B 不同可能性的效应的一种度量，该效应独立于并超越了处理 T 的效应。最后，$R(\gamma|\mu,\tau,\beta)=RSS_3-RSS_4$ 称为同时调整 $\mu,\tau,\beta$ 后 $\gamma$ 引起的缩减。它表示通过在模型的第三步中添加交互参数 $\gamma_{ij}$，可以进一步减少模型残差平方和的额外量。显然，$R(\gamma|\mu,\tau,\beta)$ 是交互作用的度量，因为第三步中的模型是一个加性模型，当且仅当没有交互作用时才成立。

表 \@ref(tab:table10-1) 给出了与**序贯分析** (sequential analysis) 对应的方差分析表。这种分析称为 I 型分析。表 \@ref(tab:table10-1) 最后四行的平方和在统计上是独立的，并且 T, B 和 T × B 均方与误差均方的比都服从非中心 $F$ 分布。确定表 \@ref(tab:table10-1) 中每个 $F$ 统计量所检验的确切假设非常有趣且具有信息价值。所检验的假设在第 \@ref(sec10-4) 节中给出。

<table>
<caption>(\#tab:table10-1)序贯分析的方差分析表（I 型分析）</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.1.png">
</td>
  </tr>
</tbody>
</table>

为了说明这一点，对于 I 型分析所需的四个模型，均需拟合到表 \@ref(tab:table9-1) 中的数据。要理解此处进行的计算，需要熟悉第 \@ref(chap6) 章的内容。然而，对于那些愿意利用统计计算软件包执行所需计算的读者来说，这种理解不是必需的。对详细内容感兴趣的读者可以考虑接下来的几页内容进行深入学习和探讨。第一步中给出的模型为 $y_{ijk}={\mu}+{\varepsilon}_{ijk}$. 该模型中 $\mu$ 的最佳估计是所有观测值的平均值 $\hat{{\mu}}=\bar{y}_{\cdot\cdot\cdot}=24.9375$，残差平方和为

$$RSS_1=\sum_{i,j,k}(y_{ijk}-\hat{\mu})^2=\sum_{i,j,k}y_{ijk}^2-n_{\cdot\cdot}\bar{y}^2_{\cdot\cdot\cdot}=10209-16(24.9375)^2=258.9375$$

该残差平方和基于 $n_{\cdot\cdot} = 16 - 1 = 15$ 个自由度。

第二步中定义的模型的正规方程组为：

$$\begin{bmatrix}16&8&8\\8&8&0\\8&0&8\end{bmatrix}\begin{bmatrix}\hat{{\mu}}\\\hat{{\tau}}_1\\\hat{{\tau}}_2\end{bmatrix}=\begin{bmatrix}399\\182\\217\end{bmatrix}$$

这些方程的一种可能解是通过使用第 \@ref(sec6-2) 节中讨论的置零限制获得的，可得 $\hat{{\tau}}_2=0,\hat{{\tau}}_1=-4.375,\hat \mu=27.125$ （回想第 \@ref(chap6) 章，不存在唯一解）。残差平方和为

$$
\begin{aligned}RSS_2=\sum_{i,j,k}(y_{ijk}-\hat{\mu}-\hat{\tau}_i)^2&=\sum_{i,j,k}y_{ijk}^2-\hat{\mu}\cdot y_{\cdot\cdot\cdot}-\hat{\tau}_1\cdot y_{1\cdot\cdot}-\hat{\tau}_{2} \cdot y_{2\cdot\cdot} \\
&=10209-27.125\cdot399-(-4.375)\cdot 182-0\cdot 217 \\
&=182.375
\end{aligned}
$$

基于 16-2=14 个自由度。因此，$R(\tau|\mu)=258.9375-182.375=76.5625$ 基于 15-14=1 个自由度。

第三步中模型的正规方程组为

$$\begin{bmatrix}16&8&8&5&5&6\\8&8&0&3&2&3\\8&0&8&2&3&3\\5&3&2&5&0&0\\5&2&3&0&5&0\\6&3&3&0&0&6\end{bmatrix}\begin{bmatrix}\hat{\mu}\\\hat{\tau}_1\\\hat{\tau}_2\\\hat{\beta}_1\\\hat{\beta}_2\\\hat{\beta}_3\end{bmatrix}=\begin{bmatrix}399\\182\\217\\112\\119\\168\end{bmatrix}$$

为获得一个解，可以令 $\tau_2 = 0$ 且 $\beta_3 = 0$（参见第 \@ref(chap6) 章）。然后，通过删除对应于 $\tau = 0$ 和 $\beta = 0$ 的行和列，可以将该方程组简化为等价的方程组。该系统简化为：

$$\begin{bmatrix}16&8&5&5\\8&8&3&2\\5&3&5&0\\5&2&0&5\end{bmatrix}\begin{bmatrix}\hat{\mu}\\\hat{\tau}_1\\\hat{\beta}_1\\\hat{\beta}_2\end{bmatrix}=\begin{bmatrix}399\\182\\112\\119\end{bmatrix}$$

该简化系统的解为 $\hat{\mu}=30.154,\hat{\tau}_1=-4.308,\hat{\beta}_1=-5.169,\hat{\beta}_2=-4.631$。该模型的残差平方和为

$$\begin{aligned}
RSS_3& =\sum_{i,j,k}y_{ijk}^2-\hat{\mu}\cdot y_{\cdot\cdot\cdot}-\hat{\tau}_1\cdot y_{1\cdot\cdot}-\hat{\beta}_1\cdot y_{\cdot1\cdot}-\hat{\beta}_2\cdot y_{\cdot2\cdot}  \\
&=10209-30.154\cdot399-(-4.308)\cdot182-(-5.169)\cdot112-(-4.631)\cdotp119 \\
&=91.631
\end{aligned}$$

基于 16-4=12 个自由度。因此 $R(\beta|\mu,\tau)=RSS_2-RSS_3=182.375-91.631=90.744$ 基于 14-12=2 个自由度。

第四步中模型的正规方程组为

$$\begin{bmatrix}16&8&8&5&5&6&3&2&3&2&3&3\\8&8&0&3&2&3&3&2&3&0&0&0\\8&0&8&2&3&3&0&0&0&2&3&3\\5&3&2&5&0&0&3&0&0&2&0&0\\5&2&3&0&5&0&0&2&0&0&3&0\\6&3&3&0&0&6&0&0&3&0&0&3\\3&3&0&3&0&0&3&0&0&0&0&0\\2&2&0&0&2&0&0&2&0&0&0&0\\3&3&0&0&0&3&0&0&3&0&0&0\\2&0&2&2&0&0&0&0&0&2&0&0\\3&0&3&0&3&0&0&0&0&0&3&0\\3&0&3&0&0&3&0&0&0&0&0&3\end{bmatrix}\left.\left[\begin{array}{c}\hat{\mu}\\\hat{\tau}_1\\\hat{\tau}_2\\\hat{\beta}_1\\\hat{\beta}_2\\\hat{\beta}_3\\\hat{\gamma}_{11}\\\hat{\gamma}_{12}\\\hat{\gamma}_{13}\\\hat{\gamma}_{21}\\\hat{\gamma}_{22}\\\hat{\gamma}_{23}\end{array}\right.\right]=\begin{bmatrix}399\\182\\217\\112\\119\\168\\60\\50\\72\\52\\69\\96\end{bmatrix}$$

为获得一个解，令 $\hat{\tau}_2=0,\hat{\beta}_3=0,\hat{\gamma}_{13}=0,\hat{\gamma}_{21}=0,\hat{\gamma}_{22}=0,\hat{\gamma}_{23}=0$ （参见第 \@ref(chap6) 章）。使用缩减技术，系统简化为

$$\begin{bmatrix}16&8&5&5&3&2\\8&8&3&2&3&2\\5&3&5&0&3&0\\5&2&0&5&0&2\\3&3&3&0&3&0\\2&2&0&2&0&2\end{bmatrix}\begin{bmatrix}\hat{\mu}\\\hat{\tau}_1\\\hat{\beta}_1\\\hat{\beta}_2\\\hat{\gamma}_{11}\\\hat{\gamma}_{12}\end{bmatrix}=\begin{bmatrix}399\\182\\112\\119\\60\\50\end{bmatrix}$$

该简化系统的解为 $\hat{\mu}=32,\hat{\tau}_1=-8,\hat{\beta}_1=-6,\hat{\beta}_2=-9,\hat{\gamma}_{11}=2,\hat{\gamma}_{12}=10$. 因此 $\hat\sigma^2=2$，因为完全模型 (full model) 的残差平方和为 $RSS_4=10,209-10,189=20$，基于 10 个自由度。而且 $R(\gamma|\mu,\tau,\beta)=RSS_3-RSS_4=91.631-20=71.631$ 基于 12-10=2 个自由度。所有前面的结果都可以总结在方差分析表中，如表 \@ref(tab:table10-2) 所示。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table10-2)I 型方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.2.png">
</td>
  </tr>
</tbody>
</table>

I 型分析中的平方和以及交互作用的检验统计量与第 \@ref(chap9) 章中使用均值模型和矩阵程序获得的结果相同；然而，这两个程序对于 T 和 B 主效应给出了不同的平方和和检验统计量。双向处理结构中的数据也可以通过先拟合 $\mu$ 后 $\beta$ 再 $\tau$，最后 $\gamma$ 来分析。表 \@ref(tab:table10-2) 中尚未给出的唯一所需的新平方和是 $R(\beta|\mu)$ 和 $R(\tau|\mu,\beta)$. 对于不均衡数据情况，T 和 B 主效应的相应 $F$ 检验通常与按表 \@ref(tab:table10-1) 和/或表 \@ref(tab:table10-2) 中给出的顺序拟合效应所获得的结果不同。

应回忆第 \@ref(chap6) 章中的参数估计 $\hat{{\mu}}=32,\hat{{\tau}}_1=-8,\hat{{\beta}}_1=-6,\hat{\beta}_2=-9,\hat{\gamma}_{11}=2,\hat{\gamma}_{12}=10$ 不是相应参数 $\mu,\tau_1,\beta_1,\beta_2,\gamma_{11},\gamma_{12}$ 的无偏估计。事实上，这些单独的参数是不可估的。在用于求解正规方程组的置零限制下，可以证明

\begin{align}
&\hat{\mu} \,\,\text{is an unbiased estimate of}\,\, \mu+\tau_2+\beta_3+\gamma_{23} \\
&\hat{\tau}_1 \,\,\text{is an unbiased estimate of}\,\, \tau_1-\tau_2+\gamma_{13}-\gamma_{23}  -\gamma_{23} \\
&\hat{\beta}_1 \,\,\text{is an unbiased estimate of}\,\, \beta_1-\beta_3+\gamma_{21}-\gamma_{23} -\gamma_{23} \\
&\hat{\beta}_2 \,\,\text{is an unbiased estimate of}\,\, \beta_2-\beta_3+\gamma_{22}-\gamma_{23} \gamma_{23}\\
&\hat{\gamma}_{11} \,\,\text{is an unbiased estimate of}\,\, \gamma_{11}-\gamma_{13}-\gamma_{21}+\gamma_{23}\\
&\hat{\gamma}_{12} \,\,\text{is an unbiased estimate of}\,\, \gamma_{12}-\gamma_{13}-\gamma_{22}+\gamma_{23}
(#eq:10-2)
\end{align}

有关可估函数及其估计的更多信息，请参阅下一节。

## 在 SAS 中使用可估函数 {#sec10-3}

SAS^®^-GLM 程序有一个选项，可用于识别模型参数的可估函数（不使用 SAS 的读者可以跳过本节）。由于使用示例来描述可估函数是最简单的，因此请再次考虑表 \@ref(tab:table9-1) 中的数据。可以使用以下语句获得该数据的 SAS-GLM 分析：


```sas
PROC GLM; 
CLASSES T B;
MODEL y = T B T*B/<selected options>;
```

许多选项可以与 MODEL 语句一起使用。其中最重要的选项之一是 `E` 选项。此选项要求 SAS-GLM 打印模型参数可估函数的一般形式。回想一下第 \@ref(chap6) 章，设计模型中参数的所有线性函数不一定都是可估的。使用 `E` 选项，SAS-GLM 会打印信息，这些信息可用于确定可估参数和不可估参数的线性组合。

SAS-GLM 给出的可估函数的一般形式如表 \@ref(tab:table10-3) 所示。这意味着模型参数的线性函数 $\boldsymbol\ell^\prime\boldsymbol \beta$ （其中 $\boldsymbol\ell^{\prime}=[\ell_1,\ell_2,\ldots,\ell_{12}]$ 和 $\boldsymbol\beta^{\prime}=[\mu,\tau_1,\tau_2,\beta_1,\beta_2,\beta_3,\gamma_{11},\gamma_{12},\gamma_{13},\gamma_{21},\gamma_{22},\gamma_{23}]$）是可估的当且仅当存在常数 $L1,L2,L4,L5,L7,L8$ 使得

$$\begin{aligned}
\boldsymbol{\ell'\beta} =
&({L1})\mu+({L2})\tau_1+({L1}-{L2})\tau_2+({L4})\beta_1+({L5})\beta_2+({L1}-{L4}-{L5})\beta_3+({L7})\gamma_{11}\\&+(\mathrm{L}8)\gamma_{12}+(\mathrm{L}2-\mathrm{L}7-\mathrm{L}8)\gamma_{13}+(\mathrm{L}4-\mathrm{L}7)\gamma_{21}+(\mathrm{L}5-\mathrm{L}8)\gamma_{22} \\
&+(L1-L2-L4-L5+L7+L8)\gamma_{23} 
\end{aligned}$$

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table10-3)可估函数的一般形式</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.3.png">
</td>
  </tr>
</tbody>
</table>

例如，从可估函数的一般形式，我们可以看到：

1. $\mu$ 是不可估的。因为要使得 $\mu$ 可估，至少需要同时满足 $L1 = 1,L2 = 0$ 以及 $L1 - L2 = 0$ 这三个条件，而这三个条件无法同时为真。
2. $\tau_1$ 是不可估的。因为要使得 $\tau_1$ 可估，至少需要同时满足 $L1 = 0、L2 = 0$ 以及 $L1 - L2 = 0$ 这三个条件，而这三个条件无法同时为真。
3. $\tau_1 - \tau_2$ 不可估。因为要使得 $\tau_1 - \tau_2$ 可估，至少需要同时满足 $L1 = 0,L2 = 1,L1 - L2 =-1,L4 = 0,L5 = 0,L1 - L4 - L5 = 0,L7 = 0,L 8 = 0$ 以及 $L2 - L7 -L 8 = 0$. 然而 $L2=1,L7=0,L8=0,L2-L7-L8$ 无法同时为真。

很明显，模型参数有许多函数是不可估的，那么模型参数有哪些有趣的函数是可估的呢？从第 \@ref(chap6) 章中，我们知道模型参数存在可估的函数。实际上，模型参数的一组可估函数基 (a basis set of estimable functions) 定义为一组线性独立的可估函数，这样的一组函数确保了任何其他的可估函数都能够表示为此集合内可估函数的线性组合。在可估函数的一般形式中，我们可以看到，我们可以自由地选择 $L$ 中的六个，即 $L1,L2,L4,L5,L7,L8$. 因此在一个基中有六个线性无关的可估函数，一个容易得到的基是依次让六个 $L$ 中的一个 等于 1，其余的 $L$ 都等于 0. 例如，取 $L1 = 1$ 和 $L2 = L4 = L5 = L7 = L8 = 0$，一般形式简化为 $\mu + \tau_2 + \beta_3 + \gamma_{23}$，这给出了一个模型参数的可估线性函数。取 $L2 = 1$ 其它的 $L$ 等于 0，一般形式简化为 $\tau_1-\tau_2+\gamma_{13}-\gamma_{23}$，这给出了令一个模型参数的可估线性函数。以如此方式继续可得到表 \@ref(tab:table10-4) 中给出的可估函数的一组基。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table10-4)可估函数的一组基</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.4.png">
</td>
  </tr>
</tbody>
</table>

请注意，在该基中的线性函数数量为六，这与第 \@ref(chap6) 章讨论的 $\boldsymbol X^\prime\boldsymbol X$ 矩阵的秩相等，并且也等于该 2×3 实验中处理组合的数量。另外，请注意表 \@ref(tab:table10-4) 中的可估函数与满足置零限制正规方程组的解所估计的函数完全相同（如方程 \@ref(eq:10-2) 所示）。这不是巧合；在使用 SAS-GLM 分析数据时，可以通过将一般形式中的每个 L 设为 1，其他 L 设为 0，从而找到由正规方程组置零解估计的函数，这总是成立的。

可估函数的基不是唯一的，表 \@ref(tab:table10-5) 给出了另一个可估函数的基，还给出了产生这个基的 $L1,L2,L4,L5,L7, L8$ 的值。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table10-5)可估函数的另一组基</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.5.png">
</td>
  </tr>
</tbody>
</table>

当在 SAS-GLM model 语句上使用 `SOLUTION` 选项时，计算机将使用置零限制打印出模型参数的最小二乘估计。`SOLUTION` 选项的结果如表 \@ref(tab:table10-6) 所示。注意，该 `SOLUTION` 选项的结果与第 \@ref(sec10-2) 节中给出的置零解一致，为方便起见，在此重复：

$$\hat{\mu}=32,\hat{\tau}_1=8,\hat{\beta}_1=-6,\hat{\beta}_2=-9,\hat{\gamma}_{11}=2,\hat{\gamma}_{12}=10$$

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;border-bottom: 0;">
<caption>(\#tab:table10-6)SAS-GLM 中使用 `SOLUTION` 选项获得的结果</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.6.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> 已发现 $\boldsymbol X^\prime \boldsymbol X$ 矩阵是奇异的，并且使用广义逆来求解正规方程。估计值后面跟着字母“B”的项不是唯一可估的</td></tr></tfoot>
</table>


[^table104]: The functions of the parameters that these estimators are really unbiased estimates of are those given in Table 10.4.

然而，如前所述，这些最小二乘估计值实际上并不直接估计它们各自对应的参数。事实上，正如先前所展示的那样，各个单独参数是不可估的。SAS-GLM 通过在表 \@ref(tab:table10-6) 中每个最小二乘估计旁边放置字母 B 来表明这一点。表 \@ref(tab:table10-4) 给出了这些估计是真正无偏估计的参数函数[^table104]。也就是说，$\hat{\mu}=32$ 是 $\mu+\tau_2+\beta_3+\gamma_{23}$ 的最佳无偏估计，$\hat{{\tau}}_1={-}8$ 是 $\tau_1-\tau_2+\gamma_{13}-\gamma_{23}$ 的最佳无偏估计，$\hat{{\tau}}_2=0$ 估计值为 0（这意味着它在估计零值上做得相当好），$\hat{{\beta}}_1=-6$ 是 $\beta_1-\beta_3+\gamma_{21}-\gamma_{23}$ 的最佳无偏估计，以此类推。表 \@ref(tab:table10-6) 中打印的标准误是估计的实际标准误估计。即，$\widehat {s.e.}( \hat{\mu})=0.8165,\widehat{s.e.}(\hat{\tau}_1)=1.1547$，以此类推。表 \@ref(tab:table10-6) 中的 $t$ 统计量检验被估计的相应函数是否等于零。例如对应于 $\hat \mu$ 的 $t=39.19$，检验 $H_0\colon\mu+\tau_2+\beta_3+\gamma_{23}=0$. 这样的检验通常不是很有趣。通过在 SAS-GLM 程序中使用 CONTRAST 或 ESTIMATE 语句，人们可以对选择考虑的参数的任何可估线性组合进行推断。线性组合不必是对比，但它必须是可估的。幸运的是，SAS-GLM 总是检查指定的线性组合是否是可估的函数。如果是，则 ESTIMATE 语句给出其最佳无偏估计、标准误估计并检验被估参数函数是否等于零的 $t$ 统计量和 $p$ 值。CONTRAST 语句允许同时检验几个可估函数，并给出 $F$ 统计量及其 $p$ 值，如第 \@ref(chap1) 章所述。在我们的例子中，正确使用 ESTIMATE 语句需要以下形式：


```sas
ESTIMATE 'label' INTERCEPT C1 T C2 C3 B C4 C5 C6 T*B C7 C8 C9 C10 C11 C12;
```

如果一个特定效应的所有系数都为零，则该效应及其系数不需要包含在 ESTIMATE 语句中。要使用 CONTRAST 语句，只需在上述形式中将单词 ESTIMATE 替换为 CONTRAST。例如，要获得表 \@ref(tab:table10-5) 中给出的可估函数的最佳估计，可以使用：

1. `ESTIMATE 'OVERALL MEAN' INTERCEPT 1 T -5 -5 B .33333 .33333 .33333 T*B .16667 .16667 .16667 .16667 .16667 .16667;`
2. `ESTIMATE 'TI -T2' T 1 —1 T*B .33333 .33333 .33333 —.33333 —.33333 -.33333;`
3. `ESTIMATE 'Bl -B3' B 1 0 -1 T*B .5 0 -.5 .5 0 -.5;`
4. `ESTIMATE 'B2 —B3' B 0 1 -1 T*B 0 .5 -.5 0 .5 -.5;`
5. `ESTIMATE 'INT1' T*B 1 0 -1 -1 0 1`;
6. `ESTIMATE 'INT2' T*B O 1 -1 0 -1 1`;

## I–IV 型假设 {#sec10-4}

本书的许多读者可能已经意识到，SAS-GLM 为用户提供了从四种类型的平方和中选择一种来检验假设的选项。本节主要涉及定义和解释相应的四种类型的假设。表 \@ref(tab:table9-1) 中的数据将用于说明这些假设。如第 \@ref(sec10-2) 节所述，I 型平方和是通过以序贯方式拟合双向效应模型获得的。每一步所获平方和是衡量该步添加特定项的重要性的指标，是通过在模型中包含该项可以减少残差平方和的量。II 型分析也是利用模型比较技术获得的。针对模型中处于相同或较低水平的每个其他效应，调整与每个效应对应的平方和。因此，对应于 T 效应的 II 型平方和为 $R(\tau|\mu,\beta)$，对应于 B 效应的 II 型平方和为 $R(\beta|\mu,\tau)$. 稍有困惑的读者应该参阅第 \@ref(sec16-3) 节，因为三向处理结构的 I 型和 II 型分析的定义将有助于澄清 I 型平方和与 II 型平方和之间的差异。表 \@ref(tab:table10-7) 显示了双向效应模型的 I 型和 II 型平方和，表 \@ref(tab:table10-8) 给出了表 \@ref(tab:table9-1) 中数据的 II 型平方和以及检验统计量。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table10-7)I 型和 II 型平方和的定义</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.7.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table10-8)II 型方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.8.png">
</td>
  </tr>
</tbody>
</table>

我们建议实验者根据均值模型中的参数而不是效应模型中的参数来思考。使用效应模型来拟合双向实验的主要优点是，许多有趣的假设能够被自动检验；这使得实验者进行相关工作变得更为便捷。效应模型的缺点是实验者很难准确理解 I 型和 II 型平方和检验的内容。均值模型的主要优点是能够准确地理解正在检验和/或估计的内容。均值模型的缺点是实验者必须编写自己的 ESTIMATE 和 CONTRAST 语句来获得感兴趣的检验统计量。**我们希望实验者对效应模型和均值模型都有所了解：效果模型可用于自动获得许多感兴趣的统计量，均值模型可用于识别感兴趣统计量所估计或检验的内容**。

[^table109]: 原文：Table 10.9 gives the hypotheses tested by a type I analysis of the data in Table 9.1 for model (10.1) in terms of the parameters in a means model.

I 型和 II 型平方和检验了哪些假设？用 I 型平方和检验的假设称为 I 型假设。表 \@ref(tab:table10-9) 给出了 I 型分析检验的假设，对应的模型是在均值模型参数框架下的模型 \@ref(eq:10-1)，数据来源于表 \@ref(tab:table9-1) [^table109] 。稍后我们将展示如何从 SAS-GLM 程序中确定这些假设。显然，实验者很少对表 \@ref(tab:table10-9) 中对应于 T 和 B 的假设感兴趣，并且随着表 \@ref(eq:9-1) 中一个或多个单元格中样本量的变化，假设也会发生变化。与 T×B 对应的假设相当于检验无交互作用假设。也就是说，在这个特定示例中，对比

$$\mu_{11}-\mu_{13}-\mu_{21}+\mu_{23}\quad\mathrm{and}\quad\mu_{12}-\mu_{13}-\mu_{22}+\mu_{23}$$
 
张成交互作用空间 (span the iteraction space). 表 \@ref(tab:table10-10) 给出了以均值模型参数表示的双向实验的 I 型假设的一般公式。注意，从一般形式中我们可以看到，T 的 I 型假设使用每个单元格内的样本量作为权重来比较行单元格均值的加权平均值。例如，对于表 \@ref(tab:table9-1) 中的数据，T 的 I 型假设可以写成

$$H_0{:}\frac{3\mu_{11}+2\mu_{12}+3\mu_{13}}8=\frac{2\mu_{21}+3\mu_{212}+3\mu_{213}}8$$

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table10-9)表 \@ref(tab:table9-1) 中数据的均值模型 I 型分析的假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.9.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table10-10)用均值模型参数表示的双向效应模型 I 型假设的一般形式</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.10.png">
</td>
  </tr>
</tbody>
</table>

对于表 \@ref(tab:table9-1) 中的数据，表 \@ref(tab:table10-11) 中给出了根据平均模型参数进行 II 型分析检验的假设。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table10-11)表 \@ref(tab:table9-1) 中数据的均值模型进行 II 型分析的假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.11.png">
</td>
  </tr>
</tbody>
</table>

对于一般模型，T 效应对应的行检验的 II 型假设为

\begin{equation}
\sum_{j=1}^b\left(n_{ij}-\frac{n_{ij}^2}{n_{\cdot j}}\right)\mu_{ij}=\sum_{i^{\prime}=1}^t\sum_{j=1}^b\frac{n_{ij}n_{i^{\prime}j}}{n_{\cdot j}}\mu_{i^{\prime}j},\quad i=1,2,\ldots,t
(#eq:10-3)
\end{equation}

表 \@ref(tab:table10-12) 和表 \@ref(tab:table10-13) 给出了根据方程 \@ref(eq:10-1) 中的效应模型参数对表 \@ref(tab:table9-1) 中的数据进行的 I 型和 II 型假设。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table10-12)表 \@ref(tab:table9-1) 中数据的效应模型进行 I 型分析的假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.12.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table10-13)表 \@ref(tab:table9-1) 中数据的效应模型进行 II 型分析的假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.13.png">
</td>
  </tr>
</tbody>
</table>

对表 \@ref(tab:table10-10)-\@ref(tab:table10-13) 的检查表明，I 型和 II 型分析检验的主效应假设可能不是很有趣。此外，拒绝或接受这些假设可能不容易解释。

计算效应模型平方和的第三种方法如下：

1. 对于处理 T 的 t 个水平，生成 t - 1 个哑变量 (dummy variables)，对于处理 B 的 b 个水平，生成 b - 1 个哑变量。创建与 $\tau$ 和 $\beta$ 对应的哑变量，以便在具有哑变量的模型中有 $\tau_t=-(\tau_1+\tau_2+\cdots+\tau_{t-1})$ 和 $\beta_t=-(\beta_1+\beta_2+\cdots+\beta_{t-1})$.
2. T 和 B 之间的交互作用由它们相应的哑变量的乘积表示。特别的，对于 $j=1,2,\cdots,b$ 有 $\gamma_{tj}=-(\gamma_{1j}+\gamma_{2j}+\cdots+\gamma_{t-1j})$ 以及对于 $i=1,2,\cdots,t$ 有 $\gamma_{ib}=-(\gamma_{i1}+\gamma_{i2}+\cdots+\gamma_{ib-1})$.
3. 拟合具有处理变量所有哑变量及其交互作用的模型，并获得残差平方和。这相当于全效应模型 (full-effects model) 的残差平方和。
4. 接下来，拟合一个模型，其中包含所有哑变量，除了与正在检验的主效应或交互作用对应的哑变量。该缩减模型的残差平方和与步骤 3 中模型的残差平方和之间的差值就是与该效应对应的平方和。

由此产生的分析称为 III 型分析。它也被称为 Yates's weighted squares of means technique. 当观测到所有处理组合时，III 型分析检验的假设与均衡数据集检验的假设相同。表 \@ref(tab:table10-14) 给出了这些 III 型假设的平均模型参数，表 \@ref(tab:table10-15) 给出了效应模型参数。这些假设通常是实验者想要的，表 \@ref(tab:table10-16) 给出了表 \@ref(tab:table9-1) 中数据的 III 型平方和。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table10-14)基于均值模型的效应模型的 III 型假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.14.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table10-15)根据效应模型参数对效应模型的 III 型假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.15.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table10-16)III 型方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.16.png">
</td>
  </tr>
</tbody>
</table>

对于表 \@ref(tab:table9-3) 中的数据，通过 Yates's method 重新参数化的模型的矩阵形式为

$$\begin{bmatrix}19\\20\\21\\24\\22\\25\\25\\25\\27\\21\\24\\21\\32\\33\end{bmatrix}=\begin{bmatrix}1&-1&1&0&-1&0\\1&-1&1&0&-1&0\\1&-1&0&1&0&-1\\1&-1&0&1&0&-1\\1&-1&0&1&0&-1\\1&-1&-1&-1&1&1\\1&-1&-1&-1&1&1\\1&-1&-1&-1&1&1\\1&1&1&0&1&0\\1&1&1&0&1&0\\1&1&1&0&1&0\\1&1&0&1&0&1\\1&1&0&1&0&1\\1&1&-1&-1&-1&-1\\1&1&-1&-1&-1&-1\\1&1&-1&-1&-1&-1\end{bmatrix}\begin{bmatrix}\mu\\\tau_1\\\beta_1\\\beta_2\\\gamma_{11}\\\gamma_{12}\end{bmatrix}+\begin{bmatrix}{\varepsilon}_{111}\\{\varepsilon}_{112}\\{\varepsilon}_{113}\\{\varepsilon}_{121}\\{\varepsilon}_{122}\\{\varepsilon}_{131}\\{\varepsilon}_{132}\\{\varepsilon}_{133}\\{\varepsilon}_{211}\\{\varepsilon}_{212}\\{\varepsilon}_{213}\\{\varepsilon}_{221}\\{\varepsilon}_{222}\\{\varepsilon}_{231}\\{\varepsilon}_{232}\\{\varepsilon}_{233}\end{bmatrix}$$

请注意，上述设计矩阵中的第五列（对应于 $\gamma_{11}$）是第二列和第三列（对应于 $\tau_{1}$ 和 $\beta_1$）的乘积，第六列（对应于 $\gamma_{12}$）是第二列和第四列的乘积（对应于 $\tau_1$ 和 $\beta_2$）。

SAS-GLM 引入了第四种生成与主效应及其交互作用对应的平方和的方法。当观测到所有处理组合时，IV 型分析检验的假设与 III 型分析检验的假设相同；然而，当未观测到某些处理组合时，III 型和 IV 型分析不一致。我们将在第 \@ref(chap14) 章讨论 IV 型假设的构建。

作为本节的总结，在分析没有缺失处理组合的双向处理结构模型中的数据时，我们提出以下建议：

1. 如果实验者想比较两种处理的效应，她应该关注 III 型分析所检验的假设。这些假设等价于在均衡或相等子类数量情况下检验的假设。
2. 如果实验者有兴趣建立一个模型来预测特定处理组合的效应，那么她可以使用 I 型和/或 II 型分析。
3. 在抽样调查实验中，每个处理组合的观测数量通常与这些组合在总体中实际发生的频率成正比。在这种情况下，实验者可能对基于 $R(\tau|\mu,\beta)$ 和 $R(\beta|\mu,\tau)$ 的假设最感兴趣，因为这些平方和检验了有关行均值和列均值的加权平均值的假设，权重与观测到的样本量成正比。如果使用 SAS-GLM，这可能需要两次 I 型分析，一次先在模型中纳入 T，另一次则先在模型中纳入 B.

## 在 SAS-GLM 中使用 I–IV 型可估函数 {#sec10-5}

在本节中，我们展示了如何使用 SAS-GLM 提供的信息来确定由不同类型的平方和检验的假设（不使用 SAS-GLM 的读者可以跳过本节）。在第 \@ref(sec10-3) 节中，我们讨论了通过使用 MODEL 语句上的 `E` 选项获得的可估函数的一般形式。如果选择了 `E1` 选项，SAS-GLM 将为每个效应打印 I 型可估函数的一般形式。当使用效应模型时，该选项对表 \@ref(tab:table9-1) 中数据的结果如表 \@ref(tab:table10-17) 所示。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table10-17)表 \@ref(tab:table9-1) 中数据的 SAS-GLM 的 I 型可估函数</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.17.png">
</td>
  </tr>
</tbody>
</table>

根据表 \@ref(tab:table10-17)，我们发现 $\boldsymbol \beta=[\mu,\tau_1,\tau_2,\beta_1,\beta_2,\beta_3,\gamma_{11},\gamma_{12},\gamma_{13},\gamma_{21},\gamma_{22},\gamma_{23}]$ 的线性组合 $\boldsymbol \ell^\prime \boldsymbol \beta$ 是 I 型可估函数当且仅当存在常数 $L2$ 使得

$$\begin{aligned}
\boldsymbol \ell^\prime \boldsymbol \beta  =&(L2)\tau_1-(L2)\tau_2+(0.125\cdot L2)\beta_1-(0.125\cdot L2)\beta_2 \\
&+(0.375\cdot L2)\gamma_{11}+(0.25\cdot L2)\gamma_{12}+(0.375\cdot L2)\gamma_{13} \\
&-(0.25\cdot L2)\gamma_{21}-(0.375\cdot L2)\gamma_{22}-(0.375\cdot L2)\gamma_{23} 
\end{aligned}$$

T 的 I 型可估函数的基可以通过选择 $L2$ 的特定值来获得，例如 $L2 = 1$ 或 $L2 = 8$. 我们可以自由地仅选择可估函数一般形式的 $L$ 中的一个，其对应于与 T 的 I 型平方和相关的 1 个自由度；剩下的所有 $L$ 都由我们对 $L2$ 的选择决定。对于我们的示例，$L1 = 0,L4 = 0.125 · L,L5 = -0.125 · L2,L7 = 0.375 · L2$ 和 $L8 = 0.25 · L2$.

取 $L2=1$，I 型可估函数的基为

$$\{\tau_1-\tau_2+(1/8)(\beta_1-\beta_2)+(1/8)(3\gamma_{11}+2\gamma_{12}+3\gamma_{13}-2\gamma_{21}-3\gamma_{22}-3\gamma_{23})\}$$

该参数函数是通过 I 型分析与零进行比较的函数。参见表 \@ref(tab:table10-11). 另一个基可以通过取 $L2 = 8$ 来构造，由下式给出

$$\{8\tau_1-8\tau_2+\beta_1-\beta_2+3\gamma_{11}+2\gamma_{12}+3\gamma_{13}-2\gamma_{21}-3\gamma_{22}-3\gamma_{23}\}$$

由于 $\mu_{ij} = \mu + \tau_i + \beta_j + \gamma_{ij}$，在均值模型参数框架下所检验的假设可以通过将效应模型表示中的 $\gamma_{ij}$ 上的系数分配给均值模型表示中的 $\mu_{ij}$ 来确定。因此，以均值模型参数表达的 T 的 I 型假设为

$$\begin{array}{l}(3\mu_{11}+2\mu_{12}+3\mu_{13}-2\mu_{21}-3\mu_{22}-3\mu_{23})/8=0\text{ or equivalently that}\\3\mu_{11}+2\mu_{12}+3\mu_{13}-2\mu_{21}-3\mu_{22}-3\mu_{23}=0\end{array}$$

这是表 \@ref(tab:table10-9) 中给出的关于 T 的假设。

如果选择 MODEL 语句上的 `E2` 选项，SAS-GLM 将打印每种效应的 II 型可估函数的一般形式。表 \@ref(tab:table9-1) 中的数据结果如表 \@ref(tab:table10-18) 所示。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table10-18)表 \@ref(tab:table9-1) 中数据的 SAS-GLM 的 II 型可估函数</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.18.png">
</td>
  </tr>
</tbody>
</table>

根据表 \@ref(tab:table10-18) 可以看出，$\boldsymbol \ell^\prime \boldsymbol \beta$ 是 II 型可估函数当且仅当存在常数 $L4,L5$ 使得

$$
\begin{aligned}\boldsymbol{\ell}^{\prime}\boldsymbol{\beta} =&(L4)\beta_1+(L5)\beta_2+(-L4-L5)\cdot\beta_3+(0.5692\cdot L4+0.0308\cdot L5)\gamma_{11}\\
&+(-0.0308\cdot L4+0.4308\cdot L5)\gamma_{12}+(-0.5385\cdot L4-0.4615\cdot L5)\gamma_{13} \\
&+(0.4308\cdot\text{L4- }0.0308\cdot L5)\gamma_{21}+(0.0308\cdot L4+0.5692\cdot L5)\gamma_{22} \\
&+(-0.4615\cdot L4-0.5385\cdot L5)\gamma_{23} \end{aligned}
$$

在这种情况下，我们可以选择两个 L 的值，即 $L4$ 和 $L5$. 因此，B 的 II 型平方和有 2 个自由度。运气好的话，或者使用表 \@ref(tab:table10-9) 中给出的一般形式，我们可以确定上述表达式中给出的十进制数有一个最小公分母 65. 选择 $L4 = 1$ 和 $L5 = 0$，然后选择 $L4 = 0$ 和 $L5 = 1$， 可以为 B 的 II 类可估函数提供一个基，在效应模型参数框架下为

$$\begin{array}{l}\{\beta_1-\beta_3+(1/65)(37\gamma_{11}-2\gamma_{12}-35\gamma_{13}+28\gamma_{21}+2\gamma_{22}-30\gamma_{23})\mathrm{~and}\\\beta_2-\beta_3+(1/65)(2\gamma_{11}+28\gamma_{12}-30\gamma_{13}-2\gamma_{21}+37\gamma_{22}-35\gamma_{23})\}\end{array}$$

对于第一个函数，令 $L4 = 65$ 且 $L5 = 0$；对于第二个函数，令 $L4 = 0$ 且 $L5 = 65$，可得到在均值模型参数框架下的基

$$\begin{gathered}
\{37{\mu}_{11}-2{\mu}_{12}-35{\mu}_{13}+28{\mu}_{21}+2{\mu}_{22}-30{\mu}_{23}\mathrm{~and} \\
2\mu_{11}+28\mu_{12}-30\mu_{13}-2\mu_{21}+37\mu_{22}-35\mu_{23}\} 
\end{gathered}$$

从表 \@ref(tab:table10-17) 和/或表 \@ref(tab:table10-18) 中 T × B 的 I 型和/或 II 型可估函数的一般形式，可以看出有两个 L 值可供选择。 通过令 $L7 = 1$ 和 $L8 = 0$，然后令 $L7 = 0$ 和 $L8 = 1$，可以得到 T × B 的 I 型可估函数的一个基，在效应模型参数框架下将产生

$$\{\gamma_{11}-\gamma_{13}-\gamma_{21}+\gamma_{23}\mathrm{~and~}\gamma_{12}-\gamma_{13}-\gamma_{22}+\gamma_{23}\}$$

在均值模型参数框架下将产生

$$\{\mu_{11}-\mu_{13}-\mu_{21}+\mu_{23}\mathrm{~and~}\mu_{12}-\mu_{13}-\gamma_{22}+\mu_{23}\}$$

读者可能希望验证，对于所有 $i\ne i^\prime$ 和 $j\ne j^\prime$，每个 2 × 2 表的差 $\mu_{ij}-\mu_{ij'}-\mu_{i'j}+\mu_{i'j'}=0$ 都可以作为这个基中这两个函数的线性组合。

根据表 \@ref(tab:table10-18)，T 的 II 型可估函数的一般形式由下式给出

$$\begin{array}{c}(L2)\tau_1-(L2)\tau_2+(0.3077\cdot L2)\gamma_{11}+(0.3077\cdot L2)\gamma_{12}+(0.3846\cdot L2)\gamma_{13}\\+(-0.3077\cdot L2)\gamma_{21}+(-0.3077\cdot L2)\gamma_{22}+(-0.3846\cdot L2)\gamma_{23}\end{array}$$

这些十进制分数的最小公分母是 13；因此，取 $L2 = 1$，T 的 II 类可估函数的基，在效应模型参数框架下为

$$\{\tau_1-\tau_2+(1/13)(4\gamma_{11}+4\gamma_{12}+5\gamma_{13}-4\gamma_{21}-4\gamma_{22}-5\gamma_{23})\}$$

在均值模型参数框架下为（取 $L2=13$）

$$\{4\mu_{11}+4\mu_{12}+5\mu_{13}-4\mu_{21}-4\mu_{22}-5\mu_{23}\}$$

最后，我们考虑 III 和 IV 型可估函数。对于表 \@ref(tab:table9-1) 中的数据，它们是相同的；如表 \@ref(tab:table10-19) 所示。T 的 III 型可估函数的基，在效应模型参数框架下为

$$\{\tau_1-\tau_2+\bar{\gamma}_{1\cdot}-\bar{\gamma}_{2\cdot}\}$$

在均值模型参数框架下为（在两种情况下都取 $L2=1$）

$$\{\bar{\mu}_{\cdot 1}-\bar{\mu}_{\cdot 2}\}$$

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table10-19)表 \@ref(tab:table9-1) 中数据的 SAS-GLM 的 III 型可估函数</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.19.png">
</td>
  </tr>
</tbody>
</table>

以类似的方式，我们可得到 B 的 III 型可估函数的基，在效应模型参数框架下为

$$\{\beta_1-\beta_3+\bar{\gamma}_{\cdot1}-\bar{\gamma}_{\cdot 3}\mathrm{~and~}\beta_2-\beta_3+\bar{\gamma}_{\cdot 2}-\bar{\gamma}_{\cdot 3}\}$$

在均值模型参数框架下为

$$\{\bar{\mu}_{\cdot1}-\bar{\mu}_{\cdot3}\mathrm{~and~}\bar{\mu}_{\cdot2}-\bar{\mu}_{\cdot3}\}$$

请注意，$H_0\colon\bar{{\mu}}_{\cdot1}-\bar{{\mu}}_{\cdot3}=0\mathrm{~and~}\bar{{\mu}}_{\cdot2}-\bar{{\mu}}_{\cdot3}=0$ 为真当且仅当 $H_0\colon\bar{\mu}_{1\cdot}=\bar{\mu}_{2\cdot}=\bar{\mu}_{3\cdot}$ 为真。读者可以验证表 \@ref(tab:table9-1) 中数据的 IV 型可估函数与 III 型可估函数是否相同。

## 总体边际均值与最小二乘均值 {#sec10-6}

双向效应模型的总体边际均值定义为，对于处理 T

$$\bar{\mu}_{i\cdot}=\mu+\tau_i+\bar{\beta}_{\cdot}+\bar{\gamma}_{i\cdot},\quad i=1,2,\ldots,t$$

对于处理 B

$$\bar{\mu}_{\cdot j}=\mu+\bar{\tau}_{\cdot}+\bar{\beta}_j+\overline{\gamma}_{\cdot j},\quad j=1,2,\ldots,b$$

这些边际均值的最佳估计分别为

$$\hat{\bar{\mu}}_{i\cdot}=\hat{\mu}+\hat{\tau}_i+\hat{\bar{\beta}}_{\cdot}+\hat{\bar{\gamma}}_i.,i=1,2,\ldots,t\quad\mathrm{and}\quad\hat{\bar{\mu}}_{\cdot j}=\hat{\mu}+\hat{\bar{\tau}}_{\cdot}+\hat{\bar{\beta}}_j+\hat{\bar{\gamma}}_{\cdot j},~j=1,2,\ldots,b$$

这些边际均值的估计通常被称为最小二乘均值。它们各自的标准误由方程 \@ref(eq:9-5) 和 \@ref(eq:9-7) 给出。为了对总体边际均值的线性组合进行推断，可以使用 \@ref(eq:9-5) 和 \@ref(eq:9-6)。如果在 SAS-GLM 程序中使用


```sas
LSMEANS T B T*B/PDIFF;
```

选项，则得到边际均值和双向单元格均值的最佳无偏估计，以及这些均值某些子集之间的成对比较。这些数据如表 \@ref(tab:table10-20) 所示。请注意，T 和 B 的最小二乘均值与表 \@ref(tab:table9-2) 中给出的 T 和 B 均值相同，双向均值与表 \@ref(tab:table9-2) 中的单元格均值相同。表 \@ref(tab:table10-20) 还提供了最小二乘均值的标准误估计和用于最小二乘均值之间成对比较的 $p$ 值。例如，$\hat{{\mu}}_{1\cdot}$ 的标准误估计为 0.51，这与第 \@ref(sec9-5) 节中计算的标准误估计一致。此外，表 \@ref(tab:table10-20) 中给出了将 $\hat{{\mu}}_{1\cdot}$ 与 $\hat{{\mu}}_{2\cdot}$ 进行比较的 $p=0.0002$，这与第 \@ref(sec9-5) 节中计算的 $t$ 统计量的 $p$ 值一致。可以对其他估计、标准误和检验统计量进行类似的比较。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table10-20)边际均值和单元格均值的最佳估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.20-1.png">
</td>
  </tr>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.20-2.png">
</td>
  </tr>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2010.20-3.png">
</td>
  </tr>
</tbody>
</table>

## 计算机分析 {#sec10-7}

**几乎所有的统计计算包都是为了处理效应模型而不是均值模型而开发的**。由于三种主要类型的假设可以检验具有不相等子类数量的均衡双向处理结构，因此必须小心确定所使用的统计软件包检验的是哪种类型。尽管这些计算软件包是为了处理效应模型而开发的，但均值模型也可以很容易地实现。使用均值模型允许用户指定处理均值之间的有意义对比。本讨论中使用的名称，即 I, II, III 和 IV 型，与 SAS 使用的名称相对应。SPSS 使用相同的表示法，并将生成 I, II 和 III 型分析。SPSS 默认使用 III 型分析。默认情况下，SAS-GLM 同时生成 I 型和 III 型分析。鼓励使用其他计算软件包的读者分析表 \@ref(tab:table9-1) 中的数据，然后将其分析结果与本章和第 \@ref(chap9) 章中的分析进行比较。这种比较应该能够揭示出他们的计算软件包生成的是哪种类型的分析。SAS-GLM 还允许如下语句来计算均值


```sas
MEANS T B T*B/<options>;
```

然而，对于不均衡的数据，Means 语句并没有给出总体边际均值的无偏估计，而是计算行均值和列均值的加权平均值的估计值。也就是说，上述 Means 选项给出了加权均值的最佳无偏估计，其定义为

$$\tilde{\mu}_{i.}=\frac1{n_{i\cdot}}\sum_{j=1}^{b}n_{ij}\mu_{ij},\mathrm{~}i=1,2,\ldots,t\quad\mathrm{and~}\quad\tilde{\mu}_{\cdot j}=\frac1{n_{\cdot j}}\sum_{i=1}^{t}n_{ij}\mu_{ij},\mathrm{~}j=1,2,\ldots,b$$

实验者可能对这些加权平均值感兴趣的唯一情况可能是抽样调查实验（见第 \@ref(sec10-4) 节末尾的第三条建议）。

对于表 \@ref(tab:table9-1) 中的数据，这些加权平均值的最佳估计值，T 加权均值为 $\bar y_{1\cdot\cdot}=22.75$ 和 $\bar y_{2\cdot\cdot}=27.125$，B 加权均值为 $\bar y_{\cdot1\cdot}=22.2,\bar y_{\cdot2\cdot}=23.8$ 和 $\bar y_{\cdot3\cdot}=28.0$.

## 结束语 {#sec10-8}

在本章中，假设所有处理组合都被观测到，使用效应模型对子类数量不等的双向处理结构进行了分析。尽管均值模型可以回答所有可能提出的问题，但在这里使用了效应模型，因为现有的大部分统计计算软件都使用效应模型。讨论了 I–III 型分析以及每种分析类型的使用条件。在几乎所有情况下，III 型分析将是首选的分析。III 型分析与第九章中给出的分析相同。将总体边际均值与加权边际均值进行对比，并考虑每种均值的适用性。

值得注意的是，如果关于表 \@ref(tab:table9-1) 中的数据建立双向加性模型，

$$y_{ijk}=\mu+\tau_i+\beta_j+\varepsilon_{ijk},\quad i=1,2,\ldots,t;\quad j=1,2,\ldots,b;\quad k=1,2,\ldots,n_{ij}$$

在包括交互项的双向模型分析中，T 和 B 的 II 型和 III 型平方和将彼此相同，并且它们都等于包含交互项的双向模型分析中 T 和 B 的 II 型平方和。
