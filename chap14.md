
# 使用效应模型分析缺失处理组合的双向处理结构 {#chap14}

> > "An approximate answer to the right problem is worth a good deal more than an exact answer to an approximate problem." - John Tukey

在第 \@ref(chap13) 章中，我们讨论了当某些处理组合未观测到时，如何使用均值模型来分析双向处理结构。在本章中，我们考虑使用效应模型来分析相同类型的情况。使用效应模型并不能回答任何使用均值模型无法回答的问题，反之亦然。虽然均值模型非常简单且易于理解，但效应模型似乎比实际情况复杂得多。我们更喜欢使用均值模型，并且在这里讨论效应模型只是因为许多统计包似乎推荐并鼓励使用效应模型。本章考虑的效应模型是

$$y_{ijk}=\mu+\tau_i+\beta_j+\gamma_{ij}+\varepsilon_{ijk},\quad i=1,2,\ldots,t;~j=1,2,\ldots,b;~k=0,1,2,\ldots,n_{ij}$$

其中 $\varepsilon_{ijk}\sim i.i.d.N(0,\sigma^2)$ 且 $n_{ij}=0$ 表示第 (i,j) 个处理组合未观测到。

## I 型和 II 型假设

对于缺少处理组合的双向处理结构的 I 型和 II 型分析可以定义为观测到所有组合的处理结构。也就是说，可以拟合连续的模型，并在将不同效应添加到模型中时确定残差平方和的减少量。 为了说明这一点，请考虑对表 \@ref(tab:table13-1) 中的数据进行 I 型分析。我们在这里只考虑结果。对实际模型拟合结果感兴趣的读者可以参考本书的上一版本。表 \@ref(tab:table14-1) 和 \@ref(tab:table14-2) 分别给出了这些数据的 I 型和 II 型方差分析表。

<table>
<caption>(\#tab:table14-1)I 型方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2014.1.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table14-2)II 型方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2014.2.png">
</td>
  </tr>
</tbody>
</table>

在第 \@ref(chap10) 章中已经指出，当数据不均衡时，即使所有的处理组合都至少被观测一次，I 型和 II 型假设也可能没有多大意义。那么很明显，在有缺失处理组合的情况下，这两种假设也不会突然有意义。事实上，为了证明这一点是正确的，表 \@ref(tab:table14-3) 和表 \@ref(tab:table14-4) 给出了 I 型和 II 型分析检验的假设（均值模型），表 \@ref(tab:table14-5) 和表 \@ref(tab:table14-6) 给出了效应模型。这些表中的条目可以通过 SAS^®^-GLM 对表 \@ref(tab:table13-1) 中的数据进行分析来确定，也可以使用表 \@ref(tab:table10-9) 和公式 \@ref(tab:table10-3) 中给出的通用公式来确定，因为这些公式也适用于缺失单元格的问题。

<table>
<caption>(\#tab:table14-3)根据均值模型对表 \@ref(tab:table13-1) 中数据进行 I 型分析的假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2014.3.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table14-4)根据均值模型对表 \@ref(tab:table13-1) 中数据进行 II 型分析的假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2014.4.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table14-5)根据效应模型对表 \@ref(tab:table13-1) 中数据进行 I 型分析的假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2014.5.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table14-6)根据效应模型对表 \@ref(tab:table13-1) 中数据进行 II 型分析的假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2014.6.png">
</td>
  </tr>
</tbody>
</table>

接下来，我们讨论 I 型和 II 型主效应假设的可能解释。I 型和 II 型假设通常没有多大意义，除非目标是建立一个简单的模型来进行预测，而不是检验有关不同处理组合效应的假设。对于模型构建，解释与第 \@ref(chap10) 章中的解释完全相同，我们在第 \@ref(chap10) 章中讨论了没有缺失处理组合的案例。如第 \@ref(chap10) 章所示，如果每个单元格中的观测数量与总体中存在的每种处理组合的实际数量成正比，那么实验者可能会对 $R(\tau|\mu)$ 和 $R(\beta|\mu)$ 感兴趣。两者都可以通过进行两个 I 型分析来获得，一个是 model 语句中首先出现 T，另一个是 model 语句中首先出现 B.

## III 型假设

当观测到所有处理组合时，III 型假设与子类数相等时检验的假设相同。当某些处理组合缺失时，无法检验此类假设，因为它们涉及没有相关信息的参数。对于表 \@ref(tab:table13-1) 中的数据，我们无法估计 $\hat \mu_{1\cdot}$ 和 $\hat \mu_{3\cdot}$，因为我们无法估计 $\mu_{12}$ 和 $\mu_{33}$. 同样，我们也无法估计 $\hat \mu_{\cdot 2}$ 和 $\hat \mu_{\cdot 3}$. 因此无法检验 $\bar{\mu}_{1\cdot}=\bar{\mu}_{2\cdot}=\bar{\mu}_{3\cdot}$ 或 $\bar{\mu}_{\cdot 1}=\bar{\mu}_{\cdot 2}=\bar{\mu}_{\cdot 3}$.

主效应的 I 型和 II 型假设都取决于每个单元格中的观测次数。只要在一个单元格中至少有一个观测结果，那么该单元格均值是可估的。因此，可估计的参数函数仅取决于观测到哪种处理组合，而不取决于观测了多少次。

III 型假设的发展使得它们不依赖于单元格的大小，而仅依赖于观测到的单元格。这与双向实验的 III 型假设的定义一致，其中观测到所有处理组合。即，假设 $\bar{\mu}_{1\cdot}=\bar{\mu}_{2\cdot}=\bar{\mu}_{3\cdot}$ 和 $\bar{\mu}_{\cdot 1}=\bar{\mu}_{\cdot 2}=\bar{\mu}_{\cdot 3}$ 不依赖于单元格大小。我们不打算讨论缺失数据情况下 III 型假设的构建。尽管目标可能看起来合理，但我们认为当存在缺失单元格时， III 型假设是最糟糕的假说，因为似乎没有合理的方法来解释它们。为了便于说明，表 \@ref(tab:table14-7) 中给出了表 \@ref(tab:table13-1) 中数据的 III 型分析。表 \@ref(tab:table14-8) 和 \@ref(tab:table14-9) 给出了通过 III 型分析测试的假设。对表 \@ref(tab:table14-8) 和 \@ref(tab:table14-9) 的研究表明，除了 T×B 交互作用之外，III 型假设没有意义。表 \@ref(tab:table14-7) - \@ref(tab:table14-9) 中的结果取自 SAS-GLM 对表 \@ref(tab:table13-1) 中数据的分析。

<table>
<caption>(\#tab:table14-7)III 型方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2014.7.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table14-8)根据均值模型对表 \@ref(tab:table13-1) 中数据进行 III 型分析的假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2014.8.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table14-9)根据效应模型对表 \@ref(tab:table13-1) 中数据进行 III 型分析的假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2014.9.png">
</td>
  </tr>
</tbody>
</table>

## IV 型假设 {#sec14-3}

在前两节中我们试图表明当存在缺失处理组合时，通过 I 型、II 型或 III 型分析检验的所谓主效应假设中没有一个是完全令人满意的，因为如果两个因素之间存在交互作用的可能性，它们很少有任何合理的解释。这种假设极难解释，因为发生在同一行或同一列中的单元格均值系数很少相同。构建 IV 型假设，使得单元格均值系数是平衡的；因此，由此产生的假设是可解释的。

为了说明这一点，让我们看看表 \@ref(tab:table13-1) 数据集中所有可能的 IV 型假设。表 \@ref(tab:table14-10) 给出了非缺失单元格的单元格均值参数。**基本上，对于双向处理结构，如果假设比较的是一个处理水平在另一个处理的一个或多个共同水平上的平均值，那么该假设就被定义为 IV 型假设**。因此，这些假设应该是边际均值假设，只不过当处理组合缺失时，人们不能对所有其他处理水平进行平均，而只能对其他因素的部分水平进行平均。

<table>
<caption>(\#tab:table14-10)表 \@ref(tab:table13-1) 中数据的单元格均值参数</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2014.10.png">
</td>
  </tr>
</tbody>
</table>

在对 $B_1$ 和 $B_3$ 求平均后比较 $T_1$ 和 $T_2$ 的 IV 型假设是：$H_0{:}(\mu_{11}+\mu_{13})/2=(\mu_{21}+\mu_{23})/2$. 另一种比较 $T_1$ 和 $T_2$ 的 IV 型假设是 $H_0{:\mu_{11}}=\mu_{21}$. 后一种假设仅在 $B_1$ 上求平均后比较比较 $T_1$ 和 $T_2$. 关于均值模型参数的 T 的所有可能 IV 型假设见表 \@ref(tab:table14-11). $B$ 的 IV 型假设是类似的。表 \@ref(tab:table14-12) 给出了 B 的所有可能的 IV 型假设。

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table14-11)表 \@ref(tab:table13-1) 中数据关于 T 的所有可能的 IV 型假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2014.11.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">*: </span> <sup></sup> 通过 SAS-GLM IV 型分析自动检验的假设</td></tr></tfoot>
</table>

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table14-12)表 \@ref(tab:table13-1) 中数据关于 B 的所有可能的 IV 型假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2014.12.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">*: </span> <sup></sup> 通过 SAS-GLM IV 型分析自动检验的假设</td></tr></tfoot>
</table>

SAS-GLM 会自动生成通常可以解释的 IV 型假设，但在解释这些假设之前，首先需要检查 IV 型可估函数，以明确 SAS-GLM 生成了哪些假设。也就是说，与没有缺失处理组合的情况不同，不存在适合所有数据集的唯一解释。事实上，在分析之前对处理重新标记 (relabeling) 可能会导致生成不同的 IV 型假设，因此在 IV 型分析中会有不同的平方和以及 $F$ 值。因此，获得的 IV 型分析并不是数据的唯一特征，而是取决于处理是如何标记的 (labeled). 显然，这并不是非常理想的，但也是不可避免的。SAS-GLM 通过在打印的自由度上放置星号，并指出 "Other Type IV Testable Hypotheses exist which may yield different SS." 来表示这种情况已经发生。对于表 \@ref(tab:table13-1) 中的数据，SAS-GLM 分析自动检验的 T 的 IV 型假设等价于同时检验

$$\mu_{11}=\mu_{31}\quad\mathrm{and}\quad\frac{\mu_{21}+\mu_{22}}2=\frac{\mu_{31}+\mu_{32}}2$$

因此，T 的 IV 型假设同时比较了 $B_1$ 上 $T_1$ 和 $T_3$ 的效应以及 $B_1$ 和 $B_2$ 上 $T_2$ 和 $T_3$ 的平均效应。这两个函数在表 \@ref(tab:table14-11) 中用星号标识。我们注意到，$B_3$ 根本不涉及此特定集合。由 SAS-GLM 得到的 IV 型方差分析表如表 \@ref(tab:table14-13) 所示。

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table14-13)IV 型方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2014.13.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">*: </span> <sup></sup> 存在其他可检验的 IV 型假设，这些假设可能产生不同的 $SS$</td></tr></tfoot>
</table>

表 \@ref(tab:table14-11) 中给出了未经 SAS-GLM 检验的 IV 型假设，但对于实验者来说可能与自动检验的假设一样有趣。为了检验这种有趣的 IV 型假设，人们可以（并且应该）使用 Estimate 或 Contrast 语句。例如，要在使用效应模型时检验表 \@ref(tab:table14-11) 中的所有 IV 型假设，我们将在 SAS-GLM 分析中使用以下 Estimate 语句


```sas
ESTIMATE 'T1 VS T2 AVE OVER Bl AND B3' T 1 —1 0 T*B .5 .5 —.5 0 —.5 0 0;
ESTIMATE 'T2 VS T3 AVE OVER Bl AND B2' T 0 1 —1 T*B 0 0 .5 .5 0 —.5 —.5;
ESTIMATE 'T1 VS T2 AT Bl' T 1 —1 0 T*B 1 0 —1 0 0 0 0;
ESTIMATE 'TI VS T3 AT Bl' T 1 0 —1 T*B 1 0 0 0 0 —1 0;
ESTIMATE 'T2 VS T3 AT B1' T 0 1 —1 T*B 0 0 1 0 0 —1 0;
ESTIMATE 'T2 VS T3 AT B2' T 0 1 —1 T*B 0 0 0 1 0 0 —1;
ESTIMATE 'T1 VS T2 AT B3' T 1 —1 0 T*B 0 1 0 0 —1 0 0;
```

上述 Estimate 语句的结果如表 \@ref(tab:table14-14) 所示。

<table>
<caption>(\#tab:table14-14)`Estimate` 选项的结果</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2014.14.png">
</td>
  </tr>
</tbody>
</table>

如果读者没有可用的 SAS-GLM，则可以使用与第 \@ref(sec13-2) 节中给出的类似的均值模型和程序来检验表 \@ref(tab:table14-11) 和 \@ref(tab:table14-12) 中的所有假设。

## 总体边际均值和最小二乘均值 {#sec14-4}

总体边际均值和最小二乘均值的定义方式与第 \@ref(sec9-5) 节和第 \@ref(sec10-6) 节中的定义方式相同。然而，如果没有观测到特定处理与其他处理因素的所有可能性，则相应的总体边际平均值是不可估的。在这种情况下，双向单元格均值表（例如 $\bar \mu_{ij}$）可用于将每个观测到的处理组合与所有其他观测到的处理组合进行比较。如果数据集非常稀疏，则很少有总体边际均值可估计。

对于表 \@ref(tab:table13-1) 中的数据，$\bar\mu_{2\cdot}$ 和$\bar\mu_{1\cdot}$ 是唯一可估的总体边际均值。它们的最佳估计值分别为 $\hat{\bar\mu}_{2\cdot}=8.333$ 和 $\hat{\bar\mu}_{1\cdot}=4.00$. 一般来说，$\sum_{i,j}c_{ij}\mu_{ij}$ 的最佳估计为 $\sum_{i,j}c_{ij}\hat{\mu}_{ij}$，且标准误估计为 ${\hat{\sigma}}\sqrt{\Sigma_{i,j}(c_{ij}^2/n_{ij})}$，其中关于所有非空单元格求和。$\sum_{i,j}c_{ij}\mu_{ij}$ 的一个 $(1-\alpha)100\%$ 置信区间为

$$\sum_{i,j}c_{ij}\hat{\mu}_{ij}\pm t_{\alpha/2,v}\hat{\sigma}\sqrt{\frac{c_{ij}^2}{n_{ij}}}$$

用于检验 $\sum_{i,j}c_{ij}\mu_{ij}=0$ 自由度为 $v$ 的 $t$ 统计量为

$$t=\frac{\sum_{i,j}c_{ij}\hat{{\mu}}_{ij}}{\hat{{\sigma}}\sqrt{\sum_{i,j}\frac{c_{ij}^2}{n_{ij}}}}$$

在这两种情况下，$n$ 是对应于误差均方 $\hat\sigma^2$ 的自由度。

对于希望进行多重比较的研究人员，我们建议在比较所有处理组合的 $F$ 值显著时，使用上述 $t$ 检验给出的观察到的 $p$ 值。如果该 $F$ 值不显著，那么应该在所有感兴趣的比较中使用 Bonferroni's method。也就是说，如果获得的 $p$值小于 $\alpha/p$，其中 $p$ 是计划比较的总次数，则声明线性组合与零显著不同。对于数据窥探和计划外比较，应使用 Scheffé's procedure.

## 计算机分析 {#sec14-5}

读者应该使用他或她自己的统计软件包来分析本章以及第 \@ref(chap15) 章和第 \@ref(chap17) 章中给出的例子。将由此获得的分析结果与本书中给出的结果进行比较，将使读者对她习惯使用的软件包所检验的各种假设有宝贵的见解。**据我们所知，没有任何一个包能够充分或完整地处理缺失处理组合的数据分析**。有几个软件包能够很好地处理不均衡数据，前提是没有缺失处理组合。

## 结束语 {#sec14-6}

总之，对于存在缺失处理组合的数据，要做出可接受的分析需要深思熟虑。实验人员或统计学家不能仅仅在计算机程序上运行数据，然后从该程序中选取数字来撰写论文。不幸的是，这种做法已经被不知情的实验人员和数据分析师采用，并且还在被采用。**我们希望学习过本章的人再也不会这样做**。对于那些愿意付出必要努力来正确分析数据的人，建议使用第 \@ref(chap13) 章中讨论的均值模型。

第 \@ref(chap15) 章将讨论一个更现实的例子。

## 练习 {#sec14-7}
