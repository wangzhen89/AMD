# 案例研究：具有许多缺失处理组合的三向处理结构 {#chap17}

> "Statistics are no substitute for judgment." - Henry Clay

在本章中，我们详细分析了当缺少许多处理组合时的三向处理结构。

## 营养评分示例 {#sec17-1}

一位家庭经济学家进行了一项抽样调查实验，研究低收入经济水平的母亲对营养的了解程度，并判断旨在增加她们营养知识的培训计划的效果。在训练项目前后对这些母亲进行了测试，并测量了她们测试成绩的变化。这些变化如表 \@ref(tab:table17-1) 所示。受试的母亲根据三个因素进行分类：age, race 以及是否收到食品券 (food stamps).

<table>
<caption>(\#tab:table17-1)培训前后营养测试成绩的变化</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2017.1.png">
</td>
  </tr>
</tbody>
</table>

## SAS-GLM 分析 {#sec17-2}

表 \@ref(tab:table17-2) 给出了使用 SAS 命令从 SAS^®^-GLM 获得的表 \@ref(tab:table17-1) 中数据的 IV 型方差分析表。


```sas
PROC GLM;  
CLASSES GROUP AGE RACE;  
MODEL GAIN=GROUP|AGE|RACE/SOLUTION E4;  
LSMEANS GROUP|AGE|RACE/PDIFF STDERR;  
RUN;
```

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table17-2)IV 型方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2017.2.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">*: </span> <sup></sup> 存在其他可检验的 IV 型假设，这些假设可能会产生不同的 SS</td></tr></tfoot>
</table>

表 \@ref(tab:table17-2) 显示表 \@ref(tab:table17-1) 中数据的 $\sigma^2$ 估计为 $\hat\sigma^2 = 2627.4724/92 = 28.56$，自由度为 92. 表 \@ref(tab:table17-2) 还显示三因素交互作用假设的自由度为零，这表明这些数据中没有可用于估计三因素交互作用的对比。这并不意味着 group, age 和 race 之间不存在三因素交互作用，只是三因素交互作用效应中没有可检验的假设。

IV 型 $F$ 值似乎表明任何主效应和双向交互作用都没有显著差异。这似乎有点奇怪，特别是因为对表 \@ref(tab:table17-1) 中的数据进行目视检查显示，未收到食品券的组中的 GAIN 存在大量负值，而收到食品券的组中的 GAIN 值主要是正的。因此，人们预计食品券会产生影响，至少在这两个亚组中。还要注意，比较具有数据的单元格对应的 15 个单元格均值的 $F$ 值为 2.67，自由度为 14 和 92. 其对应的 $p$ 值为 0.0026，这似乎也表明这 15 个处理组合均值之间存在显著差异。那么，为什么组的主效应检验（$F = 2.65,p = 0.1068$）不显著呢？当然，答案可能在于 GROUP 的 IV 型 $F$ 值实际正在检验的假设。

<table>
<caption>(\#tab:table17-3)SAS-GLM 检验的 IV 型假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2017.3.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table17-3) 显示了 SAS-GLM IV 型分析所检验的假设。对表 \@ref(tab:table17-3) 的检查显示，组的 IV 型假设将食品券组与无食品券组进行比较，对它们共有的六个单元格平均值进行平均，即对应于 (age = 1, white)、(age = 2, white)、(age = 3, black)、(age = 3, Hispanic)、(age = 3, white) 和 (age = 4, black) 的单元格。这似乎是一个非常合理的假设，因为它对尽可能多的相似类别进行了平均。为了进一步探讨这个检验为什么不显著，让我们看一下对这两组进行比较的相应的 IV 型对比的标准误。标准误估计为

$$\begin{aligned}
\widehat{s.e.}\text{ (type IV group contras}t)& =\hat{\sigma}\sqrt{\sum_i\sum_j\sum_k\left(\frac{c_{ijk}^2}{n_{ijk}}\right)}  \\
&=5.344\sqrt{\left(\frac12+\frac18+\frac14+\frac11+\frac1{31}+\frac11\right)+\left(\frac13+\frac18+\frac16+\frac12+\frac1{20}+\frac11\right)} \\
&=12.047
\end{aligned}$$

请注意，标准误的大小更多地取决于单元格内的小样本量，而不是大样本量。由于 group 的 IV 型对比涉及三个单元格，其中只有一个观测值，因此观测到的标准误必大于

$$5.344\sqrt{\frac11+\frac11+\frac11}=9.256$$

无论其他单元格中的观测数量如何。这说明了数据分析师必须了解的数据分析的另一个方面。如果相应的假设涉及样本量较小的单元格，某些检验可能没有太大的功效。即使观测到所有单元格，这也可能是正确的。也就是说，如果某些单元格样本量较小，即使具有均衡的处理结构，其功效也可能较低。

似乎在 SAS-GLM 的 IV 型分析中，几乎没有自动检验的假设能引起实验者对这些数据的特别兴趣。我们甚至不考虑 I-III 型假设，因为在有缺失单元格的情况下，它们通常没有什么意义。

在这种混乱的实验中，获取有用信息的最安全、最简单的方法是查看三向最小二乘均值，以及它们之间的成对比较。表 \@ref(tab:table17-4) 给出了三向最小二乘均值及其标准误估计，表 \@ref(tab:table17-5) 给出了与最小二乘均值之间成对比较对应的 $p$ 值。

<table>
<caption>(\#tab:table17-4)三向最小二乘均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2017.4.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table17-5)所有最小二乘均值的成对比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2017.5.png">
</td>
  </tr>
</tbody>
</table>

假设我们希望检验 IV 型 group 假设 $\mu_{\mathrm{Y3W}}=\mu_{\mathrm{N3W}}$。也就是说，具有最大样本量的两个单元格均值如何相互比较？这两个单元格对应于表 \@ref(tab:table17-5) 中的 LSMeans 6 和 13，比较这两个均值的相应 $p$ 值为 $p = 0.0004$. 因此，在 age = 3 和 race = white 的母亲亚组，由于食品券的发放，存在显著差异。接下来，考虑 age = 2, whites. 相应的单元格均值为 LSMeans 3 和 10，相应的 $p$ 值为 $p=0.0424$，这也表明由于食品券的发放而产生显著差异。最后，考虑与 age = 3 的 black 对应的群体。这两个均值对应于 LSMeans 4 和 11，并且相应的 $p$ 值为 $p＝0.0017$. 因此，很明显，如果只考虑具有合理样本量的单元格，那么我们将发现由于食品券的发放而产生显著差异。

## 一个完整的分析 {#sec17-3}

由于我们无法检验三因素交互作用，因此不知道是否存在三因素交互作用，我们接下来在第三种处理因素的每个水平上进行双向分析。让我们假设实验者对 race × group 的效应或差异最感兴趣。因此，我们在每个 age 水平上检查了这些双向组合，如表 \@ref(tab:table17-6) 所示。

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table17-6)表 \@ref(tab:table17-1) 中数据每个 age 水平观测到的 race × age 组合</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2017.6.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> “X”表示该单元格在实验中至少观测到一次</td></tr></tfoot>
</table>

对于表 \@ref(tab:table17-6) 中的 age  = 1，我们观察到以下情况：1）不可能检验该 age 组中的 race × group 交互作用，因为不存在测量双因素交互作用的对比；2） 唯一可以检验的比较组的 IV 型假设是 $\mu_{\mathrm{Y1W}}=\mu_{\mathrm{N1W}}$ ；3）age = 1 组中唯一涉及可检验 race 的 IV 型假设是 $\mu_{\mathrm{Y1B}}=\mu_{\mathrm{Y1W}}$；以及 4）假设 $\mu_{\mathrm{Y1B}}=\mu_{\mathrm{N1W}}$ 也可以在 age = 1 组中进行检验，但该假设可能是次要的，因为它涉及 race 和 group 中的不同水平。根据表 \@ref(tab:table17-5)，检验 $\mu_{\mathrm{Y1W}}=\mu_{\mathrm{N1W}}$ 的 $p$ 值为 0.5628。与 3 和 4 中的假设相对应的 $p$ 值分别为 0.5141 和 0.1975.

表 \@ref(tab:table17-7) 和 \@ref(tab:table17-8) 分别给出了可检验的假设、各自检验统计量的 $p$ 值以及 age 组（age = 2 和 age = 4）的（主观）重要性评级。

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table17-7)表 \@ref(tab:table17-1) 中数据关于 age = 2 的可检验假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2017.7.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> $p$ 值来自表 \@ref(tab:table17-5)</td></tr></tfoot>
</table>

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table17-8)表 \@ref(tab:table17-1) 中数据关于 age = 4 的可检验假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2017.8.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> $p$ 值来自表 \@ref(tab:table17-5)</td></tr></tfoot>
</table>

最后，我们检查了 age = 3 的组。在该组中，观测到了所有 race × group 组合；因此，可以检验 race × group 交互作用。如果这种交互作用不显著，我们可以检验 age = 3 的母亲亚组中 race 和 group 的主效应均值。这三个假设由

$$\begin{aligned}
&H_{01}\colon\mu_\text{N3B}{ - \mu _\text{N3W }{ - \mu _\text{Y3B}} }+\mu_\text{Y3W}{ = 0,\,\text{and}} 
\,\,\mu_\text{N3H}{ - \frac 1 2 }(\mu_\text{N3B}+\mu_\text{N3W})-\mu_\text{Y3H }+\frac12(\mu_\text{Y3B}+\mu_\text{Y3W}){=}0 \\
&H_{02}\colon\bar{\mu}_{\text{YЗ}\cdot}=\bar{\mu}_{\text{NЗ}\cdot} \\
&H_{03}\colon\bar{\mu}_{\cdot\text{ЗВ}} = \bar { \mu }_{\cdot\text{ЗН}} = \bar { \mu }_{\cdot\text{ЗW}}
\end{aligned}$$

关于 $H_{01}$ 的检验统计量实际上是由表 \@ref(tab:table17-2) 中所示的原始 SAS-GLM 分析给出的，而其他两个假设的检验则不是。假设 $H_{01}$ 等价于 race × group IV 型 $F$ 值检验的假设，这一点可以通过检查表 \@ref(tab:table17-3) 中 race × group 的可估计函数看出。这三个假设的检验可以通过使用 SAS-GLM 中的 Contrast 语句（如下所示）或手工计算获得。在这种情况下，手工检验可能更容易。当使用 SAS-GLM 时，使用均值模型也更容易。使用此模型，只需要输入与 group × age × race 效应相对应的系数。效应模型所需的 SAS-GLM Contrast 语句如表 \@ref(tab:table17-9) 所示。

<table>
<caption>(\#tab:table17-9)用于检验 $H_{01}$ - $H_{03}$ 的 Contrast 语句</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2017.9.png">
</td>
  </tr>
</tbody>
</table>

上述三个假设的检验统计量可以通过在 \@ref(sec17-2) 节中给出的 SAS 命令中包含以下 Contrast 语句来获得。表 \@ref(tab:table17-9) 中给出了所需的语句。这些选项的结果如表 \@ref(tab:table17-10) 所示。如果实验者也对种族年龄组合的影响感兴趣，可以根据组因子的每个值来分析这些组合。

<table>
<caption>(\#tab:table17-10)SAS Contrast 语句的结果</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2017.10.png">
</td>
  </tr>
</tbody>
</table>

如果实验者感兴趣于 race × age 组合的效应，则可以在 group 因素的每个值上分析这些组合。如果实验者想对每种 race 检查 group × age 组合的效应，也存在类似的情况。这两种分析都可以通过手工进行、通过使用 Contrast 语句，或者在可能的情况下，如本节开头所示，对于 age 因子的不同水平使用三向最小二乘均值。

高阶交叉分类处理结构的分析可以采用与本章所示类似的方法进行。

## 结束语 {#sec17-4}

本章对具有大量缺失处理组合的三向处理结构进行了分析。进行了 SAS-GLM 分析并给出了解释。还提出了 SAS-GLM 分析未回答的问题，并说明了回答这些问题的技术。

## 练习 {#sec17-5}
