
# 案例研究：子类数不等的均衡双向处理结构 {#chap12}

在本章中，我们分析了在随机完全区组设计中进行的双向处理结构产生的一组数据。实验本应是均衡的，但由于不可预见的情况，一些区组缺少一些处理组合。然而，我们仍然假设每种处理组合至少观察一次。第 \@ref(chap13)-\@ref(chap15) 章讨论了从未观测到某些处理组合的情况。

## 脂肪-表面活性剂示例 {#sec12-1}

一位烘焙科学家希望研究三种不同脂肪分别与三种不同表面活性剂 (surfactant) 混合对面包体积的影响，这些面包由九种不同处理组合（即每种脂肪与表面活性剂的搭配）所制备的面团烘焙而成。实验中采用了四个相同类型但来源各异的面粉 (flour) 作为区组因子，这意味着每种面粉都要使用全部九种处理组合来制作面包。遗憾的是，其中一批酵母失效了，因此使用这批酵母制作出的 10 个面包的数据无法纳入分析范围。庆幸的是，所有九种脂肪与表面活性剂的组合至少有一次观测结果。具体数据见表 \@ref(tab:table12-1).

<table>
<caption>(\#tab:table12-1)烘焙实验中的特定体积</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2012.1.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table12-1) 中的数据采用 SAS^®^-GLM 程序进行分析。由于所有处理组合都至少观测一次，并且在每个区组 （面粉） 中观测到一些相同的处理组合，III 型平方和检验了有趣且易于解释的假设。III 型分析检验的假设与没有缺失数据的情况下检验的假设相同。因此，我们可以预测这些假设是什么，并且我们不需要包括 `E3` 选项来识别正在检验的假设。正规方程是使用零和限制来求解的，因此 `SOLUTION` 选项的结果将不令人感兴趣，因此该选项不用于 SAS-GLM 分析。所有边际均值都是可估的，并且它们的估计值是通过使用 `LSMeans` 选项获得的。

使用表 \@ref(tab:table12-2) 中给出的 SAS 命令对数据进行分析。表 \@ref(tab:table12-3) 给出了表 \@ref(tab:table12-2) 中命令的方差分析结果。脂肪与表面活性剂交互作用的 III 型 $F$值为 $F=8.52$，这在 $p=0.0011$ 水平上显著。因此，应在每个脂肪水平内比较表面活性剂，并在每个表面活性剂水平内比较脂肪。因此，表 \@ref(tab:table12-4) 给出了脂肪表面活性剂最小二乘均值，表 \@ref(tab:table12-5) 给出了这些双向最小二乘均值之间的成对比较。

<table>
<caption>(\#tab:table12-2)表 \@ref(tab:table12-1) 中数据的 SAS 分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2012.2.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table12-3)模型方差分析以及主效应和交互作用的检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2012.3.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table12-4)脂肪 × 表面活性剂最小二乘均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2012.4.png">
</td>
  </tr>
</tbody>
</table>

图 \@ref(fig:figure12-1) 给出了双向最小二乘均值的图形；位于同一圆圈内的样本均值没有显著差异。使用的 $p$ 值如表 \@ref(tab:table12-5) 所示。

<div class="figure" style="text-align: center">
<img src="figure/figure 12.1.png" alt="最小二乘均值图。位于同一圆圈内的均值没有显著差异" width="758" />
<p class="caption">(\#fig:figure12-1)最小二乘均值图。位于同一圆圈内的均值没有显著差异</p>
</div>

<table>
<caption>(\#tab:table12-5)脂肪 × 表面活性剂最小二乘均值之间的成对比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2012.5.png">
</td>
  </tr>
</tbody>
</table>

从图 \@ref(fig:figure12-1) 中，我们可以得出以下观察结果
1. 脂肪 3 与表面活性剂 3 的组合给出的响应显著高于所有其它处理组合给出的响应。
2. 脂肪 3 通常给出的响应显著高于脂肪 1 给出的响应。
3. 与脂肪 1 一起使用时，表面活性剂水平没有差异。

## 结束语 {#sec12-2}

在下一章中，我们将讨论某些处理组合从未被观测到的情况。这种情况下，在选择合适分析方法时需要格外注意。本章中，我们探讨了随机完全区组设计下，当某些处理组合在某些区组中缺失时，如何对均衡双向处理结构进行分析。所描述的分析方法仅适用于每个处理组合至少被观测一次的情况。这里使用 SAS-GLM 进行了数据分析。正如第 \@ref(chap22) 章讨论的那样，若将面粉来源视为随机效应，则能更恰当地分析这份数据。
