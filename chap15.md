# 案例研究：缺失处理组合的双向处理结构 {#chap15}

## 案例研究 {#sec15-1}

在第 \@ref(chap13) 章和第 \@ref(chap14) 章中，我们讨论了在缺失处理组合时，完全随机设计结构中双向处理结构的分析。在本章中，我们将说明当未观测到某些处理组合时，如何在随机完全区组设计中分析双向处理结构。考虑表 \@ref(tab:table15-1) 中的数据，该数据是从第 \@ref(chap12) 章中描述的实验中获得的，但在任何区组中都有一些处理组合未被观测到。图 \@ref(fig:figure15-1) 显示了至少观测一次的处理组合。任何涉及处理组合（脂肪 1、表面活性剂 3）或（脂肪 2、表面活性剂 2）的假设都无法进行检验，除非做出额外的假设。在本次讨论中，令 $FS_{ij}$ 代表将脂肪 $i$ 和表面活性剂 $j$ 分配给随机选择的实验单元时预期的响应。

<table>
<caption>(\#tab:table15-1)第 \@ref(chap5) 章烘焙实验的特定体积</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2015.1.png">
</td>
  </tr>
</tbody>
</table>

<div class="figure" style="text-align: center">
<img src="figure/figure 15.1.png" alt="烘焙实验中观测到的处理组合" width="737" />
<p class="caption">(\#fig:figure15-1)烘焙实验中观测到的处理组合</p>
</div>

为了获得该实验的误差平方和，我们可以使用任何可用的统计包在随机区组设计结构中拟合效应模型或均值模型。拟合效应模型的模型是


```sas
MODEL SPVOL = BLK FAT SURF FAT*SURF;
```

拟合双向脂肪 × 表面活性剂组合的均值模型的模型为


```sas
MODEL SPVOL = BLK FAT*SURF / NOINT;
```

拟合这两个模型后，我们发现误差平方和等于 2.0941，自由度为 11. 因此

$$\hat{\sigma}^2=2.0941/11=0.1904$$

如果观测到所有处理组合，则交互作用假设将有 4 个自由度。由于从未观测到两种处理组合，因此交互作用假设仅保留 2 个自由度。交互空间中两个独立的对比是 $FS_{11}-FS_{12}-FS_{31}+FS_{32}$ 和 $FS_{21}-FS_{23}-FS_{31}+FS_{33}$. SAS^®^-GLM IV 型分析检验这两个对比是否同时为零。检验两个对比同时为零的 $F$ 统计量值为

$$F=(5.4002/2)/0.1904=14.18$$

具有 2 和 11 个自由度。

表 \@ref(tab:table15-2) 给出了所有可能的脂肪的 IV 型假设。4 中的假设通过 SAS-GLM IV 型的脂肪分析自动检验。如果愿意，可以使用 Contrast 语句来检验假设 1-3，也可以使用 Estimate 语句来检验假设 2 和 3。所进行的统计检验的结果也显示在表 \@ref(tab:table15-2) 中。

<table>
<caption>(\#tab:table15-2)关于脂肪的 IV 型假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2015.2.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table15-3) 列出了表面活性剂所有可能的 IV 型假设。通过 SAS-GLM IV 型分析自动检验表面活性剂假设 3 和假设 5 中的最后一个等式。给出的 $F$ 值为 $F = 6.34$，自由度为 2 和 11。所有五个假设都可以使用 Contrast 语句进行检验，除了假设 3 之外的所有假设都可以使用 Estimate 语句进行检验。

<table>
<caption>(\#tab:table15-3)关于表面活性剂的 IV 型假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2015.3.png">
</td>
  </tr>
</tbody>
</table>

由于这些数据中存在显著的交互作用，因此最好通过检查观测到的七种处理组合的最小二乘均值来比较所有观测到的处理组合。表 \@ref(tab:table15-4) 中给出了使用成对 $t$ 检验将它们相互比较的最小二乘均值和 $p$ 值。我们注意到，当设计结构是完全随机设计时，总体单元格均值的最佳估计是每个单元格中观测值的平均值（这在随机区组设计中是不正确的）。使用计算包或使用第 \@ref(chap6) 章的方法可以很容易地获得最佳估计。根据这些估计，我们可以构建图 \@ref(fig:figure15-2) ，其中没有显著差异的均值包含在同一个圆圈中。

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table15-4)成对比较的最小二乘均值、$t$ 统计量和 $p$ 值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2015.4.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> 为了确保整体保护水平，应仅使用与预先计划的比较相关的概率</td></tr></tfoot>
</table>

<div class="figure" style="text-align: center">
<img src="figure/figure 15.2.png" alt="脂肪 × 表面活性剂最小二乘法均值。位于同一圆圈内的均值没有显著差异" width="882" />
<p class="caption">(\#fig:figure15-2)脂肪 × 表面活性剂最小二乘法均值。位于同一圆圈内的均值没有显著差异</p>
</div>

## 结束语 {#sec15-2}

本章说明了在未观测到某些处理组合的情况下，随机完全区组设计中双向处理结构实验的分析。
