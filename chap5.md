
# 多水平设计：裂区、裂条、重复测量及组合 {#chap5}

咨询统计学家并不总是有机会设计他们必须帮助构建适当分析的实验。相反，统计学家必须首先确定研究人员使用的设计实验的类型。识别过程中的第一步也是最重要的一步是确定是否使用了不止一种尺寸的实验单元，如果使用了，则识别每种尺寸的试验单元。正如在本节中显而易见的那样，每种尺寸的实验单元都将具有相关的设计结构和处理结构。在识别出不同尺寸的实验单元后，可以通过组合用于描述与每个尺寸的实验单位相对应的设计结构和处理结构的模型来构建用于进行适当分析的模型。

## 实验单元的尺寸识别——四种基本设计结构 {#sec5-1}

涉及一种以上尺寸的实验单元的设计结构被称为**多水平设计结构** (multilevel designs structures )，并且包括**裂区** (split-plot) 类型的设计结构、**裂条区** (strip-plot) 设计结构、**重复测量** (repeated measures) 设计结构、**分层** (hierarchical) 或**嵌套** (nested) 类型的设计构造以及涉及上述各种组合的设计结构。重复测量和裂区类型的设计结构是相似的，尽管用于进行分析的假设可能不同。裂区和裂条区设计结构是从农业科学发展而来的，但在包括工程和制造业在内的许多其他学科中都有使用，第 \@ref(chap24) 章和第 \@ref(chap25) 章讨论了它们的分析。重复测量设计在社会科学和生物科学中广泛使用，但在大多数领域都适用，因为随着时间的推移，评估处理效应是有意义的，其分析见第 \@ref(chap26) - \@ref(chap28) 章。嵌套处理结构不同于重复测量和裂区设计使用的处理结构，嵌套处理结构在第 \@ref(chap30) 章中进行了描述。

有四种基本的设计结构，而最复杂的设计结构是这些结构的组合。四种基本设计结构是**完全随机设计结构** (completely randomized design 
structure)、**随机完全区组设计结构** (randomized complete block design structure)、**裂区设计结构** (split-plot design structure) 和**裂条区设计结构** (strip-plot design structure). 这些基本设计结构中的每一个都有自己的分析，具有计算所需误差平方和的独特过程。本节介绍了基本设计结构以及计算误差平方和所需的过程。讨论中使用了一个例子，其中处理结构是双向的，一个因素有两个水平，而另一个因素有三个水平，设计结构涉及 $18$ 个实验单元。

该实验包括评估烘烤后纸杯蛋糕的体积，其中有三种配方和两种烹饪温度；因此，处理结构由三种配方和两种烹饪温度的六种组合组成。希望每种处理组合重复 $3$ 次；因此需要烘烤 $18$ 个纸杯蛋糕。图 \@ref(fig:figure5-1) 中的图用于演示使用完全随机设计结构的过程。该过程是将处理结构中的六种处理组合完全随机分配给设计结构中的 $18$ 个实验单元（纸杯蛋糕）。箭头表示每个处理组合分配给三个纸杯蛋糕。通常，该过程是一次制作一个纸杯蛋糕，因此纸杯蛋糕混合和烘焙的顺序对应于实验单元。完全随机设计结构是通过混合一批给定的配方，填写一个纸杯蛋糕表格，然后在设定为特定温度的烤箱中烘烤纸杯蛋糕来完成的。每个纸杯蛋糕都必须由自己的批次制成，也就是说，每批次蛋糕配方只有一个纸杯蛋糕，每个纸杯蛋糕必须在设定的特定温度下在烤箱中自行烘焙。
该过程需要使用 $18$ 批蛋糕混合物和 $18$ 次使用一个或多个烤箱来烘焙纸杯蛋糕。与完全随机化设计结构相关的误差是根据相同处理的纸杯蛋糕之间的变异性计算的。根据温度组合，每个配方中有三个纸杯蛋糕，因此三个纸杯蛋糕体积的变化提供了两个测量纸杯蛋糕误差的自由度。人们应该检验六种处理组合方差的相等性（见第 \@ref(chap2) 章），如果可能的话，将六种处理组的误差估计值合并在一起，以提供 $12$ 个误差自由度，用于测量相同处理的纸杯蛋糕之间的变异性。可用于描述纸杯蛋糕体积的模型为

$$y_{ijk}=\mu+\tau_{i}+\beta_{j}+(\tau\beta)_{ij}+\varepsilon_{ijk},\quad i=1,2,j=1,2,3,\mathrm{and}\,k=1,2,3$$

其中 $y_{ijk}$ 表示用第 $j$ 个配方制作并在第 $i$ 个温度下烘焙的第 $k$ 个纸杯蛋糕的体积，$\mu$ 表示总平均值，$\tau_i$ 表示第 $i$ 个温度的效应，$\beta_j$ 表示第 $j$ 个配方的效应，$(\tau\beta)_{ij}$ 是配方于温度的交互作用，以及 $\varepsilon_{ijk}$ 表示与批次、批次内的纸杯蛋糕以及不同烤箱烘烤的差异性。表 \@ref(tab:table5-1) 包含了完全随机化设计结构中双向处理结构模型的方差分析表，其中根据实验单元或相同处理的纸杯蛋糕的变异性计算了 $12$ 个自由度。与处理结构相关的自由度为五 个，与设计结构相关的自由度为 $12$ 个，所有自由度均分配给误差项。表 \@ref(tab:table5-1) 中标有 $EMS$ 的列给出了方差分析表各行的预期均方形式。函数 $\phi^2(\tau),\phi^2(\beta),\phi^2(\tau\beta)$ 分别表示温度主效应均值、配方主效应均值和交互作用效应的二次函数。当相应的效应不存在时，这些函数是非负的且等于零。在本章的其余部分中可以使用类似的解释。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.1.png" alt="完全随机设计结构中双向处理结构的随机化" width="879" />
<p class="caption">(\#fig:figure5-1)完全随机设计结构中双向处理结构的随机化</p>
</div>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-1)完全随机设计结构中双向处理结构的方差表分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.1.png">
</td>
  </tr>
</tbody>
</table>

第二种基本设计结构是随机完全区组设计，图 \@ref(fig:figure5-2) 中的图表显示了将处理从处理结构分配到设计结构中的实验单元的过程。假设研究人员每天可以制作和烘烤六个纸杯蛋糕，那么实验必须持续三天，才能实现每种处理组合的三次重复。设计结构的实验单元分为三组或六块。接下来，处理结构中的 6 个处理被随机分配到每个块内的 6 个实验单元，如图 \@ref(fig:figure5-2) 中的箭头所示。正如第 4 章所讨论的，随机完全区组设计结构的误差平方和是通过计算设计结构处理结构的交互作用或按处理与区组的交互作用计算来获得的。该设计涉及六种处理和三个区组，因此处理与区组交互作用提供了 (3 - 1)(6 - 1) = 10 个与实验误差相关的自由度。这 10 个自由度合并了区组与温度、区组与配方和区组与温度与配方的交互作用相关的自由度。可用于描述随机完全区组设计结构中这种双向处理结构的纸杯蛋糕体积的模型是

$$y_{ijk}=\mu+\tau_i+\beta_j+(\tau\beta)_{ij}+d_k+\varepsilon_{ijk},\quad i=1,2,~j=1,2,3,~\mathrm{and~}~k=1,2,3$$

其中 $dzk$ 表示第 $k$ 天的效应，即区组因子，$\varepsilon_{ijk}$ 表示与批次、批次内纸杯蛋糕相关的变异性，以及同一天内不同烤箱烘烤的变异性。表 \@ref(tab:table5-2) 包含随机完全区组设计结构中双向处理结构模型的方差分析表，其中有 10 个自由度是根据区组与处理组合的交互作用测量的相同处理的实验单位的变异性计算得出的。该分析具有与完全随机设计结构相同的与处理结构相关的 5 个自由度，但现在，设计结构的 12 个自由度分布在天数（块）和误差之间。表 5.2 中的术语 $\sigma^2_{day}$ 表示 $d_k$ 的方差，$k = 1, 2, 3$. 在本章的其余部分中可以做出类似的解释。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.2.png" alt="随机完全区组设计结构中双向处理结构的随机化方案" width="998" />
<p class="caption">(\#fig:figure5-2)随机完全区组设计结构中双向处理结构的随机化方案</p>
</div>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-2)随机完全区组设计结构中双向处理结构的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.2.png">
</td>
  </tr>
</tbody>
</table>

第三种基本设计结构是裂区设计结构。这里 18 个纸杯蛋糕被分成 6 个大小为 3 的区组，如图 \@ref(fig:figure5-3) 所示。由于处理结构中有六种处理组合，所有六种处理不能发生在一个区组内，因此这是一个不完全区组设计结构。三个纸杯蛋糕，每个配方各一个，将在同一烤箱中在指定温度下烘烤（前两个基本设计结构仅在烤箱中烘烤一个纸杯蛋糕）。图 \@ref(fig:figure5-3) 中的图表显示处理结构已分为两部分，一部分指定为纸杯蛋糕部分，另一部分指定为烤箱部分。分配给每个烤箱的尺寸为三的区组或三个纸杯蛋糕构成温度水平的实验单元。随机化程序的第一部分是将每个温度随机分配给三个烤箱或区组，如图 \@ref(fig:figure5-3) 中从温度水平到烤箱的箭头所示。烤箱内的纸杯蛋糕是配方水平的实验单位，随机化程序是将配方水平随机分配给每个烤箱内的一个纸杯蛋糕。实验单元有两种尺寸，有两种设计和处理结构。烤箱实验单元的处理和设计结构由一个单向处理结构（两个温度水平）组成，该结构为具有六个烤箱的完全随机设计结构。对于单个纸杯蛋糕，处理和设计结构由随机完全区组设计结构中的单向处理结构（配方的三个水平）组成，其中烤箱代表类似实验单元的区组。因此烤箱是温度的实验单位。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.3.png" alt="采用完全随机整区设计结构的裂区随机化方案" width="1004" />
<p class="caption">(\#fig:figure5-3)采用完全随机整区设计结构的裂区随机化方案</p>
</div>

这种设计结构是嵌套或分层结构，因为纸杯蛋糕嵌套在烤箱内。因此，裂区设计也是一种分层设计结构。由于实验单元有两种尺寸，因此称为多水平设计。烤箱是较大尺寸的实验单元，通常被称为**整区** (whole-plot) 。纸杯蛋糕是较小尺寸的实验单元，通常称为**子区** (subplot) 或**裂区** (split-plot). 

分析这种裂区设计的第一步是忽略单个纸杯蛋糕或配方，只考虑两个温度和六个烤箱。图 \@ref(fig:figure5-4) 的显示表明，与实验单元烤箱尺寸对应的设计是完全随机设计结构中的单向处理结构。描述在每个烤箱上测量的响应的模型是

$$y_{ik}^*=\mu+\tau_i+e_{ik},\quad i=1,2,\mathrm{~and~}k=1,2,3$$

其中 $y_{ik}^*$ 表示在分配给第 $i$ 个温度的第 $k$ 个烤箱上测量的响应，$\varepsilon_{ik}$ 表示与烤箱相关的误差项。表 \@ref(tab:table5-3) 包含烤箱模型的方差分析表。对于每种温度水平，都有三个相同处理的烤箱，因此每个温度都有两个自由度，可以测量烤箱在相同处理时的变化。如果烤箱在两个温度之间的方差相等，则可以将这两个方差合并，以提供具有四个自由度的烤箱的误差项。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.4.png" alt="用于裂区设计的烤箱设计和处理结构" width="994" />
<p class="caption">(\#fig:figure5-4)用于裂区设计的烤箱设计和处理结构</p>
</div>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-3)以烤箱为实验单位忽略配方水平的完全随机设计结构中单向处理结构的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.3.png">
</td>
  </tr>
</tbody>
</table>

分析的下一步是忽略温度水平，从而提供一种设计，该设计是随机完全区组设计结构（六个块）中的单向处理结构（三个配方）。随机化过程如图 \@ref(fig:figure5-5) 所示。可以用来描述纸杯蛋糕体积的模型是

$$y_{ijk}=\mu+\beta_j+o_{ik}+\varepsilon_{ijk}^*,\quad i=1,2,\quad j=1,2,3,\mathrm{~and~}\quad k=1,2,3$$

其中 $o_{ik}$ 表示区组效应或烤箱效应，$\varepsilon_{ijk}^*$ 表示与纸杯蛋糕相关的误差。方差分析表在表 \@ref(tab:table5-4) 中，其中残余平方和由配方与烤箱交互作用组成。如果所有烤箱都进行同样的处理，则残差平方和将提供纸杯蛋糕与纸杯蛋糕之间变异性的估计。然而，一些烤箱经受一种温度而另一些烤箱经受另一种温度；因此，配方与烤箱交互作用包括温度与配方交互作用。表 \@ref(tab:table5-4) 中的预期均方并不完全由该模型确定，这由表中的“？”表示。为了获得纸杯蛋糕之间的变异性，需要在设计上进行一次额外的缩简 (reduction). 只考虑那些温度为1的烤箱，如图 \@ref(fig:figure5-6) 中的黑线所示。缩简设计是指随机完全区组设计结构中的单向处理结构，其中三个区组中的处理相同（温度1），方差分析表如表 \@ref(tab:table5-5) 所示。配方与烤箱交互作用提供了烤箱中经过相同处理的纸杯蛋糕变异的测量。该过程通过在第二温度下考虑其他三个烤箱而继续，其中配方与烤箱交互作用提供了额外的四个自由度，用于测量纸杯蛋糕在烤箱内进行相同处理时的变化。当两个温度下蛋糕的方差相等时（参见第 2 章的检验），可以将两个源合并在一起，以提供误差平方和，用于测量蛋糕在 8 个自由度下的变异情况。纸杯蛋糕误差平方和是通过计算温度内的烤箱与配方交互作用来计算的。表 \@ref(tab:table5-6) 给出了最终完整的裂区方差分析表，其中配方与温度交互作用从纸杯蛋糕残差中分离出来。描述具有完全随机设计整区设计结构的裂区设计结构中双向处理结构的数据的模型是

$$y_{ijk}=\mu+\tau_i+\beta_j+(\tau\beta)_{ij}+o_{ik}+\varepsilon_{ijk},\quad i=1,2,~j=1,2,3,\mathrm{~and~}~k=1,2,3$$

其中 $o_{ik}$ 表示温度范围内的烤箱变异，$\varepsilon_{ijk}$ 表示烤箱内纸杯蛋糕的变异。表 \@ref(tab:table5-6) 中仍有 5 个与处理结构相关的自由度，但与设计结构相关的 12 个自由度分布在两个误差项之间，其中有 4 个与烤箱误差分量相关的自由度和 8 个与纸杯蛋糕误差分量相关的自由度。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.5.png" alt="纸杯蛋糕实验单元的随机完全区组设计结构" width="986" />
<p class="caption">(\#fig:figure5-5)纸杯蛋糕实验单元的随机完全区组设计结构</p>
</div>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-4)忽略温度水平的随机完全区组设计结构中单向处理结构的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.4.png">
</td>
  </tr>
</tbody>
</table>

<div class="figure" style="text-align: center">
<img src="figure/figure 5.6.png" alt="仅为温度1烤箱运行的部分设计结构" width="995" />
<p class="caption">(\#fig:figure5-6)仅为温度1烤箱运行的部分设计结构</p>
</div>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-5)随机完全区组设计结构中单向处理结构温度1数据的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.5.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-6)完全随机整区设计结构裂区设计中双向处理结构的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.6.png">
</td>
  </tr>
</tbody>
</table>

作为裂区设计结构的变体，假设研究人员一天内只能运行两次烤箱；因此，这些烤箱被分为三组，每组两个，研究需要三天时间才能完成。图 \@ref(fig:figure5-7) 中的图表描述了每天向烤箱随机分配温度水平，然后随机分配每个烤箱内的配方水平。烤箱设计结构是随机完全区组设计，因此烤箱误差项是通过天与温度交互作用计算的。在烤箱或整区上施加区组结构不会改变纸杯蛋糕或子区的设计结构。仍有六个尺寸为三的区组，所以分析中的纸杯蛋糕部分没有改变。描述具有随机完全区组全区设计结构的裂区设计结构中双向处理结构的纸杯蛋糕体积的模型为

$$y_{iik}=\mu+\tau_i+\beta_j+(\tau\beta)_{ij}+d_k+o_{ik}+\varepsilon_{iik},\quad i=1,2,~j=1,2,3,~\mathrm{and~}~k=1,2,3$$

其中 $d_k$ 表示天的效应，$o_{ik}$ 表示一天内的烤箱变异，$\varepsilon_{ijk}$ 表示一天内烤箱内纸杯蛋糕的变异。表 \@ref(tab:table5-7) 给出了与上述模型相对应的方差分析表，其中包括天、误差（烤箱）和误差（纸杯蛋糕）的行。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.7.png" alt="具有随机完全区组全区（烤箱）设计结构的裂区设计结构图" width="989" />
<p class="caption">(\#fig:figure5-7)具有随机完全区组全区（烤箱）设计结构的裂区设计结构图</p>
</div>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-7)随机完全区组整区设计结构裂区设计中双向处理结构的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.7.png">
</td>
  </tr>
</tbody>
</table>

如上所述，裂区设计结构是不完全区组设计。图 \@ref(fig:figure5-8) 显示了每天内将处理组合分配给烤箱或区组的情况。每个区组内只能出现三种处理组合，并且事实上，在一个区组内只能出现具有相同温度水平的那些处理组合。最终的设计是**部分均衡不完全区组** (partially balanced incomplete block)，其中一些处理组合总是在一个区组内一起出现，而一些处理组合永远不会在一个区组内一起出现。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.8.png" alt="裂区设计结构表示为不完全区组设计结构" width="916" />
<p class="caption">(\#fig:figure5-8)裂区设计结构表示为不完全区组设计结构</p>
</div>

第四种基本设计结构是裂条区. 裂条区设计结构的构建方法是首先将实验单元排列成具有行和列的矩形，如图 \@ref(fig:figure5-9) 所示。在此实验中，使用其中一种配方混合一批蛋糕面团，并从该批次中提取两个纸杯蛋糕。批次是由配方制成的实体，因此批次是配方水平的实验单位。其中一个纸杯蛋糕将在温度 1 下烘烤，另一个纸杯蛋糕将在温度 2 下烘烤。因此，每个烤箱将包含三个纸杯蛋糕，每个配方各有一个。烤箱是分配温度水平的实体，并且是温度水平的实验单元。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.9.png" alt="纸杯蛋糕实验裂条区设计结构图" width="984" />
<p class="caption">(\#fig:figure5-9)纸杯蛋糕实验裂条区设计结构图</p>
</div>

图 \@ref(fig:figure5-9) 中矩形的列对应于蛋糕面团的批次，三种配方被随机分配到列中；也就是说，一列中的两个纸杯蛋糕都来自同一批次，该批次是配方水平的实验单位。行由三个纸杯蛋糕组成，与烤箱相对应，烤箱是温度水平的实验单元。该设计的分析可以通过考虑每种尺寸的实验单元的设计和处理结构来构建。首先，忽略配方，只考虑矩形的行，那么烤箱设计是随机完全区组设计结构中的单向处理结构。如表 \@ref(tab:table5-8) 所示，与烤箱对应的误差项是通过矩形（或天）与温度交互作用来计算的，这误差（烤箱）提供了两个自由度。接下来，忽略温度，只考虑矩形的列，那么批次设计是随机完全区组设计结构中的单向处理结构。如表 \@ref(tab:table5-9) 所示，与批次对应的误差项是通过矩形（或天）与配方的交互作用来计算的，为误差（批次）提供了四个自由度。最后，配方与温度交互作用是没有行效应和列效应的对比 (contrasts that are free of the row effects and free of the column effects)，使纸杯蛋糕成为交互作用比较的实验单元。配方、温度和矩形之间的交互作用提供了纸杯蛋糕误差项，在这种情况下，它有四个自由度。该设计涉及三种尺寸的实验单元，可用于描述该结构数据的模型是

$$y_{ijk}=\mu+\tau_i+\beta_j+(\tau\beta)_{ij}+d_k+o_{ik}+b_{jk}+\varepsilon_{ijk},\quad i=1,2,\quad j=1,2,3,\mathrm{~and~}\quad k=1,2,3$$

其中 $d_k$ 表示矩形或天的效应，$o_{ik}$ 表示矩形内的烤箱效应，$b_{jk}$ 表示矩形内的批次效应，$\varepsilon_{ijk}$ 表示批次、烤箱和矩形内的纸杯蛋糕效应。裂条区模型对应的方差分析表如表 \@ref(tab:table5-10) 所示，其中有 3 个误差项。仍有 5 个与处理结构相关的自由度，但与设计结构相关的 12 个自由度分布为两个矩形、两个误差（烤箱）、四个误差（批次）和四个误差（纸杯蛋糕）。裂条区设计结构是一种多水平设计，但它不是分层设计结构，因为行不嵌套在列内，列也不嵌套在行内。行和列嵌套在一个矩形内，但这就是嵌套停止的地方。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-8)裂条区设计结构中处理结构温度部分的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.8.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-9)裂条区设计结构中处理结构配方部分的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.9.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-10)裂条区设计结构双向处理结构方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.10.png">
</td>
  </tr>
</tbody>
</table>

正如第 \@ref(chap20) 章所讨论的，与效应相对应的预期均方将指示在评估与该效应相关联的假设时使用哪个误差项。一些预期均方涉及带“*”上标的方差，以表明除下标所示之外，还有其他因素影响变异性。

在方差分析表中写入误差项的过程中，只要研究涉及一个以上的实验单位，将与误差相关的实验单位的尺寸包括在括号中，这是本书其余部分中使用的惯例。这一惯例使研究人员能够容易地确定研究中的变异来源。

更复杂的设计结构通常在进行研究所需的材料和时间方面更节省。对于纸杯蛋糕烘焙实验，完全随机和随机完全区组设计结构需要制作18批纸杯蛋糕面团，并且需要使用18次烤箱进行18次烘焙。裂区设计结构涉及18批蛋糕面团（每个纸杯蛋糕一个），但只有六次烘焙。因此，裂区设计只使用了完全随机设计所需的三分之一的时间来烘焙纸杯蛋糕。裂条区需要九批蛋糕面团和六次烘焙。裂条区设计结构只使用了一半批次的蛋糕面团和三分之一的烘焙时间。使用更复杂的设计结构通常节省时间和资源，并且是执行实验更方便的方式，前提是可以解决如上所述的必要的误差项。

随着在设计上施加更多的结构，与最小实验单元尺寸对应的误差相关的自由度会减小。然而，在前面的例子中，与完全随机设计相关的误差包括由于一个批次中的烤箱、批次和纸杯蛋糕引起的变异。裂区设计有两个误差项，其中一部分误差被指定为烤箱之间的变异性，另一部分是纸杯蛋糕之间的变异性（也包括批次间的可变性）。裂条区设计有三个误差项，其中研究中的变异性分为烤箱间的变异性、批次间的变异性和批次内纸杯蛋糕间的变异性。使用更复杂的设计结构为纸杯蛋糕误差项提供了更少的自由度，但该误差项被细化为裂条区设计结构中纸杯蛋糕与纸杯蛋糕之间的变异性，其中涉及完全随机和随机化完全区组中烤箱和批次变异性以及裂区的批次变异性。因此，自由度较小的事实并不一定意味着处理因素之间比较的功效较小，因为重要方差分量的大小也可以减小，从而增大功效。与更简单的设计结构的误差项相比，更复杂的设计具有更少可变性来源的误差项。表 \@ref(tab:table5-11) 总结了与纸杯蛋糕设计的每个误差项相关的变异性来源。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-11)与纸杯蛋糕示例相关的各种设计的每个误差项的变异性来源</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.11.png">
</td>
  </tr>
</tbody>
</table>

以下部分提供了几个多级水平的示例，每个多水平设计可以表示为上述四种基本设计结构中的一种或几种的组合。能够识别与四种基本设计结构中的一个相对应的设计成分，为确定研究中每个尺寸的实验单元的误差源提供了一种方法。能够描述模型中的所有误差项是很重要的，因为在使用软件提取适当的分析时需要这些描述。当在模型的残差或最小实验单位水平上存在重复测量或不等方差时，需要残差平方和的表达式，以便软件对适当的方差和协变量进行建模。重复测量实验的建模将需要知道如何计算残差平方和。大多数作者通过减法计算残差平方和，**但当需要对残差方差和协方差进行建模时，减法方法是不够的**。 

多水平设计有两个重要特点。首先，处理结构至少由一组双向处理组合组成。将多水平设计与第 \@ref(chap4) 章中的设计区分开来的第二个特征是，在一个实验中使用了多个尺寸的实验单元。在设计结构中，每个水平都有一个尺寸的实验单元。每种尺寸的实验单元都有自己的设计结构和处理结构，可以通过组合每种尺寸实验单元的模型来构建模型。由于存在不止一种尺寸的实验单元，因此在分析中使用了不止一个误差项；也就是说，实验中每个实验单元的尺寸都有一个误差项，该误差项也反映在模型中。

本章提供了几个例子来证明使用第 \@ref(sec5-2) 节中讨论的四个基本设计结构所需的原则，以正确识别研究中使用的设计实验。一旦实验者能够使用这些原理来识别本章中讨论的设计实验，她将能够使用它们来识别其他设计的特征。

多水平设计可以用许多不同的方式进行结构化。接下来的一系列例子用于演示识别不同尺寸实验单元的过程，然后将这些信息用于构建一个适当的模型，作为分析的基础。每个例子都包括一个方差分析表，其中列出了变异源、自由度和预期均方（关于计算预期均方的讨论，请参见第 \@ref(chap18)。设计结构与四个基本设计结构相关，从而可以确定误差项的形式。在对数据进行计算机分析之前，列出适当的变化源和相应的分析自由度是很重要的，因为这可以很好地检查是否使用了适当的模型指定来描述数据。

## 分层设计：一种多水平的设计结构 {#sec5-2}

社会科学中经常使用分层设计，其中个体的组形成较大尺寸的实验单元，而组内的个体是较小尺寸的实验单位。例如，一项评估五年级学生数学教学方法的研究涉及从学校系统中选择五年级学生的六个班级，并随机将两种方法中的每一种分配给其中三个班级。班级是教学方法的实验单元。确定教学方法是否对男女学生产生不同的影响是很有意义的。学生是学生性别的实验单位。个体嵌套于班级中的学生性别中，班级又嵌套于教学方法中。本研究采用嵌套设计结构和双向处理结构，采用两种教学方法，由两种性别的学生进行。可以用来描述学生在接受其中一种教学方法后的数学测试成绩的模型是

$$y_{ijkm}=\mu+\tau_i+c_{ij}+\beta_k+(\tau\beta)_{jk}+\varepsilon_{ijkm},\quad i=1,2,j=1,2,3,k=1,2,m=1,2,\ldots,n_{ijk}$$

其中 $y_{ijkm}$ 为第 $i$ 个方法授课的第 $j$ 个班级第 $k$ 个性别的第 $m$ 个学生的成绩，$\mu$ 表示平均分，$\tau_i$ 表示教学方法效应，$c_{ij}$ 表示第 $i$ 个方法授课的第 $j$ 个班级的效应， $\beta_k$ 表示第 $k$ 个性别效应，$(\tau\beta)_{jk}$ 表示通过性别互动的教学方法，$\varepsilon_{ijkm}$ 表示通过教学方法教授的教室的性别内的学生效应。上述模型的方差分析表如表 \@ref(tab:table5-12) 所示，其中班级和学生被假定为随机效应（见第 \@ref(chap18) 章）。系数 $k_1$ 和 $k_2$ 取决于教室内每个性别的学生人数（有关预期均方的评估，请参阅第 \@ref(chap18) 章）。该模型有两个误差项，一项针对班级，一项针对学生，这与裂区设计结构中的结构类似。因此，分层设计与裂区基本设计结构相同。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-12)分层或裂区设计结构的教学方法研究方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.12.png">
</td>
  </tr>
</tbody>
</table>

## 裂区设计结构：两水平设计结构 {#sec5-3}

裂区设计主要用于农业、工业和生物研究，但也可以有效地用于大多数其他研究领域。裂区设计结构是具有两个或多个水平的多水平设计。裂区设计涉及两种或多种处理结构，具有不完全区组设计结构和至少两个不同尺寸的实验单元。裂区设计与重复测量设计的区别在于，处理结构中每个因素的水平可以随机应用于各种规模的实验单元。相反，重复测量设计涉及这样一个步骤：处理结构中至少一个因素（通常是时间）的水平不能随机分配给各个实验单元。以下示例演示了裂区设计的用途，并为识别适当的设计提供指南。

### 例 5.1：烹饪豆子——最简单的裂区或两水平设计结构 {#sec5-3-1}

一位实验者想要研究五种豆子对三种烹饪方法的响应。感兴趣的因变量是豆子煮熟后的嫩度和风味。实验者有一个由 15 个同质行组成的区域。他将五种单向处理结构中的每一种随机分配到三行，从而在完全随机设计结构中生成单向处理结构。品种被分配到行中，如图 \@ref(fig:figure5-10) 所示；因此，行是与品种相关的实验单位。收获时节，每行的豆子都会被放入一个盒子中。对于对一行豆子（或一个盒子）进行的某些测量，该行实验单元的模型是

$$y_{ij}=\mu_i+r_{ij},\quad i=1,2,3,4,5,\quad j=1,2,3$$

其中 $\mu_i$ 表示第 $i$ 个品种的平均值，$r_{ij}$ 表示与分配给第 $j$ 行的第 $i$ 个品种相关的误差。表 \@ref(tab:table5-13) 给出了用于比较品种平均响应的行模型的分析或方差表。每个品种分配了三行，因此，品种内相同处理的行之间的变异性为误差（行）提供了两个自由度。如果五个品种内的行方差相等，则将方差合并在一起以提供与误差（行）相关的 10 个自由度。行设计是完全随机设计结构。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.10.png" alt="用于将品种分配到烹饪豆实验的行的随机化方案" width="831" />
<p class="caption">(\#fig:figure5-10)用于将品种分配到烹饪豆实验的行的随机化方案</p>
</div>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-13)行分析的方差分析表，以比较烹饪豆子例子中的品种</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.13.png">
</td>
  </tr>
</tbody>
</table>

接下来，实验者想要检查烹饪方法。有几种可能的方法来进行这部分实验，这里讨论了两种。首先，实验者可以为种植了给定品种的三行中的每一行分配不同的烹饪方法，如图 \@ref(fig:figure5-11) 所示。图\@ref(fig:figure5-11) 中的排列产生了一个完全随机设计结构中的双向处理结构，其中行是实验单元。然而，每个品种与烹饪方法的组合只有一次重复或行。因此，没有相同处理的行，这意味着并没有实验误差或行方差的度量。方差分析表的结果如表 \@ref(tab:table5-14) 所示。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.11.png" alt="将烹饪方法分配到每个品种内的行的随机化方案，以烹饪豆子为例" width="995" />
<p class="caption">(\#fig:figure5-11)将烹饪方法分配到每个品种内的行的随机化方案，以烹饪豆子为例</p>
</div>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-14)行分析的方差分析表，以烹饪豆子为例</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.14.png">
</td>
  </tr>
</tbody>
</table>

具有零误差自由度的设计不是很理想 [尽管一些分析可以使用 Milliken & Johnson, 1989; Millikan & Graybill,1970; Johnson & Grayboll, 1972 中讨论的双向非重复实验技术 (two-way non-replicated experiment techniques) 进行]。另一种分配烹饪方法的方式避免了零自由度问题，但确实使实验设计和分析更加复杂。

另一种方法是将每盒豆子（每行一盒）分成三批，然后将每种烹饪方法随机分配给一行中三批中的一批。由于烹饪方法被分配给一个批次，因此烹饪方法的实验单元是一个批次。因此，该实验有两种尺寸的实验单元；行（大尺寸）是品种的实验单元，批（小尺寸）是烹饪处理的实验单元。这种分配如图 \@ref(fig:figure5-12) 所示。

批次实验单元的处理和设计结构是随机完全区组设计结构中的单向处理结构，其中行（或盒子）是区组。表 \@ref(tab:table5-15) 给出了此设计批次部分的方差分析表。

行的平方和包括变化的平方和加上表 \@ref(tab:table5-13) 中误差（行）的平方和。如果所有行都经相同处理，则烹饪方法与行的交互作用将提供批次误差项，但一些行种植了品种1，另一些行种植了品种2，依此类推。因此，烹饪方法与行的交互作用包括了烹饪方法与品种的交互作用。所以，只考虑种植了品种1的那些行。品种1中的烹饪方法与行的交互作用提供了(3-1)*(3-1)=4个自由度，用于测量在同一行中经相同处理的批次的变异。批次误差平方和是通过在五个品种中合并烹饪方法与品种的交互作用平方和来获得的，从而为误差（批次）提供20个自由度。描述烹饪豆子示例数据的模型如下

$$y_{ijk}=\mu_{ik}+r_{ij}+\varepsilon_{ijk},\quad i=1,2,3,4,5,j=1,2,3,k=1,2,3$$

其中，$\mu_{ik}$ 表示通过第 $k$ 种方法烹饪的第 $i$ 个品种的平均值，$r_{ij}$ 是分配给第 $i$ 个品种的第 $j$ 行的随机效应，假设其分布为 $N(0,\sigma^2_{\mathrm{row}})$， $\varepsilon_{ijk}$ 表示用第 $k$ 种方法烹饪的第 $i$ 个品种的第 $j$ 行的批次的随机效应，假设其分布为 $N(0,\sigma^2_{\mathrm{batch}})$，还假设 $r_{ij}$ 和 $\varepsilon_{ijk}$ 是独立的随机变量。

在效应模型中，均值 $\mu_{ik}$ 可表示为

$$\mu_{ik}=\mu+\upsilon_i+\omega_k+(\upsilon\omega)_{ik}$$

其中 $\mu$ 表示总体平均值，$\upsilon_i$ 表示第 $i$ 个品种的效应，$\omega_k$ 表示第 $k$ 种烹饪方法的效应，$(\upsilon\omega)_{ik}$ 表示品种与烹饪方法的交互作用。

$$\begin{aligned}y_{ijk}=&\mu+\upsilon_i+r_{ij}\,&&\text{\}row part of model}\\&+\omega_k+(\upsilon\omega)_{ik}+\varepsilon_{ijk}\,&&\text{\}batch part of model}\end{aligned}$$

其中模型的行部分也是模型的批次部分的区组结构。裂区设计可以分两步进行分析，先进行行分析，然后进行批次分析。当数据集均衡时，在将单个模型拟合到数据或进行两步分析过程的情况下，将获得相同的结果。表 \@ref(tab:table5-16) 给出了烹饪豆子例子的方差分析表，它划分了每个实验单元尺寸、行尺寸和批次尺寸的分析。行分析也是分析的批次部分的区组结构，如箭头所示。

豆子实验的整区或行设计结构是完全随机设计结构；因此，这被称为最简单的裂区或两水平设计结构。通常的整区设计结构通常涉及随机完全区组设计结构，如下例所示。

### 示例5.2：磨小麦——通常的裂区或两水平设计结构 {#sec5-3-2}

谷物碾磨实验包括在将小麦粒碾磨或研磨成小麦粉后评估各种小麦品种的特性。实验包括将研磨辊之间的间隙（称为辊缝 (roll     gap)）设置为一个值，然后研磨每种品种的一批（按随机顺序）。接下来，将辊缝更改为另一个值，并研磨新的批次。假设研究人员希望评估三个辊缝和五个品种。因此，五个品种的三个辊缝的双向处理结构的一次重复需要15次轧机运行。一个工作日可以完成一次重复或15次运行，因此需要四天时间才能获得四次重复。环境条件如湿度可能会对研磨过程产生影响，并且这些条件可能每天都不同。因此，为了帮助控制这些条件，将天用作区组因子，其中在一天内获得15种处理组合的一次复制。随机化过程是将三个辊隙的顺序随机分配给每天运行的一个辊缝，然后在每个辊缝内，将品种的顺序随机指定给要研磨的批次。图 \@ref(fig:figure5-13) 显示了随机化过程，其中一天对应于一个区组，一天内的一组五次运行 (run) 是整区实验单元，单次运行是子区实验单元。并不是所有显示品种分配给运行的箭头都会显示出来，因为完整的一组线会使展示混乱不堪。整区设计结构是一个随机完全区组，有四个尺寸为三的区组。整区模型可以表示为

$$y_{ij}=\mu+R_i+d_j+e_{ij},\quad i=1,2,3,j=1,2,3,4,\mathrm{~where~}d_j\thicksim i.i.d.\mathrm{~N}(0,\sigma_{\mathrm{day}}^2)\mathrm{~and~}e_{ij}\thicksim i.i.d.\mathrm{~N}(0,\sigma_{\mathrm{run}}^2)$$

其中 $R_i$ 表示辊缝效应，$d_j$ 表示随机天效应 (random day effect )，$e_{ij}$ 表示随机辊缝运行效应或整区误差。整区方差分析表如表 \@ref(tab:table5-17) 所示，其中运行或整区误差是通过天与辊缝的交互作用计算的；也就是说，随机完全区组设计结构的误差是处理结构与设计结构的交互作用。分析的下一步是确定批此误差的来源。这是通过仅考虑具有一个辊缝（例如1.5mm）的品种来实现的。表 \@ref(tab:table5-18) 中显示了用于比较1.5 mm辊缝下品种的方差分析表。误差项是通过天与品种（或一天内的运行，因为给定的某天只进行一次运行）的交互作用来计算的，这提供了12个自由度来测量在一次运行中相同处理的批次如何变化。这一过程是针对其他两个辊缝设置执行的，每个设置都为批次误差提供12个自由度。如果这三个方差相等，那么它们可以合并到具有36个自由度的批次误差项中。批次误差项可以表示为品种×天（辊缝），也可以表示为“跨辊缝水平的品种与天的交互作用” (variety by day interaction pooled across the levels of roll gap).

包含运行部分和批次部分的模型可以表示为

$$y_{ijk}=\mu+R_i+d_j+e_{ij}+V_k+(RV)_{ik}+\varepsilon_{ijk},\quad i=1,2,3,\quad j=1,2,3,4,\quad k=1,2,\ldots,5,$$

其中

$$d_j\thicksim i.i.d.\mathrm{~}N(0,\sigma_{\mathrm{day}}^2)\text{, }e_{ij}\thicksim i.i.d.\mathrm{~}N(0,\sigma_{\mathrm{run}}^2)\text{, and }\varepsilon_{ijk}\thicksim i.i.d.\mathrm{~}N(0,\sigma_{\mathrm{batch}}^2)$$

在上述模型中，$R_i$ 表示第 $i$ 个辊缝效应，$V_k$ 表示第 $k$ 个品种效应，$(RV)_{ik}$ 表示辊隙与品种交互作用，$\varepsilon_{ijk}$ 表示随机批次效应。$\mu+R_i+d_j+e_{ij}$ 是模型的整区或运行部分，也是模型的子区（或模型的批次部分）的区组结构。模型的批次部分是 $V_k+(RV)_{ik}+\varepsilon_{ijk}$.

包括运行分析和批次分析的最终方差分析表如表 \@ref(tab:table5-19) 所示。括号和箭头表示，模型的运行部分或整区部分是模型批次部分的区组结构。如前所述，**通常的裂区设计是具有随机完全区组整区设计结构的裂区设计**。

### 示例5.3：烘焙面包——具有不完全块设计结构的裂区 {#sec5-3-3}

一位烘焙科学家设计了一项研究，以评估温度对由两种不同配方制成并在三种不同温度下烘焙的面包体积的影响。图 \@ref(fig:figure5-14) 显示了本实验的随机化过程。一天只有两个烤箱可用，而且当天只能使用一次。每个烤箱可以容纳每个配方中的一条面包，但在一天中只能观察到三种温度中的两种。研究人员希望对每个温度进行四次重复，这意味着该研究需要在不同的六天进行。表 \@ref(tab:table5-20) 包含六天的温度分配。图 \@ref(fig:figure5-14) 显示了将温度随机分配到烤箱的情况，其中天是区组，烤箱是温度水平的实验单位。烤箱设计是一种不完全区组设计结构中的单向处理结构（温度的水平）。可以用来描述每个烤箱数据的模型是

$$y_{ij}=\mu+T_i+d_j+o_{ij},\quad(i,j)\in\{(1,1),(2,1),(1,2),(3,2),(2,3),(3,3),(1,4),(3,4),(2,5),(3,5),(1,6),(2,6)\}$$

其中

$$d_j\thicksim i.i.d.N(0,\sigma_{\mathrm{day}}^2)\quad\mathrm{and}\quad o_{ij}\thicksim i.i.d.N(0,\sigma_{\mathrm{oven}}^2)$$

下标 $(i,j$ 在上述指标集中取值，方差分析表见表 \@ref(tab:table5-21). 预期均方是为III型平方和计算的（见第 \@ref(chap10) 章）。烤箱误差项是根据设计结构与处理结构交互作用    （或天与温度的交互作用）来计算的。如果观察到表 \@ref(tab:table5-20) 中的所有单元格，则烤箱误差将有(3-1)(6-1)=10个自由度，但其中六个单元格为空，为烤箱误差提供10-6=4个自由度。面包设计是一种随机完全区组设计结构中的单向处理结构（配方的水平），其中烤箱是区组。并不是所有区组都是一样的，所以考虑那些温度相同的烤箱。温度为160°C的处理结构部分如表 \@ref(tab:table5-22) 所示。表 \@ref(tab:table5-22) 中的设计是具有四个区组或烤箱或天的随机完全区组设计结构中的单向处理结构（两个配方）。根据160°C的温度数据，面包误差项有三个自由度，对应于配方与天的交互作用。160°C下面包分析的方差分析表如表 \@ref(tab:table5-23) 所示。可以为每个温度计算面包误差，如果方差相等，则将三个方差合并，为误差（面包）提供九个自由度。可以用来表示所有数据的模型是

$$\begin{aligned}
&y_{ijk}=\mu+T_i+d_j+o_{ij}+R_k+(TR)_{ik}+\varepsilon_{ijk} \\
&(i,j)\in\{(1,1),(2,1),(1,2),(3,2),(2,3),(3,3),(1,4),(3,4),(2,5),(3,5),(1,6),(2,6)\},k=1,2
\end{aligned}$$

其中

$$d_j\thicksim N(0,\sigma_\mathrm{day}^2),\quad o_{ij}\thicksim i.i.d.\mathrm{~N}(0,\sigma_\mathrm{oven}^2)\quad\mathrm{and}\quad\varepsilon_{ijk}\thicksim i.i.d.\mathrm{~N}(0,\sigma_\mathrm{loaf}^2)$$

该模型的最终方差分析表见表 \@ref(tab:table5-24)，其中预期均方是根据III型平方和计算的。由于不完全区组整区设计结构，因此使用了III型平方和。该模型的烤箱部分为 $\mu+T_i+d_j+o_{ij}$，这也是该模型面包部分的阻塞结构。模型的面包部分是 $R_k+(TR)_{ik}+\varepsilon_{ijk}$. 这个例子表明，**任何类型的设计结构都可以用于实验单元的每个水平或每种尺寸**。

### 示例5.4：展示柜中的肉——复杂裂区或四水平设计 {#sec5-3-4}

一位肉类科学家想研究三个水平的温度（T）、两个水平的包装类型（P）、四个水平的光照类型（L）和四个水平的光照强度（I）对在肉类冷却器中储存七天的肉的颜色的影响。实验有六个冷却器，三种温度（1、3和5°C）分别随机分配给两个冷却器，如图 \@ref(fig:figure5-15) 所示。

每个冷却器在一个4×4的网格上分为16个隔间（图 \@ref(fig:figure5-16)）。光强度通过它们在冷却器表面上方的距离来调节，因此，列中的所有分区 (partitions) 都被分配相同的光强度。四种类型的光随机分配到每列内的四个分区。最后，将两种包装分配给牛排，并将两种包装放入每个分区中。图 \@ref(fig:figure5-16) 显示了这种冷却器的安排方式。

在构建适当的分析之前，必须首先正确识别实验单元的不同尺寸或实验水平。温度水平的实验单位是冷却器。冷却器设计是一种完全随机设计结构中的单向处理结构（T 或温度水平）。如果对每个冷却器进行一次测量，则响应可以通过以下方式建模

$$y_{ij}=\mu+T_i+c_{ij},\quad i=1,2,3,j=1,2,\mathrm{~and~}c_{ij}\thicksim i.i.d.N(0,\sigma_{\mathrm{cooler}}^2)$$

表 \@ref(tab:table5-25) 给出了冷却器模型的方差分析表，其中冷却器误差项是根据两个冷却器在一定温度内的变异、跨三种温度计算得出的。因此，误差（冷却器）有三个自由度。

光照强度水平的实验单位是冷却器中四个隔间的柱 (columns). 柱设计由单向处理结构（I或强度水平）组成，采用随机完全区组设计结构，具有六个区组或冷却器。如果所有冷却器都经相同处理，那么柱误差项将根据强度与冷却器交互作用来计算。但有三种不同的温度，因此将分析限制在分配给1°C的两个冷却器上。分配给1°C的两个冷却器如图 \@ref(fig:figure5-17) 所示。在这一点上，设计是一个单向处理结构（I的四个水平）的随机完全区组设计结构（两个冷却器）。如果对这两个冷却器的每个柱进行测量，则可用于描述响应的模型为

$$y_{1jk}=\mu+I_k+c_{1j}+d_{1jk},\quad j=1,2,\,k=1,2,3,4,\,c_{1j}\sim i.i.d. N(0,\sigma_{\mathrm{cooler}}^2)\,\mathrm{~and~}\,d_{1jk}\sim i.i.d.N(0,\sigma_{\mathrm{colum}n}^2)$$

柱模型的方差分析表如图 \@ref(fig:figure5-17) 所示，其中柱误差项计算为强度与冷却器交互作用，提供三个自由度。对于另外两个温度重复该过程，每个温度都提供具有三个自由度的柱误差项。如果这三个方差相等，则可以将它们合并到具有九个自由度的误差（柱）中。误差（柱）可以用强度×冷却器（温度）表示，读作“跨温度水平的强度与冷却器交互作用” (intensity by cooler interaction pooled across the levels of temperature).

照明类型水平的实验单位是柱的分区。分区设计是一种单向处理结构（L的水平），采用随机完全区组设计结构，有24个区组，由六个冷却器的四个柱组成。如果所有柱都经相同处理，则分区误差项将被计算为列交互的照明类型。但所有列的处理方式并不相同，因为有三种温度和四种强度级别。将结构限制为仅涉及温度为 1°C 且强度为 I~1~ 的那些列，如图 \@ref(fig:figure5-18) 所示。与图 \@ref(fig:figure5-18) 相关的设计是随机完整块设计结构（两列，但每列来自不同的冷却器）中的单向处理结构（四级照明类型）。如果对图 \@ref(fig:figure5-18) 中的每个分区进行一次测量，则可用于描述响应的模型为

$$\begin{aligned}y_{1j1m}&=\mu+L_m+d_{1j1}^*+p_{1j1m},\quad j=1,2,m=1,2,3,4,d_{1j1}^*\thicksim i.i.d.N(0,\sigma_{\mathrm{column}}^{2^*})\mathrm{~and}\\p_{1j1m}&\sim i.i.d.N(0,\sigma_{\mathrm{partition}}^2)\end{aligned}$$

其中 $d_{1j1}^*$ 表示冷却器和冷却器作用柱的组合。该分区模型的方差分析表如图 \@ref(fig:figure5-18) 所示，其中分区误差项计算为照明与冷却器交互作用，提供三个自由度。
这个过程需要对12种温度和强度的组合进行，每种组合为误差（分区）提供三个自由度。如果这12个方差相等，则可以将它们合并为一个项，为误差（分区）提供三个自由度提供36个自由度。误差（分区）项可以表示为光照×强度×冷却器（温度）。

最后，包装水平的实验单位是半分区 (half-partitions)（或牛排）。半分区设计是一种单向处理结构（包装的水平或 P），采用随机完全区组设计结构，具有96个分区，由六个冷却器的四个柱中的四个分区组成。如果所有分区都经相同处理，那么半分区误差项将通过包装与分区交互作用来计算。但并不是所有的分区都经相同处理，因为它们被分配到三个温度水平、四个强度水平和四个照明水平。选择分配给温度为1°C、强度为I~1~和照明为L~1~的分区，如图 \@ref(fig:figure5-19) 所示。如果在图 \@ref(fig:figure5-19) 中对每块牛排或半分区进行一次测量，则可用于描述响应的模型为

$$\begin{aligned}y_{1j11n}&=\mu+P_n+p_{1j1n}^*+\varepsilon_{1j11n},\quad j=1,2,n=1,2,p_{1j1m}^*\thicksim i.i.d.N(0,\sigma_{\mathrm{partition}}^{2*})\quad\mathrm{and}\\\varepsilon_{1j11n}&\sim i.i.d.N(0,\sigma_{\frac12\text{partition}} ^ 2 ) \end{aligned}$$

该半分区模型的方差分析表如图 \@ref(fig:figure5-19) 所示，其中半分区误差项计算为包装与冷却器交互作用，提供一个自由度。这个过程需要对 48 种温度、强度和光照的组合进行。这 48 个平方和合并在一起为误差（半分区）提供 48 个自由度。误差（半分区）项可以表示为包装×光照×强度×冷却器（温度）。

上述讨论提供了研究中的四个水平的模型：冷却器模型、柱模型、分区模型和半分区模型。这些模型可以组合成一个模型，其中因素之间的交互作用被加在一起。由于基本设计结构是每个水平的裂区，因此每个水平的处理结构以及与其上一水平处理结构的交互作用都包含在模型的该水平部分中。可以用来描述数据的模型是

$$\begin{aligned}
\mathcal{Y}_{ijkmn}=& \mu+T_i+\mathcal{C}_{ij}  &&\{\text{cooler part of the model}\\
&+I_k+(TI)_{ik}+d_{ijk} &&\{\text{column part of the model}\\
&+L_m+(TL)_{im}+(IL)_{km}+(TIL)_{ikm}+p_{ijkm} &&\{\text{partition part of the model}\\
&+P_n+(TP)_{in}+(IP)_{kn}+(TIP)_{ikn}+(LP)_{mn}+(TLP)_{imn}+(ILP)_{km}+(TILP)_{ikmn}+\varepsilon_{ijkmn} &&\{\text{half-partition part of the model}\\
&i=1,2,3,\quad j=1,2,\quad k=1,2,3,4,\quad m=1,2,3,4,\quad n=1,2,\\
&c_{ij}\thicksim i.i.d.N(0,\sigma_{\mathrm{cooler}}^2),\quad d_{ijk}\thicksim i.i.d.N(0,\sigma_{\mathrm{column}}^2),\quad  p_{ijkn}\thicksim i.i.d.N(0,\sigma_{\mathrm{partition}}^2),\quad \boldsymbol{\varepsilon}_{ijkmn}\thicksim i.i.d.N(0,\sigma_{\frac12\text{partition}} ^ 2 )
\end{aligned}$$

表 \@ref(tab:table5-26) 给出了上述模型的方差分析表，其中较大的实验单元分析也是下一个较小实验单元的区组结构。这是一个裂-裂-裂区 (split-split-split-plot) 实验模型，由四个水平（或实验单元尺寸）组成，涉及四个误差项，每个水平（或每种实验单元尺寸）对应一个。这也是一种分层设计结构，其中半分区嵌套在分区内，分区嵌套在列内，列嵌套在冷却器内。

## 裂条图设计结构：一种无层次的多水平设计 {#sec5-4}

构造裂条图设计结构的过程是将实验单元排列成矩形。一组因子的水平随机分配给每个矩形的行，另一组因子的水平随机分配给每个矩形的列。因此，行是与第一组因子相关的实验单元，列是与第二组因子相关的实验单元。但因此，一行和一列的单元格或交叉点是与两组因素之间的交互作用比较相关的实验单元。通过一个例子来演示裂条图设计结构的一些用途。

### 示例5.5：制作奶酪 {#sec5-4-1}

一位乳制品科学家设计了一项奶酪制造研究，涉及两种脂肪水平、三种奶酪类型、两种储存温度和两种储存湿度。该实验是一个分两步进行的过程，第一步是用两种脂肪水平为每种类型的奶酪制作一批奶酪。每批奶酪要制成四包一磅重的奶酪。第二步是将奶酪储存在一系列环境条件下，用四周的时间使其固化。该研究的储存部分包括将每批奶酪中的一包放入一个指定储存温度和湿度的室 (chamber) 中。这个过程相当于将24包奶酪排列成四行六列的矩形。该列对应一批奶酪，该行对应一个环境室，如图 \@ref(fig:figure5-20) 所示。这位奶制品科学家有四个实验室可用于研究，但他希望有四次重复。因此，通过在四个不同的月份进行实验获得了四次重复。可以通过评估每种尺寸的实验单元的处理和设计结构来构建模型并分析。室设计如图 \@ref(fig:figure5-21) 所示，其中仅显示了四个月中每个月的温度和湿度水平。室设计是一种双向处理结构（温度水平与湿度水平），采用随机完全区组设计结构，每个月为一个区组。可用于描述在每个室中进行的测量的模型为

$$\begin{aligned}y_{ijk}=&\mu+T_j+H_k+(TH)_{jk}+m_i+c_{ijk},\quad i=1,2,3,4,j=1,2,k=1,2\\&m_i\sim i.i.d.N(0,\sigma_{\mathrm{month}}^2),\text{ and }c_{ijk}\thicksim i.i.d.N(0,\sigma_{\mathrm{chamber}}^2)\end{aligned}$$

表 \@ref(tab:table5-27) 包含室模型的方差分析表，其中室误差项是处理结构与设计结构的交互作用。在处理结构中有四种处理组合，在设计结构中有四个区组；因此，室误差项基于九个自由度。误差（室）项是温度与月交互作用、湿度与月交互作用和温度与湿度与月交互作用项的总和。如果在 SAS-Mixed 程序中使用三向交互项，则三个交互项将合并为一个误差项。室的方差分量用星号表示，因为它也包括由于干酪块引起的变化。月份的方差分量用星号表示，因为它还包括由于批次和奶酪块引起的方差。
