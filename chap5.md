
# 多水平设计：裂区、裂条区、重复测量及其组合 {#chap5}

统计学家在咨询工作中并不总是有机会设计他们必须帮助构建适当分析的实验。相反，统计学家首先需要确定研究人员所采用的设计实验的类型。识别过程中的第一步也是最重要的一步是确定是否使用了不止一种**尺寸** (size) 的实验单元，如果使用了，则需要识别每种尺寸的实验单元。正如在本节中显而易见的那样，每种尺寸的实验单元都将具有相关的设计结构和处理结构。在识别出不同尺寸的实验单元后，可以通过组合用于描述与每个尺寸的实验单元对应的设计结构和处理结构的模型，从而构建用于进行适当分析的模型。

## 识别实验单元的尺寸——四种基本设计结构 {#sec5-1}

涉及一种以上尺寸的实验单元的设计结构被称为**多水平设计结构** (multilevel designs structures)，并且包括**裂区** (split-plot) 类型的设计结构、**裂条区** (strip-plot) 设计结构、**重复测量** (repeated measures) 设计结构、**分层** (hierarchical) 或**嵌套** (nested) 类型的设计结构以及涉及上述各种组合的设计结构。重复测量和裂区类型的设计结构是相似的，尽管用于进行分析的假设可能不同。裂区和裂条区设计结构是从农业科学发展而来的，但在包括工程和制造业在内的许多其他学科中都有使用，第 \@ref(chap24) 章和第 \@ref(chap25) 章讨论了它们的分析。重复测量设计在社会科学和生物科学中广泛使用，但在大多数领域都适用，因为随着时间的推移，评估处理效应是有意义的，其分析见第 \@ref(chap26) - \@ref(chap28) 章。嵌套处理结构不同于重复测量和裂区设计使用的处理结构，嵌套处理结构在第 \@ref(chap30) 章中进行了描述。

有四种基本的设计结构，而最复杂的设计结构是这些结构的组合。四种基本设计结构是**完全随机设计结构** (completely randomized design 
structure)、**随机完全区组设计结构** (randomized complete block design structure)、**裂区设计结构** (split-plot design structure) 和**裂条区设计结构** (strip-plot design structure). 这些基本设计结构中的每一个都有自己的分析，计算所需误差平方和的过程也是独特的。本节介绍了基本设计结构以及计算误差平方和所需的过程。讨论中使用了一个例子，其中处理结构是双向的，一个因素有两个水平，另一个因素有三个水平，设计结构涉及 18 个实验单元。

该实验包括评估纸杯蛋糕 (cupcake) 烘焙后的体积，其中有三种配方和两种烹饪温度；因此，处理结构由三种配方和两种烹饪温度的六种组合组成。希望每种处理组合重复 3 次；因此需要烘焙 18 个纸杯蛋糕。图 \@ref(fig:figure5-1) 中的图用于演示使用完全随机设计结构的过程。该过程是将处理结构中的六种处理组合完全随机分配给设计结构中的 18 个实验单元（纸杯蛋糕）。箭头表示每个处理组合分配给三个纸杯蛋糕。通常，该过程是一次制作一个纸杯蛋糕，因此纸杯蛋糕混合和烘焙的顺序与实验单元相对应。完全随机设计结构是通过混合一批给定的配方，填充一个纸杯蛋糕形状，然后在设置为特定温度的烤箱 (oven) 中烘焙纸杯蛋糕来完成的。每个纸杯蛋糕都必须用单独一批材料制作，也就是说，每个批次的蛋糕配方只能制作一个纸杯蛋糕，并且每个纸杯蛋糕都必须在设置为特定温度的烤箱中单独烘烤。该过程需要使用 18 批蛋糕混合物和 18 次使用一个或多个烤箱来烘焙纸杯蛋糕。与完全随机设计结构相关的误差是根据相同处理的纸杯蛋糕之间的变异性计算的。根据温度组合，每个配方中有三个纸杯蛋糕，因此三个纸杯蛋糕体积的变异提供了两个测量纸杯蛋糕误差的自由度。人们应该检验六种处理组合方差的相等性（见第 \@ref(chap2) 章），如果可能的话，将六种处理组的误差估计值合并在一起，以提供 $12$ 个误差自由度，用于测量相同处理的纸杯蛋糕之间的变异性。可用于描述纸杯蛋糕体积的模型为

$$y_{ijk}=\mu+\tau_{i}+\beta_{j}+(\tau\beta)_{ij}+\varepsilon_{ijk},\quad i=1,2,j=1,2,3,\mathrm{and}\,k=1,2,3$$

其中 $y_{ijk}$ 表示用第 $j$ 个配方制作并在第 $i$ 个温度下烘焙的第 $k$ 个纸杯蛋糕的体积，$\mu$ 表示总平均值，$\tau_i$ 表示第 $i$ 个温度的效应，$\beta_j$ 表示第 $j$ 个配方的效应，$(\tau\beta)_{ij}$ 是温度与配方交互作用，以及 $\varepsilon_{ijk}$ 表示与批次、批次内的纸杯蛋糕以及不同烤箱烘焙的差异相关的变异性。表 \@ref(tab:table5-1) 包含了完全随机设计结构中双向处理结构模型的方差分析表，其中根据实验单元或相同处理的纸杯蛋糕的变异性计算了 12 个自由度。与处理结构相关的自由度为 5，与设计结构相关的自由度为 12 个，所有自由度均分配给误差项。表 \@ref(tab:table5-1) 中标有 $EMS$ 的列给出了方差分析表各行的预期均方形式。函数 $\phi^2(\tau),\phi^2(\beta),\phi^2(\tau\beta)$ 分别表示温度主效应均值、配方主效应均值和交互作用效应的二次函数。当相应的效应不存在时，这些函数是非负的且等于零。在本章的其余部分中可以使用类似的解释。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.1.png" alt="完全随机设计结构中双向处理结构的随机化" width="879" />
<p class="caption">(\#fig:figure5-1)完全随机设计结构中双向处理结构的随机化</p>
</div>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-1)完全随机设计结构中双向处理结构的方差分析表</caption>
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

第二种基本设计结构是随机完全区组设计，图 \@ref(fig:figure5-2) 中的图表显示了将处理从处理结构分配到设计结构中的实验单元的过程。假设研究人员每天可以制作和烘焙六个纸杯蛋糕，那么实验必须持续三天，才能实现每种处理组合的三次重复。设计结构的实验单元分为三个尺寸为六的区组。接下来，处理结构中的六个处理被随机分配到每个块内的六个实验单元，如图 \@ref(fig:figure5-2) 中的箭头所示。正如第 \@ref(chap4) 章所讨论的，随机完全区组设计结构的误差平方和是通过设计结构与处理结构交互作用来计算的，或通过处理与区组的交互作用来计算的。该设计涉及六种处理和三个区组，因此处理与区组交互作用提供了 (3 - 1)(6 - 1) = 10 个与实验误差相关的自由度。这 10 个自由度合并了区组与温度、区组与配方和区组与温度与配方的交互作用相关的自由度。可用于描述随机完全区组设计结构中这种双向处理结构的纸杯蛋糕体积的模型是

$$y_{ijk}=\mu+\tau_i+\beta_j+(\tau\beta)_{ij}+d_k+\varepsilon_{ijk},\quad i=1,2,~j=1,2,3,~\mathrm{and~}~k=1,2,3$$

其中 $d_k$ 表示第 $k$ 天的效应，即区组因素，$\varepsilon_{ijk}$ 表示与批次、批次内纸杯蛋糕相关的变异性，以及同一天内不同烤箱烘焙的变异性。表 \@ref(tab:table5-2) 包含随机完全区组设计结构中双向处理结构模型的方差分析表，其中有 10 个自由度是根据区组与处理组合的交互作用测量的相同处理的实验单元的变异性计算得出的。该分析具有与完全随机设计结构相同的与处理结构相关的 5 个自由度，但现在，设计结构的 12 个自由度分布在天（区组）和误差之间。表 \@ref(tab:table5-2) 中的 $\sigma^2_{day}$ 项表示 $d_k$ 的方差，$k = 1, 2, 3$. 在本章的其余部分中可以做出类似的解释。

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

第三种基本设计结构是裂区设计结构。这里 18 个纸杯蛋糕被分成六个尺寸为三的区组，如图 \@ref(fig:figure5-3) 所示。由于处理结构中有六种处理组合，所有六种处理不能发生在一个区组内，因此这是一个不完全区组设计结构。三个纸杯蛋糕，每个配方各一个，将在同一烤箱中在指定温度下烘焙（前两个基本设计结构仅在一个烤箱中烘焙一个纸杯蛋糕）。图 \@ref(fig:figure5-3) 中的图表显示处理结构分为两部分，一部分指定为纸杯蛋糕部分，另一部分指定为烤箱部分。分配给每个烤箱的尺寸为三的区组或三个纸杯蛋糕构成温度水平的实验单元。随机化程序的第一部分是将每个温度随机分配给三个烤箱或区组，如图 \@ref(fig:figure5-3) 中从温度水平到烤箱的箭头所示。烤箱内的纸杯蛋糕是配方水平的实验单元，随机化程序是将配方水平随机分配给每个烤箱内的一个纸杯蛋糕。实验单元有两种尺寸，有两种设计和处理结构。烤箱实验单元的处理和设计结构由单向处理结构（两个温度水平）和完全随机设计结构（六个烤箱）组成。对于单个纸杯蛋糕，处理和设计结构由随机完全区组设计结构中的单向处理结构（配方的三个水平）组成，其中烤箱代表类似实验单元的区组。因此烤箱是温度的实验单元。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.3.png" alt="采用完全随机整区设计结构的裂区随机化方案" width="1004" />
<p class="caption">(\#fig:figure5-3)采用完全随机整区设计结构的裂区随机化方案</p>
</div>

这种设计结构是嵌套或分层结构，因为纸杯蛋糕嵌套在烤箱内。因此，**裂区设计也是一种分层设计结构**。由于实验单元有两种尺寸，因此称为多水平设计。烤箱是较大尺寸的实验单元，通常被称为**整区** (whole-plot) 。纸杯蛋糕是较小尺寸的实验单元，通常称为**子区** (subplot) 或裂区。 

分析这种裂区设计的第一步是忽略单个纸杯蛋糕或配方，只考虑两个温度和六个烤箱。图 \@ref(fig:figure5-4) 的显示表明，与实验单元烤箱尺寸对应的设计是完全随机设计结构中的单向处理结构。描述在每个烤箱上测量的响应的模型是

$$y_{ik}^*=\mu+\tau_i+e_{ik},\quad i=1,2,\mathrm{~and~}k=1,2,3$$

其中 $y_{ik}^*$ 表示在分配给第 $i$ 个温度的第 $k$ 个烤箱上测量的响应，$\varepsilon_{ik}$ 表示与烤箱相关的误差项。表 \@ref(tab:table5-3) 包含烤箱模型的方差分析表。对于每种温度水平，都有三个相同处理的烤箱，因此每个温度都有两个自由度，可以测量烤箱在相同处理时的变异。如果烤箱在两个温度之间的方差相等，则可以将这两个方差合并，以提供具有四个自由度的烤箱的误差项。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.4.png" alt="用于裂区设计的烤箱设计和处理结构" width="994" />
<p class="caption">(\#fig:figure5-4)用于裂区设计的烤箱设计和处理结构</p>
</div>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-3)以烤箱为实验单元忽略配方水平的完全随机设计结构中单向处理结构的方差分析表</caption>
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

分析的下一步是忽略温度水平，从而提供一种设计，该设计是随机完全区组设计结构（六个区组）中的单向处理结构（三个配方）。随机化过程如图 \@ref(fig:figure5-5) 所示。可以用来描述纸杯蛋糕体积的模型是

$$y_{ijk}=\mu+\beta_j+o_{ik}+\varepsilon_{ijk}^*,\quad i=1,2,\quad j=1,2,3,\mathrm{~and~}\quad k=1,2,3$$

其中 $o_{ik}$ 表示区组效应或烤箱效应，$\varepsilon_{ijk}^*$ 表示与纸杯蛋糕相关的误差。方差分析表在表 \@ref(tab:table5-4) 中，其中残差平方和由配方与烤箱交互作用组成。如果所有烤箱都进行同样的处理，则残差平方和将提供纸杯蛋糕与纸杯蛋糕之间变异性的估计。然而，一些烤箱经受一种温度而另一些烤箱经受另一种温度；因此，配方与烤箱交互作用包括温度与配方交互作用。表 \@ref(tab:table5-4) 中的预期均方并不完全由该模型确定，这由表中的“?”表示。为了获得纸杯蛋糕之间的变异性，需要在设计上进行一次额外的缩减 (reduction). 只考虑那些温度为 1 的烤箱，如图 \@ref(fig:figure5-6) 中的黑线所示。缩减设计是指随机完全区组设计结构中的单向处理结构，其中三个区组中的处理相同（温度 1），方差分析表如表 \@ref(tab:table5-5) 所示。配方与烤箱交互作用提供了烤箱中经相同处理的纸杯蛋糕变异的测量。该过程继续考虑第二温度下的其他三个烤箱，其中配方与烤箱交互作用提供了额外的四个自由度，用于测量纸杯蛋糕在烤箱内经相同处理时的变异。当两个温度下蛋糕的方差相等时（参见第 \@ref(chap2) 章的检验），可以将两个来源合并在一起，提供八个自由度的误差平方和以测量蛋糕的变异情况。纸杯蛋糕的误差平方和是通过计算每个温度下烤箱与配方的交互作用，并将这些交互作用在不同温度下进行合并来计算的 (computing the oven by recipe interaction within a temperature pooled across temperatures). 表 \@ref(tab:table5-6) 给出了最终完整的裂区方差分析表，其中配方与温度交互作用从纸杯蛋糕残差中分离出来。描述在裂区设计结构中具有完全随机设计整区设计结构的双向处理结构数据的模型是

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

作为裂区设计结构的变体，假设研究人员一天内只能运行两次烤箱；因此，这些烤箱被分为三组，每组两个，研究需要三天时间才能完成。图 \@ref(fig:figure5-7) 中的图表描述了每天向烤箱随机分配温度水平，然后随机分配每个烤箱内的配方水平。烤箱设计结构是随机完全区组设计，因此烤箱误差项通过天与温度交互作用来计算。在烤箱或整区上施加区组结构不会改变纸杯蛋糕或子区的设计结构。仍有六个尺寸为三的区组，所以分析中的纸杯蛋糕部分没有改变。描述在裂区设计结构中具有完全随机区组整区设计结构的双向处理结构纸杯蛋糕体积的模型为

$$y_{iik}=\mu+\tau_i+\beta_j+(\tau\beta)_{ij}+d_k+o_{ik}+\varepsilon_{iik},\quad i=1,2,~j=1,2,3,~\mathrm{and~}~k=1,2,3$$

其中 $d_k$ 表示天的效应，$o_{ik}$ 表示一天内的烤箱变异，$\varepsilon_{ijk}$ 表示一天内烤箱内纸杯蛋糕的变异。表 \@ref(tab:table5-7) 给出了与上述模型相对应的方差分析表，其中包括天、error (oven) 和 error (cupcake) 的行。

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

第四种基本设计结构是裂条区。 裂条区设计结构的构建方法是首先将实验单元排列成具有行和列的矩形，如图 \@ref(fig:figure5-9) 所示。在此实验中，使用其中一种配方混合一批蛋糕面团，并从该批次中提取两个纸杯蛋糕。由于这批面团是由一个配方制成的，因此批次是配方水平的实验单元。其中一个纸杯蛋糕将在温度 1 下烘焙，另一个纸杯蛋糕将在温度 2 下烘焙。因此，每个烤箱将包含三个纸杯蛋糕，每个配方各有一个。烤箱是分配温度水平的实体，并且是温度水平的实验单元。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.9.png" alt="纸杯蛋糕实验裂条区设计结构图" width="984" />
<p class="caption">(\#fig:figure5-9)纸杯蛋糕实验裂条区设计结构图</p>
</div>

图 \@ref(fig:figure5-9) 中矩形的列对应于蛋糕面团的批次，三种配方被随机分配到列中；也就是说，一列中的两个纸杯蛋糕都来自同一批次，该批次是配方水平的实验单元。行由三个纸杯蛋糕组成，与烤箱相对应，烤箱是温度水平的实验单元。该设计的分析可以通过考虑每种尺寸的实验单元的设计和处理结构来构建。首先，忽略配方，只考虑矩形的行，那么烤箱设计是随机完全区组设计结构中的单向处理结构。如表 \@ref(tab:table5-8) 所示，与烤箱对应的误差项通过矩形（或天）与温度交互作用来计算，这为 error (oven) 提供了两个自由度。接下来，忽略温度，只考虑矩形的列，那么批次设计是随机完全区组设计结构中的单向处理结构。如表 \@ref(tab:table5-9) 所示，与批次对应的误差项通过矩形（或天）与配方交互作用来计算，为 error (batch) 提供了四个自由度。最后，配方与温度交互作用是通过排除行效应和列效应的对比来计算的 (contrasts that are free of the row effects and free of the column effects)，使纸杯蛋糕成为交互作用比较的实验单元。配方、温度和矩形之间的交互作用提供了纸杯蛋糕误差项，在这种情况下，它有四个自由度。该设计涉及三种尺寸的实验单元，可用于描述该结构数据的模型是

$$y_{ijk}=\mu+\tau_i+\beta_j+(\tau\beta)_{ij}+d_k+o_{ik}+b_{jk}+\varepsilon_{ijk},\quad i=1,2,\quad j=1,2,3,\mathrm{~and~}\quad k=1,2,3$$

其中 $d_k$ 表示矩形或天的效应，$o_{ik}$ 表示矩形内的烤箱效应，$b_{jk}$ 表示矩形内的批次效应，$\varepsilon_{ijk}$ 表示批次、烤箱和矩形内的纸杯蛋糕效应。裂条区模型对应的方差分析表如表 \@ref(tab:table5-10) 所示，其中包含了三个误差项。处理结构仍然与五个自由度相关联，而与设计结构相关的12个自由度则分配为：两个用于矩形、两个用于 error (oven)、四个用于 error (batch) 以及四个用于 error (cupcake). **裂条区设计结构是一种多水平设计，但它不是分层设计结构**，因为行不嵌套在列内，列也不嵌套在行内。行和列是嵌套在矩形（或天）内的，但嵌套结构在这里就停止了。

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

正如第 \@ref(chap20) 章所讨论的，与效应相对应的预期均方将指示在评估与该效应相关的假设时使用哪个误差项。一些预期均方涉及带“*”上标的方差，以表明除下标所示之外，还有其他因素影响变异性。

在向方差分析表中写入误差项的过程中，只要设计涉及一个以上尺寸的实验单元，将与误差相关的实验单元的尺寸包括在括号中，这是本书其余部分中使用的惯例。这一惯例使研究人员能够轻松地识别研究的变异来源。

更复杂的设计结构在材料和时间方面往往更加节约。对于纸杯蛋糕烘焙实验，完全随机和随机完全区组设计结构需要制作 18 批纸杯蛋糕面团，并且需要使用 18 次烤箱进行 18 次烘焙。裂区设计结构涉及 18 批蛋糕面团（每个纸杯蛋糕一批），但只需要烘焙六次。因此，与完全随机设计相比，裂区设计只使用了三分之一的时间来烘焙纸杯蛋糕。裂条区需要九批蛋糕面团和六次烘焙。裂条区设计结构只使用了一半的蛋糕面团批次和三分之一的烘焙时间。如上所述，**只要能够解决必要的误差项，使用更复杂的设计结构通常可以节省时间和资源，并且是进行实验的更方便的方式**。

随着在设计上施加更多的结构，与最小实验单元尺寸对应的误差自由度会减少。然而，在前面的例子中，与完全随机设计相关的误差包括烤箱、批次和纸杯蛋糕引起的变异。裂区设计有两个误差项，其中一部分误差是由于烤箱之间的变异性，另一部分是由于纸杯蛋糕之间的变异性（也包括批次间的变异性）。裂条区设计有三个误差项，其中研究中的变异性分为烤箱间的变异性、批次间的变异性和批次内纸杯蛋糕间的变异性。使用更复杂的设计结构为纸杯蛋糕误差项提供了更少的自由度，但该误差项在裂条区设计结构中被细化为纸杯蛋糕与纸杯蛋糕之间的变异性，其中涉及在完全随机和随机完全区组设计中的烤箱和批次的变异性，以及裂区设计中的批次的变异性。因此，尽管自由度较少，但并不一定意味着处理因素之间比较的功效较小，因为重要方差分量的幅度也可能减小，从而增大功效。与更简单的设计结构的误差项相比，更复杂的设计具有更少变异性来源的误差项。表 \@ref(tab:table5-11) 总结了与纸杯蛋糕设计的每个误差项相关的变异性来源。

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

以下各节提供了几个多水平设计的示例，每个多水平设计可以表示为上述四种基本设计结构中的一种或几种的组合。能够识别对应于四种基本设计结构之一的设计部分，为确定研究中每个实验单元尺寸的误差来源提供了一种方法。能够描述模型中的所有误差项是很重要的，因为在使用软件进行适当的分析时需要这些描述。当在模型的残差或最小实验单元水平上存在重复测量或不等方差时，需要残差平方和的表达式，以便软件对适当的方差和协变量进行建模。重复测量实验的建模将需要知道如何计算残差平方和。大多数作者通过减法计算残差平方和，**但当需要对残差方差和协方差建模时，减法方法是不够的**。 

多水平设计有两个重要特点。首先，处理结构至少由一组双向处理组合组成。将多水平设计与第 \@ref(chap4) 章中的设计区分开来的第二个特征是，在一个实验中使用了多个尺寸的实验单元。在设计结构中，每个水平都有一种尺寸的实验单元。每种尺寸的实验单元都有自己的设计结构和处理结构，并且可以通过组合每种尺寸实验单元的模型来构建模型。由于存在不止一种尺寸的实验单元，因此在分析中使用了多个误差项；也就是说，实验中每种尺寸的实验单元都有一个误差项，这也反映在模型中。

本章提供了几个例子来演示使用第 \@ref(sec5-2) 节中讨论的四个基本设计结构所需的原则，以正确识别研究中使用的设计实验。一旦实验者能够使用这些原则来识别本章中讨论的设计实验，她就能够使用它们来识别其他设计的特征。

多水平设计可以用许多不同的方式进行结构化。接下来的一系列例子用于演示识别不同尺寸实验单元的过程，然后将这些信息用于构建一个适当的模型，作为分析的基础。每个例子都包括一个方差分析表，其中列出了变异来源、自由度和预期均方（关于计算预期均方的讨论，请参见第 \@ref(chap18)。设计结构与四个基本设计结构相关，从而可以确定误差项的形式。在对数据进行计算机分析之前，列出适当的变异来源和相应的自由度进行分析是很重要的，因为这可以很好地检查是否使用了适当的模型指定来描述数据。

## 分层设计：一种多水平的设计结构 {#sec5-2}

社会科学中经常使用分层设计，其中个体的组形成较大尺寸的实验单元，而组内的个体是较小尺寸的实验单元。例如，一项评估五年级学生数学教学方法的研究涉及从学校系统中选择五年级学生的六个班级，并随机将两种方法中的每一种分配给其中三个班级。班级是教学方法的实验单元。确定教学方法是否对男女学生产生不同的影响是很有意义的。学生是学生性别的实验单元。个体嵌套于班级中的学生性别中，班级又嵌套于教学方法中。这项研究涉及嵌套设计结构和双向处理结构，其中有两种教学方法和两种性别的学生。可以使用以下模型来描述学生在接受其中一种教学方法后在数学测试中的得分

$$y_{ijkm}=\mu+\tau_i+c_{ij}+\beta_k+(\tau\beta)_{jk}+\varepsilon_{ijkm},\quad i=1,2,j=1,2,3,k=1,2,m=1,2,\ldots,n_{ijk}$$

其中 $y_{ijkm}$ 为第 $i$ 个方法授课的第 $j$ 个班级第 $k$ 个性别的第 $m$ 个学生的成绩，$\mu$ 表示平均分，$\tau_i$ 表示教学方法效应，$c_{ij}$ 表示第 $i$ 个方法授课的第 $j$ 个班级的效应， $\beta_k$ 表示第 $k$ 个性别效应，$(\tau\beta)_{jk}$ 表示教学方法与性别交互作用，$\varepsilon_{ijkm}$ 表示采用某种教学方法的班级内某个性别的学生效应。上述模型的方差分析表如表 \@ref(tab:table5-12) 所示，其中班级和学生被假定为随机效应（见第 \@ref(chap18) 章）。系数 $k_1$ 和 $k_2$ 取决于每个班级内每个性别的学生人数（有关预期均方的评估，请参阅第 \@ref(chap18) 章）。该模型有两个误差项，一项针对班级，一项针对学生，这与裂区设计结构中的结构类似。因此，分层设计与裂区基本设计结构相同。

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

裂区设计主要用于农业、工业和生物研究，但也可以有效地用于大多数其他研究领域。裂区设计结构是具有两个或多个水平的多水平设计。裂区设计涉及两种或多种处理结构，具有不完全区组设计结构和至少两种不同尺寸的实验单元。裂区设计与重复测量设计的区别在于，处理结构中每个因素的水平可以随机应用于各种尺寸的实验单元。相反，重复测量设计涉及这样一个步骤：处理结构中至少一个因素（通常是时间）的水平不能随机分配给各个实验单元。以下示例演示了裂区设计的用途，并为识别适当的设计提供指南。

### 示例 5.1：烹饪大豆——最简单的裂区或两水平设计结构 {#sec5-3-1}

一位实验者想要研究五种大豆对三种烹饪方法的响应。感兴趣的因变量是大豆煮熟后的嫩度和风味。实验者有一个由 15 个同质行 (rows) 组成的区域。他将五种单向处理结构中的每一种随机分配到三行，从而在完全随机设计结构中生成单向处理结构。品种被分配到行中，如图 \@ref(fig:figure5-10) 所示；因此，行是与品种相关的实验单元。在收获时节，每行的大豆都会被放入一个盒子 (box) 中。对一行大豆（或一个盒子）进行的某些测量，行实验单元的模型是

$$y_{ij}=\mu_i+r_{ij},\quad i=1,2,3,4,5,\quad j=1,2,3$$

其中 $\mu_i$ 表示第 $i$ 个品种的平均值，$r_{ij}$ 表示与分配给第 $j$ 行的第 $i$ 个品种相关的误差。表 \@ref(tab:table5-13) 给出了用于比较品种平均响应的行模型的分析或方差表。每个品种分配了三行，因此，品种内相同处理的行之间的变异性为 error (row) 提供了两个自由度。如果五个品种内的行方差相等，则将方差合并在一起以提供与 error (row) 相关的 10 个自由度。行设计是完全随机设计结构。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.10.png" alt="用于将品种分配到烹饪豆实验的行的随机化方案" width="831" />
<p class="caption">(\#fig:figure5-10)用于将品种分配到烹饪豆实验的行的随机化方案</p>
</div>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-13)行分析的方差分析表，以比较烹饪大豆例子中的品种</caption>
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

接下来，实验者想要检查烹饪方法。有几种可能的方法来进行这部分实验，这里讨论了两种。首先，实验者可以为种植了给定品种的三行中的每一行分配不同的烹饪方法，如图 \@ref(fig:figure5-11) 所示。图 \@ref(fig:figure5-11) 中的安排产生了一个完全随机设计结构中的双向处理结构，其中行是实验单元。然而，每个品种与烹饪方法的组合只有一次重复或行。因此，没有相同处理的行，这意味着并没有实验误差或行方差的度量。方差分析表的结果如表 \@ref(tab:table5-14) 所示。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.11.png" alt="将烹饪方法分配到每个品种内的行的随机化方案，以烹饪大豆为例" width="995" />
<p class="caption">(\#fig:figure5-11)将烹饪方法分配到每个品种内的行的随机化方案，以烹饪大豆为例</p>
</div>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-14)行分析的方差分析表，以烹饪大豆为例</caption>
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

另一种方法是将每盒大豆（每行一盒）分成三批，然后将每种烹饪方法随机分配给一行中三批中的一批。由于烹饪方法被分配给一个批次 (batch)，因此烹饪方法的实验单元是一个批次。因此，该实验有两种尺寸的实验单元；行（大尺寸）是品种的实验单元，批（小尺寸）是烹饪处理的实验单元。这种分配如图 \@ref(fig:figure5-12) 所示。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.12.png" alt="将烹饪方法分配给批次的随机化过程的图表（仅显示第1行和第15行的线）" width="974" />
<p class="caption">(\#fig:figure5-12)将烹饪方法分配给批次的随机化过程的图表（仅显示第1行和第15行的线）</p>
</div>

批次实验单元的处理和设计结构是随机完全区组设计结构中的单向处理结构，其中行（或盒子）是区组。表 \@ref(tab:table5-15) 给出了此设计批次部分的方差分析表。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-15)用于烹饪大豆示例的批次分析（忽略品种）的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.15.png">
</td>
  </tr>
</tbody>
</table>

行的平方和包括品种的平方和加上表 \@ref(tab:table5-13) 中 error (row) 的平方和。如果所有行都经相同处理，则烹饪方法与行的交互作用将提供批次误差项，但一些行种植了品种 1，另一些行种植了品种 2，依此类推。因此，烹饪方法与行的交互作用包括了烹饪方法与品种的交互作用。所以，只考虑种植了品种 1 的那些行。品种 1 中的烹饪方法与行的交互作用提供了 (3-1)*(3-1)=4 个自由度，用于测量在同一行中经相同处理的批次的变异。批次误差平方和是通过在五个品种中合并烹饪方法与品种交互作用平方和而获得的 (obtained by pooling the row by cooking method interaction sum of squares within a variety across the five varieties)，从而为 error (batch) 提供 20 个自由度。描述烹饪大豆示例数据的模型如下

$$y_{ijk}=\mu_{ik}+r_{ij}+\varepsilon_{ijk},\quad i=1,2,3,4,5,j=1,2,3,k=1,2,3$$

其中，$\mu_{ik}$ 表示通过第 $k$ 种方法烹饪的第 $i$ 个品种的平均值，$r_{ij}$ 是分配给第 $i$ 个品种的第 $j$ 行的随机效应，假设其分布为 $N(0,\sigma^2_{\mathrm{row}})$， $\varepsilon_{ijk}$ 表示用第 $k$ 种方法烹饪的第 $i$ 个品种的第 $j$ 行的批次的随机效应，假设其分布为 $N(0,\sigma^2_{\mathrm{batch}})$，还假设 $r_{ij}$ 和 $\varepsilon_{ijk}$ 是独立的随机变量。

在效应模型中，均值 $\mu_{ik}$ 可表示为

$$\mu_{ik}=\mu+\upsilon_i+\omega_k+(\upsilon\omega)_{ik}$$

其中 $\mu$ 表示总体平均值，$\upsilon_i$ 表示第 $i$ 个品种的效应，$\omega_k$ 表示第 $k$ 种烹饪方法的效应，$(\upsilon\omega)_{ik}$ 表示品种与烹饪方法的交互作用。

$$\begin{aligned}y_{ijk}=&\mu+\upsilon_i+r_{ij}\,&&\text{\}row part of model}\\&+\omega_k+(\upsilon\omega)_{ik}+\varepsilon_{ijk}\,&&\text{\}batch part of model}\end{aligned}$$

其中模型的行部分也是模型的批次部分的区组结构。裂区设计可以分两步进行分析，先进行行分析，然后进行批次分析。当数据集均衡时，在将单个模型拟合到数据或进行两步分析过程的情况下，将获得相同的结果。表 \@ref(tab:table5-16) 给出了烹饪大豆例子的方差分析表，它划分了每个实验单元尺寸、行尺寸和批次尺寸的分析。行分析也是分析的批次部分的区组结构，如箭头所示。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-16)烹饪大豆实验的方差分析表，显示了每个实验单元尺寸的分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.16.png">
</td>
  </tr>
</tbody>
</table>

大豆实验的整区或行设计结构是完全随机设计结构；因此，这被称为最简单的裂区或两水平设计结构。通常的整区设计结构通常涉及随机完全区组设计结构，如下例所示。

### 示例 5.2：磨小麦——通常的裂区或两水平设计结构 {#sec5-3-2}

谷物碾磨实验包括在将小麦粒碾磨或研磨成小麦粉后评估各种小麦品种 (variety) 的特性。实验包括将研磨辊之间的间隙[称为辊缝 (roll     gap)]设置为一个值，然后研磨每种品种的一批（按随机顺序）。接下来，将辊缝更改为另一个值，并研磨新的批次。假设研究人员希望评估三个辊缝和五个品种。因此，五个品种的三个辊缝的双向处理结构的一次重复需要 15 次轧机的运行 (runs)。一个工作日 (day) 可以完成一次重复或 15 次运行，因此需要四天时间才能获得四次重复。环境条件如湿度可能会对研磨过程产生影响，并且这些条件可能每天都不同。因此，为了帮助控制这些条件，将天用作区组因素，其中在一天内获得 15 种处理组合的一次重复。随机化过程是将三个辊隙的顺序随机分配给每天运行的一个辊缝，然后在每个辊缝内，将品种的顺序随机指定给要研磨的批次。图 \@ref(fig:figure5-13) 显示了随机化过程，其中一天对应于一个区组，一天内的一组五次运行是整区实验单元，单次运行是子区实验单元。未显示所有将品种分配给运行的箭头，因为完整的一组线会使展示混乱不堪。整区设计结构是一个随机完全区组，有四个尺寸为三的区组。整区模型可以表示为

$$y_{ij}=\mu+R_i+d_j+e_{ij},\quad i=1,2,3,j=1,2,3,4,\mathrm{~where~}d_j\thicksim i.i.d.\mathrm{~N}(0,\sigma_{\mathrm{day}}^2)\mathrm{~and~}e_{ij}\thicksim i.i.d.\mathrm{~N}(0,\sigma_{\mathrm{run}}^2)$$

其中 $R_i$ 表示辊缝效应，$d_j$ 表示随机天效应 (random day effect)，$e_{ij}$ 表示随机辊缝运行效应或整区误差。整区方差分析表如表 \@ref(tab:table5-17) 所示，其中运行或整区误差是通过天与辊缝交互作用计算的；也就是说，随机完全区组设计结构的误差是处理结构与设计结构的交互作用。分析的下一步是确定批次误差的来源。这是通过仅考虑具有一个辊缝（例如 1.5mm）的品种来实现的。表 \@ref(tab:table5-18) 中显示了用于比较 1.5mm 辊缝下品种的方差分析表。误差项是通过天与品种（或一天内的运行，因为给定的某天只进行一次运行）的交互作用来计算的，这提供了 12 个自由度来测量在一次运行中相同处理的批次如何变异。这一过程是针对其他两个辊缝设置执行的，每个设置都为批次误差提供 12 个自由度。如果这三个方差相等，那么它们可以合并到具有 36 个自由度的批次误差项中。批次误差项可以表示为 variety × day (roll gap)，也可以表示为“跨辊缝水平的品种与天交互作用” (variety by day interaction pooled across the levels of roll gap).

<div class="figure" style="text-align: center">
<img src="figure/figure 5.13.png" alt="磨粉实验随机化过程的显示" width="998" />
<p class="caption">(\#fig:figure5-13)磨粉实验随机化过程的显示</p>
</div>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-17)磨粉实验的整区分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.17.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-18)磨粉实验中辊缝为1.5mm时品种比较的方差分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.18.png">
</td>
  </tr>
</tbody>
</table>

包含运行部分和批次部分的模型可以表示为

$$y_{ijk}=\mu+R_i+d_j+e_{ij}+V_k+(RV)_{ik}+\varepsilon_{ijk},\quad i=1,2,3,\quad j=1,2,3,4,\quad k=1,2,\ldots,5,$$

其中

$$d_j\thicksim i.i.d.\mathrm{~}N(0,\sigma_{\mathrm{day}}^2)\text{, }e_{ij}\thicksim i.i.d.\mathrm{~}N(0,\sigma_{\mathrm{run}}^2)\text{, and }\varepsilon_{ijk}\thicksim i.i.d.\mathrm{~}N(0,\sigma_{\mathrm{batch}}^2)$$

在上述模型中，$R_i$ 表示第 $i$ 个辊缝效应，$V_k$ 表示第 $k$ 个品种效应，$(RV)_{ik}$ 表示辊隙与品种交互作用，$\varepsilon_{ijk}$ 表示随机批次效应。$\mu+R_i+d_j+e_{ij}$ 是模型的整区或运行部分，也是模型的子区（或模型的批次部分）的区组结构。模型的批次部分是 $V_k+(RV)_{ik}+\varepsilon_{ijk}$.

包括运行分析和批次分析的最终方差分析表如表 \@ref(tab:table5-19) 所示。括号和箭头表示，模型的运行部分或整区部分是模型批次部分的区组结构。如前所述，**通常的裂区设计是具有随机完全区组整区设计结构的裂区设计**。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-19)磨粉实验的方差分析表，显示了每个实验单元尺寸的分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.19.png">
</td>
  </tr>
</tbody>
</table>

### 示例 5.3：烘焙面包——具有不完全块设计结构的裂区 {#sec5-3-3}

一位烘焙科学家设计了一项研究，以评估温度对由两种不同配方制成并在三种不同温度下烘焙的面包 (bread) 体积的影响。图 \@ref(fig:figure5-14) 显示了本实验的随机化过程。一天 (day) 只有两个烤箱 (oven) 可用，而且当天只能使用一次。每个烤箱可以容纳每种配方下的一条面包，但在一天中只能观察到三种温度中的两种。研究人员希望对每个温度进行四次重复，这意味着该研究需要在不同的六天进行。表 \@ref(tab:table5-20) 包含六天的温度分配。图 \@ref(fig:figure5-14) 显示了将温度随机分配到烤箱的情况，其中天是区组，烤箱是温度水平的实验单元。烤箱设计是一种不完全区组设计结构中的单向处理结构（温度的水平）。可以用来描述每个烤箱数据的模型是

$$y_{ij}=\mu+T_i+d_j+o_{ij},\quad(i,j)\in\{(1,1),(2,1),(1,2),(3,2),(2,3),(3,3),(1,4),(3,4),(2,5),(3,5),(1,6),(2,6)\}$$

其中

$$d_j\thicksim i.i.d.N(0,\sigma_{\mathrm{day}}^2)\quad\mathrm{and}\quad o_{ij}\thicksim i.i.d.N(0,\sigma_{\mathrm{oven}}^2)$$

下标 $(i,j)$ 在上述指标集中取值，方差分析表见表 \@ref(tab:table5-21). 预期均方是根据III型平方和计算的（见第 \@ref(chap10) 章）。烤箱误差项是根据设计结构与处理结构交互作用    （或天与温度的交互作用）来计算的。如果观察到表 \@ref(tab:table5-20) 中的所有单元格，则烤箱误差将有 (3-1)(6-1)=10 个自由度，但其中六个单元格为空，从而为烤箱误差提供 10-6=4 个自由度。面包设计是一种随机完全区组设计结构中的单向处理结构（配方的水平），其中烤箱是区组。并不是所有区组都是一样的，所以考虑那些温度相同的烤箱。温度为 160°C 的处理结构部分如表 \@ref(tab:table5-22) 所示。表 \@ref(tab:table5-22) 中的设计是具有四个区组或烤箱或天的随机完全区组设计结构中的单向处理结构（两个配方）。根据 160°C 的温度数据，面包误差项有三个自由度，对应于配方与天的交互作用。160°C 下面包分析的方差分析表如表 \@ref(tab:table5-23) 所示。可以为每个温度计算面包误差，如果方差相等，则将三个方差合并，为 error (bread) 提供九个自由度。可以用来表示所有数据的模型是

$$\begin{aligned}
&y_{ijk}=\mu+T_i+d_j+o_{ij}+R_k+(TR)_{ik}+\varepsilon_{ijk} \\
&(i,j)\in\{(1,1),(2,1),(1,2),(3,2),(2,3),(3,3),(1,4),(3,4),(2,5),(3,5),(1,6),(2,6)\},k=1,2
\end{aligned}$$

其中

$$d_j\thicksim N(0,\sigma_\mathrm{day}^2),\quad o_{ij}\thicksim i.i.d.\mathrm{~N}(0,\sigma_\mathrm{oven}^2)\quad\mathrm{and}\quad\varepsilon_{ijk}\thicksim i.i.d.\mathrm{~N}(0,\sigma_\mathrm{loaf}^2)$$

<div class="figure" style="text-align: center">
<img src="figure/figure 5.14.png" alt="一天内烤箱的温度分配和烤箱内面包的配方分配，使用不完全区组整区设计结构" width="940" />
<p class="caption">(\#fig:figure5-14)一天内烤箱的温度分配和烤箱内面包的配方分配，使用不完全区组整区设计结构</p>
</div>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;border-bottom: 0;">
<caption>(\#tab:table5-20)不完全区组整区设计结构中，将温度分配到天</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.20.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> X 表示当天未观测温度</td></tr></tfoot>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-21)不完全区组设计结构的烤箱分析（忽略配方）的方差分析表，预期均方(EMS)使用III型平方和计算</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.21.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-22)160°C 烘焙的配方设计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.22.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-23)不完全区组设计结构的面包分析的方差分析表（忽略温度）</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.23.png">
</td>
  </tr>
</tbody>
</table>

该模型的最终方差分析表见表 \@ref(tab:table5-24)，其中预期均方是根据III型平方和计算的。由于不完全区组整区设计结构，因此使用了III型平方和。该模型的烤箱部分为 $\mu+T_i+d_j+o_{ij}$，这也是该模型面包部分的区组结构。模型的面包部分是 $R_k+(TR)_{ik}+\varepsilon_{ijk}$. 这个例子表明，**任何类型的设计结构都可以用于实验单元的每个水平或每种尺寸**。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-24)采用不完全区组整区设计结构的面包烘焙研究的完整方差分析表，预期均方(EMS)根据III型平方和计算（见第 \@ref(chap10) 章）</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.24.png">
</td>
  </tr>
</tbody>
</table>

### 示例 5.4：展示柜中的肉——复杂裂区或四水平设计 {#sec5-3-4}

一位肉类科学家想研究三个水平的温度（temperature, T）、两个水平的包装类型（packaging, P）、四个水平的光照类型（lightning, L）和四个水平的光照强度（intensity, I）对在肉类冷却器 (cooler) 中储存七天的肉的颜色的影响。实验有六个冷却器，三种温度（1、3和5°C）分别随机分配给两个冷却器，如图 \@ref(fig:figure5-15) 所示。

每个冷却器在一个 4×4 的网格上分为 16 个分区 (partitions)（图 \@ref(fig:figure5-16)）。光照强度通过它们在冷却器表面上方的距离来调节，因此，一列 (column) 中的所有分区都被分配相同的光照强度。四种类型的光照随机分配到每列内的四个分区。最后，将两种包装分配给牛排，并将两种包装放入每个分区中。图 \@ref(fig:figure5-16) 显示了这种冷却器的安排方式。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.15.png" alt="展示柜研究中肉类冷却器的温度分配" width="680" />
<p class="caption">(\#fig:figure5-15)展示柜研究中肉类冷却器的温度分配</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 5.16.png" alt="在展示柜研究中，为每个冷却器中的肉分配列的光照强度、分区的光照类型和半分区的包装类型" width="710" />
<p class="caption">(\#fig:figure5-16)在展示柜研究中，为每个冷却器中的肉分配列的光照强度、分区的光照类型和半分区的包装类型</p>
</div>

在构建适当的分析之前，必须首先正确识别实验单元的不同尺寸或实验水平。温度水平的实验单元是冷却器。冷却器设计是一种完全随机设计结构中的单向处理结构（T 或温度水平）。如果对每个冷却器进行一次测量，则响应可以通过以下方式建模

$$y_{ij}=\mu+T_i+c_{ij},\quad i=1,2,3,j=1,2,\mathrm{~and~}c_{ij}\thicksim i.i.d.N(0,\sigma_{\mathrm{cooler}}^2)$$

表 \@ref(tab:table5-25) 给出了冷却器模型的方差分析表，其中冷却器误差项是根据两个冷却器在一个温度内的变异、跨三种温度计算得出的 (computed from the variation of the two coolers within a temperature, pooled across the three temperatures). 因此，error (cooler) 有三个自由度。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-25)冷却器实验单元部分的模型的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.25.png">
</td>
  </tr>
</tbody>
</table>

光照强度水平的实验单元是冷却器中四个分区所在的列。列设计由单向处理结构（I 或强度水平）组成，采用随机完全区组设计结构，具有六个区组或冷却器。如果所有冷却器都经相同处理，那么列误差项将根据强度与冷却器交互作用来计算。但有三种不同的温度，因此将分析限制在分配给 1°C 的两个冷却器上。分配给 1°C 的两个冷却器如图 \@ref(fig:figure5-17) 所示。这时，设计是一个单向处理结构（I 的四个水平）的随机完全区组设计结构（两个冷却器）。如果对这两个冷却器的每个列进行测量，则可用于描述响应的模型为

$$y_{1jk}=\mu+I_k+c_{1j}+d_{1jk},\quad j=1,2,\,k=1,2,3,4,\,c_{1j}\sim i.i.d. N(0,\sigma_{\mathrm{cooler}}^2)\,\mathrm{~and~}\,d_{1jk}\sim i.i.d.N(0,\sigma_{\mathrm{colum}n}^2)$$

列模型的方差分析表如图 \@ref(fig:figure5-17) 所示，其中列误差项是通过强度与冷却器交互作用来计算的，提供三个自由度。对于另外两个温度重复该过程，每个温度都提供具有三个自由度的列误差项。如果这三个方差相等，则可以将它们合并到具有九个自由度的 error (column) 中。error (column) 可以用 intensity × cooler (temperature) 表示，读作“跨温度水平的强度与冷却器交互作用” (intensity by cooler interaction pooled across the levels of temperature).

<div class="figure" style="text-align: center">
<img src="figure/figure 5.17.png" alt="展示柜研究中，比较温度 1 下的强度水平以开始列分析" width="849" />
<p class="caption">(\#fig:figure5-17)展示柜研究中，比较温度 1 下的强度水平以开始列分析</p>
</div>

光照类型水平的实验单元是列的分区。分区设计是一种单向处理结构（L 的水平），采用随机完全区组设计结构，有 24 个区组，由六个冷却器的四个列组成。如果所有列都经相同处理，则分区误差项将通过光照类型与列交互作用来计算。但所有列的处理方式并不相同，因为有三种温度和四种强度水平。将结构限制为仅涉及温度为 1°C 且强度为 I~1~ 的那些列，如图 \@ref(fig:figure5-18) 所示。与图 \@ref(fig:figure5-18) 相关的设计是随机完全区组设计结构（两个列，但每个列来自不同的冷却器）中的单向处理结构（光照类型的四个水平）。如果对图 \@ref(fig:figure5-18) 中的每个分区进行一次测量，则可用于描述响应的模型为

$$\begin{aligned}y_{1j1m}&=\mu+L_m+d_{1j1}^*+p_{1j1m},\quad j=1,2,m=1,2,3,4,d_{1j1}^*\thicksim i.i.d.N(0,\sigma_{\mathrm{column}}^{2^*})\mathrm{~and}\\p_{1j1m}&\sim i.i.d.N(0,\sigma_{\mathrm{partition}}^2)\end{aligned}$$


其中 $d_{1j1}^*$ 表示冷却器和该冷却器效果的列的组合 (the combination of cooler and column of that cooler effect). 该分区模型的方差分析表如图 \@ref(fig:figure5-18) 所示，其中分区误差项通过光照类型与冷却器交互作用来计算，提供三个自由度。这个过程需要对 12 种温度和强度的组合进行，每种组合为 error (partition) 提供三个自由度。如果这 12 个方差相等，则可以将它们合并为一个项，为 error (partition) 提供三个自由度提供 36 个自由度。error (partition) 项可以表示为 light × intensity × cooler (temperature).

<div class="figure" style="text-align: center">
<img src="figure/figure 5.18.png" alt="展示柜研究中，比较温度 1 和强度 1 下的光照水平以开始分区分析" width="811" />
<p class="caption">(\#fig:figure5-18)展示柜研究中，比较温度 1 和强度 1 下的光照水平以开始分区分析</p>
</div>

最后，包装水平的实验单元是半分区 (half-partitions)（或牛排）。半分区设计是一种单向处理结构（包装的水平或 P），采用随机完全区组设计结构，具有 96 个分区，由六个冷却器的四个列中的四个分区组成。如果所有分区都经相同处理，那么半分区误差项将通过包装与分区交互作用来计算。但并不是所有的分区都经相同处理，因为它们被分配到三个温度水平、四个强度水平和四个光照水平。选择分配给温度为 1°C、强度为 I~1~ 和光照为 L~1~ 的分区，如图 \@ref(fig:figure5-19) 所示。如果在图 \@ref(fig:figure5-19) 中对每块牛排或半分区进行一次测量，则可用于描述响应的模型为

$$\begin{aligned}y_{1j11n}&=\mu+P_n+p_{1j1n}^*+\varepsilon_{1j11n},\quad j=1,2,n=1,2,p_{1j1m}^*\thicksim i.i.d.N(0,\sigma_{\mathrm{partition}}^{2*})\quad\mathrm{and}\\\varepsilon_{1j11n}&\sim i.i.d.N(0,\sigma_{\frac12\text{partition}} ^ 2 ) \end{aligned}$$

该半分区模型的方差分析表如图 \@ref(fig:figure5-19) 所示，其中半分区误差项通过包装与冷却器交互作用来计算，提供一个自由度。这个过程需要对 48 种温度、强度和光照的组合进行。这 48 个平方和合并在一起为 error (half-partition) 提供 48 个自由度。error (half-partition) 项可以表示为 packaging × light × intensity × cooler (temperature).

<div class="figure" style="text-align: center">
<img src="figure/figure 5.19.png" alt="展示柜研究中，比较温度 1、强度 1 和光照 1 下的包装水平以开始半分区分析" width="716" />
<p class="caption">(\#fig:figure5-19)展示柜研究中，比较温度 1、强度 1 和光照 1 下的包装水平以开始半分区分析</p>
</div>

上述讨论提供了研究中的四个水平的模型：冷却器模型、列模型、分区模型和半分区模型。这些模型可以组合成一个模型，其中因素之间的交互作用被加在一起。由于基本设计结构是每个水平的裂区，因此每个水平的处理结构以及与其上一水平处理结构的交互作用都包含在模型的该水平部分中。可以用来描述数据的模型是

$$\begin{aligned}
\mathcal{Y}_{ijkmn}=& \mu+T_i+\mathcal{C}_{ij}  &&\{\text{cooler part of the model}\\
&+I_k+(TI)_{ik}+d_{ijk} &&\{\text{column part of the model}\\
&+L_m+(TL)_{im}+(IL)_{km}+(TIL)_{ikm}+p_{ijkm} &&\{\text{partition part of the model}\\
&+P_n+(TP)_{in}+(IP)_{kn}+(TIP)_{ikn}+(LP)_{mn}+(TLP)_{imn}+(ILP)_{km}+(TILP)_{ikmn}+\varepsilon_{ijkmn} &&\{\text{half-partition part of the model}\\
&i=1,2,3,\quad j=1,2,\quad k=1,2,3,4,\quad m=1,2,3,4,\quad n=1,2,\\
&c_{ij}\thicksim i.i.d.N(0,\sigma_{\mathrm{cooler}}^2),\quad d_{ijk}\thicksim i.i.d.N(0,\sigma_{\mathrm{column}}^2),\quad  p_{ijkn}\thicksim i.i.d.N(0,\sigma_{\mathrm{partition}}^2),\quad \boldsymbol{\varepsilon}_{ijkmn}\thicksim i.i.d.N(0,\sigma_{\frac12\text{partition}} ^ 2 )
\end{aligned}$$

表 \@ref(tab:table5-26) 给出了上述模型的方差分析表，其中较大的实验单元分析也是下一个较小实验单元的区组结构。这是一个裂-裂-裂区 (split-split-split-plot) 实验模型，由四个水平（或实验单元尺寸）组成，涉及四个误差项，每个水平（或每种实验单元尺寸）对应一个。这也是一种分层设计结构，其中半分区嵌套在分区内，分区嵌套在列内，列嵌套在冷却器内。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;border-bottom: 0;">
<caption>(\#tab:table5-26)展示柜研究中复杂裂区设计完整的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.26.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> 括号和箭头指示哪些效应形成下一个较小尺寸的实验单元的区组</td></tr></tfoot>
</table>

## 裂条区设计结构：一种无层次的多水平设计 {#sec5-4}

构造裂条区设计结构的过程是将实验单元排列成矩形。一组因素的水平随机分配给每个矩形的行，另一组因素的水平随机分配给每个矩形的列。因此，行是与第一组因素相关的实验单元，列是与第二组因素相关的实验单元。但因此，一行和一列的单元格或交叉点是与两组因素之间的交互作用比较相关的实验单元。通过一个例子来演示裂条区设计结构的一些用途。

### 示例 5.5：制作奶酪 {#sec5-4-1}

一位乳制品科学家设计了一项奶酪制造研究，涉及两种脂肪水平、三种奶酪类型、两种储存温度和两种储存湿度。该实验是一个分两步进行的过程，第一步是用两种脂肪水平为每种类型的奶酪制作一批 (batch) 奶酪。每批奶酪要制成四包一磅重的奶酪块 (cheese block). 第二步是将奶酪储存在一系列环境条件下，用四周的时间使其固化。该研究的储存部分包括将每批奶酪中的一包放入一个指定储存温度和湿度的室 (chamber) 中。这个过程相当于将 24 包奶酪排列成四行六列的矩形。该列对应一批奶酪，该行对应一个环境室，如图 \@ref(fig:figure5-20) 所示。这位奶制品科学家有四个实验室可用于研究，但他希望有四次重复。因此，通过在四个不同的月份进行实验获得了四次重复。可以通过评估每种尺寸的实验单元的处理和设计结构来构建模型并分析。室设计如图 \@ref(fig:figure5-21) 所示，其中仅显示了四个月中每个月的温度和湿度水平。室设计是一种双向处理结构（温度水平与湿度水平），采用随机完全区组设计结构，每个月为一个区组。可用于描述在每个室中进行的测量的模型为

$$\begin{aligned}y_{ijk}=&\mu+T_j+H_k+(TH)_{jk}+m_i+c_{ijk},\quad i=1,2,3,4,j=1,2,k=1,2\\&m_i\sim i.i.d.N(0,\sigma_{\mathrm{month}}^2),\text{ and }c_{ijk}\thicksim i.i.d.N(0,\sigma_{\mathrm{chamber}}^2)\end{aligned}$$

表 \@ref(tab:table5-27) 包含室模型的方差分析表，其中室误差项是处理结构与设计结构的交互作用。在处理结构中有四种处理组合，在设计结构中有四个区组；因此，室误差项基于九个自由度。error (chamber) 项是温度与月交互作用、湿度与月交互作用和温度与湿度与月交互作用项的总和。如果在 SAS-Mixed 程序中使用三向交互项，则三个交互项将合并为一个误差项。室的方差分量用星号表示，因为它也包括由于奶酪块引起的变异。月份的方差分量用星号表示，因为它还包括由于批次和奶酪区组引起的变异。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.20.png" alt="采用裂条区设计结构的奶酪制作实验的随机化方案" width="1006" />
<p class="caption">(\#fig:figure5-20)采用裂条区设计结构的奶酪制作实验的随机化方案</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 5.21.png" alt="奶酪制作实验中室设计的随机化方案" width="855" />
<p class="caption">(\#fig:figure5-21)奶酪制作实验中室设计的随机化方案</p>
</div>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-27)奶酪制作实验中室实验单元设计的方差分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.27.png">
</td>
  </tr>
</tbody>
</table>

批次设计如图 \@ref(fig:figure5-22) 所示，其中仅考虑了奶酪类型和脂肪百分比。批次设计是一种双向处理结构（奶酪类型水平与脂肪水平），采用随机完全区组设计结构，其中每个月为一个区组。可用于描述对每个批次进行测量的模型为

$$\begin{aligned}y_{imn}=&\mu+F_m+C_n+(FC)_{mn}+m_i+b_{imn},\quad i=1,2,3,4,m=1,2,n=1,2,3,\\&m_i\thicksim i.i.d.N(0,\sigma_{\mathrm{month}}^2)\mathrm{~and~}b_{imn}\thicksim i.i.d.N(0,\sigma_{\mathrm{batch}}^2)\end{aligned}$$

表 \@ref(tab:table5-28) 包含室模型的方差分析表，其中室误差项是处理结构与设计结构的交互作用。处理结构中有六种处理组合，设计结构中有四个区组，因此批次误差项基于 15 个自由度。error (batch) 项是月份与脂肪、月份与奶酪类型以及脂肪与奶酪类型与月份交互作用项的总和。如果在 SAS-Mixed 程序中使用三向交互项，则三个交互项将合并为一个误差项。批次的方差分量用星号表示，因为它还包括由于奶酪块引起的变异。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.22.png" alt="奶酪制作实验的批设计的随机化方案。为了简化绘图，有些线没有包括在内" width="924" />
<p class="caption">(\#fig:figure5-22)奶酪制作实验的批设计的随机化方案。为了简化绘图，有些线没有包括在内</p>
</div>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-28)奶酪制作实验中批次实验单元设计的方差分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.28.png">
</td>
  </tr>
</tbody>
</table>

1 磅奶酪块是分配给室的因素和分配给批次的因素之间交互作用的实验单元。模型的奶酪块部分是

$$\begin{aligned}&(TF)_{jm}+(TC)_{jn}+(TFC)_{jmn}+(HF)_{kn}+(HC)_{kn}+(HFC)_{kmn}+(THFC)_{jkmn}+\varepsilon_{ijkmn}\\&\varepsilon_{ijkmn}\thicksim i.i.d.N(0,\sigma_{\mathrm{block}}^2)\end{aligned}$$

该裂条区设计的模型是通过将三个实验单元的模型合并为一个模型而获得的：

$$\begin{aligned}y_{ijk}&=\mu+m_i&&\{\text{blocking part of the model}\\&+T_j+H_k+(TH)_{jk}+c_{ijk}&&\{\text{chamber part of the model}\\&+F_m+C_n+(FC)_{mn}+b_{imn}&&\{\text{batch part of the model}\\&+(TF)_{jm}+(TC)_{jn}+(TFC)_{jmn}+(HF)_{km}+(HC)_{kn}+(HFC)_{kmn}+(THFC)_{jkmn}+\varepsilon_{ijkmn}&&\{\text{cheese block part of the model}\\
&i=1,2,3,4,\quad j=1,2,\quad k=1,2,\quad m=1,2,\quad n=1,2,3 \\
&m_i\thicksim i.i.d.N(0,\sigma_{\mathrm{month}}^2),\quad c_{ijk}\thicksim i.i.d.N(0,\sigma_{\mathrm{chamber}}^2),\quad b_{imn}\thicksim i.i.d.N(0,\sigma_{\mathrm{batch}}^2), \\
&\mathrm{and~}\varepsilon_{ijkmn}\sim i.i.d.N(0,\sigma_{\mathrm{block}}^2)\end{aligned}$$

表 \@ref(tab:table5-29) 给出了奶酪实验的方差分析表，其中各行按实验单元的尺寸进行分离。奶酪块误差项通过室处理结构与批次处理结构与设计结构的交互作用来计算，提供 (4 - 1)(6 - 1)(4 - 1) = 45 个自由度。该设计涉及三种尺寸的实验单元：室尺寸、批次尺寸和奶酪块尺寸。奶酪块嵌套在批次内，奶酪块嵌套在室内，但批次不嵌套在室内，室也不嵌套在批次内。因此，**这种多层设计结构不是分层设计结构**。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-29)奶酪制作实验完整的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.29.png">
</td>
  </tr>
</tbody>
</table>

这个例子说明了与给定尺寸的实验单元相关的设计结构可以涉及任何所需的处理结构。室和批次实验单元的处理结构都是双向析因安排。
我们可以采用这样的双向析因安排：将一个对照作为其中一个实验单元的处理结构。**这种先确定实验单元，然后指定每个实验单元的设计和处理结构的过程，为复杂实验提供了确定适当的设计和相应模型的一般方法**。

## 重复测量设计 {#sec5-5}

重复测量设计在许多研究领域中得到了有效的应用。这些设计包括至少两个因素的处理结构，不完全区组设计结构和至少两种实验单元尺寸。重复测量设计与裂区设计具有相同的设计结构类型，但重复测量设计与裂区类型设计的不同之处在于，一个或多个因素的水平不能随机分配给相应的实验单元。通常，时间是这样一个因素：其水平不能随机分配给个体的时间间隔。涉及时间的重复测量设计是用于纵向研究的设计。因此，重复测量设计涉及一个或多个步骤，其中不可能将某些因素的水平随机分配给其实验单元，而裂区型设计则可以在每一步使用随机化。下面的例子演示了重复测量设计的一些结构，并为正确识别设计实验提供了指导。

### 示例 5.6：马足——基本重复测量设计 {#sec5-5-1}

这个特殊的实验引出了使用两种不同尺寸的实验单元的重复测量设计的想法。然而，如果实验者不小心，他或她可能会无意中错过有两种不同尺寸的实验单元的事实，从而导致错误的分析。

兽医有两种技术可用于融合马足 (horse feet) 骨折后的关节，她希望确定哪种技术更好。 这个实验包括取一些马，在每匹马身上折断一个关节，用两种技术中的一种进行修复，四个月后确定融合关节的强度。 她还想确定同样的技术对前足和后足是否同样有效。 由于马匹稀缺且价格昂贵，她计划在每匹马的前足上折断关节，等待其愈合，然后再在后足上折断关节，或者反过来，先折断后足后折断前足。这个愈合过程也在设计中引入了时间因素。因此，处理结构是由两种融合技术（F）、两个位置（P）和两个愈合时间（T）产生的 2^3^ 析因安排。设计结构是一个不完全区组设计，每匹马都是一个区组，每个区块（马）有两个观测值。由于区组不完全，一些处理结构信息将与区组（或马）效应混淆 (Cochran and Cox, 1957). 有多种方法可以将处理组合分配给马的双足。假设有四匹马。一个过程是为每匹马指定一种融合技术，其中前足和后足都接受相同水平的融合处理，如表 \@ref(tab:table5-30) 所示。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;border-bottom: 0;">
<caption>(\#tab:table5-30)马足实验处理组合的首次分配</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.30.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> 两种融合技术为 $F_1$ 和 $F_2$，两个时间为 $T_1$ 和 $T_2$，两个位置为 $P_1$ 和 $P_2$</td></tr></tfoot>
</table>

设 $y_{ijkm}$ 表示通过第 $i$ 种融合技术、第 $j$ 个位置、第 $k$ 个时间和第 $m$ 匹马获得的观测值。用于描述数据的模型是

$$y_{ijkm}=\mu_{ijk}+h_m+\varepsilon_{ijkm},\quad i=1,2,\quad j=1,2,k=1,2,m=1,2,\ldots,4$$

其中 $\mu_{ijk}$ 表示第 $i$ 种融合技术、第 $j$ 个位置、第 $k$ 个时间的平均响应，$h_m$ 表示第 $m$ 匹马的效应，$\varepsilon_{ijkm}$ 表示在马足上进行的测量的响应误差。在处理组合均值 $\mu_{ijk}$ 中可以进行两种类型的比较。一种类型的比较是马内 (intra-horse or within-horse) 比较，而另一种类型的比较是马间 (inter-horse or between-horse) 比较。析因效应可以定义为

$$\begin{aligned}
\text{Mean}& =\bar{\mu}_{\cdot\cdot\cdot}  \\
F& =\bar{\mu}_{1\cdot\cdot}-\bar{\mu}_{2\cdot\cdot},P=\bar{\mu}_{\cdot1\cdot}-\bar{\mu}_{\cdot2\cdot},T=\bar{\mu}_{\cdot\cdot1}-\bar{\mu}_{\cdot\cdot2}  \\
F\times P& \begin{aligned}&=\bar{\mu}_{11\cdot}-\bar{\mu}_{12\cdot}-\bar{\mu}_{21\cdot}+\bar{\mu}_{22\cdot}\end{aligned}  \\
F\times T& =\bar{\mu}_{1\cdot1}-\bar{\mu}_{1\cdot2}-\bar{\mu}_{2\cdot1}+\bar{\mu}_{2\cdot2}  \\
P\times T& =\bar{\mu}_{\cdot11}-\bar{\mu}_{\cdot12}-\bar{\mu}_{\cdot21}+\bar{\mu}_{\cdot22}  \\
F\times P\times T& =\mu_{111}-\mu_{112}-\mu_{121}+\mu_{122}-\mu_{211}+\mu_{212}+\mu_{221}-\mu_{222} 
\end{aligned}$$

$\mu_{ijk}$ 的最佳估计为 $\bar y_{ijk\cdot}$. 主效应 $P$ 的最佳估计可以表示为

$$\begin{aligned}\hat{P}&=\bar{y}_{\cdot1\cdot\cdot}-\bar{y}_{\cdot2\cdot\cdot}\\&=\frac{1}{4}[({y}_{1111}-{y}_{1221})+({y}_{1122}-{y}_{1212})+({y}_{2113}-{y}_{2223})+({y}_{2124}-{y}_{2214})]\end{aligned}$$

这是马内比较。将上述模型的右侧替换为 $P$ 中的 ${y}_{ijkm}$，可以很容易地证明马内比较的合理性，这给出

$$\begin{aligned}\hat{P}&=\frac14[(\mu_{111}+h_1+\varepsilon_{1111}-\mu_{122}-h_1-\varepsilon_{1221})+(\mu_{112}+h_2+\varepsilon_{1122}-\mu_{121}-h_2-\varepsilon_{1212})\\&+(\mu_{211}+h_3+\varepsilon_{2113}-\mu_{22}-h_3-\varepsilon_{223})+(\mu_{212}+h_4+\varepsilon_{2124}-\mu_{221}-h_4-\varepsilon_{2214})]\\&=\bar{\mu}_{\cdot1}-\bar{\mu}_{\cdot2}+\frac{1}{4}[(\varepsilon_{1111}-\varepsilon_{1221})+(\varepsilon_{1122}-\varepsilon_{1212})+(\varepsilon_{2113}-\varepsilon_{222})+(\varepsilon_{2124}-{\varepsilon}_{2214})]\end{aligned}$$

请注意到马效应 $h_m$ 从表达式中消去。$\hat P$ 的方差取决于 $\varepsilon_{ijkm}$ 的方差，而不是 $h_m$ 的方差。 $\hat P$ 的方差为 $\operatorname{Var}(P) = \sigma^2_{\varepsilon} /2$.

$F$ 的最佳估计可以表示为

$$\begin{aligned}\hat{F}&=\bar{y}_{1\cdot\cdot\cdot}-\bar{y}_{2\cdot\cdot\cdot}\\&=\frac14[(y_{1111}+y_{1221})+(y_{1122}+y_{1212})-(y_{2113}+y_{2223})-(y_{2124}+y_{2214})]\end{aligned}$$

此估计是马间的比较；也就是说，它是马 1 和马 2 平均值与马 3 和马 4 平均值的比较，因此它取决于马效应。通过将模型右侧的 $F$ 表示为

$$\begin{aligned}
\hat{F}& =\frac14[(\mu_{111}+h_1+\varepsilon_{1111}+\mu_{122}+h_1+\varepsilon_{1221})+(\mu_{112}+h_2+\varepsilon_{1122}+\mu_{121}+h_2+\varepsilon_{1212})  \\
&=[\bar{\mu}_{\cdot1}-\bar{\mu}_{\cdot2}+\frac14[(\varepsilon_{1111}+\varepsilon_{1221})+(\varepsilon_{1122}+\varepsilon_{1212})-(\varepsilon_{2113}+\varepsilon_{2223})-(\varepsilon_{2124}+\varepsilon_{2214})] \\
&+\frac12[h_1+h_2-h_3-h_4]
\end{aligned}$$

其中涉及 $h_m$，从而此估计依赖于马效应。$\hat F$ 的方差取决于 $h_m$ 的方差和 $\varepsilon_{ijkm}$ 的方差，$\mathrm{Var}(\hat{F})=\frac{1}{2}\left[\sigma_{\varepsilon}^{2}+2\sigma_{\mathrm{horse}}^{2}\right]$。同样可以证明，$P,T,F\times T,F\times P$ 是马内效应，而 $F,P\times T,F\times P\times T$ 是马间效应。由于马间效应涉及 $h_m$，因此称它们与马效应混淆 (confounded)；即 $F,P\times T,F\times P\times T$ 与马效应混淆。

该设计由不完全区组设计结构中的三向处理结构组成，其中每匹马是一个区组。由于一些比较是马间比较，另一些是马内（足间）比较，因此该实验涉及两种尺寸的实验单元。足是较小尺寸的实验单元，而马是较大尺寸的实验单元。模型中的 $h_m$ 项代表马误差（即，经相同处理的马之间的差异引起的变异），而 $\varepsilon_{ijkm}$ 项表示足误差（即，同一匹马经相同处理的足之间的差异引起的变异）。当设计涉及不同尺寸的实验单元的实验时，最好选择那些最重要的效应，以便它们涉及较小实验单元之间的比较，并让那些最不重要的效应涉及较大实验单元之间的比较。然而，这类安排并不总是可行的。例如，如果马足实验涉及两种类型的马（例如赛马和工作马），那么马的类型就只能作为马间比较。在只有三个因素的马足实验中，实验者最感兴趣的是比较两种融合技术。表 \@ref(tab:table5-30) 中给出的设计使得融合效应与马效应混淆，导致融合技术均值比较的精度低于预期。表 \@ref(tab:table5-31) 中给出的设计实现了融合技术均值的比较，且不与马效应混淆，因此达到了使用较小实验单元（马足）的变异性来比较最重要效应的目标。

第一种安排的模型也可以用来表示第二种将处理组合分配给马或区组的数据。使用与第一种安排相同的技术，可以证明 $F,P,T,F\times P\times T$ 效应是马内比较，$F\times T,F\times P,P\times T$ 效应是马间比较，这些效应与马效应混淆。上述两种设计都没有产生足够的观测值来为估计两个误差项提供任何自由度。为了获得两种类型的误差方差的一定自由度，可以使用八匹马重复表 \@ref(tab:table5-31) 中的设计，其中两匹马被随机分配到每组处理组合。分析仍然由两部分组成，马间分析和马内（足间）分析。将进行八次马内比较（每匹马一次），可划分为对 $F,T,P,F \times P \times T$ 效应的估计，以及对马内误差方差的估计，该误差方差用 error (feet) 表示。将进行七此马间比较，可划分为 $F\times T,F\times P , P\times T$ 效应的估计，以及马间误差方差的估计，该误差方差用 error (horse) 表示。方差分析结果表显示在表 \@ref(tab:table5-32) 中。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-31)马足实验处理组合的第二次分配</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.31.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-32)马间和马内方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.32.png">
</td>
  </tr>
</tbody>
</table>

这个设计的实验属于重复测量设计类型，因为每匹马都有两次测量（重复）；也就是说，对每匹马测量前足和后足，并且不能将位置水平随机分配给其实验单元。

### 示例 5.7：舒适度研究——重复测量设计 {#sec5-5-2}

一位实验者想研究六种环境条件对人的舒适度的影响。他有六个环境室，每个环境室都可以设置不同的环境条件。该实验包括将一个人放在一个室内，然后在暴露 1、2 和 3 小时后测量他的舒适度。研究中有 36 名个体，其中 6 名个体被随机分配到每种环境中。研究人员可以在一天内从六种环境条件中获得数据，因此天被用作一个区组因素，每天将每种环境条件分配给六个室中的一个。暴露时间是研究中的一个重要因素，因此处理结构是双向的，环境水平与三个暴露时间交叉。受试者被随机分配一个从 1 到 36 的数字，其中前 6 人在第一天参与实验，以此类推。环境和人员的分配如图 \@ref(fig:figure5-23) 所示。图 \@ref(fig:figure5-23) 中的每个矩形代表一个人，T~1~、T~2~和T~3~分别代表暴露 1、2 和 3 小时后的舒适度。环境的实验单元是一个人或一个室，而时间的实验单元则是该人员“内”的1小时间隔。实际上，人被“分成”三部分，但这种设计结构被称为重复测量设计而不是裂区设计，因为暴露时间水平不能随机分配给一个人的三个 1 小时暴露时间。该设计的结构与随机完全区组整区设计结构中具有双向处理结构的通常裂区相同。对于这种重复测量设计结构，较大尺寸的实验单元（人，或整区）设计由随机完全区组设计结构（六天）中的单向处理结构（六个环境水平）组成。较小尺寸的实验单元（1 小时时间间隔，或裂区）设计是随机完全区组设计结构中的单向处理结构（三次暴露），其中每个人或室是一个区组。可以用来描述该实验数据的模型是

$$y_{ijk}=\mu+E_i+d_j+c_{ij}+T_k+(ET)_{ik}+\varepsilon_{ijk},\quad i=1,2,\ldots,6,\quad j=1,2,\ldots,6,\quad k=1,2,3$$

其中 $d_j$ 表示第 $j$ 天的效应，$c_{ij}$ 表示第 $j$ 天用于第 $i$ 个环境的室（或人）的效应，$\varepsilon_{ijk}$ 是时间间隔效应。如果对每人进行一次观测，则所得数据集将来自随机完全区组设计结构中的单向处理结构。室（或人）的误差项是通过天与环境交互作用计算的。如果比较第一种环境条件下的时间，则所得设计是随机完全区组设计结构中的单向处理结构，其中人（或室，或天）是区组因素。时间间隔误差项是根据天与时间交互作用计算的。可以为每个环境计算时间间隔的误差项，如果方差相等，则可以将其合并到时间间隔误差中。表 \@ref(tab:table5-33) 中给出了该重复测量设计的方差分析表。这种基本的重复测量结构与通常的裂区设计相同。
裂区模型假设 $d_j,c_{ij},\varepsilon_{ijk}$ 是独立的随机变量，方差为 $\sigma^2_\text{day},\sigma^2_\text{chamber},\sigma^2_\varepsilon$. 重复测量设计也假设 $d_j,c_{ij},\varepsilon_{ijk}$ 是独立的随机变量，但由于 $\varepsilon_{ijk}$ 是在同一个人身上测量的，它们可以是相关的，相关结构可能需要建模，协方差结构的建模将在第 \@ref(chap26) 章和第 \@ref(chap27) 章讨论。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.23.png" alt="舒适度研究的数据安排" width="927" />
<p class="caption">(\#fig:figure5-23)舒适度研究的数据安排</p>
</div>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-33)舒适度研究的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.33.png">
</td>
  </tr>
</tbody>
</table>

将一个人分配到每个环境室的另一种方法是将六人一组分配到一个环境室并在一天内完成研究。这个分配过程似乎是一个好主意，但生成的数据集没有提供任何这样的测量：经相同处理时，室是如何变异的。室误差项不存在是因为六个人组成了环境的实验单元，因此每个环境只有一个独立的观测。如果没有室误差项，就无法评估环境对舒适度的效应。

### 示例 5.8：交叉或转换设计

一种将处理应用于个体（或对象）的比较方法，就是让个体自身作为对照，将两种或多种处理应用于同一个体。如果有两种处理方法，那么过程是先将处理 A 应用于第一个时间段并测量响应，让处理 A 的效应减弱 (diminish) 或洗脱 (wash out)，然后在第二个时间段将处理 B 应用于个体并观察对处理 B 的响应。随机化过程是构建处理序列 (sequence)（这里是先 A 后 B 和先 B 后 A），然后随机将个体或动物分配给这些序列。这种方法还可以用于植物、地块或其他对象，在将实验单元暴露于处理后，该实验单元应至少能够从第一次处理的效应中恢复一部分。

在这种比较两种处理方法的方法中，动物有两种处理分配序列——先 A 后 B 和先 B 后 A。这两个序列通常用 AB 序列和 BA 序列表示。处理结构是具有两个水平（A 和 B）的单向处理集，但由于处理按顺序应用于每个实验单元，因此生成的序列成为另一种类型的处理。因此，本实验的设计涉及双向处理结构，处理与序列交叉。

这种交叉设计是一种具有两种尺寸的实验单元的重复测量设计。将处理序列分配给个体，使个体成为更大的实验单元。观察处理的时间段或时间是较小尺寸的实验单元。较大尺寸实验单元的设计是完全随机设计结构中具有两个水平的序列（两个可能的序列）的单向处理结构。尽管可以使用任何设计结构，但最常用的是完全随机设计。较小尺寸实验单元（时间间隔）的设计是完全随机区组设计结构中的单向处理结构（两个时间段或两次处理），个体是区组。数据可按表 \@ref(tab:table5-34) 所示进行排列，可用于描述数据的模型为

$$y_{ijk}=\mu_{ik}+s_{ij}+\varepsilon_{ijk},\quad\mathrm{~where~}(i,j,k)\in I_\text{D}{ , s _ { i j }}\thicksim i.i.d.\text{ N}(0,\sigma_{\mathrm{subject}}^2),\mathrm{~and~}\varepsilon_{ijk}\thicksim i.i.d.\text{ N}(0,\sigma_{\mathrm{time}}^2)$$

其中 $I_D$ 是实验中实际使用的 $(i,j,k)$ 三元组的指标集，$\mu_{ik}$ 表示第 $i$个序列中第$k$次处理的效应，$s_{ij}$表示分配给第$i$个序列的第$j$个受试者的随机效应，$\varepsilon_{ijk}$表示在应用第$k$次处理的第$i$个序列的时间段中测量的随机误差。第 \@ref(chap29) 章介绍了一些交叉设计的分析。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-34)两阶段两处理交叉设计的数据安排</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.34.png">
</td>
  </tr>
</tbody>
</table>

## 涉及嵌套因素的设计 {#sec5-6}

在给定的设计中，可能在设计结构、处理结构或两个结构中都具有**嵌套效应** (nested effects). 嵌套最常见于实验的设计结构中，其中较小尺寸的实验单元嵌套在较大尺寸的实验中。如果每个较大实验单元中的较小实验单元都不同，那么一个尺寸的实验单元就嵌套在另一个较大的尺寸中。当设计结构由几种不同尺寸的实验单元组成，并且存在一种排序：最小尺寸嵌套在次小尺寸内，直到次大尺寸嵌套在最大尺寸的实验单元内，这种设计也称为分层设计结构。本章前几节讨论的裂区设计和重复测量设计是这样的设计结构的良好示例：其中较小尺寸实验单元嵌套在较大尺寸实验单元内。当一个因素的水平仅与第二个因素的一个水平一起出现时，处理结构中就会出现嵌套。在这种情况下，第一个因素的水平被认为是嵌套在第二个因素的水平内。表 \@ref(tab:table5-35) 说明了示例 [4.2](#sec4-3-2) 第二部分的设计结构中的嵌套，其中房屋嵌套在方块内（“X”表示房屋包含在相应的方块内）。每个方块是较大的实验单元或房屋区组，而房屋是较小的实验单元。由于第一个方块的房屋与第二个方块的房屋不同，因此房屋嵌套在方块内。这种嵌套效应通常在模型中用 $s_k + h_{m(k)}$ 表示，其中 $s_k$ 表示第 $k$ 个方块的效应，$h_{m(k)}$ 表示嵌套在第 $k$ 个方块中的第 $m$ 个房屋的效应。方块和方块内房屋的平方和分别用 SSSQUARE 和 SSHOUSES(SQUARES) 表示。如果模型中包含房屋而不是方块，则只有一个由房屋（SSHOUSES）引起的平方和，可以将其划分为 SSHOUSES = SSSQUARES + SSHOUSES(SQUARES)。此过程可以再进行一步，其中房屋的侧面嵌套房屋内（而房屋又嵌套在方块内）。这样，这三种尺寸的实验单元就形成了一个分层设计结构。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-35)显示嵌套在方块内的房屋的设计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.35.png">
</td>
  </tr>
</tbody>
</table>

### 示例 5.9：动物遗传学 {#sec5-6-1}

一位动物科学家想要研究基因对羔羊生长速度的影响。她有 4 只公羊和 12 只母羊。育种结构如表 \@ref(tab:table5-36) 所示（“X”表示交配）。在这个例子中，每只公羊都与三只母羊交配，而这三只母羊对于每只公羊来说都是不同的。因此，母羊的水平被称为嵌套因素，母羊嵌套在公羊中。如果动物科学家对这个公羊组合和这个母羊组合的效应感兴趣，那么嵌套在公羊水平内的母羊水平就构成了一个嵌套处理结构。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-36)育种结构显示母羊嵌套在公羊中</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.36.png">
</td>
  </tr>
</tbody>
</table>

当处理结构中出现嵌套时，处理结构必须至少包含两个因素。在这种情况下，嵌套因素的每个水平只与其他因素的一个或几个水平一起出现一次。接下来的两个示例演示了处理结构中的嵌套。

### 示例 5.10：大豆的生育期组 {#sec5-6-2}

农作物大豆对于处理结构中的嵌套是一个很好的例子。大豆品种按生育期组 (maturity groups) 进行分类。在美国中西部地区种植着生育期组 4、5 和 6 的大豆品种，但某一生育期组的品种可能比另一生育期组的品种更适合特定地区。一项研究设计用于评估 8 个大豆品种，其中 2 个来自成熟期组 4，4 个来自生育期组 5，2 个来自生育期组 6。表 \@ref(tab:table5-37) 展示了处理结构，表明品种水平嵌套在生育期水平内。一个可以描述来自这个嵌套处理结构的数据的模型，是一个具有四个区组的随机完全区组设计结构

$$\begin{aligned}&y_{ijk}=\mu+M_i+V(M)_{j(i)}+b_k+\varepsilon_{ijk},\quad i=1,2,3,j=1,\ldots,n_i,k=1,2,3,4\\&b_k\sim i.i.d.N(0,\sigma_{\mathrm{row}}^2),\varepsilon_{ijk}\sim i.i.d.N(0,\sigma_{\varepsilon}^2),n_1=n_3=2,\mathrm{~and~}n_2=4.\end{aligned}$$

其中 $M_i$ 表示第 $i$ 个生育期组的效应，$V(M)_{j(i)}$ 表示嵌套在第 $i$ 个生育期组内的第 $j$ 个品种的效应，$b_k$ 表示第 $k$ 个区组的效应，$\varepsilon_{ijk}$ 表示实验单元误差。表 \@ref(tab:table5-38) 给出了上述模型的方差分析表，其中误差项是通过品种与区组交互作用来计算的。品种与区组交互作用可以划分为：生育期组与区组交互作用+嵌套在生育期组内的品种与区组的交互作用。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-37)品种嵌套于生育期组的大豆研究的处理结构</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.37.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-38)随机完全区组设计结构中大豆品种嵌套于生育期组的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.38.png">
</td>
  </tr>
</tbody>
</table>

上述设计的一个问题是，由于生育期组 4 的大豆比生育期组 5 的大豆更早收获，而生育期组 5 的大豆又比生育期组 6 的大豆更早收获，因此地块 (plots) 需要在不同的时间进行收获。随机完全区组设计的一种替代方案是使用裂区设计，其中整区由生育期组形成，子区由生育期组的品种形成。裂区设计的优点是可以同时收获整区内的所有品种。随机化过程是将生育期组水平随机分配给每个整区内的地块集 (sets of plots)，然后将各自的品种随机分配给生育期组内的地块。这种裂区设计结构的模型是

$$\begin{aligned}y_{ijk}&=\mu+M_i+b_k+w_{ik}+V(M)_{j(i)}+\varepsilon_{ijk},\quad i=1,2,3,j=1,\ldots,n_i,k=1,2,3,4\\b_k&\sim i.i.d.N(0,\sigma_{\mathrm{row}}^2),v_{ik}\sim i.i.d.N(0,\sigma_{\mathrm{wplot}}^2),\mathrm{and}\,\,\varepsilon_{ijk}\sim i.i.d.N(0,\sigma_{\varepsilon}^2)\end{aligned}$$

其中 $w_{ik}$ 表示整区效应，$\varepsilon_{ijk}$ 表示子区效应。整图设计是随机完全区组设计结构中的单向处理结构（三个生育期水平），其中整区误差是根据区组与生育期交互作用计算的。子区设计由三个独立的单向处理结构组成，采用随机完全区组设计结构，每个生育期组一个。由于每个生育期组中有不同数量的品种，因此整区的尺寸不同。这种双向嵌套处理结构在裂区设计结构中的方差分析表如表 \@ref(tab:table5-39) 所示。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-39)裂区设计结构中大豆品种嵌套于生育期组的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.39.png">
</td>
  </tr>
</tbody>
</table>

### 示例 5.11：飞机引擎 {#sec5-6-3}

某飞机公司希望用三种飞机来评估七种引擎的性能。由于某些机械特性，每种飞机只能使用某些引擎。表 \@ref(tab:table5-40) 给出了可能的引擎-飞机配置（标记为“X”）。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-40)嵌套处理结构中可观测的引擎-飞机配置</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.40.png">
</td>
  </tr>
</tbody>
</table>

如表 \@ref(tab:table5-40) 所示，每种引擎类型只能出现在一种飞机类型上，因此引擎水平嵌套在飞机水平内。飞机公司为七种处理组合（引擎类型与飞机类型）中的每一种制造了三架飞机。
数据收集过程是将 21 种引擎型号-飞机型号的配置随机排序，让试飞员按随机顺序驾驶飞机并测量性能特征。设 $y_{ijk}$ 表示由第 $i$ 种飞机型号和第 $j$ 种引擎型号制成的第 $k$ 架飞机的性能测量。可用于描述性能测量的均值模型为

$$\begin{aligned}y_{ijk}&=\mu_{ij}+p_{ijk},\quad\text{for }(i,j)\in\Theta\text{ and }p_{ijk}\thicksim i.i.d.N(0,\sigma_{\text{plane}}^2)\end{aligned}$$

或效应模型

$$y_{ijk}=\mu+A_i+E_{j(i)}+p_{ijk},\quad\mathrm{for~}(i,j)\in\Theta,$$

其中

$$\Theta=\{(1,A),(1,B),(1,C),(2,D),(2,E),(3,F),(3,G)\}$$

表 \@ref(tab:table5-41) 显示了效应模型的方差分析表。使用同一引擎类型和飞机类型的配置制造三架飞机，可计算这三架飞机性能得分的变异，并将这七种配置的性能得分变异合并在一起以计算飞机误差项 (The plane  error term is computed from the variation of the performance scores of the three planes 
made with the same engine type and aircraft type configuration pooled across the seven 
configurations). 就均值模型而言，飞机平方和用于检验飞机均值的相等性，即检验原假设

$$\bar{\mu}_{1\cdot}=\bar{\mu}_{2\cdot}=\bar{\mu}_{3\cdot}$$

或

$$\frac{{\mu}_{1A}+{\mu}_{1B}+{\mu}_{1C}}3=\frac{{\mu}_{2D}+{\mu}_{2E}}2=\frac{{\mu}_{3F}+{\mu}_{3G}}2$$

引擎（嵌套在飞机内）平方和用于检验原假设

$${\mu}_{1A}={\mu}_{1B}={\mu}_{1C},\quad{\mu}_{2D}={\mu}_{2E},\quad\mathrm{and}\quad{\mu}_{3F}={\mu}_{3G}$$

当**处理结构**中出现嵌套时必须非常小心，因为有这样一种执行分析的趋势：就好像存在嵌套或分层**设计结构**。单个飞机是本研究的唯一实验单元，即只有一个尺寸的实验单元和一个误差项。当设计结构中出现嵌套时，存在不止一个尺寸的实验单元，因此模型中存在不止一个误差项。下一个示例说明了设计结构中的嵌套。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-41)完全随机设计结构中双向嵌套处理结构的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.41.png">
</td>
  </tr>
</tbody>
</table>

### 示例 5.12：简单的舒适度实验 {#sec5-6-4}

为了研究温度（三个水平，分别为 18°C、21°C 和 24°C）和性别（两个水平，男性 M 和女性 F）对人员舒适度的影响，进行了一项舒适度实验。测量人的舒适度的方法有很多种，因此，在此讨论中，假设对每个人员只进行一次舒适度测量。三个温度分别随机分配给九个可用环境室 (environmental chambers) 中的三个。室 (chambers) 是温度水平的实验单元，室设计是完全随机设计结构中的单向处理结构。

18 名男性和 18 名女性被随机分配到各个室，因此每个室都分配了两名男性和两名女性。人员性别的实验单元是人员，人员设计是随机区组设计中的单向处理结构，其中室是区组。在每个区组的每个性别水平上都有两次重复。室关于温度的分配和人员关于室的分配如图 \@ref(fig:figure5-24) 所示。

将这些人员置于环境条件 3 小时后，测量他们的舒适度。可用于描述这些数据的均值模型为

$$y_{ijkm}=\mu_{ik}+c_{j(i)}+p_{m(ijk)}\quad i=1,2,3,~j=1,2,3,~k=1,2,~m=1,2$$

其中

$$c_{j(i)}\thicksim i.i.d.N(0,\sigma_{\mathrm{chamber}}^2)\mathrm{~and~}p_{m(ijk)}\thicksim i.i.d.N(0,\sigma_{\mathrm{person}}^2)$$

效应模型为

$$y_{ijkm}=\mu+T_i+c_{j(i)}+S_k+(TS)_{ik}+p_{m(ijk)}$$

其中 $\mu_{ik}$ 表示第 $i$ 个温度和第 $k$ 个性别的均值，$c_{j(i)}$ 是分配给第 $i$ 个温度的第 $j$ 个腔室的随机效应，$p_{m(ijk)}$ 表示分配给第 $j$ 个温度的第 $k$ 个腔室的第 $k$ 个性别的第 $m$ 个人的随机效应，$T_i,S_k,(TS)_{ik}$ 分别表示温度、性别和温度与性别交互作用的效应。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.24.png" alt="将人员和温度分配给室以进行简单的舒适性实验" width="733" />
<p class="caption">(\#fig:figure5-24)将人员和温度分配给室以进行简单的舒适性实验</p>
</div>

这个实验涉及两个层次的嵌套。首先，室嵌套在温度内。其次，同性别的人嵌套在室内。本实验的方差分析表如表 \@ref(tab:table5-42) 所示。性别效应和温度与性别交互作用是人员间比较 (between-person comparisons)，因此使用人员误差项来进行这些均值之间的比较。室误差项用于在温度均值之间进行比较。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-42)简单的舒适度实验嵌套设计结构的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.42.png">
</td>
  </tr>
</tbody>
</table>

### 示例 5.13：重复测量的多地点研究 {#sec5-6-5}

典型的农业研究涉及在多个地点 (locations) 评估处理，且通常会随着时间的推移测量响应。本研究使用三个地点对四个苜蓿品种进行评估。苜蓿是一种在生长季节收获三次或三次以上的作物，因此数据集包含重复测量。在每个地点，该设计是一个随机完全区组设计结构中的单向处理结构，具有三次重复。图 \@ref(fig:figure5-25) 显示了每个地点品种关于地块 (plots) 的分配情况。图 \@ref(fig:figure5-26) 是每个地点的实验单元的图形。由于一个地点的所有地块都是同时收割 (cut) 或收获 (harvest) 的，因此同一地点的所有地块在收割之间都处于相同的环境条件下。因此，重复测量的是给定地点的整个实验，而不是单个地块。图 \@ref(fig:figure5-26) 中的一列表示一个品种被应用和收获三次的实体 (entity)，是品种的实验单元。地块的一层 (layer) 代表受环境条件影响的实体，并且是收割水平的实验单元。切割水平在所有给定地点的地块上形成裂条区结构。单元格表示在一个收割时间收获的品种单元。

<div class="figure" style="text-align: center">
<img src="figure/figure 5.25.png" alt="将苜蓿品种分配给地块以进行多地点实验" width="558" />
<p class="caption">(\#fig:figure5-25)将苜蓿品种分配给地块以进行多地点实验</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 5.26.png" alt="在一个地点显示实验单元，用于多地点研究" width="988" />
<p class="caption">(\#fig:figure5-26)在一个地点显示实验单元，用于多地点研究</p>
</div>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table5-43)重复测量多地点研究的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%205.43.png">
</td>
  </tr>
</tbody>
</table>

可以通过忽略层来构建列的分析，这在每个地点的随机完全区组设计结构中提供了单向处理结构。因此，列误差是通过品种与区组交互作用来计算的，其中合并了所有地点 (pooled 
across locations)。可以通过忽略列来构建层的分析，这在随机完全区组设计结构中提供了一种单向处理结构（收割水平），其中地点是区组。层误差项是通过收割与地点交互作用来计算的。最后，单元格是收割与品种交互作用的单元，并且单元格误差项是这样计算的：收割与品种与区组交互作用，其中合并了所有地点，再加上嵌套于地点内的收割与区组交互作用。可用于描述每个单元格响应的模型是

$$y_{ijkm}=\mu+l_i+b_{j(i)}+V_k+(lV)_{ik}+(bV)_{jk(i)}+C_m+(lC)_{im}+(VC)_{km}+(lVC)_{ikm}+\varepsilon_{ijkm}$$

其中

$$\begin{aligned}&l_i\thicksim N(0,\sigma_{\mathrm{Loc}}^2),\quad b_{j(i)}\thicksim N(0,\sigma_{blk(\mathrm{Loc})}^2),\quad(lV)_{ik}\thicksim N(0,\sigma_{V*\mathrm{Loc}}^2),\quad(bV)_{jk(i)}\thicksim N(0,\sigma_{V*Blk(\mathrm{Loc})}^2),\\&(lC)_{im}\thicksim i.i.d.N(0,\sigma_{\mathrm{C^*Loc}}^2),\quad(lVC)_{ikn}\thicksim i.i.d.N(0,\sigma_{V*C *\text{Loc}}^2),\quad\varepsilon_{ijkn}\thicksim i.i.d.N(0,\sigma_{\mathrm{cell}}^2)\end{aligned}$$

该模型的假设是，将地点和嵌套在地点内的区组视为随机效应（参见第 \@ref(chap18) 章）。因此，所有涉及地点或嵌套在地点内的区组的交互作用也都是随机效应。在每个地点的收割水平被视为重复测量：最先是第一次收割，其次是第二次收割，最后是第三次收割。由于切割水平不能随机分配给地点内的层，因此提供了重复测量。该模型的假设是各层之间相互独立且方差相等，但由于存在重复测量，因此使用更复杂的协方差结构可能更合适。有关重复测量设计的协方差结构建模的讨论，参见第 \@ref(chap26) 章和第 \@ref(chap27) 章。该模型的方差分析表如表 \@ref(tab:table5-43) 所示，其中期望均方是从关于随机效应的模型假设中获得的（参见第 \@ref(chap18) 章）。

## 结束语 {#sec5-7}

在本章中，考虑了几个不同尺寸的实验单元的设计结构。这些设计也称为多水平设计，一些多水平设计称为分层设计。讨论的设计类型包括裂区设计、裂条区设计、嵌套设计、重复测量设计以及它们的变体和组合。本章讨论的设计可以组合成一些非常复杂的设计。在这种情况下，确定一个合适的模型来描述数据的关键在于识别实验单元的尺寸并为每个尺寸的实验单元确定处理结构和设计结构。几乎总是这样：给定尺寸的实验单元的设计结构将是四种基本设计结构之一，并且误差平方和将使用特定基本设计结构的过程进行计算。本章的重点在于识别这样的设计以及何时和如何使用它们。涉及多于一种尺寸的实验单元的设计分析将在第 \@ref(chap24) 章至第 \@ref(chap30) 章中介绍，其中将详细讨论有关模型的假设。有关多水平设计的进一步讨论和示例，请参阅 Milliken (2003a, b) 和 Milliken et al. (1998).

## 练习 {#sec5-8}

