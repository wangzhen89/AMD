
# 分析裂区型设计的方法 {#chap24}

> "Numerical quantities focus on expected values, graphical summaries on unexpected values." - John Tukey

## 介绍 {#sec24-1}

裂区型设计 (split-plot type design) 涉及一种具有多个实验单元尺寸的设计结构，其中较小尺寸的实验单元嵌套在较大尺寸的实验单元内。第 \@ref(chap5) 章介绍了裂区设计结构的一些示例，包括分层设计结构类 (class of hierarchal design structures). 裂区型设计结构的设计和分析存在两个主要问题。第一个问题包括选择和/或识别设计结构中使用的不同尺寸的实验单元，然后将处理从处理结构分配到设计结构中的不同实验单元尺寸。成功识别不同尺寸的实验单元对于指定可以描述结果数据的适当模型至关重要。第二个问题是构建适当的模型来描述处理和设计结构的相关特征。能够识别与每种尺寸的实验单元相关的变异源是很重要的。这些变异源用于计算各自的误差项，所述误差项用于计算均值估计的标准误估计 (estimates of the standard errors of estimated means) 以及用于均值之间的成对比较。由于这些设计结构涉及一个以上尺寸的实验单元，因此固定效应参数的标准误估计及其之间的比较涉及一个或多个变异源。裂区型设计模型的一个非常重要的特征是，它们是第 \@ref(chap26) 章中讨论的重复测量模型构建的基本模型。第 \@ref(chap5) 章介绍了几个概念的例子。

第 \@ref(sec24-1) 节解释了具有两种尺寸实验单元的裂区或分层设计结构的设计和分析，第 \@ref(sec24-2) 节描述了与固定效应相关的标准误的确定和估计。第 \@ref(sec24-3) 节讨论了在一般裂区设计结构中确定适当标准误及其对固定效应参数估计的一般方法。均值对比的标准误的计算在第 \@ref(sec24-4) 节中讨论。第 \@ref(sec24-6) 节介绍了裂区设计结构的四个示例，其中每个示例都演示了分析此类设计的一些显著特征。第 \@ref(sec24-7) 节讨论了裂区设计结构样本量的确定和功效的计算。本章给出了使用 SAS^®^-Mixed 和 JMP^®^ 的分析，其中 JMP 分析如第 \@ref(sec24-8) 节所示。

构建裂区设计模型的关键概念是识别实验单元的不同尺寸，然后确定相应的设计结构和处理结构。整体模型是通过整合为每种尺寸的实验单元开发的模型而构建的。第 \@ref(chap5) 章中给出了模型构建的几个例子，但没有说明模型的基本假设。这些假设是，表示各种实验单元误差项的分量都是独立分布的，均值为零，并具有相关的方差 (associated variance)（更一般的假设请参见第 \@ref(chap26) 章）。在理想条件下，误差项呈正态分布。分析的目的是利用模型假设来获得总体参数的估计，并对其进行推断。在以下示例中，将使用矩法和 REML 方法来演示估计固定效应标准误所需的计算。在实践中，REML 是可在大多数情况下推荐使用的方法。

### 示例 24.1：面包配方和烘焙温度 {#sec24-1-1}

烤面包的过程包括根据配方 (recipe) 的规格混合一批面包面团，将面团放入烤盘（容器）中，让面包发酵，然后将装有面包面团的烤盘放入烤箱 (oven) 中，在特定的温度 (temperature) 和时间组合下进行烘烤。每个烤箱都足够大，可以同时放入四个装有面包面团的烤盘。我们设计了一项实验，以评估四种不同的面包配方在三种不同温度下的烘烤效果，其中测量的响应是最终面包的体积。过程是，根据这四种配方分别制作面团，并将每种配方的一个面包放入设定为特定温度的单个烤箱中。这些批次在烤箱中放置指定的时间长度，然后将面包冷却至室温后再测量其体积。表 \@ref(tab:table24-1) 中的数据是在三天内重复该过程时，由四种配方和三种温度制成的面包的体积。在设计结构中，天 (day) 被视为一个区组因子。

<table>
<caption>(\#tab:table24-1)示例 [24.1](#sec24-1-1) 各种配方和温度下的面包体积 (cm^3^)</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.1.png">
</td>
  </tr>
</tbody>
</table>

图 \@ref(fig:figure24-1) 中的图表演示了为烤箱分配温度的过程，因此烤箱是温度水平的实验单元（请注意，配方不包括在该过程的此步骤中）。与烤箱尺寸实验单元相关的设计是随机完全区组设计结构（三天）中的单向处理结构（三个温度水平）。可用于描述每个烤箱内四个面包的平均面包体积（即每个烤箱一次观测）的模型是

$$y_{ik}^o=\mu_i^o+d_k^o+o_{ik}^o$$

其中 $y_{ik}^o$ 表示观测到的平均面包体积，$\mu_{i}^o$ 表示第 i 个温度水平的平均面包体积，$d_{k}^o$ 表示第 k 天的随机效应，$o_{ik}^o$ 表示第 i 个温度和第 k 天的随机烤箱效应。假设 $d_{k}^o\thicksim i.i.d.N(0,\sigma_d^2),o_{ik}^o\thicksim i.i.d.N(0,\sigma_o^2)$ 并且所有 $d_{k}^o,o_{ik}^o$ 独立。与烤箱模型相关的方差分析显示在表 \@ref(tab:table24-2) 中，其中烤箱误差项计算为温度与天的交互均方。

<div class="figure" style="text-align: center">
<img src="figure/figure 24.1.png" alt="显示每天将温度随机化到烤箱的示意图" width="850" />
<p class="caption">(\#fig:figure24-1)显示每天将温度随机化到烤箱的示意图</p>
</div>

<table>
<caption>(\#tab:table24-2)烤箱水平分析的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.2.png">
</td>
  </tr>
</tbody>
</table>

图 \@ref(fig:figure24-2) 中的图表显示了将配方分配到每个烤箱内的位置。一天内的烤箱是四个配方的一个区组。每个配方提供一条面包。面包设计是一个单向处理结构（四个配方），采用随机完全区组设计结构，共有九个区组（三天中每天三个区组或烤箱）。

<div class="figure" style="text-align: center">
<img src="figure/figure 24.2.png" alt="显示每天将配方随机化到一个烤箱内的位置的示意图" width="852" />
<p class="caption">(\#fig:figure24-2)显示每天将配方随机化到一个烤箱内的位置的示意图</p>
</div>

如果所有烤箱处于相同温度，则数据结构将是随机完全区组设计结构中的单向处理结构。但并非所有烤箱都经相同处理，三个烤箱设置为 325°F，三个设置为 340°F，三个设置为 355°F。接下来通过考虑 325°F 的所有数据来简化设计，如图 \@ref(fig:figure24-3) 所示。所得数据来自具有三个区组的随机完全区组设计结构中的单向处理结构。该数据的误差项是配方与天的交互作用。可用于描述 325°F 烤箱数据的模型是

$$y_{jk}^+=\mu_j^++o_k^++\varepsilon_{jk}^+$$

其中 $y_{jk}^+$ 表示第 k 天在第 j 个配方的 325°F 烤箱中观测到的面包体积，$\mu_{j}^+$ 表示第 j 个配方的平均面包体积，$o_{ik}^o$ 表示第 k 天使用的烤箱的随机效应，以及 $\varepsilon_{jk}^+$ 表示第 j 个配方和第 k 天的随机面包效应。假设 $o_k^+\thicksim i.i.d.N(0,\sigma_0^2),\varepsilon_{jk}^+\thicksim i.i.d.N(0,\sigma_\varepsilon^2)$ 且所有 $o_k^+,\varepsilon_{jk}^+$ 独立分布。与面包模型相关的方差分析显示在表 \@ref(tab:table24-3) 中，其中面包误差项根据天与配方的交互均方来估计。

<div class="figure" style="text-align: center">
<img src="figure/figure 24.3.png" alt="325°F 烤箱设计的一部分" width="852" />
<p class="caption">(\#fig:figure24-3)325°F 烤箱设计的一部分</p>
</div>

<table>
<caption>(\#tab:table24-3)325°F 下面包体积数据的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.3.png">
</td>
  </tr>
</tbody>
</table>

测量面包间变异性的误差项是通过合并三天内不同温度下的配方与天交互作用计算得出的，得到 Error(loaf)=Recipe × Day(Temperature). 将这两个模型结合起来得到模型

$$y_{ijk}=\mu_{ij}+d_k+o_{ik}+\varepsilon_{ijk},\quad i=1,2,3,\quad j=1,2,3,4,\quad k=1,2,3$$

其中 ${\mu}_{ij}={\mu}+T_i+R_j+(TR)_{ij}$，$T_i$ 表示第 i 个温度的效应，$R_i$ 表示第 j 个配方的效应，$(TR)_{ij}$ 表示温度与配方交互作用，$d_k$ 表示第 k 天的效应，$o_{ik}$ 表示第 k 天的第 i 个烤箱的效应，以及 $\varepsilon_{ijk}$ 表示误差项。在理想条件下 $d_k\thicksim i.i.d.N(0,\sigma_{\mathrm{day}}^2),o_{ik}\thicksim i.i.d.N(0,\sigma_{\mathrm{oven}}^2),\varepsilon_{ijk}\thicksim i.i.d.N(0,\sigma_{\mathrm{loaf}}^2)$，且所有 $d_k,o_{ik},\varepsilon_{ijk}$  独立分布。烤箱是整区 (whole plot) 或较大尺寸的实验单元，面包是子区 (subplot) 或裂区或较小尺寸的实验单元。

该模型可以用实验单元的尺寸表示为

$$\begin{aligned}y_{ijk}=&\mu_{ij}+d_k+T_i+o_{ik} &&\}\quad\text{whole-plot or oven part of the model}
\\&+R_j+(TR)_{ij}+\varepsilon_{ijk}&&\}\quad\text{subplot or loaf part of the model}\end{aligned}$$

结合表 \@ref(tab:table24-2) 和 \@ref(tab:table24-3) 中的方差分析表给出表 \@ref(tab:table24-4) 中该模型的方差分析表。预期望均方决定了计算固定效应检验统计量的适当分母。Error(oven) 用作检验温度主效应的误差，Error(loaf) 用作检验配方主效应和配方与温度交互效应的误差。表 \@ref(tab:table24-5) 中列出了面包体积数据的 SAS-Mixed 代码和方差分析结果表。有迹象表明存在显著的温度与配方交互作用 (p = 0.0021)，因此应遵循 temperature × recipe 双向均值之间的进一步比较（参见第 \@ref(tab:table24-2) 节）

<table>
<caption>(\#tab:table24-4)示例 [24.1](#sec24-1-1) 面包体积数据的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.4.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table24-5)示例 [24.1](#sec24-1-1) 面包体积数据的方差分析表和 SAS-Mixed 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.5.png">
</td>
  </tr>
</tbody>
</table>

### 示例 24.2：在不同肥力条件下生长的小麦品种 {#sec24-1-2}

图 \@ref(fig:figure24-4) 中的数据是在四种不同的肥力条件 (fertility regimes) 下（A1, A2, A3 和 A4）种植的两种小麦（B1 和 B2）的产量（单位为磅）。田地被分成两个区组 (blocks)，每个区组有四个整区 (whole plots). 每个肥料 (fertilizer) 水平被随机分配到每个区组内的一个整区。因此，整区设计由单向处理结构（四个肥料水平）和随机完全区组设计结构（两个区组）组成。每个区组包含四个整区实验单元，这些单位被分成两部分（称为子区，subplots）。每种小麦被随机分配到每个整区内的一个子区。子区设计由单向处理结构（两个品种）和随机完全区组设计组成，共有八个区组，每个区组包含两个子区实验单元。可用于描述这些数据的模型是

$$y_{ijk}=\mu_{ij}+b_k+w_{ik}+\varepsilon_{ijk},\quad i=1,2,3,4,\quad j=1,2,k=1,2$$

其中，$\mu_{ij}$ 表示第 i 个肥料水平和第 j 个品种的预期响应（产量），$y_{ijk}$ 表示第 k 个肥料水平和第 j 个品种的观测产量（响应），$b_k$ 表示假设分布为 $i.i.d.N(0,\sigma_{\mathrm{block}}^2)$ 的区组效应，$w_{ik}$ 表示假设分布为 $i.i.d.~N(0,~\sigma_{wp}^2)$ 的整区误差，以及 $\varepsilon_{ijk}$ 表示假设分布为 $i.i.d.~N(0,~\sigma_{\varepsilon}^2)$ 的子区误差。还假设所有 $b_k,w_{ik},\varepsilon_{ijk}$ 独立分布。平均响应可以使用效应模型表示为 $\mu_{ij}=\mu+F_i+V_j+(FV)_{ij}$. 该效应模型的方差分析表见表 \@ref(tab:table24-6). 三个固定效应比较的 $F$ 检验的分母由期望均方确定；也就是说，Error(whole plot) 用于检验肥料主效应，Error(subplot) 用于检验通过品种主效应和肥料与品种交互效应。表 \@ref(tab:table24-7) 给出了 SAS Mixed 代码和使用 III 型平方和的数值结果，表 \@ref(tab:table24-8) 给出了使用 REML 选项的结果。这两种分析的固定效应的 $F$ 检验是相同的，因为数据集是均衡的，并且两个方差分量的估计都大于零。

<div class="figure" style="text-align: center">
<img src="figure/figure 24.4.png" alt="品种与肥力条件的裂区示例数据" width="708" />
<p class="caption">(\#fig:figure24-4)品种与肥力条件的裂区示例数据</p>
</div>

<table>
<caption>(\#tab:table24-6)示例 [24.2](#sec24-1-2) 小麦产量数据的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.6.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table24-7)示例 [24.2](#sec24-1-2) 小麦产量数据使用 III 型平方和的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.7.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table24-8)示例 [24.2](#sec24-1-2) 小麦产量数据使用 REML 的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.8.png">
</td>
  </tr>
</tbody>
</table>

## 模型定义和参数估计 {#sec24-2}

考虑随机完全区组整区设计结构的裂区设计，其中具有 r 个区组、整区因子 A 有 a 个水平以及子区因子 C 有 c 个水平，其一般模型可以表达为：

$$y_{ijk}=\mu+\alpha_i+b_k+w_{ik}+\gamma_j+(\alpha\gamma)_{ij}+\varepsilon_{ijk},\quad i=1,2,\ldots,a,\quad j=1,2,\ldots,c,\quad k=1,2,\ldots,r$$

其中，$y_{ijk}$ 是观测到的响应，$b_k$ 表示第 k 个区组效应，假设其分布为 $N(0,\sigma_B^2)$，$w_{ik}$ 表示整区误差，假设其分布为 $N(0,\sigma_w^2)$，$\varepsilon_{ijk}$ 表示子区误差，假设其分布为 $N(0,\sigma_\varepsilon^2)$. 还假设所有的 $b_k,w_{ik},\varepsilon_{ijk}$ 都是独立分布的。值得注意的是，最重要的假设是所有的 $b_k,w_{ik},\varepsilon_{ijk}$ 都是独立分布的。幸运的是，由于固定效应因子被随机分配到其适当尺寸的实验单位，这一假设可以通过随机化过程得到保证。该模型中的固定效应包括总均值 $\mu$，整区因子（A）的效应 $\alpha_i$、子区因子（C）的效应 $\gamma_j$，以及整区因子水平和子区因子水平之间的交互效应 $(\alpha\gamma)_{ij}$. 均值模型可以用效应模型表示为 $\mu_{ij}=\mu+\alpha_i+\gamma_j+(\alpha\gamma)_{ij}$. 表 \@ref(tab:table24-9) 给出了这个一般模型的方差分析表，包括变异来源、自由度和期望均方。

<table>
<caption>(\#tab:table24-9)\@ref(sec24-1) 节一般裂区模型的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.9.png">
</td>
  </tr>
</tbody>
</table>

整区误差计算为 Block × A 交互作用，子区误差是通过合并 A 的所有水平下 Block × C 交互作用计算得出的，记作 Block × C(A). 用于获得方差分量矩法估计的方程组为

$$\begin{aligned}
\text{MSBlock}& =\tilde{{\sigma}}_\varepsilon^2+c\tilde{{\sigma}}_w^2+ac\tilde{{\sigma}}_B^2  \\
\text{MSError(whole plot)}& =\tilde{\sigma}_\varepsilon^2+c\tilde{\sigma}_w^2  \\
MSError(subplot)&=\tilde{\sigma}_\varepsilon^2
\end{aligned}$$

这些方程的矩法解是

$$\begin{aligned}
\tilde{{\sigma}}_\varepsilon^2& =MSError(subplot)  \\
\tilde{{\sigma}}_w^2& =\frac{MSError(wholeplot)-MSError(subplot)}c  \\
\tilde{{\sigma}}_B^2& =\frac{MSBlock-MSError(wholeplot)}{ac}
\end{aligned}$$

方差分量的矩法估计为

$$\begin{aligned}
&\hat{{\sigma}}_\varepsilon^2 =\tilde{\sigma}_\varepsilon^2  \\
&\hat{{\sigma}_w^2} =\begin{cases}\tilde{{\sigma}}_w^2&\operatorname{~if~}\tilde{{\sigma}}_w^2>0\\0&\operatorname{~if~}\tilde{{\sigma}}_w^2\leq0&\end{cases}  \\
&\hat{{\sigma}}_{B}^{2} =\begin{cases}\tilde{{\sigma}}_B^2&\mathrm{~if~}\tilde{{\sigma}}_B^2>0\\0&\mathrm{~if~}\tilde{{\sigma}}_B^2\leq0&\end{cases} 
\end{aligned}$$

$\mu_{ij},\bar{\mu}_{i\cdot},\bar{\mu}_{\cdot j}$ 的估计分别为 $\bar{{y}}_{ij\cdot},\bar{{y}}_{i{\cdot\cdot}},\bar{{y}}_{\cdot j\cdot}$. A 水平之间的比较是整区间比较 (between whole plot comparisons)，用于检验 A 水平相等性的适当 $F$ 统计量是 $F_A=MSA/MSError(whole plot)$. C 水平之间的比较和 A × C 交互作用的比较是整区内比较 (within whole plot comparisons) 或整区内的子区间比较 (between subplot comparisons within a whole plot)，并且适当的 $F$ 统计是 $F_C=MSC/MSError(subplot)$ 和 $F_{A×C}=MSA×C/MSError(subplot)$. 这些 $F$ 统计量是通过查看表 \@ref(tab:table24-9) 中的期望均方来构建的。

一旦计算出 $F$ 检验以确定均值之间是否存在显著差异，下一步就是进行多重比较以确定差异发生的位置。以下部分介绍了计算裂区设计的各种均值差异的标准误的方法。

## 均值间比较的标准误 {#sec24-3}

使用处理均值或处理组合均值的对比来研究处理效应，特别是当方差分析表明一种或多种固定效果与零显著不同时。样本均值对比的标准误是必要的，以确定均值对比是否等于零，或者构造关于均值对比的置信区间。通常，对比涉及两个均值的比较。因此，本节讨论了两个均值的比较，而一般的对比将在第 \@ref(sec24-4) 节中讨论。

为了演示确定适当标准误的方法，使用了一般裂区设计模型，其中整区设计是随机完全区组设计中的单向处理结构，子区处理结构涉及单向处理结构。这种情况下相应的效应模型可以表示为

$$y_{ijk}=\mu+\alpha_i+b_k+iv_{ik}+\gamma_j+(\alpha\gamma)_{ij}+\varepsilon_{ijk},\quad i=1,2,\ldots,a,\quad j=1,2,\ldots,c,\quad k=1,2,\ldots,r$$

均值模型可以表示为

$$y_{ijk}=\mu_{ij}+b_k+w_{ik}+\varepsilon_{ijk},\quad i=1,2,\ldots,a,~j=1,2,\ldots,c,~k=1,2,\ldots,r$$

其中各项如第 \@ref(sec24-1) 节所述。

根据 A 水平和 C 水平之间是否存在任何交互作用，我们可能会感兴趣于四种类型的比较。如果没有交互作用，则有兴趣将 A 的水平相互比较，并将 C 的水平相互比较。

[^expressingy]: 原文：The process of determining the appropriate standard error involves expressing $\bar{{y}}_{\cdot j\cdot}$ in terms of the quantities in the model obtained by summing over i and k.

为比较 C 的水平，需要比较两个 $\bar{\mu}_{\cdot j}$，这是由 $\bar{{y}}_{\cdot j\cdot}$ 估计的。确定适当标准误的过程涉及用模型的一些量来表示 $\bar{{y}}_{\cdot j\cdot}$，这些量是通过对 i 和 k 求和得到的[^expressingy]。C 的第 j 个主效应均值的模型是 $\bar{y}_{{\cdot}j{\cdot}}=\bar{{\mu}}_{{\cdot}j}+\bar{b}_{{\cdot}}+\bar{{w}}_{{\cdot}\cdot}+\bar{{\varepsilon}}_{{\cdot}j{\cdot}}$. 考虑差值 $\bar{{\mu}}_{\cdot1}-\bar{{\mu}}_{\cdot2}$， 其估计为 $\bar{y}_{\cdot 1\cdot}-\bar{y}_{\cdot 2\cdot}$，可以用关于 C 的均值模型表示为 $\bar{y}_{\cdot1\cdot}-\bar{y}_{\cdot2\cdot}=\bar{\mu}_{\cdot1}-\bar{\mu}_{\cdot2}+\bar{\varepsilon}_{\cdot1\cdot}-\bar{\varepsilon}_{\cdot2\cdot}$，因为涉及 $\bar{b}_{{\cdot}}+\bar{{w}}_{{\cdot}\cdot}$ 的项抵消了；也就是说，比较 $\bar{y}_{\cdot1\cdot}-\bar{y}_{\cdot2\cdot}$ 不取决于整区误差，也不取决于区组误差。$\bar{y}_{\cdot1\cdot}-\bar{y}_{\cdot2\cdot}$ 的方差可以证明等于

$$\mathrm{Var}(\bar{y}_{\cdot1\cdot}-\bar{y}_{\cdot2\cdot})=\mathrm{Var}(\bar{\varepsilon}_{\cdot1\cdot}-\bar{\varepsilon}_{\cdot2\cdot})=\frac{2\sigma_\varepsilon^2}{ar}$$


其中均值 $\bar{\varepsilon}_{\cdot1\cdot}$ 的方差为 $\mathrm{Var}(\bar{{\varepsilon}}_{\cdot1\cdot})=\sigma_\varepsilon^2/ar$ 其中 ar 为均值中的观测次数。类似地，对任何 $j\ne j^\prime$，$\mathrm{Var}(\bar{y}_{\cdot j\cdot}-\bar{y}_{\cdot j^\prime\cdot})=2\sigma_\varepsilon^2/ar$. $\bar{y}_{\cdot j\cdot}-\bar{y}_{\cdot j^\prime\cdot}$ 的标准误估计为

$$\widehat{s.e.}(\bar{y}_{.j.}-\bar{y}_{.j'.})=\sqrt{\frac{2\hat{\sigma}_{\varepsilon}^2}{ar}}=\sqrt{\frac{2MSError(subplot)}{ar}}\quad\text{for all }j\neq j'$$

基于 a(c-1)(r-1) 自由度。如果要进行多重比较（见第 \@ref(chap3) 章），a(c-1)(r-1) 是在确定所需多重比较程序的分位数时需要使用的自由度。为了简单起见，计算 LSD 值，但 LSD 值可能不是用于给定情况的合适方法。用于比较两个子区处理均值的 LSD 值为

$$\mathrm{LSD}_\alpha=[t_{\alpha/2,a(c-1)(r-1)}]\widehat{s.e.}(\bar{y}_{\cdot j\cdot}-\bar{y}_{\cdot j^{\prime}\cdot})$$

为了比较 A 的水平，需要比较由 $\bar{{y}}_{i\cdot \cdot}$ 估计的 $\bar{{\mu}}_{i\cdot}$. 量 $\bar{{\mu}}_{i\cdot}$ 可以用一般模型表示，通过对 j 和 k 求和，得到 $\bar{y}_{i\cdot\cdot}=\bar{\mu}_{i\cdot}+\bar{b}_{\cdot}+\bar{w}_{i\cdot}+\bar{\varepsilon}_{ij\cdot\cdot}$

对比 $\bar{{\mu}}_{1\cdot}-\bar{{\mu}}_{2\cdot}$ 的估计为 $\bar{{y}}_{1\cdot \cdot}-\bar{{y}}_{2\cdot \cdot}$，可以用关于 $\bar{{y}}_{i\cdot \cdot}$ 的模型表示为 $\bar{y}_{1\cdot\cdot}-\bar{y}_{2\cdot\cdot}=\bar{\mu}_{1\cdot}-\bar{\mu}_{2\cdot}+\bar{w}_{1\cdot}-\bar{w}_{2\cdot}+\bar{\varepsilon}_{1\cdot\cdot}-\bar{\varepsilon}_{2\cdot\cdot}$. 这种比较取决于整区和子区方差分量。$\bar{{\mu}}_{1\cdot}-\bar{{\mu}}_{2\cdot}$ 的方差为

$$\begin{aligned}
\mathrm{Var}(\bar{y}_{1\cdot\cdot}-\bar{y}_{2\cdot\cdot}) &=\operatorname{Var}(\bar{w}_{1\cdot}-\bar{w}_{2\cdot}+\bar{\varepsilon}_{1\cdot\cdot}-\bar{\varepsilon}_{2\cdot\cdot}) \\
&=\frac{2\sigma_w^2}r+\frac{2\sigma_\varepsilon^2}{rc} \\
&=\frac{2(\sigma_\varepsilon^2+c\sigma_w^2)}{rc} 
\end{aligned}$$

$\bar{{\mu}}_{1\cdot}-\bar{{\mu}}_{2\cdot}$ 的标准误估计为

$$\widehat{s.e.}(\bar{y}_{1\cdot\cdot}-\bar{y}_{2\cdot\cdot})=\sqrt{\frac{2(\hat{\sigma}_\varepsilon^2+c\hat{\sigma}_w^2)}{rc}}=\sqrt{\frac{2MSError(wholeplot)}{rc}}$$

基于 (r-1)(a-1) 个自由度。比较两个整区处理均值的 LSD 值为 $\mathrm{LSD}_\alpha=[t_{\alpha/2,(a-1)(r-1)}]\widehat{\mathrm{s.e.}}(\bar{y}_{1\cdot\cdot}-\bar{y}_{2\cdot\cdot})$

当存在显著的 A × C 交互作用时，比较通常必须基于双向单元格均值的集合。在研究这些单元格均值时，必须考虑两种不同类型的比较。**第一种类型出现在比较同一整区处理（A）水平下的两个子区处理（C）均值时**，例如 $\mu_{11}-\mu_{12}$. $\mu_{11}-\mu_{12}$ 的最佳估计为 $\bar{y}_{11\cdot}-\bar{y}_{12\cdot}$. 通过在一般模型中关于 k 求和，$\bar{y}_{ij\cdot}$ 可表示为 ${\bar{y}}_{ij{\cdot}}={\mu}_{ij}+\bar{b}_\cdot+\bar{w}_{i{\cdot}}+{\bar{\varepsilon}}_{ij{\cdot}}$，以及 $\bar{y}_{11\cdot}-\bar{y}_{12\cdot}$ 的方差为 $\mathrm{Var}(\bar{y}_{11\cdot}-\bar{y}_{12\cdot})=\mathrm{Var}(\bar{w}_{1{\cdot}}+\bar{{\varepsilon}}_{11{\cdot}}-\bar{w}_{1{\cdot}}\bar{{\varepsilon}}_{12\cdot})=2{\sigma}_{{\varepsilon}}^2/r$. 因此，**在整区处理的相同水平上，子区处理间比较的方差仅取决于子区误差**。$\bar{y}_{11\cdot}-\bar{y}_{12\cdot}$ 标准误估计为

$$\widehat{s.e.}(\bar{y}_{11\cdot}-\bar{y}_{12\cdot})=\sqrt{\frac{2\hat{\sigma}_\varepsilon^2}r}=\sqrt{\frac{2MSError(subplot)}r}$$

相应的 LSD 值为 $\mathrm{LSD}_{\alpha}=[t_{\alpha/2,a(r-1)(c-1)}]\widehat{s.e.}(\bar{y}_{1\cdot}-\bar{y}_{1\cdot2})$. 该 LSD 值可用于在整区处理的相同水平上比较任何一对子区处理。

**第二种类型的比较发生在相同或不同子区处理水平下比较两个整区处理**，例如 $\mu_{11}-\mu_{21}$ 或 $\mu_{11}-\mu_{22}$. 这两种类型的比较具有相同的标准误。$\mu_{11}-\mu_{21}$ 的最佳估计为 $\bar{y}_{11\cdot}-\bar{y}_{21\cdot}$，用一般模型可表示为 $\bar{y}_{11\cdot}-\bar{y}_{21\cdot}={\mu}_{11}-{\mu}_{21}+\bar{{w}}_{1\cdot}-\bar{{w}}_{2{\cdot}}+\bar{{\varepsilon}}_{11{\cdot}}-\bar{{\varepsilon}}_{21{\cdot}}$. 那么

$$\begin{aligned}
\mathrm{Var}(\bar{y}_{11\cdot}-\bar{y}_{21\cdot}) =&\operatorname{Var}(\bar{w}_{1\cdot}-\bar{w}_{2\cdot}+\bar{\varepsilon}_{11\cdot}-\bar{\varepsilon}_{21\cdot}) \\
=&\frac{2\sigma_w^2}r+\frac{2\sigma_\varepsilon^2}r \\
=&\frac{2(\sigma_\varepsilon^2+\sigma_w^2)}r 
\end{aligned}$$

**这种比较取决于整区和子区方差分量**。

$\sigma_\varepsilon^2+\sigma_w^2$ 的一个无偏估计为

$$\widehat{\sigma_\varepsilon^2+\sigma_w^2}=\frac{MSError(wholeplot)+(c-1)MSError(subplot)}c$$

有关 $\widehat{\sigma_\varepsilon^2+\sigma_w^2}$ 的抽样分布不是卡方分布，而是卡方分布的线性组合。有关 $\widehat{\sigma_\varepsilon^2+\sigma_w^2}$ 的自由度可使用 Satterthwaite 近似估计为

$$\hat{v}=\frac{\left(\widehat{\sigma_\varepsilon^2+\sigma_w^2}\right)^2}{\frac{\begin{bmatrix}MSE(wholeplot)/c\end{bmatrix}^2}{(r-1)(a-1)}+\frac{\left[\frac{c-1}cMSE(subplot)\right]^2}{a(r-1)(c-1)}}$$

在相同或不同子区处理下比较两个整区处理的近似 LSD 值为

$$\mathrm{LSD}_\alpha=(t_{\alpha/2,\hat{v}})\sqrt{\frac{2(\widehat{\sigma_\varepsilon^2+\sigma_w^2})}r}$$

对于示例 [24.1](#sec24-1-1)（表 \@ref(tab:table24-1)）中的数据，表 \@ref(tab:table24-10) 中给出了温度均值、配方均值和温度与配方均值。四种类型比较的标准误估计计算如下：

<table>
<caption>(\#tab:table24-10)表 \@ref(tab:table24-1) 中数据的配方、温度以及配方与温度组合的均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.10.png">
</td>
  </tr>
</tbody>
</table>

1) **比较配方的主效应均值**：标准误估计为

$$\widehat{s.e.}(\bar{y}_{\cdot1\cdot}-\bar{y}_{\cdot2\cdot})=\sqrt{\frac{2(657.87)}{3(3)}}=12.09$$

并且该标准误估计基于 18 个自由度。注意到 $t_{0.025,18}=2.101$，因此用于比较配方主效应均值的 5% LSD 值为 $\text{LSD}_{0.05}=2.101(12.09)=25.40$.

2) 比较温度主效应均值：标准误估计为

$$\widehat{s.e.}(\bar{y}_{1\cdot\cdot}-\bar{y}_{2\cdot\cdot})=\sqrt{\frac{2\left(4096.42\right)}{3(4)}}=26.13$$

并且该标准误估计基于 4 个自由度。那么 $t_{0.025,4}=2.776$ 用于比较温度主效应均值的 5% LSD 值为 $\text{LSD}_{0.05}=2.776(26.13) = 72.55$.

3) 在相同温度水平下比较配方均值：标准误估计为

$$\widehat{s.e.}(\bar{y}_{11\cdot}-\bar{y}_{12\cdot})=\sqrt{\frac{2(657.87)}3}=20.94$$

并且该估计的标准误差基于 18 个自由度。因此，在相同温度水平下比较两个配方均值的 5% LSD 值为 $LSD _{0.05}=2.101(20.94)=43.99$.

4) 比较相同或不同配方水平下的温度均值：标准误估计为

$$\widehat{s.e.}(\bar{y}_{11\cdot}-\bar{y}_{21\cdot})=\sqrt{\frac{2\left[\frac144096.42+\frac{(4-1)}4657.87\right]}3}=\sqrt{\frac{2\left[1517.51\right]}3}=31.80$$

基于

$$\hat{v}=\frac{\left(1517.51\right)^2}{\frac{\left(\frac144096.42\right)^2}4+\frac{\left(\frac34657.87\right)^2}{18}}=8.35$$

个自由度。注意到 $t_{0.025,8.35}=2.289$，因此，比较处于不同温度水平的两均值的 5% LSD 值为 $LSD _{0.05}=2.289(31.80) = 72.79$.

接下来描述了为裂区型设计计算标准误的一般过程。

## 计算均值差标准误的一般方法 {#sec24-4}

将上述技术应用于涉及多种实验单位尺寸的比较的标准误计算并不总是那么直接。本节描述了一种计算标准误和近似自由度的一般方法，该方法可应用于更复杂的情况。普通读者可以跳过这一节。该方法通过模型示例进行描述。

$$y_{ijk}=\mu+\alpha_i+b_k+w_{ik}+\gamma_j+(\alpha\gamma)_{ij}+\varepsilon_{ijk},\quad i=1,2,\ldots,a,\quad j=1,2,\ldots,c,\quad k=1,2,\ldots,r$$

基本上，该技术包括将给定的均值比较表示为分量的总和，其中每个分量是仅涉及一个实验单元尺寸的均值的比较。然后，比较的分量是独立分布的，比较的方差是通过将各分量的方差相加得到的。例如，${\mu_{11}}-{\mu_{21}}$ 可表示为

$$\mu_{11}-\mu_{21}=(\bar{\mu}_{1\cdot}-\bar{\mu}_{2\cdot})+[(\mu_{11}-\bar{\mu}_{1\cdot})-(\mu_{21}-\bar{\mu}_{2\cdot})]$$

分量 $(\bar{\mu}_{1\cdot}-\bar{\mu}_{2\cdot})$ 是整区比较，分量 $[(\mu_{11}-\bar{\mu}_{1\cdot})-(\mu_{21}-\bar{\mu}_{2\cdot})]$ 是子区比较。这些分量的估计为 $\hat{\bar{{\mu}}}_{1\cdot}-\hat{\bar{{\mu}}}_{2\cdot}=\bar{{y}}_{1\cdot\cdot}-\bar{{y}}_{2\cdot\cdot}$ 和 $(\hat{\mu}_{11}-\hat{\bar{\mu}}_{1\cdot})-(\hat{{\mu}}_{21}-\hat{\bar{\mu}}_{2\cdot})=(\bar{y}_{11\cdot}-\bar{y}_{1\cdot\cdot})-(\bar{y}_{21\cdot}-\bar{y}_{2\cdot\cdot})$. ${\mu_{11}}-{\mu_{21}}$ 的估计为

$$\begin{aligned}
\hat{{\mu}}_{11}-\hat{{\mu}}_{21}& =\hat{\bar{\mu}}_{1\cdot}-\hat{\bar{\mu}}_{2\cdot}+(\hat{\mu}_{11}-\hat{\bar{\mu}}_{1\cdot})-(\hat{\mu}_{21}-\hat{\bar{\mu}}_{2\cdot})  \\
&=\bar{y}_{11\cdot}-\bar{y}_{21\cdot}
\end{aligned}$$

由于从模型的整区部分计算出的比较与从子区部分计算出的比较是独立分布的，因此 $\hat{{\mu}}_{11}-\hat{{\mu}}_{21}$ 的方差是

$$\mathrm{Var}(\hat{\mu}_{11}-\hat{\mu}_{21})=\mathrm{Var}(\bar{y}_{1\cdot\cdot}-\bar{y}_{2\cdot\cdot})+\mathrm{Var}[(\bar{y}_{11\cdot}-\bar{y}_{1\cdot\cdot})-(\bar{y}_{21\cdot}-\bar{y}_{2\cdot\cdot})]$$

量 $\bar{y}_{1\cdot\cdot}-\bar{y}_{2\cdot\cdot}$ 可用一般模型表示为 $\bar{y}_{1\cdot\cdot}-\bar{y}_{2\cdot\cdot}=\bar{\mu}_{1\cdot}-\bar{\mu}_{2\cdot}+\bar{w}_{1\cdot}-\bar{w}_{2\cdot}+\bar{\varepsilon}_{1\cdot\cdot}-\bar{\varepsilon}_{2\cdot\cdot}$，其方差为

$$\begin{aligned}
\mathrm{Var}(\bar{y}_{1\cdot\cdot}-\bar{y}_{2\cdot\cdot}) &=\operatorname{Var}(\bar{w}_{1\cdot}-\bar{w}_{2\cdot}+\bar{\varepsilon}_{1\cdot\cdot}-\bar{\varepsilon}_{2\cdot\cdot}) \\
&=\frac2{rc}(\sigma_\varepsilon^2+c\sigma_w^2) 
\end{aligned}$$

$\mathrm{Var}(\bar{y}_{1\cdot\cdot}-\bar{y}_{2\cdot\cdot})$ 的一个估计为 $(2/rc) MSError(WholePlot)$.

量 $(\bar{y}_{11\cdot}-\bar{y}_{1\cdot\cdot})-(\bar{y}_{21\cdot}-\bar{y}_{2\cdot\cdot})$ 可用一般模型表示为

$$\begin{aligned}(\bar{y}_{11}.-\bar{y}_{1\cdot\cdot})-(\bar{y}_{21\cdot}-\bar{y}_{2\cdot\cdot})=&(\mu_{11}+\bar{b}_{\cdot}+\bar{w}_{1\cdot}+\bar{\varepsilon}_{11\cdot})-(\bar{\mu}_{1\cdot}+\bar{b}_{\cdot}+\bar{w}_{1\cdot}+\bar{\varepsilon}_{1\cdot\cdot})\\&-(\mu_{21}+\bar{b}_{\cdot}+\bar{w}_{2\cdot}+\bar{\varepsilon}_{21\cdot})+(\bar{\mu}_{2\cdot}+\bar{b}_{\cdot}+\bar{w}_{2\cdot}+\bar{\varepsilon}_{2\cdot\cdot})\\
=&[(\mu_{11}-\bar{\mu}_{1\cdot})-(\mu_{21}-\bar{\mu}_{2\cdot})]+[(\bar{\varepsilon}_{11\cdot}-\bar{\varepsilon}_{1\cdot\cdot})-(\bar{\varepsilon}_{21\cdot}-\bar{\varepsilon}_{2\cdot\cdot})]
\end{aligned}$$

方差为

$$\begin{aligned}
\mathrm{Var}[(\bar{y}_{11\cdot}-\bar{y}_{1\cdot\cdot})-(\bar{y}_{21\cdot}-\bar{y}_{2\cdot\cdot})]& =\mathrm{Var}[(\bar{{\varepsilon}}_{11{\cdot}}-\bar{{\varepsilon}}_{1{\cdot\cdot}})-(\bar{{\varepsilon}}_{21{\cdot}}-\bar{{\varepsilon}}_{2{\cdot\cdot}})]  \\
&=\frac{2(c-1)}{cr}\sigma_\varepsilon^2
\end{aligned}$$

$\mathrm{Var}[(\bar{y}_{11\cdot}-\bar{y}_{1\cdot\cdot})-(\bar{y}_{21\cdot}-\bar{y}_{2\cdot\cdot})]$ 的一个估计为

$$\widehat{\mathrm{Var}}[(\bar{y}_{11\cdot}-\bar{y}_{1\cdot\cdot})-(\bar{y}_{21\cdot}-\bar{y}_{2\cdot\cdot})]=\frac{2(c-1)}{rc}MSError(subplot)$$

将整区分量方差与子区分量方差相结合产生

$$\mathrm{Var}[\bar{y}_{11\cdot}-\bar{y}_{21\cdot}]=\left(\frac2{rc}\right)(\sigma_\varepsilon^2+c\sigma_w^2)+\frac{2(c-1)}{rc}\sigma_\varepsilon^2$$

$\bar{y}_{11\cdot}-\bar{y}_{21\cdot}$ 的标准误估计为

$$s.e.[\bar{y}_{11\cdot}-\bar{y}_{21\cdot}]=\sqrt{\left(\frac2{rc}\right)MSE(wholeplot)+\frac{2(c-1)}{rc}MSE(subplot)}$$

有关 $s.e.[\bar{y}_{11\cdot}-\bar{y}_{21\cdot}]$ 的自由度可用 Satterthwaite 近似获得

$$\begin{aligned}\hat{v}&=\frac{\left\{\left(\frac2{rc}\right)[MSE(wholeplot)]+\frac{2(c-1)}{rc}[MSE(subplot)]\right\}^2}{\frac{\left\{\left(\frac2{rc}\right)[MSE(wholeplot)]\right\}^2}{(a-1)(r-1)}+\frac{\left\{\frac{2(c-1)}{rc}[MSE(subplot)]\right\}^2}{a(r-1)(c-1)}}\end{aligned}$$

用一些简单的代数可以证明上面的 $s.e.[\bar{y}_{11\cdot}-\bar{y}_{21\cdot}]$ 与第 \@ref(sec24-2) 节中得到的结果相同。

### 通过一般对比进行比较 {#sec24-5}

可以为第 \@ref(sec24-2) 节中讨论的每个比较构建均值的一般对比。C 的水平或 $\bar\mu_{\cdot j}$ 之间的对比为

$$\theta=d_1\bar{\mu}_{\cdot1}+d_2\bar{\mu}_{\cdot2}+\cdots+d_c\bar{\mu}_{\cdot c}\mathrm{~where~}\sum_{j=1}^cd_j=0$$

该对比的一个估计为 $\hat\theta=d_1\bar{y}_{\cdot1\cdot}+d_2\bar{y}_{\cdot2\cdot}+\cdots+d_c\bar{y}_{\cdot c\cdot}$，方差为

$$\mathrm{Var}(\hat{\theta})=\frac{\sigma_\varepsilon^2}{ar}\sum_{j=1}^cd_j^2$$

$\hat\theta$ 的标准误估计为

$$\widehat{s.e.}(\hat{\theta})=\sqrt{\frac{MSError(subplot)}{ar}\sum_{j=1}^cd_j^2}$$

A 的水平或 ${\bar{\mu}}_{i\cdot}$ 的对比为 ${\tau}={h}_1{\bar{\mu}}_{1\cdot}+{h}_2{\bar{\mu}}_{2\cdot}+\cdots+{h}_a{\bar{\mu}}_{a\cdot}$，其中 $\sum_{i=1}^a{h}_i=0$；它的估计为 $\hat{\tau}={h}_1{\bar{y}}_{1\cdot\cdot}+{h}_2{\bar{y}}_{2\cdot\cdot}+\cdots+{h}_a{\bar{y}}_{a\cdot\cdot}$，方差为

$$\mathrm{Var}(\hat{\tau})=\left(\frac{\sigma_\varepsilon^2+c\sigma_w^2}{rc}\right)\sum_{i=1}^ah_i^2$$

$\hat{\tau}$ 的标准误估计为

$$\widehat{s.e.}(\hat{\tau})=\sqrt{\left(\frac{MSError(wholeplot)}{rc}\right)\sum_{i=1}^ah_i^2}$$

在整区处理相同水平上子区处理的对比，即，$\mu_{i1},\mu_{i2},\ldots,\mu_{ic}$ 的一个对比为 $\delta_i=s_1\mu_{i1}+s_2\mu_{i2}+\cdots+s_c\mu_{ic}$ 其中 $\Sigma_{j=1}^cs_j=0$. 它的估计为 $\hat\delta_i=s_1\bar y_{i1\cdot}+s_2\bar y_{i2\cdot}+\cdots+s_c\bar y_{ic\cdot}$，方差为

$$\mathrm{Var}(\hat{\delta}_i)=\frac{\sigma_\varepsilon^2}r\sum_{j=1}^cs_j^2$$

$\hat{\delta}_i$ 的标准误估计为

$$\widehat{s.e.}(\hat{\delta}_i)=\sqrt{\frac{MSError(subplot)}r\sum_{j=1}^cs_j^2}$$

在子区处理相同水平上整区处理的对比，即，${\mu}_{1{j}},{\mu}_{2{j}},\ldots,{\mu}_{{a}{j}}$ 的一个对比为 ${\lambda}_j={u}_1{\mu}_{1{j}}+{u}_2{\mu}_{2{j}}+\cdots+{u}_a{\mu}_{{a}{j}}$ 其中 $ {\Sigma}_{i=1}^a{u}_i=0$. 它的估计为 $\hat\lambda_j=u_1\bar y_{1j\cdot}+u_2\bar y_{2j\cdot}+\cdots+u_a\bar y_{aj\cdot}$，方差为

$$\mathrm{Var}(\hat{\lambda}_j)=\frac{(\sigma_\varepsilon^2+\sigma_w^2)}r\sum_{i=1}^au_i^2$$

$\hat{\lambda}_j$ 的标准误估计为

$$\widehat{s.e.}(\hat{\lambda}_j)=\sqrt{\frac{(\hat{\sigma}_\varepsilon^2+\hat{\sigma}_w^2)}r\sum_{i=1}^au_i^2}$$

其中

$$\hat{\sigma}_\varepsilon^2+\hat{\sigma}_w^2=\frac{MSError(wholeplot)+(c-1)MSError(subplot)}c$$

有关 $\hat{\sigma}_\varepsilon^2+\hat{\sigma}_w^2$ 的自由度是通过 Satterthwaite 近似（见第 \@ref(sec24-2) 节）获得的

$$\hat{v}=\frac{(\hat{\sigma}_\varepsilon^2+\hat{\sigma}_w^2)^2}{\frac{\left[\frac1cMSE(wholeplot)\right]^2}{(r-1)(a-1)}+\frac{\left[\frac{c-1}cMSE(subplot)\right]^2}{a(r-1)(c-1)}}$$

接下来考虑 $\mu_{ij}$ 的任何线性组合，如 $\sum_{i=1}^a\sum_{j=1}^cg_{ij}\mu_{ij}$. 根据第 \@ref(chap8) 章，我们知道如果对于每个 j 有 ${\sum}_{i=1}^{{a}}{g}_{ij}={0}$ 且对于每个 i 有 ${\sum}_{j=1}^{{c}}{g}_{ij}={0}$，则这样的对比是交互作用对比。可以证明

$$\begin{aligned}
\mathrm{Var}\left(\sum_{i=1}^a\sum_{j=1}^cg_{ij}\bar{y}_{ij\cdot}\right)& =\mathrm{Var}\left(\sum_{i=1}^a\sum_{j=1}^cg_{{ij}}\left(\bar{b}_\cdot+\bar{w}_{i\cdot}+\bar{\varepsilon}_{{ij\cdot}}\right)\right)=\mathrm{Var}\left[\sum_{i=1}^a\sum_{j=1}^c\left(g_{{ij}}\bar{b}_\cdot+g_{_{ij}}\bar{w}_{{i\cdot}}+g_{{ij}}\bar{\varepsilon}_{{ij\cdot}}\right)\right]  \\
&=\mathrm{Var}\left[\sum_{i=1}^a\sum_{j=1}^c(g_{ij}\bar{b}_\cdot)\right]+\mathrm{Var}\left[\sum_{i=1}^a\sum_{j=1}^c(g_{ij}\bar{w}_{i\cdot})\right]+\mathrm{Var}\left[\sum_{i=1}^a\sum_{j=1}^c(g_{ij}\bar{\varepsilon}_{ij.})\right] \\
&=(g_{\cdot\cdot})^2\frac{\sigma_B^2}r+\left[\sum_{i=1}^a(g_{i\cdot})^2\right]\frac{\sigma_W^2}r+\left[\sum_{i=1}^a\sum_{j=1}^c(g_{ij})^2\right]\frac{\sigma_\varepsilon^2}r
\end{aligned}$$

如果 $\sum_{i=1}^a\sum_{j=1}^cg_{ij}\mu_{ij}$ 是交互作用对比，则

$$\mathrm{Var}\Bigg(\sum_{i=1}^a\sum_{j=1}^cg_{ij}\bar{y}_{ij\cdot}\Bigg)=\Bigg[\sum_{i=1}^a\sum_{j=1}^c(g_{ij})^2\Bigg]\frac{\sigma_\varepsilon^2}r$$

其方差估计为

$$\mathrm{Var}\left(\sum_{i=1}^a\sum_{j=1}^cg_{ij}\bar{y}_{ij\cdot}\right)=\left[\sum_{i=1}^a\sum_{j=1}^c(g_{ij})^2\right]\frac{\hat{\sigma}_\varepsilon^2}r$$

相应的自由度为 a(r-1)(c-1).

接下来考虑由 $\bar y_{i\cdot}$ 定义的 $A_{i\cdot}$ 主效应均值。作为 $\sum_{i=1}^a\sum_{j=1}^cg_{ij}\bar{y}_{ij\cdot}$ 的一个特例，我们可以通过下式得到这个主效应均值

$$g_{ij}=\begin{cases}\frac1c\quad\mathrm{if~}i=i^{\prime}\quad\mathrm{~for~}j=1,2,...,c\\0\quad\mathrm{otherwise}\end{cases}$$

那么 $(g_{\cdot\cdot})^2=1$ 以及 $\sum_{i=1}^a(g_{i\cdot})^2=(g_{i^\prime \cdot})^2=1$，以及

$$\sum_{i=1}^a\sum_{j=1}^c(g_{ij})^2=\sum_{j=1}^c(g_{i'j})^2=\sum_{j=1}^c(\frac1c)^2=\frac1c$$

因此

$$\mathrm{Var}(\bar{y}_{i\cdot\cdot})=\frac{\sigma_B^2}r+\frac{\sigma_W^2}r+\frac{\sigma_\varepsilon^2}{rc}=\frac1{rc}(\sigma_\varepsilon^2+c\sigma_W^2+c\sigma_B^2)$$

注意到 A 主效应均值的方差估计取决于区组、整区和子区方差分量。类似地，可以证明 C 主效应均值 $\bar y_{\cdot j\cdot}$ 的方差估计为

$$\mathrm{Var}(\bar{y}_{\cdot j\cdot})=\frac{\sigma_B^2}r+\frac{\sigma_W^2}r+\frac{\sigma_\varepsilon^2}{ra}=\frac1{ra}(\sigma_\varepsilon^2+a\sigma_W^2+a\sigma_B^2)$$

最后，$\mu_{ij}$ 的估计为 $\bar y_{ij\cdot}$，其方差为

$$\mathrm{Var}(\bar{y}_{ij\cdot})=\frac{\sigma_B^2}r+\frac{\sigma_W^2}r+\frac{\sigma_\varepsilon^2}r=\frac1r(\sigma_\varepsilon^2+\sigma_W^2+\sigma_B^2)$$

注意，与这三个方差分量函数的每个估计相关的自由度将需要通过 Satterthwaite 法来估计。

许多研究人员使用效应模型而不是均值模型来描述他们的数据。重要的是能够将均值模型的对比表示为其在效应模型参数中的相应对比。例如，$\delta_i$ 的表达式为

$$\begin{aligned}
{\delta}_i& =\sum_{j=1}^cs_j\mu_{ij}  \\
&=\sum_{j=1}^cs_j[\mu+\alpha_i+\gamma_j+(\alpha\gamma)_{ij}] \\
&=\sum_{j=1}^cs_j\mu+\sum_{j=1}^cs_j\alpha_i+\sum_{j=1}^cs_j\gamma_j+\sum_{j=1}^cs_j(\alpha\gamma)_{ij} \\
&=\mu\sum_{j=1}^cs_j+\alpha_i\sum_{j=1}^cs_j+\sum_{j=1}^cs_j\gamma_j+\sum_{j=1}^cs_j(\alpha\gamma)_{ij} \\
&=\sum_{j=1}^cs_j\gamma_j+\sum_{j=1}^cs_j(\alpha\gamma)_{ij}\quad\mathrm{~since~}\sum_{j=1}^cs_j=0
\end{aligned}$$

因此，这种均值的对比涉及 $\gamma_j$ 和 $(\alpha\gamma)_{ij}$. 使用 SAS Mixed 进行对比的示例将在本章后面进行说明。

上面给出的估计和估计的标准误可以用于通过 $t$ 检验来检验假设，并/或构建关于固定效应参数可估函数的置信区间。感兴趣的对比可能有多种选择。例如，如果一个因素的水平是定量的，可以使用线性、二次等对比来研究均值参数的可能趋势。其他对比可能包括将一组对照与每种处理进行比较。

## 其他示例 {#sec24-6}

### 示例 24.3：水分和肥料 {#sec24-6-1}

<table>
<caption>(\#tab:table24-11)示例 [24.3](#sec24-6-1) 每盆中的干物质量的测量</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.11.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table24-11) 中的数据来自一个实验，该实验测量了在不同水分 (moisture) 水平和不同肥料 (fertilizer) 量下种植的小麦植株的干物质量 (amount of dry matter). 实验中有 48 个不同的泥炭盆 (peat pots) 和 12 个塑料托盘 (plastic trays)；每个托盘可以放入四个盆。水分处理包括每天向托盘中添加 10, 20, 30 或 40 毫升的水，水被泥炭盆吸收。水分水平随机分配给托盘。托盘是大型实验单位或整个地块，整个地块设计是完全随机设计结构中的单向处理结构（四个水分水平）。肥料水平为每盆 2, 4, 6 或 8 毫克。四个肥料水平随机分配给每个托盘中的四个盆，以便每个肥料水平在每个托盘中出现一次。盆是最小的实验单位或裂区或子区，子区设计是随机完全区组设计结构中的单向处理结构（四个肥料水平），其中 12 个托盘是区组。在每个盆中种植小麦种子，30 天后测量每个盆中生长的小麦植株的干物质量。可用于描述分配给第 i 个水分水平和第 j 个肥料水平的第 k 个托盘中盆的干物质量的模型是

$$y_{ijk}=\mu_{ij}+t_{ik}+p_{ijk},\quad i=1,2,3,4;~j=1,2,3,4;~k=1,2,3$$

$\mu_{ij}$ 是水分水平 i 与肥料水平 j 的平均干物质，$t_{ik}$ 是分布为 $i.i.d.N(0,\sigma^2_{\text{tray}})$ 的托盘误差项，$p_{ijk}$ 是分布为 $i.i.d.N(0,\sigma^2_{\text{plot}})$ 的盆误差项。请注意，$p_{ijk}$ 等价于该模型中的残差。方差分析表如表 \@ref(tab:table24-12) 所示，其中显示了显著的水分与肥料交互作用。由于水分和肥料之间存在交互作用，因此在进行推断时使用表\@ref(tab:table24-13) 中给出的处理组合平均值。由于水分水平和肥料水平是等距的定量水平，正交多项式可以很容易地用于研究每一水分水平的肥料水平趋势和每一肥料水平的水分水平趋势。测量第 i 个水分水平下肥料的线性和二次趋势的对比（相同整区处理下的子区处理的对比）为

$$\begin{aligned}\delta_{LinF|{M}_i}&=-3\mu_{i1}-\mu_{i2}+\mu_{i3}+3\mu_{i4},\quad &&i=1,2,3,4\\\delta_{QuadF|{M}_i}&=\mu_{i1}-\mu_{i2}-\mu_{i3}+\mu_{i4},\quad &&i=1,2,3,4\end{aligned}$$

它们的估计分别为

$$\begin{aligned}\hat{\delta}_{LinF|M_i}&=-3\bar{y}_{i\cdot1}-\bar{y}_{i\cdot2}+\bar{y}_{i\cdot3}+3\bar{y}_{i\cdot4},\quad &&i=1,2,3,4,\\\hat{\delta}_{QuadF|M_i}&=\bar{y}_{i\cdot1}-\bar{y}_{i\cdot2}-\bar{y}_{i\cdot3}+\bar{y}_{i\cdot4},\quad &&i=1,2,3,4\end{aligned}$$

这些对比的方差分别为

$$\begin{aligned}
\mathrm{Var}(\hat{\delta}_{LinF|M_i})& =\frac{\sigma_p^2}3[(-3)^2+(-1)^2+1^2+3^2]  \\
&=\frac{20{\sigma}_p^2}3 \\
\mathrm{Var}(\hat{\delta}_{QuadF|M_i})& =\frac{\sigma_p^2}3[(-1)^2+(-1)^2+1^2+1^2]  \\
&=\frac{4{\sigma}_p^2}3
\end{aligned}$$

<table>
<caption>(\#tab:table24-12)水分和肥料示例干物质量的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.12.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table24-13)肥料与水分的干物质量单元格均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.13.png">
</td>
  </tr>
</tbody>
</table>

标准误估计是通过将各自方差的平方根中的 $\sigma^2_p$ 替换为 $MSError(pot)$ 来获得的。表 \@ref(tab:table24-14) 中给出了评估这些对比的 SAS-Mixed 代码。表 \@ref(tab:table24-15) 中给出了每种水分水平的肥料线性趋势和二次趋势的估计以及相应的标准误估计和 $t$ 统计量（检验趋势是否为零）。

<table>
<caption>(\#tab:table24-14)带 estimate 语句的 SAS-Mixed 代码，用于评估线性和二次趋势</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.14.png">
</td>
  </tr>
</tbody>
</table>

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table24-15)每种水分水平下肥料线性和二次趋势的估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.15.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> 线性和二次趋势的标准误估计分别为 2.239 和 1.001. 将 $t$ 值与 $t_{lpha/2,24}$ 进行比较</td></tr></tfoot>
</table>

测量每个肥料水平下水分的线性和二次趋势的对比（相同子区处理下的整区处理的对比）为

$$\begin{aligned}\lambda_{LinM|F_j}&=-3\mu_{1j}-\mu_{2j}+\mu_{3j}+3\mu_{4j},\quad &&j=1,2,3,4\\
\lambda_{QuadM|F_j}&=\mu_{1j}-\mu_{2j}-\mu_{3j}+\mu_{4j},\quad &&j=1,2,3,4\end{aligned}$$

每种水分水平下肥料水平的估计和相应的方差为

$$\begin{aligned}\hat{\lambda}_{LinM|F_j}&=-3\bar{y}_{1j\cdot}-\bar{y}_{2j\cdot}+\bar{y}_{3j.}+3\bar{y}_{4j\cdot},\quad &&j=1,2,3,4\\\hat{\lambda}_{QuadM|F_k}&=\bar{y}_{1j\cdot}-\bar{y}_{2j\cdot}-\bar{y}_{3j\cdot}+\bar{y}_{4j\cdot},\quad &&j=1,2,3,4\end{aligned}$$

以及

$$\begin{aligned}
\mathrm{Var}(\hat{\lambda}_{LinM|F_k}) &=\frac{\sigma_{\mathrm{pot}}^2+\sigma_{\mathrm{tray}}^2}3(3^2+1^2+1^2+3^2) \\
&=\frac{20(\sigma_{\mathrm{pot}}^2+\sigma_{\mathrm{tray}}^2)}3 \\
\mathrm{Var}(\hat{\lambda}_{QuadM|F_k}) &=\frac{\sigma_{\mathrm{pot}}^2+\sigma_{\mathrm{tray}}^2}3(1^2+1^2+1^2+1^2) \\
&=\frac{4(\sigma_{\mathrm{pot}}^2+\sigma_{\mathrm{tray}}^2)}3 
\end{aligned}$$

在各个方差的平方根中将标准误估计 $\sigma_{\mathrm{pot}}^2+\sigma_{\mathrm{tray}}^2$ 替换为 $\hat{\sigma}_{\mathrm{pot}}^2+\hat{\sigma}_{\mathrm{tray}}^2=\frac{MSError(tray)+(4-1)MSError(pot)}4$. Satterthwaite 近似可用于获得与 $\hat{\sigma}_{\mathrm{pot}}^2+\hat{\sigma}_{\mathrm{tray}}^2$ 相关的近似自由度。对于这个例子

$$\hat{\sigma}_{\mathrm{pot}}^2+\hat{\sigma}_{\mathrm{tray}}^2=\frac{3.406+(4-1)0.752}4=1.416$$

基于 $\hat v$ 个自由度，其中

$$\hat{v}=\frac{[3.406+(4-1)0.752]^2}{\frac{[3.406]^2}8+\frac{[(4-1)0.752]^2}{24}}=19.3$$

表 \@ref(tab:table24-15) 给出了每种水分水平下肥料线性和二次趋势的估计以及相应的标准误估计和 $t$ 统计量。表 \@ref(tab:table24-15) 显示，每种水分水平对肥料都有线性的响应，没有显著的二次趋势。表 \@ref(tab:table24-16) 显示，每种肥料水平对水分都有线性和二次响应。图 \@ref(fig:figure24-5) 和 \@ref(fig:figure24-6) 中的图表显示了每种水分水平对肥料的响应以及每种肥料水平对水分的响应。使用 Satterthwaite 近似来确定表 \@ref(tab:table24-16) 中给出的近似自由度。

<div class="figure" style="text-align: center">
<img src="figure/figure 24.5.png" alt="各水分水平下干物质量均值关于肥料水平的图形" width="682" />
<p class="caption">(\#fig:figure24-5)各水分水平下干物质量均值关于肥料水平的图形</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 24.6.png" alt="各肥料水平下干物质量均值关于水分水平的图形" width="681" />
<p class="caption">(\#fig:figure24-6)各肥料水平下干物质量均值关于水分水平的图形</p>
</div>

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table24-16)每种肥料水平水分线性和二次趋势的估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.16.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> 线性和二次趋势的标准误估计分别为 3.072 和 1.374. 将 $t$ 值与 $t_{lpha/2,19.3}$ 进行比较</td></tr></tfoot>
</table>

### 示例 24.4：具有裂区误差的回归 {#sec24-6-2}

<table>
<caption>(\#tab:table24-17)完全回归模型的 SAS 混合代码和协方差参数估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.17.png">
</td>
  </tr>
</tbody>
</table>

示例 [24.3](#sec24-6-1) 中的水分和肥料水平是定量水平，因此可能有兴趣研究是否可以构建一个回归模型来描述数据。由于数据是相关的 (correlated)，因此这与通常的回归情况不同。需要混合模型软件来进行计算。为了获得适当的分析，需要将水分作为回归模型中的连续变量，并在 random 语句中作为 class 变量。Mst 用于表示作为连续变量的水分，mstblk（等于 mst）用于表示 class 变量。表 \@ref(tab:table24-17) 包含 SAS-Mixed 代码，其中 model 语句是水分（mst）和肥料（fr）的一般三次回归模型。random 语句纳入托盘误差项，并将相关结构强加给数据。表 \@ref(tab:table24-17) 给出了方差分量的估计，表 \@ref(tab:table24-18) 给出了回归系数的解。大多数显著性水平都相当大，大多数显著性水平都相当大，因此进行了几步删除过程（SAS-Mixed 中没有自动过程），直到所有剩余变量的系数都显著不同于零。缩减模型为

$$DM_{ijk}=\beta_0+\beta_1mst_i+\beta_2fr_j+\beta_3(mst_i)(fr_j)+\beta_4(mst_i)^2(fr_j)+\beta_5(mst_i)(fr_j)^2+t_{ik}+p_{ijk}$$

确定缩简模型是否充分描述了数据是令人感兴趣的，因此构造了一个欠拟合检验 (test for lack of fit). 设 frx=fr，并将 frx 包含在 class 语句中。在缩简模型中包括 mstblk × frx，如表 \@ref(tab:table24-19) 所示。表 \@ref(tab:table24-20) 的方差分析表中与 mstblk × frx 对应的 $F$ 检验提供了缩减模型的欠拟合检验。在这种情况下，与欠拟合检验相对应的显著性水平为 0.5561，表明缩减模型充分描述了数据。表 \@ref(tab:table24-21) 中的代码显示了缩减模型并显示了协方差参数的估计，回归系数的估计如表 \@ref(tab:table24-22) 所示。与完全模型相比，缩减模型的托盘误差项略小，而残差或盆方差相似。缩减模型的预测值如图\@ref(fig:figure24-7) 所示，水分与肥料单元格均值如图 \@ref(fig:figure24-8) 所示。正如从拟合度检验的结果中所预期的那样，回归模型充分描述了单元格均值。

<table>
<caption>(\#tab:table24-18)全响应面模型参数的解</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.18.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table24-19)用于缩减回归模型欠拟合检验的 SAS-Mixed 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.19.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table24-20)在 MSTBLK × FXR 中的欠拟合检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.20.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table24-21)最终回归模型的 SAS 混合代码和协方差参数估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.21.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table24-22)缩减回归模型参数的解</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.22.png">
</td>
  </tr>
</tbody>
</table>

<div class="figure" style="text-align: center">
<img src="figure/figure 24.7.png" alt="干物质量随水分和肥料变化的预测曲面" width="689" />
<p class="caption">(\#fig:figure24-7)干物质量随水分和肥料变化的预测曲面</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 24.8.png" alt="肥料和水分组合的均值" width="678" />
<p class="caption">(\#fig:figure24-8)肥料和水分组合的均值</p>
</div>

### 示例 24.5：混乱的裂区设计 {#sec24-6-3}

<table>
<caption>(\#tab:table24-23)分配给每个区组内每个整区的品种与除草剂的组合</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.23.png">
</td>
  </tr>
</tbody>
</table>

两名研究人员设计了一项研究，评估三种大豆品种 (variates) 对四种不同类型除草剂 (herbicides) 的响应。在 1 号试验点，第一名研究人员在两次重复中使用了裂区试验设计，其中大豆品种为整区的水平，除草剂为子区的水平。在 2 号试验点，第二名研究人员在两次重复中使用了裂区试验设计，但将除草剂水平分配给了整区，将品种水平分配给子区。希望进行组合分析 (combined analysis )，但由于两种设计非常不同，因此似乎不可能进行组合分析。然而，请记住，从第 \@ref(chap5) 章中可以看出，裂区设计结构无非是一个不完全区组设计，其中整区是不完全区组设计的区组。14 个整区或不完全区组中的每一个整区的处理组合如表 \@ref(tab:table24-23) 所示，其中 Vx × y 表示含有除草剂 y 的品种 x. 前八个整区（来自 1 号试验点）是尺寸为三的不完全区组，后六个整区是尺寸为 4 的不完全区组。使用如下模型可以完成组合分析

$$y_{ijkl}=\mu_{kl}+b_i+w_{ij}+\varepsilon_{ijkl}$$

其中 $\mu_{kl}$ 是第 k 个品种和第 l 个除草剂的平均响应。$b_i\thicksim i.i.d.~N(0,\sigma_{blk}^2)$ 每个试验点内的大区效应 (large block effects)，$w_{ij}\thicksim i.i.d.~N(0,~\sigma_{wp}^2)$ 为整区或不完全区组效应，以及 ${\varepsilon_{ijkl}}\thicksim{i.i.d.N(0,\sigma_\varepsilon^2)}$ 为子区效应。数据见表 \@ref(tab:table24-24) 和 \@ref(tab:table24-25)，SAS-Mixed 代码和协方差参数估计见表 \@ref(tab:table24-26). model 语句包含了品种（V）、除草剂（H）和品种与除草剂交互作用的处理结构。random 语句包含区组（试验点内的重复）和区组内的整区（不完全区组）。固定效应检验如表 \@ref(tab:table24-27) 所示，均值如表 \@ref(tab:table24-28) 所示。大豆品种和除草剂之间存在显著的交互作用，因此需要使用单元格均值进行比较。对于每个品种内除草剂水平的比较，差异标准误估计为 3.33，自由度为 26.1. 对于每种除草剂内品种的比较，差异标准误估计为 3.26，自由度为 24.7. 

<table>
<caption>(\#tab:table24-24)示例 [24.5](#sec24-6-3) 中来自 1 号试验点的数据。其中 H 和 V 分别表示除草剂和品种，表中的条目是大豆的重量</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.24.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table24-25)示例 [24.5](#sec24-6-3) 中来自 2 号试验点的数据。其中 H 和 V 分别表示除草剂和品种，表中的条目是大豆的重量</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.25.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table24-26)示例 [24.5](#sec24-6-3) 的 SAS-Mixed 代码和协方差参数估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.26.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table24-27)示例 [24.5](#sec24-6-3) 固定效应检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.27.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table24-28)示例 [24.5](#sec24-6-3) 品种与除草剂单元格均值、品种均值和除草剂均值的估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.28.png">
</td>
  </tr>
</tbody>
</table>

### 示例 24.6：裂-裂区设计 {#sec24-6-4}

裂-裂区设计 (split-split-plot design) 有三种尺寸的实验单元，即整区 (whole plot)、子区 (subplot) 和次子区 (subsubplot)，其中整区、子区和次子区分别是最大、中等和最小尺寸的实验单元。本研究旨在评估两种饲料 (rations)（普通玉米和高油玉米）、两种温度（3°C 和 6°C）和三种包装类型（真空、二氧化碳和低氧）的组合对肉嫩度的影响。二十头阉公牛 (steers) 随机分配了两种饲料（每种饲料 10 头）。屠宰时，动物被一分为二，从每侧 (side) 提取一块腰肉。腰肉被随机分配到两种储存温度中的一种。储存 10 天后，从每块腰肉上切下三块牛排，并将包装 (packaging) 水平分配给其中一块牛排。牛排被放入展示柜中再放置 7 天。置于展示柜 7 天后，从每块牛排上获得五个芯 (cores)，并测量剪切每个芯所需的力。表 \@ref(tab:table24-29) 中的数据是五个芯的平均值。由于一些包装问题，导致一些牛排不能代表真实状态，因此它们被从数据集中删除。

<table>
<caption>(\#tab:table24-29)示例 [24.6](#sec24-6-4) 剪切力数据</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.29.png">
</td>
  </tr>
</tbody>
</table>

本研究涉及三种尺寸的实验单位。动物是被分配饲料水平的实验单元。动物的侧是温度水平的实验单元，牛排是包装水平的实验单元。可以用来描述剪切力 (shear forces) ($f_{ijkl}$)的模型是

$$f_{ijkl}=\mu_{ikl}+a_{ij}+s_{ijk}+\varepsilon_{ijkl}\quad i=1,2,j=1,2,\ldots,10,k=1,2,l=1,2,3$$

其中 $\mu_{ikl}$ 表示第 i 种饲料、第 k 种温度、第 l 种包装的平均响应，$a_{ij}\thicksim i.i.d.~N(0,\sigma_{\mathrm{animal}}^2)$ 表示动物误差，$s_{ijk}\thicksim i.i.d.~N(0,\sigma_{\mathrm{side}}^2)$ 表示侧误差，${\varepsilon_{ijkl}}\thicksim i.i.d.~N(0,\sigma_{\mathrm{steak}}^2)$ 表示牛排误差。

动物误差项可以通过忽略处理结构中的温度和包装部分来获得。动物设计是一个完全随机设计结构中的单向处理结构，其中每个饲料分配了 10 只动物。因此，动物误差项是通过在相同饲料内经相同处理的动物之间的变异来计算的，其中合并了所有饲料 (the animal error term is computed as the variation among animals treated alike within a ration pooled across rations)，记作 animal(ration). 侧误差项是通过仅考虑分配了饲料 1 的动物并忽略包装效应来获得的。在饲料 1 内的每只动物针对两种温度形成尺寸为 2（两侧）的区组，因此侧设计是随机完全区组设计结构中（饲料 1 内的 10 只动物）的单向处理结构（两种温度）。因此，侧误差项是通过在同一饲料内的温度与动物交互作用来计算的，其中合并了所有饲料 (side error term is computed as the temperature by animal interaction within a ration pooled across rations)，记作 temp × animal(ration). 残差是牛排误差项，但其形式可以通过考虑来自饲料 1 和温度为 3°C 的数据来获得。来自饲料 1 并分配了 3°C 的动物的侧是三种包装类型的区组因子。牛排设计是随机完全设计结构中（分配了饲料 1 并在 3°C 下储存的动物的 10 个侧面）的单向处理结构（三个包装水平）。包装与动物交互作用是这部分数据的误差项。残差的计算形式是 package × animal(temp ration). 如果数据集是均衡的，则动物误差、侧误差和牛排误差项的自由度分别为 18, 18, 72.

为剪切力数据拟合裂-裂区模型的 SAS-Mixed 代码在表 \@ref(tab:table24-30) 中，其中 random 语句包含 animal(ration) 和temp × animal(ration)，分别表示整区误差和子区误差。协方差参数的估计也包含在表 \@ref(tab:table24-30) 中。固定效应检验在表 \@ref(tab:table24-31) 中给出，其中存在显著的 ration × temperature 和 ration × packaging交互作用。ration × temperature 的平均值在表 \@ref(tab:table24-32) 中给出，它们之间的成对比较在表 \@ref(tab:table24-33) 中给出。在 lsmeans 语句中使用了 `adjust = Tukey` 选项，以提供 Tukey 多重比较方法的调整 $p$ 值。比较饲料内温度的标准误估计约为 0.45，而比较温度内饲料的标准误估计约为 0.679. 比较饲料内的温度是动物内比较 (within-animal comparisons)，而比较温度内的饲料是动物间比较 (between-animal comparisons)，因此具有较大的标准误估计。表 \@ref(tab:table24-34) 和 \@ref(tab:table24-35) 包含 ration × packaging 的均值和均值的成对比较。比较饲料内包装的标准误估计约为 0.49，而比较包装内饲料的标准误估计约为 0.71. 包装的比较是在动物的一侧进行的，而饲料的比较是在动物之间进行的，从而解释了两个标准误估计的大小差异。

<table>
<caption>(\#tab:table24-30)示例 [24.6](#sec24-6-4) 的 SAS-Mixed 代码和协方差参数估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.30.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table24-31)示例 [24.6](#sec24-6-4) 固定效应检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.31.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table24-32)示例 [24.6](#sec24-6-4) 饲料与温度的力均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.32.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table24-33)饲料与温度的力均值的成对比较，其中 Tukey 调整用于多重比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.33.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table24-34)示例 [24.6](#sec24-6-4) 饲料与包装的力均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.34.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table24-35)饲料与包装的力均值的成对比较，其中 Tukey 调整用于多重比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.35.png">
</td>
  </tr>
</tbody>
</table>

## 样本量和功效考虑 {#sec24-7}

特定研究所需的样本量取决于最有趣的比较。用于检测两均值之差为 $\delta$ 的基本样本量方程，其中差异的方差为 $2\sigma^2 /n$，$\hat\sigma^2$ 是基于 v 个自由度的 $\sigma^2$ 的估计，第一类错误率为 $\alpha$，第二类错误率为 $\beta$，其公式为

$$n=\frac{2\hat{\sigma}^{2}}{\delta^{2}}[t_{\alpha/2,v}+t_{\beta,v}]^{2}$$

样本量方程可用于评估检验的功效，方法是求解 $t_{\beta,v}$，然后确定功效为 $1-\beta$ 的值。$t_{\beta,v}$ 的方程为

$$t_{\beta,v}=\sqrt{\frac{n\delta^2}{2\hat{\sigma}^2}}-t_{\alpha/2,v}$$

对于示例 [24.6](#sec24-6-4)，假设感兴趣的是检测在某一温度水平下，饲料均值之间是否存在 1.5 个力单位的差异。表示第一个温度下两种饲料差异的参数是 $\bar{{\mu}}_{11\cdot}-\bar{{\mu}}_{21\cdot}$. 如果数据集是均衡的，则 $\bar{{\mu}}_{11\cdot}-\bar{{\mu}}_{21\cdot}$ 的最佳估计为 $\bar{{f}}_{1\cdot1\cdot}-\bar{{f}}_{2\cdot2\cdot}$. $\bar{{f}}_{i\cdot k}$ 的模型为 $\bar{f}_{i\cdot k\cdot}=\bar{\mu}_{ik\cdot}+\bar{a}_{i\cdot}+\bar{s}_{i\cdot k}+\bar{\varepsilon}_{i\cdot k\cdot}$ 以及 $\bar{{f}}_{1\cdot1\cdot}-\bar{{f}}_{2\cdot2\cdot}$ 的模型为 $\bar{{\mu}}_{11{\cdot}}-\bar{{\mu}}_{21{\cdot}}+\bar{{a}}_{1{\cdot}}-\bar{{a}}_{2{\cdot}}+\bar{s}_{1\cdot1}-\bar{s}_{2\cdot1}+\bar{\varepsilon}_{1\cdot1\cdot}-\bar{\varepsilon}_{2\cdot1\cdot}$

所需的样本量是每种饲料的动物数量，因为侧的数量以及每侧的牛排数量分别固定为 2 和 3. $\bar{{f}}_{1\cdot1\cdot}-\bar{{f}}_{2\cdot2\cdot}$ 的方差为

$$\begin{aligned}
\mathrm{Var}(\bar{f}_{1{\cdot}1{\cdot}}-\bar{f}_{2{\cdot}1{\cdot}})& =\frac{2\sigma_{\mathrm{animal}}^2}n+\frac{2\sigma_{\mathrm{side}}^2}n+\frac{2\sigma_{\mathrm{steak}}^2}{3n}  \\
&=\frac2{3n}(\sigma_{\mathrm{steak}}^2+3\sigma_{\mathrm{side}}^2+3\sigma_{\mathrm{animal}}^2)
\end{aligned}$$

令第一类错误和第二类错误分别为 $\alpha,\beta$，检测出两均值之差为 $\delta$ 所需的每种饲料的动物数量为

$$n=\frac23(\hat{\sigma}_{\mathrm{steak}}^2+3\hat{\sigma}_{\mathrm{side}}^2+3\hat{\sigma}_{\mathrm{animal}}^2)[t_{\alpha/2,\hat{v}}+t_{\beta,\hat{v}}]^2/\delta^2$$

其中

$$\hat{{\sigma}}_{\mathrm{steak~}}^2+3\hat{{\sigma}}_{\mathrm{side~}}^2+3\hat{{\sigma}}_{\mathrm{animal}}^2$$

是 ${{\sigma}}_{\mathrm{steak~}}^2+3{{\sigma}}_{\mathrm{side~}}^2+3{{\sigma}}_{\mathrm{animal}}^2$ 的估计，基于 $\hat v$ 个自由度。使用均衡数据集的三个误差项的期望均方为

$$\begin{aligned}
E[MSAnimal(Ration)]&=\sigma_{\text{steak}}^2+3\sigma_{\text{side}}^2+6\sigma_{\text{animal}}^2\\
E[MSTemp\times Animal(Ratiоn)]&=\sigma_\mathrm{steak}^2+3\sigma_\mathrm{side}^2 \\
E[MSPack\times Animal(Теmp* Ratiоn)]&=\sigma_{\mathrm{steak}}^2 
\end{aligned}$$

${{\sigma}}_{\mathrm{steak~}}^2+3{{\sigma}}_{\mathrm{side~}}^2+3{{\sigma}}_{\mathrm{animal}}^2$ 的一个无偏估计为

$$\hat{\sigma}_{\mathrm{steak}}^2+3\hat{\sigma}_{\mathrm{side}}^2+3\hat{\sigma}_{\mathrm{animal}}^2=\frac12MSAnimal(Ration)+\frac12MSTemp\times Animal(Ration)$$

其近似自由度可通过 Satterthwaite 法获得

$$\hat{v}=\frac{(\hat{\sigma}_{\mathrm{steak}}^2+3\hat{\sigma}_{\mathrm{side}}^2+3\hat{\sigma}_{\mathrm{animal}}^2)^2}{\frac{[\frac12\textit{MSAnimal}(Ration)]^2}{18}+\frac{[\frac12\textit{MSTemp}\times Animal(Ration)]^2}{18}}$$

当数据集不均衡时，可以使用来自 REML 的方差分量及其相应渐近协方差阵的估计来评估近似自由度。令 $f(\sigma_\varepsilon^2,\sigma_\mathrm{side}^2,\sigma_\mathrm{animal}^2)$ 表示要待估的方差分量函数。$f(\sigma_\varepsilon^2,\sigma_\mathrm{side}^2,\sigma_\mathrm{animal}^2)$ 的估计为 $f(\hat\sigma_\varepsilon^2,\hat\sigma_\mathrm{side}^2,\hat\sigma_\mathrm{animal}^2)$，以及近似方差为

$$\mathrm{Var}\Bigl[f(\hat{\sigma}_\varepsilon^2,\hat{\sigma}_\mathrm{side}^2,\hat{\sigma}_\mathrm{animal}^2)\Bigr]=\begin{bmatrix}\frac{\partial f}{\partial\sigma_\varepsilon^2}\quad\frac{\partial f}{\partial\sigma_\mathrm{side}^2}\quad\frac{\partial f}{\partial\sigma_\mathrm{animal}^2}\end{bmatrix}\mathrm{Var}\Bigl[\hat{\sigma}_\varepsilon^2,\hat{\sigma}_\mathrm{side}^2,\hat{\sigma}_\mathrm{animal}^2\Bigr]\begin{bmatrix}\frac{\partial f}{\partial\sigma_\varepsilon^2}\\\frac{\partial f}{\partial\sigma_\mathrm{side}^2}\\\frac{\partial f}{\partial\sigma_\mathrm{animal}^2}\end{bmatrix}$$

在方差分量的 REML 估计处评估 $\mathrm{Var}\Bigl[f(\hat{\sigma}_\varepsilon^2,\hat{\sigma}_\mathrm{side}^2,\hat{\sigma}_\mathrm{animal}^2)\Bigr]$，并计算

$$\begin{aligned}Z=\frac{f(\hat{\sigma}_\varepsilon^2,\hat{\sigma}_\mathrm{side}^2,\hat{\sigma}_\mathrm{animal}^2)}{\sqrt{\widehat{\text{Var}}\left[f(\hat{\sigma}_\varepsilon^2,\hat{\sigma}_\mathrm{side}^2,\hat{\sigma}_\mathrm{animal}^2)\right]}}\end{aligned}$$

则与 $f(\hat\sigma_\varepsilon^2,\hat\sigma_\mathrm{side}^2,\hat\sigma_\mathrm{animal}^2)$ 相关的 Satterthwaite 近似自由度为 $\hat v=2Z^2$.

对于示例 [24.6](#sec24-6-4)，使用分析结果并确定检测 1.5 剪切单位差异所需的动物数量，第一类和第二类错误率分别为 0.05 和 0.05. 表 \@ref(tab:table24-30) 中方差分量的估计为 $\hat{\sigma}_{\mathrm{animal}}^2=1.2434,\hat{\sigma}_{\mathrm{side}}^2=0.1992,\hat{\sigma}_{\mathrm{steak}}^2=1.7529$. ${\sigma}_\mathrm{steak}^2+3{\sigma}_\mathrm{side}^2+3{\sigma}_{\mathrm{animal}}^2$ 的估计为 6.08086. 方差分量渐近协方差阵的估计如表 \@ref(tab:table24-36) 所示，为

$$\widehat{\text{Var}}\begin{bmatrix}\hat{{\sigma}}_\mathrm{arimal}^2\\\hat{{\sigma}}_\mathrm{side}^2\\\hat{{\sigma}}_\mathrm{steak}^2\end{bmatrix}=\begin{bmatrix}0.3723&-0.06154&0.002209\\-0.06154&0.1305&-0.05370\\0.002209&-0.05370&0.1296\end{bmatrix}=\hat{\boldsymbol V}$$

<table>
<caption>(\#tab:table24-36)示例 [24.6](#sec24-6-4) 方差分量估计的渐近协方差阵</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.36.png">
</td>
  </tr>
</tbody>
</table>

则 ${\hat\sigma}_\mathrm{steak}^2+3{\hat\sigma}_\mathrm{side}^2+3{\hat\sigma}_{\mathrm{animal}}^2$ 的近似方差为

$$\widehat{\mathrm{Var}}(\hat{\sigma}_{\mathrm{steak}}^2+3\hat{\sigma}_{\mathrm{side}}^2+3\hat{\sigma}_{\mathrm{animal}}^2)=\begin{bmatrix}3&3&1\end{bmatrix}\hat{\boldsymbol V}\begin{bmatrix}3\\3\\1\end{bmatrix}=3.23805$$

以及 $Z$ 分数为

$$Z=\frac{\hat{\sigma}_\mathrm{steak}^2+3\hat{\sigma}_\mathrm{side}^2+3\hat{\sigma}_\mathrm{animal}^2}{\sqrt{\widehat{\mathrm{Var}}(\hat{\sigma}_\mathrm{steak}^2+3\hat{\sigma}_\mathrm{side}^2+3\hat{\sigma}_\mathrm{animal}^2)}}=3.379$$

因此与 $\hat{\sigma}_{\mathrm{steak}}^2+3\hat{\sigma}_{\mathrm{side}}^2+3\hat{\sigma}_{\mathrm{animal}}^2$ 相关的自由度为 $2(3.379)^2=22.8$. 样本量估计为每种饲料下的动物 $n=\frac23(6.08086)[t_{0.025,22.8}+t_{0.05,22.8}]^2/1.5^2=25.8$ 或 $n=26$.

[^solvingbeta]: 译者注：通过最左边和最右边解出 $\beta$，例如在 R 中通过命令 `pt(0.2869,22.8)` 求解。

该研究每种饲料下有 10 只动物，在同一温度下的饲料间，检验 1.5 个单位的剪切力差异的功效通过下式计算[^solvingbeta]

$$t_{\beta,\hat{v}}=\sqrt{\frac{10\delta^2}{3(\hat{\sigma}_{\mathrm{seak}}^2+3\hat{\sigma}_{\mathrm{side}}^2+3\hat{\sigma}_{\mathrm{animal}}^2)}}-t_{\alpha/2,\hat{v}}=\sqrt{\frac{10(1.5)^2}{3(6.08086)}}-2.069=0.2869$$

当 $\alpha=0.05$ 时 $\beta=0.388$，提供了 $1-\beta=0.612$ 的功效；即，每种饲料下有 10 只动物，当两种饲料均值实际相差 1.5 个单位时，有 61% 的机会确定这两种饲料有显著差异。

[^secondproblem]: 原文：The second sample size problem is to determine the number of animals required to detect a difference of 1.5 shear force units between two temperature means at the same level of packaging averaged over the two levels of ration.

第二个样本量问题是，在相同包装下对两种饲料水平进行平均，确定检测出两个温度均值之间 1.5 个剪切力单位的差异所需的动物数量[^secondproblem]。感兴趣的比较为 $\bar{{\mu}}_{\cdot11}-\bar{{\mu}}_{\cdot21}$. 评估差值方差所需的模型是 $\bar{f}_{\cdot\cdot kl}=\bar{\mu}_{\cdot kl}+\bar{a}_{\cdot\cdot}+{\bar{s}}_{\cdot\cdot k}+{\bar{\varepsilon}}_{\cdot\cdot jk}$. 差值的方差为

$$\begin{aligned}
\mathrm{Var}(\bar{f}_{\cdot\cdot11}-\bar{f}_{\cdot\cdot21}) =\text{ Var}(\bar{s}_{{\cdot\cdot}1}-\bar{s}_{\cdot\cdot2}+\bar{{\varepsilon}}_{\cdot\cdot11}-\bar{{\varepsilon}}_{\cdot\cdot21}) \\
=\frac{2{\sigma_\mathrm{side}}^2}{2n}+\frac{2{\sigma_\mathrm{steak}}^2}{2n} \\
=\frac2{2n}(\sigma_{\mathrm{side}}^2+\sigma_{\mathrm{steak}}^2) 
\end{aligned}$$

每种饲料的动物数量为 $n=\frac22(\hat{\sigma}_{\mathrm{steak}}^2+\sigma_{\mathrm{side}}^2)[t_{0.025,49.9}+t_{0.05,49.9}]^2/1.5^2=11.8$，其中使用上述程序将自由度确定为 49.9. 

存在多种可能重要的比较，因此需要为每个比较确定样本量，所需动物数量或整区数量可以基于每个被认为重要的比较所需动物数量的最大值来确定。

## 使用 JMP 进行计算：示例 24.7 {#sec24-8}

<div class="figure" style="text-align: center">
<img src="figure/figure 24.9.png" alt="示例 [24.7](#sec24-8) 前 32 个观测数据" width="614" />
<p class="caption">(\#fig:figure24-9)示例 [24.7](#sec24-8) 前 32 个观测数据</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 24.10.png" alt="示例 [24.7](#sec24-8) 后 16 个观测数据" width="774" />
<p class="caption">(\#fig:figure24-10)示例 [24.7](#sec24-8) 后 16 个观测数据</p>
</div>

图 \@ref(fig:figure24-9) 和 \@ref(fig:figure24-10) 中的数据集来自 JMP 表，其中包含名义变量 block, variety, fert 和 rate。实验设计包括四个整区的四个区组，其中两种小麦品种 (varieties) 和两个肥料 (fertilizer) 水平的组合被随机分配给整区。每个整区被分成三个子区，三个播撒率 (seeding rate) 水平被随机分配给三个子区。整区设计是随机完全区组设计结构（四个区组）中的双向处理结构（品种与肥料水平）。由于不可控的昆虫破坏，存在一些缺失的数据点，但以下对结构的评估是在假设没有缺失数据点的情况下进行的。整区设计的方差分析表由表 \@ref(tab:table24-37) 给出。整区误差是根据处理结构与设计结构交互作用计算得出的。在这种情况下，存在一个具有四个处理组合的双向处理结构，因此整区误差由 (4-1)(4-1) = 9 个自由度组成。接下来，比较品种一和肥料 1 的播撒率水平。所得设计是在随机完全区组设计结构中的单向处理结构，方差分析表在表 \@ref(tab:table24-38) 中给出。子区误差自由度是根据区组与速率交互作用计算得出的，提供 6 个自由度。品种和肥料有四种组合，每种组合为子区误差提供 6 个自由度。合并这四组 6 个自由度为子区误差提供 (4)(6) = 24 个自由度。可用于描述该数据的模型是

$$y_{ijkl}=\mu_{ijl}+b_k+zv_{ijk}+\varepsilon_{ijkl}\quad i=1,2,j=1,2,k=1,2,3,4,l=1,2,3$$

其中  $b_k\thicksim i.i.d.N(0,\sigma_{blk}^2),w_{ijk}\thicksim i.i.d.N(0,\sigma_{ivp}^2),\varepsilon_{ijkl}\thicksim i.i.d.N(0,\sigma_{sp}^2)$ 以及 $\mu_{ijl}$ 表示品种 i, 肥料 j 和播撒率 l 的平均值。

<table>
<caption>(\#tab:table24-37)示例 [24.7](#sec24-8) 的整区分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.37.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table24-38)比较品种 1 和肥料 1 的播种率的方差分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.38.png">
</td>
  </tr>
</tbody>
</table>

[^missingdata]: 原文：Missing data cause the variances of comparisons to be different, which in turn can provide multiple comparisons that are not completely ordered as demonstrated above.

包含两个误差项和期望均方的方差分析表在表 \@ref(tab:table24-39) 中显示。JMP fit model 界面如图 \@ref(fig:figure24-11) 所示。通过使用 attributes 按钮，将 block 和 block × variety × fert 项定义为 random 效应。默认的估计方法是 REML，但可以通过更改 method 按钮来使用 III 型平方和。单击 run model 按钮以获取分析。图 \@ref(fig:figure24-12) 包含方差分量的估计和固定效应的分析。与 SAS-Mixed 不同，关于方差分量的置信区间是使用 Waldner 区间而不是 Satterthwaite 近似区间计算的。所有涉及率的效应都是显著的，包括三向交互作用，因此需要解决 variety × fert × rate 三向均值。点击 effects 按钮可提供每个效应的最小二乘均值表，如图 \@ref(fig:figure24-13) - \@ref(fig:figure24-15) 所示。least square means 选项如图 \@ref(fig:figure24-16) 所示。tables 选项是默认选项。plots 选项提供平均值的折线图。contrasts 选项可以构建感兴趣的最小二乘均值的对比。Tukey HSD 选项为多重比较提供 Tukey 调整，slice 选项为交互中每个效应的每个水平提供等均值检验。图 \@ref(fig:figure24-17) 中的结果是 Tukey HSD 选项结果的一部分，提供了线条或字母。由于研究中有缺失数据，这些字母可能会导致一些矛盾。例如，均值 2, 2, 10 与均值 2, 2, 6 显著不同（见 A 和 B），但 2, 2, 10 与 1, 2, 4 没有显著不同（见 C），即使 2, 2, 6 的均值大于 1, 2, 4 的均值。缺失数据会导致比较的方差不同，这反过来又会提供不完全排序的多重比较，如上所示[^missingdata]。

<table>
<caption>(\#tab:table24-39)示例 [24.7](#sec24-8) 的方差分析表（假设没有缺失数据）</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2024.39.png">
</td>
  </tr>
</tbody>
</table>

<div class="figure" style="text-align: center">
<img src="figure/figure 24.11.png" alt="示例 [24.7](#sec24-8) fit model 界面" width="975" />
<p class="caption">(\#fig:figure24-11)示例 [24.7](#sec24-8) fit model 界面</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 24.12.png" alt="示例 [24.7](#sec24-8) 方差分量估计 (REML) 和固定效应检验" width="967" />
<p class="caption">(\#fig:figure24-12)示例 [24.7](#sec24-8) 方差分量估计 (REML) 和固定效应检验</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 24.13.png" alt="示例 [24.7](#sec24-8) 品种、肥料及其交互作用的最小二乘均值" width="530" />
<p class="caption">(\#fig:figure24-13)示例 [24.7](#sec24-8) 品种、肥料及其交互作用的最小二乘均值</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 24.14.png" alt="率及其与品种和肥料交互作用的最小二乘均值" width="540" />
<p class="caption">(\#fig:figure24-14)率及其与品种和肥料交互作用的最小二乘均值</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 24.15.png" alt="品种与肥料与率的交互作用的最小二乘均值" width="540" />
<p class="caption">(\#fig:figure24-15)品种与肥料与率的交互作用的最小二乘均值</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 24.16.png" alt="效应的 least square means 选项" width="532" />
<p class="caption">(\#fig:figure24-16)效应的 least square means 选项</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 24.17.png" alt="三项交互作用的最小二乘均值，通过 Tukey's method，字母指示了分组" width="620" />
<p class="caption">(\#fig:figure24-17)三项交互作用的最小二乘均值，通过 Tukey's method，字母指示了分组</p>
</div>

图 \@ref(fig:figure24-18) 包含将 Tukey HSD 选项应用于 Variety × Fert 均值的完整结果集。该表包含差异估计、每个差异的标准误估计以及每个差异 95% 置信区间的上限和下限。当置信区间不包括零时，结果显示为红色。所有这些置信区间都包括零，因此它们以黑色显示。可以通过在图 \@ref(fig:figure24-19) 中显示的对比表上为每个均值选择系数来评估最小二乘平均值的对比。

<div class="figure" style="text-align: center">
<img src="figure/figure 24.18.png" alt="示例 [24.7](#sec24-8) 品种与肥料均值的多重比较，Tukey's method 的计算" width="640" />
<p class="caption">(\#fig:figure24-18)示例 [24.7](#sec24-8) 品种与肥料均值的多重比较，Tukey's method 的计算</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 24.19.png" alt="最小二乘均值的 contrast 选项窗口" width="498" />
<p class="caption">(\#fig:figure24-19)最小二乘均值的 contrast 选项窗口</p>
</div>

此示例显示了可用于使用 JMP 软件分析裂区型设计的一些选项。本章中描述的所有数据集和模型都可以使用 JMP 的 fit model 界面进行分析。

## 结束语 {#sec24-9}

在本章中，详细讨论了裂区设计，并使用 SAS Mixed 和 JMP 进行了分析。所讨论设计的方差分析表涉及多个误差项，并且描述了能够构造每个误差项的适当表示的方法。各种估计的标准误估计可能涉及多个误差项。提出了计算实验者可能感兴趣的各种类型比较的标准误估计的方法。标准误估计也可用于进行多重比较。当因子的水平是定量的时，描述了估计线性、二次和类似趋势的方法。演示了样本量和功效方法。详细讨论了七个示例。

## 练习 {#sec24-10}
