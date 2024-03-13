
# 嵌套设计的分析 {#chap30}

> "All models are wrong, but some are useful." - George E. P. Box

嵌套效应 (nesting effects) 在设计的实验中可以只发生在设计结构或处理结构中，也可以同时发生。为了在设计结构中进行嵌套，必须有多个尺寸的实验单元，其中一个小的实验单元嵌套在一个大的实验单元中。裂区、重复测量和分层设计是设计结构中存在嵌套的示例。当存在两个或多个因素时，可以在处理结构中嵌套。这些因素可能都是固定效应，都是随机效应，或两者的混合。因此，在处理结构中嵌套的实验可以被建模为固定、随机或混合模型。第 \@ref(chap5) 章介绍了设计结构中嵌套因子的概念。本章给出了一些例子，展示了模型构建、参数估计和假设检验。在设计结构中具有嵌套的嵌套设计通常被称为分层设计 (hierarchical designs).

## 定义，假设和模型 {#sec30-1}

在处理结构中，如果因子 B 的每个水平仅与因子 A 的一个水平一起出现，则因子 B 的水平将嵌套在因子 A 的水平内。以下示例演示了处理结构中的嵌套。

### 示例 30.1：公司和杀虫剂 {#sec30-1-1}

四家化学公司生产特定的杀虫剂 (insecticides).  A 公司生产三种这样的产品，B 公司和 C 公司各生产两种这样的产品，D 公司生产四种这样的产品。没有一家公司生产的产品与另一家公司的产品完全相同。处理结构是双向的，一是公司，二是产品。这种处理结构如表 \@ref(tab:table30-1) 所示，其中每个产品的水平在每个公司的水平内仅出现一次。因此，产品水平嵌套在公司水平内。

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table30-1)公司和杀虫剂示例的处理结构</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2030.1.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> X 表示特定产品来自相应的公司</td></tr></tfoot>
</table>

处理结构中两个因素的水平都是固定效应。为了进行实验，我们将一盒装有活蓝草植物 (live bluegrass plants) 的土壤和 400 只蚊子分别放入 33 个玻璃容器中。然后将三个玻璃容器随机分配给每种产品。用该产品处理玻璃容器，4 小时后统计活蚊子的数量。

这个实验的设计结构是完全随机的。可以用来描述每个容器中活蚊子数量的模型是

\begin{equation}
y_{ijk}=\mu+\gamma_i+\rho_{j(i)}+\varepsilon_{ijk},\quad i=1,2,3,4,\quad j=2,3,\mathrm{~or~}4,\quad k=1,2,3
(#eq:30-1)
\end{equation}

其中 $y_{ijk}$ 是第 i 个公司的第 j 个产品的第 k 次重复中观察到的蚊子数量，$\mu$ 是总体平均值，$\gamma_i$ 是第 i 家公司的效应，$\rho_{j(i)}$ 是第 i 公司的第 j 个产品的效应，${\varepsilon}_{ijk}\sim i.i.d.N(0,{\sigma}_\varepsilon^2)$ 表示与测量 $y_{ijk}$ 相关的误差。该模型只有一个尺寸的实验单元，因此只有一个误差项。实验的设计是完全随机设计结构中的双向嵌套处理结构（两个因素都是固定效应）。要估计的模型参数是 $\sigma^2_\varepsilon$ 以及 $\mu,\gamma_i,\rho_{j(i)}$ 的可估函数。

本例数据见表 \@ref(tab:table30-2)，分析见第 \@ref(sec30-2) 节和第 \@ref(sec30-3) 节。

<table>
<caption>(\#tab:table30-2)示例 [30.1](#sec30-1-1) 数据</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2030.2.png">
</td>
  </tr>
</tbody>
</table>

### 示例 30.2：舒适度实验回顾 {#sec30-1-2}

示例 [5.7](#sec5-5-2) 中的舒适度实验是设计结构中嵌套的一个例子。大实验单位是环境室，小实验单位是人。模型是

\begin{equation}
y_{ijkm}=\mu_{ik}+c_{j(i)}+p_{m(ijk)}\quad i=65,70,75,~j=1,2,\ldots,9,~k=M,F,~m=1,2
(#eq:30-2)
\end{equation}

其中 $\mu_{ik}$ 表示第 i 个温度和第 k 个性别的均值，$c_{j(i)}$ 表示分配给第 i 个温度的第 j 个室的随机效应且假设 $c_{j(i)}\thicksim N(0,\sigma_{\text{Сhamber}} ^ 2 )$，$p_{m(ijk)}$ 表示分配到第 i 个温度的第 j 个室的第 k 个性别的第 m 个人的随机效应且假设 $p_{m(ijk)}\thicksim N(0,\sigma_{\text{Person}} ^ 2 )$. 项 $c_{j(i)}$ 和 $p_{m(ijk)}$ 表示在第 j 个室（室嵌套在温度内）和第 i 个温度的第 j 个房间内第 k 个性别的人 m（人嵌套在室内）上的观测温度 i 的误差。由于实验单元有两种尺寸，因此分析有两个水平和两个误差项。需要估计的参数是 $\mu_{ik},i=65,70,75,k=M,F,\sigma_{\text{Chamber'}} ^ 2,\sigma_{\text{Person}} ^ 2$. 该示例的分析将在第 \@ref(sec30-2) 节和第 \@ref(sec30-3) 节中讨论。

### 示例 30.3：咖啡价格示例回顾 {#sec30-1-3}

第 \@ref(chap18) 章中的咖啡价格示例是采用三向处理结构的抽样调查设计，其中商店水平嵌套在城市水平内，而城市水平又嵌套在州水平内（所有三个因素都是随机效应）。第 \@ref(chap18) 章中用于描述咖啡价格变异的模型是

$$y_{ijk}=\mu+s_i+c_{j(i)}+a_{k(ij)},\quad i=1,2,\ldots,r,j=1,2,\ldots,t_i,\quad k=1,2,\ldots,n_{ij}$$

其中 $\mu$ 为美国咖啡价格均值，$s_i\thicksim i.i.d.~N(0\text{,}\sigma_{\mathrm{State}}^2),c_{j(i)}\thicksim i.i.d.~N(0\text{,}\sigma_{\mathrm{City}}^2),a_{k(ij)}\thicksim i.i.d.~N(0\text{,}\sigma_{\mathrm{Store}}^2)$ 分别表示随机州、城市和商店效应。感兴趣的参数是方差分量 $\sigma_{\text{State}}^2,\sigma_{\text{City}}^2,\sigma_{\text{Store}}^2$ 和总均值 $\mu$. 咖啡价格示例是多阶段抽样实验的示例，这是嵌套设计的非常常见的应用。该示例的分析将在接下来的两节中讨论。表 \@ref(tab:table30-3) 包含一项小型研究的咖啡价格，该研究用于演示参数估计、置信区间估计和假设检验。

<table>
<caption>(\#tab:table30-3)示例 [30.3](#sec30-1-3) 以美分/磅为单位的咖啡价格</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2030.3.png">
</td>
  </tr>
</tbody>
</table>

上述示例演示了可能发生在处理结构、设计结构或两者中的嵌套。处理结构可以涉及随机和/或固定效应，并且设计结构可以涉及多种尺寸的实验单元。因此，此类设计的分析涉及使用固定效应模型、随机效应模型和混合效应模型技术，如下两节所述。

## 参数估计 {#sec30-2}

涉及嵌套的模型属于已经讨论过的模型类别之一。也就是说，如果存在固定效应，则需要估计均值；如果存在随机效应，则需要估计方差分量；如果有多种尺寸的实验单元，则分析涉及多个误差项。因此，设计可能非常简单，也可能非常复杂。本章中的示例用于演示前面讨论的技术在涉及嵌套的设计分析中的应用。

### 示例 30.1：继续 {#sec30-2-1}

模型 \@ref(eq:30-1) 的可估计函数涉及 $\mu + \gamma_i + \rho_{j(i)}$ 的线性组合。例如，每个 i 的 $\rho_{j(i)}$ 的对比是可估的，其中 $\rho_{1(i)}-\rho_{2(i)}$ 是同一公司内产品的对比，在这种情况下，对比是对第 i 个公司的产品 1 和 2 进行比较。$\mu + \gamma_i + \bar \rho_{\cdot(i)}$ 的对比用于比较各公司产品的均值，其中 $\gamma_1 + \bar\rho_{\cdot(1)} - \gamma_2 + \bar\rho_{\cdot(2)}$ 将公司 1 与公司 2 进行比较。$\gamma_1 + \rho_{1(1)} - \gamma_2 + \rho_{1(2)}$ 将公司 1 的产品 1 与公司 2 的产品 1 进行比较。这些对比的估计分别为

$$\begin{aligned}
\widehat{\rho_{1(i)}-\rho_{2(i)}}&=\bar y_{i1\cdot}-\bar y_{i2\cdot}\\
\widehat{\gamma_1 + \bar\rho_{\cdot(1)} - \gamma_2 + \bar\rho_{\cdot(2)}}-\bar y_{i2\cdot}&=\bar y_{1\cdot\cdot}-\bar y_{2\cdot\cdot}\\
\widehat{\gamma_1 + \rho_{1(1)} - \gamma_2 + \rho_{1(2)}}&=\bar y_{11\cdot}-\bar y_{21\cdot}
\end{aligned}$$

模型 \@ref(eq:30-1) 也可以表示为

$$y_{ijk}=\mu_{j(i)}+\varepsilon_{ijk},\quad i=1,2,3,4,\quad j=2,3,\mathrm{~or~}4,\quad k=1,2,3$$

其中 $\mu_{j(i)}$ 为公司 i 的产品 j. $\mu_{j(i)}$ 的估计为 $\bar y_{ij\cdot}$。任何包含 $\mu_{j(i)}$ 对比的线性组合都是可估的，但在解释所选对比时必须小心。$\sigma^2_\varepsilon$ 的估计来自于将 j(i) 组合内观测值的方差合并为

$$\hat{\sigma}_\varepsilon^2=\frac1{N-q}\sum_{i=1}^m\sum_{j=1}^{p_i}\sum_{k=1}^{n_{j(i)}}(y_{ijk}-\bar{y}_{ij\cdot})^2$$

其中 $n_{j(i)}$ 表示 j(i) 单元格内的观测数，q 为 $n_{j(i)}>0$ 的个数，m 为公司数，$p_i$ 为来自第 i 个公司的产品数，$n_{j(i)}$ 为分配给第 i 个公司的第 j 个产品的玻璃容器的数量，$N=\sum_{j(i)}n_{j(i)}$ 为观测总数。

对于这个例子

$$\begin{aligned}{\hat{\sigma}}_\varepsilon^2=\frac{1}{33-11}\sum_{i=1}^4\sum_{j=1}^{p_i}\sum_{k=1}^3(y_{ijk}-\bar{y}_{ij\cdot})^2=60.818\end{aligned}$$

样本均值的估计见表 \@ref(tab:table30-4).

<table>
<caption>(\#tab:table30-4)示例 [30.1](#sec30-1-1) 均值估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2030.4.png">
</td>
  </tr>
</tbody>
</table>

$\mu_{j(i)}$ 之间可以进行多重比较。可以检验关于 $\mu_{j(i)}$ 的其他假设，例如比较公司或公司内的产品。例如，可以通过考虑比较公司 B 和公司 D

$$\bar{\mu}_{\cdot(2)}=\bar{\mu}_{\cdot(4)}\quad\mathrm{~where~}\bar{\mu}_{\cdot(2)}=\frac{\mu_{1(2)}+\mu_{2(2)}}2\quad\mathrm{and~}\quad\bar{\mu}_{\cdot(4)}=\frac{\mu_{1(4)}+\mu_{2(4)}+\mu_{3(4)}+\mu_{4(4)}}4$$

也就是说，比较将在 B 公司的两种产品的均值与 D 公司的四种产品的均值之间进行。研究人员需要决定这种比较是否有意义。第 \@ref(sec30-3) 节讨论了这种嵌套处理结构的假设检验。

如果 j(i) 单元格中的观测数量不等，则可以使用分析不均衡模型的技术来获得 $\mu_{j(i)}$ 的估计。总体边际均值的估计提供了 $\hat\mu_{j(i)}$ 的估计，而 $\hat\sigma^2_\varepsilon$ 的估计是通过合并处理组合之间的方差而获得的，该方差可以通过方差分析获得。

### 示例 30.2：继续 {#sec30-2-2}

表 \@ref(tab:table30-5) 给出了舒适度研究的数据。嵌套出现在设计结构中，人嵌套在室内，而室嵌套在温度内。带期望均方的方差分析表展示在表 \@ref(tab:table30-6) 中。两个方差分量的矩法估计（也是 REML 和 MINQUE0 估计）为 $\hat{{\sigma}}_{\mathrm{Person}}^2=1.65,\hat{{\sigma}}_{\mathrm{Chamber}}^2=2.36$，ML 估计为 $\hat{{\sigma}}_{\mathrm{Person}}^2=1.47,\hat{{\sigma}}_{\mathrm{Chamber}}^2=1.48$。用于裂区实验的方法（第 \@ref(chap24) 章中描述）可用于比较 $mu_{ij}$。如果设计不均衡，使用 REML 进行混合模型分析来估计方差分量是合适的。

<table>
<caption>(\#tab:table30-5)示例 [30.2](#sec30-2-2) 的数据，其中值为舒适度得分，其中 1 = 冷，8 = 舒适，15 = 热</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2030.5.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table30-6)示例 [30.2](#sec30-2-2) 舒适度研究的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2030.6.png">
</td>
  </tr>
</tbody>
</table>

### 示例 30.3：继续 {#sec30-2-3}

咖啡价格研究涉及随机效应模型，其感兴趣的参数是 $\mu,\sigma_{\text{State}}^2,\sigma_{\text{City}}^2,\sigma_{\text{Store}}^2$. ML 将提供参数的分析和估计，就像 MINQUE0 技术一样。采用 I 型平方和的矩法分析通常用于此类多水平抽样设计。 REML 法的结果是优选的。表 \@ref(tab:table30-7) 包含 I 型平方和、它们的期望均方以及用于一般研究的矩法估计。表 \@ref(tab:table30-8) 包含表 \@ref(tab:table30-3) 中数据的 I 型分析。表 \@ref(tab:table30-9) 使用 I 型平方和, REML, MIVQUE0 和 ML 的矩法来估计方差分量。方差分量的 95% 置信区间是使用具有指定自由度的卡方分布计算的。Wald 法用于矩法提供的区间，但使用自由度为 $2(Z\text{-value})^2$ 的卡方分布重新计算它们。表 \@ref(tab:table30-10) 给出了方差分量的每种估计方法对平均咖啡价格的估计。

<table>
<caption>(\#tab:table30-7)示例 [30.2](#sec30-2-2) I 型平方和以及矩法估计的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2030.7.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table30-8)咖啡价格数据的 I 型平方和、均方、期望均方和检验统计量</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2030.8.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table30-9)咖啡数据方差分量的四种方法的估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2030.9.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table30-10)咖啡价格均值方差分量的四种方法的估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2030.10.png">
</td>
  </tr>
</tbody>
</table>

## 假设检验和置信区间的构建 {#sec30-3}

对于嵌套在一个或两个结构中的均衡设计，可以使用均衡随机效应模型、均衡混合效应模型和均衡裂区重复测量模型的方法来构建置信区间并进行参数之间的比较。对于不均衡嵌套设计，估计方差分量的 REML 法提供了一种使用 Satterthwaite 近似构建置信区间的好方法。在这种情况下，近似自由度计算为 $2(Z\text{-value})^2$. 使用 1 型平方和估计、REML 估计、ML 估计和 MINQUE0 估计为咖啡价格数据构建的置信区间显示在表 \@ref(tab:table30-9) 中。矩法可以有效地用于提供有关方差分量的假设检验。表 \@ref(tab:table30-8) 使用预期均方来检验 $H_0\colon{\sigma}_\mathrm{State}^2=0\mathrm{vs~}H_a\colon{\sigma}_\mathrm{State}^2>0$ 和 $H_0{:{\sigma_{\mathrm{City}}^2}}=0\operatorname{vs}H_a{:{\sigma_{\mathrm{City}}^2}}>0$，分别得到 2.40 和 11.58 的 F 统计量。

当嵌套在处理结构中且涉及固定效应因素时，采用分析固定效应模型的方法。主要变化是可以用方差分析表测试的假设类型，如示例 [30.1](#sec30-1-1) 的数据所示

### 示例 30.1：继续 {#sec30-3-1}

杀虫剂数据模型的编写有两种方法

\begin{equation}
y_{ijk}=\mu_{j(i)}+\varepsilon_{ijk}\quad i=1,2,3,4,\quad j=2,3,\mathrm{~or~}4,\quad k=1,2,3
(#eq:30-3)
\end{equation}

或

\begin{equation}
y_{ijk}=\mu+\gamma_i+\rho_{j(i)}+\varepsilon_{ijk}\quad i=1,2,3,4,\quad j=2,3,\mathrm{~or~}4,\quad k=1,2,3
(#eq:30-4)
\end{equation}

表 \@ref(tab:table30-11) 包含模型 \@ref(eq:30-3) 的方差分析表。因产品导致的平方和检验的假设是

$$\mu_{1(1)}=\mu_{2(1)}=\mu_{3(1)}=\mu_{1(2)}=\mu_{2(2)}=\mu_{1(3)}=\mu_{2(3)}=\mu_{1(4)}=\mu_{2(4)}=\mu_{3(4)}=\mu_{4(4)}$$
<table>
<caption>(\#tab:table30-11)示例 30.1 基于模型 \@ref(eq:30-3) 的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2030.11.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table30-12) 包含模型 \@ref(eq:30-4) 的方差分析表。表 \@ref(tab:table30-11) 中产品的平方和已划分为公司的平方和以及公司内嵌套的产品的平方和的和。由于双向处理结构是嵌套的（产品嵌套在公司内），因此没有衡量产品水平和公司水平之间的交互作用。公司的平方和用于检验假设 $\bar{\mu}_{\cdot(1)}=\bar{\mu}_{\cdot(2)}=\bar{\mu}_{\cdot(3)}=\bar{\mu}_{\cdot(4)}$，其中 $$\bar{{\mu}}_{{\cdot}(i)}=(1/m_i)\sum_{j=1}^{m_i}{\mu}_{j(i)}={\mu}+\gamma_i+{\bar{\rho}}_{{\cdot}(i)}$$

<table>
<caption>(\#tab:table30-12)示例 30.1 基于模型 \@ref(eq:30-4) 的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2030.12.png">
</td>
  </tr>
</tbody>
</table>

公司内部的产品的平方和检验了以下假设（以均值模型参数表示）

$$\mu_{1(1)}=\mu_{2(1)}=\mu_{3(1)},\quad\mu_{1(2)}=\mu_{2(2)}=\mu_{1(3)}=\mu_{2(3)},\quad\mathrm{and}\quad\mu_{1(4)}=\mu_{2(4)}=\mu_{3(4)}=\mu_{4(4)}$$

或以效应模型参数表示

$$\rho_{1(1)}=\rho_{2(1)}=\rho_{3(1)},\quad\rho_{1(2)}=\rho_{2(2)}=\rho_{1(3)}=\rho_{2(3)},\quad\mathrm{and}\quad\rho_{1(4)}=\rho_{2(4)}=\rho_{3(4)}=\rho_{4(4)}$$

任何适当的多重比较程序（见第 \@ref(chap3) 章）都可以用于以下情况：1）比较 $\bar{{\mu}}_{\cdot (i)}$，2）比较 $\bar{{\mu}}_{j (i)}$、3）比较每个 i 的 $\bar{{\mu}}_{j (i)}$，或 4）比较上述任何组合。对于不均衡模型，III 型或 IV 型平方和适用于检验各种假设，并且总体边际均值的估计提供了 $\bar{{\mu}}_{\cdot (i)}$ 和 $\bar{{\mu}}_{j (i)}$ 的估计。

对于其他嵌套模型，方差分析表中的预期均方可以用作构建关于方差分量的适当 $F$ 检验的准则，并且混合模型 “固定效应的 III 型检验” 可以用于检验关于固定效应的各种假设。第 \@ref(chap22) 章和第 \@ref(chap23) 章中的方法将在这种情况下使用。

## 使用 JMP 进行分析 {#sec30-4}

JMP^®^ 软件用于提供示例 30.2 和 30.3 的分析。图 \@ref(fig:figure30-1) 包含数据集的部分列表，其中 chamber, temperature, gender 和 person 被定义为名义变量。图 \@ref(fig:figure30-2) 中的 fit model 界面包括 temperature, gender, temperature × gender 作为固定效应，chamber(temperature) 作为随机效应。person 误差项为模型的残差。图 \@ref(fig:figure30-3) 包含方差分量的 REML 估计，以及室误差的 Wald 置信区间，残差或人变异的卡方置信区间，固定效应的 III 型检验。Tukey 多重比较的最小二乘均值如图 \@ref(fig:figure30-4) 所示。图 \@ref(fig:figure30-5) 中的数据界面包含示例 30.3 的咖啡价格的部分列表。图 \@ref(fig:figure30-6) 中的 fit model 界面将 state 和 city(state) 作为随机效应。store 间的变异是通过残差来衡量的。图 \@ref(fig:figure30-7) 包含方差分量的估计以及总均值的估计。还给出了 state 和 city(state) 方差分量的 Wald 置信区间，为残差或 store 间变异提供了卡方置信区间。总均值的标准误估计为 3.359，而使用 REML 的 SAS^®^-Mixed 的总均值标准误估计为 3.2935（见表 \@ref(tab:table30-10)）。

<div class="figure" style="text-align: center">
<img src="figure/figure 30.1.png" alt="示例 30.2 舒适度数据的 JMP 界面" width="761" />
<p class="caption">(\#fig:figure30-1)示例 30.2 舒适度数据的 JMP 界面</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 30.2.png" alt="示例 30.2 舒适度数据的 JMP fit model 界面" width="838" />
<p class="caption">(\#fig:figure30-2)示例 30.2 舒适度数据的 JMP fit model 界面</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 30.3.png" alt="示例 30.2 的随机和固定效应的 JMP AOV 表" width="850" />
<p class="caption">(\#fig:figure30-3)示例 30.2 的随机和固定效应的 JMP AOV 表</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 30.4.png" alt="示例 30.2 中使用 Tukey 多重比较的最小二乘均值" width="476" />
<p class="caption">(\#fig:figure30-4)示例 30.2 中使用 Tukey 多重比较的最小二乘均值</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 30.5.png" alt="示例 30.3 咖啡数据的 JMP 数据界面" width="650" />
<p class="caption">(\#fig:figure30-5)示例 30.3 咖啡数据的 JMP 数据界面</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 30.6.png" alt="示例 30.3 JMP fit model 界面" width="850" />
<p class="caption">(\#fig:figure30-6)示例 30.3 JMP fit model 界面</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 30.7.png" alt="示例 30.3 JMP 方差分量的 REML 估计" width="768" />
<p class="caption">(\#fig:figure30-7)示例 30.3 JMP 方差分量的 REML 估计</p>
</div>

## 结束语 {#sec30-5}

本章定义了嵌套的概念，并通过实例证明嵌套可以发生在处理结构、设计结构或两者中。嵌套模型可以包括固定效应、随机效应或混合效应模型，对实例的分析阐述了所涉及的一些问题。
