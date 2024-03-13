
# 分析条区型设计的方法 {#chap25}

> "Statistics are the heart of democracy." - Simeon Strunsky

**条区** (strip-plot) 设计结构常用于实验单元以矩形排列的情况，例如田间、布料上或一摞笼子中。在第 \@ref(chap5) 章中，条区设计被描述为基本设计结构之一，其分析细节在本章中进行讨论。当一个因素的水平应用于每个矩形的行 (row)，而另一个因素的水平应用于每个矩形的列 (column) 时，就形成了条区设计结构。因此，当处理方法不允许对实验单元进行单独处理，而是按组（如整行或整列的单元）进行处理时，这种设计非常有用。条区设计结构也适用于在一系列步骤中处理实验单元的实验，其中在每个步骤中一起处理一组实验单元。本章讨论了六个示例，以演示复杂设计结构的分析方法。第一个涉及灌溉和氮的水平的示例用于展示最简单条区设计的基本分析。其余大部分示例涉及条区设计与裂区设计结构的不同组合。最后一个示例涉及条-条区 (strip-strip-plot) 设计结构。

## 条区设计和模型的描述 {#sec25-1}

条区设计结构类似于裂区设计结构，但实验单元构建方式不同。条区设计至少包含双向处理结构，其中基本实验单元以矩形集合的形式安排。每个矩形具有 a 行，其中 a 是第一因素（A）的水平数；同时具有 c 列，其中 c 是第二因素（C）的水平数。第一因素 A 的各个水平随机分配至矩形的各行，使得同一行内的所有实验单元接受相同的 A 水平，并一起进行处理。第二因素 C 的各个水平则随机分配至矩形的各列，使得同一列内的所有实验单元接受相同的 C 水平，并一起进行处理。图 \@ref(fig:figure25-1) 展示了将 A 的各个水平分配给行以及将 C 的各个水平分配给列的示意图，在该图中共有 r 个矩形，代表被分区的重复 (blocked replicates). 这一分配过程生成了一种具有三种尺寸实验单元的设计。行是因子 A 的实验单元，列是因子 C 的实验单元，而单元格则是 A × C 交互作用的实验单元。

<div class="figure" style="text-align: center">
<img src="figure/figure 25.1.png" alt="两个因子的水平分配给排列成矩形的实验单元的条（行和列）" width="629" />
<p class="caption">(\#fig:figure25-1)两个因子的水平分配给排列成矩形的实验单元的条（行和列）</p>
</div>

可视化适当分析的一个简单方法是首先忽略因子 C 的各个水平，只观察那些将因子 A 的水平分配给各行的矩形，如图 \@ref(fig:figure25-2) 所示。行设计由随机完全区组设计结构中的单向处理结构组成。可用于描述行设计数据的模型是

$$y_{ik}=\mu_i^r+b_k^r+e_{ik}^r,\quad i=1,2,\ldots,a,\quad k=1,2,\ldots,r$$

其中区组效应 $b_k^r\thicksim i.i.d.N(0,\sigma_{\mathrm{block}}^2)$，行效应 $e_{ik}^r\thicksim i.i.d.~N(0,\sigma_{\mathrm{row}}^2)$. 

<div class="figure" style="text-align: center">
<img src="figure/figure 25.2.png" alt="将 A 的水平分配给每个矩形的行" width="469" />
<p class="caption">(\#fig:figure25-2)将 A 的水平分配给每个矩形的行</p>
</div>

表 \@ref(tab:table25-1) 给出了行设计的方差分析表，其中行误差项计算为 rectangle × A 交互作用。

<table>
<caption>(\#tab:table25-1)条区设计结构行部分的方差分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.1.png">
</td>
  </tr>
</tbody>
</table>

接下来，忽略因子 A 的各个水平，仅关注那些将因子 C 的水平分配给各列的矩形，如图 \@ref(fig:figure25-3) 所示。列设计包含了在一个随机完全区组设计结构中的单向处理结构。可用于描述列设计数据的模型是

$$y_{jk}=\mu_j^c+b_k^c+e_{jk}^c\quad i=1,2,\ldots,c,k=1,2,\ldots,r$$

<div class="figure" style="text-align: center">
<img src="figure/figure 25.3.png" alt="将 C 的水平分配给每个矩形的列" width="470" />
<p class="caption">(\#fig:figure25-3)将 C 的水平分配给每个矩形的列</p>
</div>

[^acinteraction]: 原文：Thus the size of the experimental unit on which the A × C interaction is measured is the intersection of a row and a column which corresponds to a cell within the rectangle.

其中区组效应 $b_k^c\thicksim i.i.d.N(0,\sigma_{\mathrm{block}}^2)$，行效应 $e_{ik}^c\thicksim i.i.d.~N(0,\sigma_{\mathrm{column}}^2)$. 列设计的方差分析表如表 \@ref(tab:table25-2) 所示，其中列误差项计算为  rectangle × C 交互作用。因子 A 水平的比较是行间比较 (between-row comparisons)，因子 C 水平的比较则是列间比较 (between-column comparisons)，但 A × C 交互作用的测量是行内比较与列内比较。因此，测量 A × C 交互作用的实验单元的尺寸是与矩形内的单元相对应的行和列的交点[^acinteraction]。

<table>
<caption>(\#tab:table25-2)条区设计结构列部分的方差分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.2.png">
</td>
  </tr>
</tbody>
</table>

描述处理结构中具有两个因素的条区设计结构的模型是

$$y_{ijk}=\mu+b_k+\alpha_i+r_{ik}+\gamma_j+c_{jk}+(\alpha\gamma)_{ij}+\varepsilon_{ijk},\quad i=1,2,\ldots,a,\quad j=1,2,\ldots,c,\quad k=1,2,\ldots,r$$

其中 $\mu$ 表示总均值效应，$\alpha_i$ 表示 A 的第 i 个水平的效应，$\gamma_j$ 表示 C 的第 j 个水平的效应，$(\alpha\gamma)_{ij}$ 表示因子 A 和 C 的交互效应，$b_k$ 表示第 k 个区组的效应并假定 $b_k\thicksim i.i.d.N(0,\sigma_{\mathrm{block}}^2)$，$r_{ik}$ 表示第 k 个区组内第 i 行的效应并假定 $r_{ik}\thicksim i.i.d.N(0,\sigma_{\mathrm{row}}^2)$，$c_{jk}$ 表示第 k 个区组内第 j 列的效应并假定 $c_{jk}\thicksim i.i.d.N(0,\sigma_{\mathrm{column}}^2)$，以及 $\varepsilon_{ijk}$ 表示第 k 个区组第 ij 个单元格的随机效应并假定 ${\varepsilon}_{ijk}\thicksim i.i.d.~N(0,{\sigma}_{\mathrm{cell}}^2)$. 此外，假设所有随机项都是独立分布的。

<table>
<caption>(\#tab:table25-3)双向处理结构条区设计结构的方差表分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.3.png">
</td>
  </tr>
</tbody>
</table>

行误差被计算为 A 与区组交互作用，列误差项被计算为 C 与区组交互作用，单元格误差项被计算为 A 与 C 与区组的三向交互作用。方差分析表如表 \@ref(tab:table25-3) 所示，它是通过将表 \@ref(tab:table25-1) 和 \@ref(tab:table25-2) 中的条目与 AC交互作用和单元格误差项相结合而构建的。条区模型可以表示为均值模型，如下所示

$$y_{ijk}=\mu_{ij}+b_k+r_{ik}+c_{jk}+\varepsilon_{ijk},\quad i=1,2,\ldots,a,j=1,2,\ldots,c,\quad k=1,2,\ldots,r$$

其中 $\mu_{ij}=\mu+\alpha_i+\gamma_j+(\alpha\gamma)_{ij}$.

对于均衡的数据集，单元格均值的估计为 $\hat{\mu}_{ij}=\bar{y}_{ij\cdot}$，涉及 $\mu_{ij}$ 的各种对比、总和或均值是通过取 $\hat \mu_{ij}$ 的相应对比、总和或均值来估计的。方差分量的矩法解是通过将表 \@ref(tab:table25-3) 中的均方与预期均方相等而生成的方程组获得的

$$\begin{aligned}
\tilde{\sigma}_{\mathrm{cell}}^2& =\text{MSError(cell)}  \\
\tilde{{\sigma}}_{\mathrm{row}}^2& =\frac{MSError(row)-MSError(cell)}c  \\
\tilde{{\sigma}}_{\mathrm{column}}^2& =\frac{MSError(column)-MSError(cell)}a  \\
\tilde{{\sigma}}_{\mathrm{block}}^2& =\frac{MSBlock-MSError(column)-MSError(row)+MSError(cell)}{ac} 
\end{aligned}$$

方差分量的估计为

$$\begin{aligned}
\hat{{\sigma}}_{\mathrm{cell}}^2&= \tilde{{\sigma}}_{\mathrm{cell}}^2  \\
\hat{\sigma}_{\mathrm{row}}^2&= \begin{cases}\tilde{{\sigma}}_\mathrm{row}^2&&\mathrm{if~}\tilde{{\sigma}}_\mathrm{row}^2>0\\0&&\mathrm{if~}\tilde{{\sigma}}_\mathrm{row}^2\leq0\end{cases}  \\
\hat{{\sigma}}_{\mathrm{column}}^2&= \begin{cases}\tilde{{\sigma}}_\mathrm{column}^2&&\mathrm{if~}\tilde{{\sigma}}_\mathrm{column}^2>0\\0&&\mathrm{if~}\tilde{{\sigma}}_\mathrm{column}^2\leq0\end{cases}  \\
\hat{{\sigma}}_{\mathrm{block}}^2&= \begin{cases}{\tilde{\sigma}}_\mathrm{block}^2&&\mathrm{if~}{\tilde{\sigma}}_\mathrm{block}^2>0\\0&&\mathrm{if~}{\tilde{\sigma}}_\mathrm{block}^2\leq0\end{cases} 
\end{aligned}$$

在对该模型的参数进行推断时必须小心，因为该模型涉及三个误差项，每个实验单元尺寸对应一个误差项。下一节将讨论对假设的必要检验和在均值之间进行比较的适当标准误。

## 推断技术 {#sec25-2}

表 \@ref(tab:table25-3) 中的期望均方列提供了构建检验统计量所需的信息，以推断因子的效应。为检验无交互假设 $H_{0AC}\colon\mu_{ij}-\mu_{ij}-\mu_{ij}+\mu_{ij}=0\text{ for all }i\neq i'\text{ and }j\neq j'\text{ vs }H_{aAC}:\text{ (not }H_{0AC})$，将 A × C 均方除以 MSError(cell). 为检验因子 A 的等主效应均值假设 $H_{0A}\colon\bar{{\mu}}_{1\cdot}=\bar{{\mu}}_{2\cdot}=\cdots=\bar{{\mu}}_{a\cdot}\mathrm{~vs~}H_{aA}{:}(\text{nоt }H_{0A})$，将因子 A 的均方除以 MSError(row). 为检验因子 C 的等主效应均值假设 $H_{0C}\colon\bar{{\mu}}_{\cdot 1}=\bar{{\mu}}_{\cdot2}=\cdots=\bar{{\mu}}_{\cdot c}\mathrm{~vs~}H_{aC}{:}(\text{nоt }H_{0C})$，将因子 C 的均方除以 MSError(column).

至于第 \@ref(chap24) 章中讨论的裂区设计，其中有四个具有不同标准误的均值比较，而在条区设计中则存在六种具有不同标准误的均值比较。表 \@ref(tab:table25-4) 给出了比较、比较的估计（基于均衡数据集）和比较方差。使用第 \@ref(chap24) 章中描述的方法来获得单元格均值的每个对比的特定方差。为了演示这一点，两个均值差的方差，其中 A 和 C 的水平对于两个均值都不同，例如 $\mu_{12}-\mu_{34}$，需要首先获得 $\hat\mu_{ij}=\bar y_{ij\cdot}$ 的模型表示，该表示为 $\hat{{\mu}}_{ij}={\mu}_{ij}+\bar{{b}}_\cdot+\bar{{r}}_{i{\cdot}}+\bar{{c}}_{j{\cdot}}+\bar{{\varepsilon}}_{i{j}}$. 那么 $\hat\mu_{12}-\hat\mu_{34}$ 的模型为

$$\hat{\mu}_{12}-\hat{\mu}_{34}=\bar{y}_{12\cdot}-\bar{y}_{34\cdot}=\mu_{12}-\mu_{34}+\bar{r}_{1\cdot}-\bar{r}_{3\cdot}+\bar{c}_{2\cdot}-\bar{c}_{4\cdot}+\bar{\varepsilon}_{12\cdot}-\bar{\varepsilon}_{34\cdot}$$

<table>
<caption>(\#tab:table25-4)裂区设计的比较、估计和方差估计，其中 $i \ne m$ 和 $j \ne n$</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.4.png">
</td>
  </tr>
</tbody>
</table>

$\mu_{12}-\mu_{34}$ 的估计为 $\bar y_{12\cdot}-\bar y_{34\cdot}$，$\hat\mu_{12}-\hat\mu_{34}$ 的方差为

$$\begin{aligned}
\mathrm{Var}(\hat{\mu}_{12}-\hat{\mu}_{34})& =\mathrm{Var}(\bar{r}_1.-\bar{r}_3.+\bar{c}_2.-\bar{c}_4.+\bar{\varepsilon}_{12.}-\bar{\varepsilon}_{34.})  \\
&=\frac{2\sigma_\mathrm{row}^2}r+\frac{2\sigma_\mathrm{column}^2}r+\frac{2\sigma_\mathrm{cell}^2}r \\
&=\frac{2(\sigma_{\mathrm{cell}}^2+\sigma_{\mathrm{column}}^2+\sigma_{\mathrm{row}}^2)}{r}
\end{aligned}$$

$\hat\mu_{12}-\hat\mu_{34}$ 方差估计为

$$\widehat{\mathrm{Var}}(\hat\mu_{12}-\hat\mu_{34})
=\frac{2(\hat\sigma_{\mathrm{cell}}^2+\hat\sigma_{\mathrm{column}}^2+\hat\sigma_{\mathrm{row}}^2)}{r}$$

当 $\tilde{\sigma}_{\mathrm{row}}^2>0$ 以及 $\tilde{\sigma}_{\mathrm{column}}^2>0$ 时，该方差估计为

$$\widehat{\mathrm{Var}}(\hat{\mu}_{12}-\hat{\mu}_{34})=\frac2{acb}\left[aMSError(row)+cMSError(column)+(ac-a-c)MSError(cell)\right]$$

表 \@ref(tab:table25-4) 给出了六种类型的比较、它们的估计及其相应的方差。表 \@ref(tab:table25-5) 给出了比较的标准误估计，其中自由度对应于通过使用 Satterthwaite 近似计算的单个均方或均方组合的自由度，如下

$$\begin{aligned}
\hat{{\omega}}_{1}&=\frac{\left[MSError(column)+(c-1)MSError(cell)\right]^2}{\left[MSError(column)\right]^2}+\frac{\left[(c-1)MSError(cell)\right]^2}{(a-1)(c-1)(r-1)} \\
\hat{{\omega}} _2&=\frac{[MSError(row)+(a-1)MSError(cell)]^2}{\frac{[MSError(row)]^2}{(a-1)(r-1)}+\frac{[(a-1)MSError(cell)]^2}{(a-1)(c-1)(r-1)}} \\
\hat{{\omega}}_{3} &=\frac{[aMSError(row)+cMSError(column)+(ac-a-c)MSError(cell)]^2}{\frac{[aMSError(row)]^2}{(a-1)(r-1)}+\frac{[cMSError(column)]^2}{(c-1)(r-1)}+\frac{[(ac-a-c)MSError(cell)]^2}{(a-1)(c-1)(r-1)}} 
\end{aligned}$$

<table>
<caption>(\#tab:table25-5)裂区设计比较的标准误估计，附带相关自由度，其中 $i \ne m$ 和 $j \ne n$</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.5.png">
</td>
  </tr>
</tbody>
</table>

当有缺失的数据时，根据缺失数据的模式，会有更多不同类型的标准误。

## 示例：氮与灌溉 {#sec25-3}

<div class="figure" style="text-align: center">
<img src="figure/figure 25.4.png" alt="氮肥和灌溉条区的田间布局、处理分配和数据" width="466" />
<p class="caption">(\#fig:figure25-4)氮肥和灌溉条区的田间布局、处理分配和数据</p>
</div>

一项实验被开展以研究两种灌溉方法 (irrigation) 及三个不同氮肥 (nitrogen) 水平对小麦产量的影响和关系。实验采用四个区组或矩形布局，并按照随机化方案分配处理，其产量结果见图 \@ref(fig:figure25-4) 所示。这些数据的 SAS^®^-Mixed 代码以及 III 型平方和分析的方差分析表见表 \@ref(tab:table25-6)。比较氮肥水平、比较灌溉方法水平以及评估氮肥与灌溉交互作用的 $F$ 值分别为 60.13, 52.18 和 33.12. 这六种比较类型的标准误估计及其相应的自由度或近似自由度在表 \@ref(tab:table25-7) 中给出。最后三项比较涉及误差项的组合，因此自由度进行了近似计算。表 \@ref(tab:table25-8) 给出了用于提供方差分量的 REML 估计的 SAS-Mixed 代码，同时还提供了方差分量的估计。固定效应检验的结果显示在表 \@ref(tab:table25-9) 中。表 \@ref(tab:table25-9) 中的结果与表 \@ref(tab:table25-6) 中 III 型分析的结果一致。由于数据集是均衡的，并且方差分量的所有估计均大于零，这些检验结果是一致的。表 \@ref(tab:table25-7) 中用于比较氮肥水平、灌溉水平以及其他四种比较类型的 SAS-Mixed 代码及 estimate 语句在表 \@ref(tab:table25-10) 中给出。再次强调，这里有六种不同类型的比较，且每种都有各自的标准误估计。

<table>
<caption>(\#tab:table25-6)为氮肥和灌溉示例提供 III 型平方和分析的 SAS-Mixed 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.6.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table25-7)氮肥和灌溉数据的六个比较、标准误估计和自由度，其中 $i \ne m$ 和 $j \ne n$</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.7.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table25-8)为氮肥和灌溉示例提供 REML 分析的 SAS-Mixed 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.8.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table25-9)氮肥和灌溉示例 REML 分析的固定效应检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.9.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table25-10)在 SAS-Mixed REML 分析中使用 Estimate 语句来评估六种类型的比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.10.png">
</td>
  </tr>
</tbody>
</table>

**涉及多个尺寸实验单元设计的推断问题之一是，单元格均值的各种比较涉及不同的标准误估计。不可能通过提供单一的 LSD 值来总结这些类型的数据，就像完全随机和随机完全区组设计结构的情况一样**。

## 示例：含裂区的条区 1 {#sec25-4}

<div class="figure" style="text-align: center">
<img src="figure/figure 25.5.png" alt="田间区组之一的布局，氮肥和灌溉为条区，品种为裂区" width="351" />
<p class="caption">(\#fig:figure25-5)田间区组之一的布局，氮肥和灌溉为条区，品种为裂区</p>
</div>

第 \@ref(sec25-3) 节中的设计可以扩展到每个单元格内有三种小麦品种 (varieties) 的情况，如图 25.5 所示。图 \@ref(fig:figure25-5) 中没有显示随机化过程，但实际上，灌溉水平被随机分配到矩形的列，氮肥水平被随机分配到矩形的行，而品种则被随机分配到由每个矩形内部行与列相交形成的单元格内的三个区 (plots) 中。如果对品种水平进行平均，则氮肥与灌溉的分析结果与表 \@ref(tab:table25-3) 或表 \@ref(tab:table25-6) 中的结果相同。品种水平是每个单元格内的子区 (subplots)，品种的任何比较都是单元格内或子区间比较 (within-cell or between-subplot comparisons). 可用于描述该数据的模型是

$$y_{ijkl}=\mu_{ijk}+b_l+r_{il}+c_{jl}+d_{ijl}+\varepsilon_{ijkl},\quad i=1,2,j=1,2,3,k=1,2,3,l=1,2,3,4$$

其中 $\mu_{ijk}$ 表示第 i 种灌溉，第 j 种氮肥和第 k 种品种的均值，$b_l$ 表示第 l 个区组（矩形）的效应并假设 $b_l\sim i.i.d.N(0,\sigma_{\mathrm{block}}^2)$，$c_{jl}$ 表示第 l 个区组内第 i 行的效应并假设 $c_{jl}\sim i.i.d.N(0,\sigma_{\mathrm{column}}^2)$，$d_{ijl}$ 与第 l 个区组内第 ij 个单元格相关的随机效应 $d_{ijl}\sim i.i.d.N(0,\sigma_{\mathrm{cell}}^2)$，以及 $\varepsilon_{ijkl}$ 与分配给第 l 个区组内第 ij 个单元格的第 k 种品种相关的子区随机效应并假设 ${\varepsilon}_{ijkl}\thicksim i.i.d.N(0,{\sigma}_{\mathrm{subplot}}^2)$. 表 \@ref(tab:table25-11) 给出了该模型的方差分析，表 \@ref(tab:table25-12) 给出了可用于拟合该模型的 SAS-Mixed 代码。行误差项是 block × nit，列误差项是 block × irr，单元格误差是 block × nit × irr. 残差是子区误差项，是通过 variety 与 block 交互作用来计算的，其中合并了所有 nit 和 irr 的水平 (the variety by block interaction pooled across the levels of nit and irr).

<table>
<caption>(\#tab:table25-11)含裂区的条区设计结构方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.11.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table25-12)为第 \@ref(sec25-4) 节中的条-裂区 (strip-split-plot) 数据集拟合模型的 SAS-Mixed 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.12.png">
</td>
  </tr>
</tbody>
</table>

## 示例：含裂区的条区 2 {#sec25-5}

<div class="figure" style="text-align: center">
<img src="figure/figure 25.6.png" alt="田间区组之一的布局，氮肥和灌溉为条区，品种为氮肥水平的裂区，品种在灌溉水平下构成了条区" width="332" />
<p class="caption">(\#fig:figure25-6)田间区组之一的布局，氮肥和灌溉为条区，品种为氮肥水平的裂区，品种在灌溉水平下构成了条区</p>
</div>

[^splitplotstripplot]: 原文：The levels of variety are subplots for the levels of nitrogen, but the levels of variety and the levels of nitrogen are strip-plots with the levels of irrigation.

另一种（不同于第 \@ref(sec25-4) 节所述的）将品种水平纳入研究的方法在图 \@ref(fig:figure25-6) 中展示。在每个氮肥水平下，品种水平跨越了灌溉的两个水平。品种水平是氮肥水平的子区，但品种水平和氮肥水平在灌溉水平下构成了条区[^splitplotstripplot]。通过在每个单元格内合并品种水平，得到的灌溉与氮肥的分析结果与第 \@ref(sec25-3) 节中的相同。设计中氮肥与品种部分可通过在每个氮肥品种组合内合并灌溉水平来获得。由此产生的设计是一个裂区设计，其中氮肥水平构成整区因子，品种水平构成子区因子。裂区分析的方差表在表 \@ref(tab:table25-13) 中展示。整区误差与行误差相同，计算为区组与氮肥交互作用。子区或子行 (subrow) 或 1/3 行误差计算为 block × variety 交互作用，其中合并了所有灌溉水平。接下来考虑仅来自氮肥水平 1 的数据。在这种情况下，得到的设计是一个涉及灌溉和品种水平的条区设计。对应的条区方差分析表在表 \@ref(tab:table25-14) 中列出。我们感兴趣的误差项是品种与灌溉交点处的子单元格 (subcell) 误差，这被计算为 block × irrigation × variety 交互作用，具有 6 个自由度。这些误差在氮肥水平之间进行合并，产生了 18 个自由度。nitrogen × irrigation × variety 交互作用也是一个子单元格比较。因此，block × nitrogen × irrigation × variety 交互作用提供了额外的 12 个子单元格误差自由度，使子单元格误差项共有 30 个自由度。表 \@ref(tab:table25-15) 中的方差分析表展示了变异源及其相应的自由度。

<table>
<caption>(\#tab:table25-13)第 \@ref(sec25-5) 节中示例分析的裂区部分</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.13.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table25-14)第 \@ref(sec25-5) 节中氮肥第一水平的灌溉和品种数据</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.14.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table25-15)第 \@ref(sec25-5) 节条-裂区 (strip-split-plot) 设计的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.15.png">
</td>
  </tr>
</tbody>
</table>

## 示例：含裂区的条区 3 {#sec25-6}

[^varietiesformastripplot]: 原文：The varieties form a strip plot with the combination of levels of nitrogen and herbicide.

[^seventytwodf]: 原文：The subplot error is obtained by computing the block × seeding rate interaction for each combination of herbicide, nitrogen, and variety and then pooling these interactions across the levels of herbicide, nitrogen, and variety.

<div class="figure" style="text-align: center">
<img src="figure/figure 25.7.png" alt="田间区组之一的图形：行分配了除草剂和氮肥的不同水平，列分配了品种的不同水平，单元格分配了播撒率的不同水平" width="467" />
<p class="caption">(\#fig:figure25-7)田间区组之一的图形：行分配了除草剂和氮肥的不同水平，列分配了品种的不同水平，单元格分配了播撒率的不同水平</p>
</div>

图 \@ref(fig:figure25-7) 显示，三个品种被分配到矩形的列，除草剂 (herbicide) 和氮的组合被分配到矩形的行，播撒率 (seeding rates) 的水平被分配到单元格内的子区。假设设计由四个这样的区组组成。品种形成氮和除草剂水平组合的条区[^varietiesformastripplot]，播撒率是单元格内的子区。行的双向处理结构是该设计的独特之处。该结构指出，**任何类型的处理结构都可以与设计结构的任何特征相关联**。方差分析表列于表 \@ref(tab:table25-16) 中，其中行误差是通过合并 block × herbicide, block c nitrogen 和 block × herbicide × nitrogen 交互作用来计算的，提供九个自由度。列误差项由 block × variety 交互作用计算。单元格误差项是通过合并 block × variety × herbicide, block × variety × nitrogen 和 block × variety × herbicide × nitrogen 交互作用来计算的，这些为单元格误差项提供了 18 个自由度。子区误差是通过计算每个 herbicide, nitrogen 和 variaty 组合的 block × seeding 交互作用得到的，然后跨 herbicide, nitrogen 和 variety 的水平合并这些交互作用[^seventytwodf]，为残差提供 72 个自由度。

<table>
<caption>(\#tab:table25-16)第 \@ref(sec25-5) 节条-条区设计的方差分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.16.png">
</td>
  </tr>
</tbody>
</table>

## 示例：含裂区的条区 4 {#sec25-7}

<div class="figure" style="text-align: center">
<img src="figure/figure 25.8.png" alt="田间区组之一的布局，将品种和播种率作为整区，将条区中的除草剂和氮肥水平作为子区" width="420" />
<p class="caption">(\#fig:figure25-8)田间区组之一的布局，将品种和播种率作为整区，将条区中的除草剂和氮肥水平作为子区</p>
</div>

这个实验涉及研究两种品种、两种播撒率、两种除草剂和两种氮肥对玉米产量的影响。图 \@ref(fig:figure25-8) 显示了实验单元的四个区组中的一个是如何处理的。品种和播撒率的水平被分配给每个区组内的四个大方块 (squares). 在每个方块内，除草剂和氮肥的水平形成一个条区。条区部分是大方块的子区。如果忽略氮肥和除草剂的水平，品种和播撒率的设计是一个随机完全区组设计结构中的双向处理结构。误差项通过处理结构与区组交互作用计算，其中合并了 block × variety, block × seeding rate 和 block × variety × seeding rate 交互作用得出，共提供 9 个自由度。可以在 SAS-Mixed 中使用 variety × seeding × block 获得大区组 (big block) 误差项，该程序将自动进行合并。表 \@ref(tab:table25-17) 给出了大区组部分的方差分析表。对于品种 1 和播撒率 1，除草剂和氮肥的水平提供了一个条区设计，其方差分析总结在表 \@ref(tab:table25-18) 中。误差项是通过合并品种和播撒率的四种组合得出的，为 Error(row/block)（区组内行的误差）, Error(column/block)（区组内列的误差）和 Error(cell/block)（区组内单元格的误差）提供了 12 个自由度。SAS-Mixed 代码和完整的方差分析表显示在表 \@ref(tab:table25-19) 中。误差项的计算方法如下：Error(big squares) = Variety × Seeding × Block, Error(row/ block) = Herbicide × Block(Variety Seeding), Error(column/block) = Nitrogen × Block(Variety Seeding), and Error(cell/block) = Nitrogen × Herbicide × Block(Variety Seeding).

<table>
<caption>(\#tab:table25-17)第 \@ref(sec25-7) 节中示例的大区组方差分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.17.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table25-18)第 \@ref(sec25-7) 节中示例，来自品种 1 和播撒率 1 的除草剂和氮肥数据的条区方差分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.18.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table25-19)第 \@ref(sec25-7) 节示例中裂区中的条区的 SAS Mixed 代码和方差分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.19.png">
</td>
  </tr>
</tbody>
</table>

## 条-条区的设计与分析：基于 JMP7 {#sec25-8}

<div class="figure" style="text-align: center">
<img src="figure/figure 25.9.png" alt="应用温度、清洁和冲洗水平的过程" width="470" />
<p class="caption">(\#fig:figure25-9)应用温度、清洁和冲洗水平的过程</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 25.10.png" alt="三维矩形展示：清洁水平的分配给行，冲洗水平分配给列，温度水平分配给层。数字对应于图 \@ref(fig:figure25-9) 中的硅片编号" width="454" />
<p class="caption">(\#fig:figure25-10)三维矩形展示：清洁水平的分配给行，冲洗水平分配给列，温度水平分配给层。数字对应于图 \@ref(fig:figure25-9) 中的硅片编号</p>
</div>

当实验单元被组织成具有行、列和层 (tiers) 的三维矩形时，就会出现条形-条形-图设计结构。半导体行业的某个实验被用来演示条-条区 (split-split plot) 设计结构。该实验包括评估三个因素，每个因素在各自独立的工艺步骤中出现；也就是说，该实验涉及三个步骤。第一步是在硅片 (silicon wafer)（直径 20 厘米，厚度 0.2 厘米）表面添加一层氧化物。通过将硅片放入设定特定温度的炉子中，对硅片表面施加氧化层。研究了两种温度水平。第二步是对硅片表面进行抛光 (polishing) 或清洁 (cleaning)，以消除可能产生的任何凸起，评估了两种清洗方法。最后，第三步涉及冲洗 (washing off) 表面残留颗粒，这些颗粒可能会对正在构建的电路造成干扰。实验中包含了两种冲洗方法。在每个步骤中，均有四片硅片一同接受相应因素的处理。图 \@ref(fig:figure25-9) 的示意图展示了八片硅片如何经过这三个步骤以及在每个步骤中硅片是如何分组的。图 \@ref(fig:figure25-10) 是一个三维展示图，其中清洁的水平被分配到行，清洁后冲洗的水平被分配到列，温度的水平被分配到层。该研究重复了四次，产生了四个区组。清洁水平的实验单元是区组内的行，并且行误差项计算为 clean × block 交互作用。冲洗水平的实验单元是区组内的列，并且列误差项计算为 wash × block 交互作用。如果对温度水平进行求和，所得结构就是以清洁和冲洗作为两个因素的条区设计。清洁与冲洗交互作用的实验单元是来自两种温度下的两片硅片，而这两片硅片的误差计算为 clean × wash × block 的交互作用。温度水平的实验单元是三维矩形的层。层误差项计算为 temperature × block 交互作用。如果对清洁水平进行求和，所得结构就是以温度和冲洗作为两个因素的条区设计。温度与冲洗交互作用的实验单元是来自两个清洁水平的硅片，其误差项计算为 temperature × wash × block 交互作用。如果对冲洗水平进行求和，所得结构就是以温度和清洁作为两因素的条区。温度与清洁交互作用的实验单元是来自两个冲洗水平的两个硅片，并且两个冲洗硅片的误差计算为 temperature × clean × block 交互作用。测量温度、清洁与冲洗交互作用的实验单元是单个硅片，硅片误差计算为 temperature × clean × wash × block     交互作用。测量的响应是经三个步骤后在每个硅片上的九个位置（以埃为单位）测量的所得氧化物层的平均厚度。可用于描述厚度数据的模型是

$$y_{ijkl}=\mu_{ijk}+b_l+r_{il}+c_{jl}+t_{kl}+w_{ijl}^{2W}+w_{ikl}^{2C}+w_{jkl}^{2T}+\varepsilon_{ijkl},\quad i=1,2,j=1,2,k=1,2,l=1,2,3,4$$

$\mu_{ijk}$ 表示经水平 i 温度、水平 j 清洁和水平 k 冲洗的平均厚度，以及

$$\begin{aligned}&b_{l}\sim i.i.d.N(0,\sigma_{\mathrm{block}}^{2}),\quad r_{il}\sim i.i.d.N(0,\sigma_{\mathrm{row}}^{2}),\quad c_{jl}\sim i.i.d.N(0,\sigma_{\mathrm{column}}^{2})\\&t_{kl}\sim i.i.d.N(0,\sigma_{\mathrm{tier}}^{2}),\quad w_{ijl}^{2W}\sim i.i.d.N(0,\sigma_{2\text{washwafer}}^{2})\\&w_{ikl}^{2\text{С}} \sim i . i . d . N ( 0 , \sigma _ { 2\text{cleanwafer}}^{2}),\quad w_{jkl}^{2T}\sim i.i.d.N(0,\sigma_{2\text{tempwafer}r}^{2})\end{aligned}$$

以及

$$\varepsilon_{ijkl}\sim i.i.d.~N(0,\sigma_{\mathrm{wafer}}^2)$$

方差分析在表 \@ref(tab:table25-20) 中给出，其中 w, 2Tw, 2Cw 和 2Ww 分别表示硅片、两个温度硅片、两个清洁硅片和两个冲洗硅片。每个因子效应都有自己的一组实验单元和相应的误差项。

<table>
<caption>(\#tab:table25-20)第 \@ref(sec25-8) 节条-条区设计的方差分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2025.20.png">
</td>
  </tr>
</tbody>
</table>

图 \@ref(fig:figure25-11) 包含了数据集，而图 \@ref(fig:figure25-12) 是 JMP model specification 界面。通过使用 Attributes 按钮，将所有涉及 block 的交互作用定义为随机效应。图 \@ref(fig:figure25-13) 展示了固定效应的 REML 方差成分估计以及 III 型检验结果。存在显著的 temperature × clean × wash 交互作用，因此需要评估这些单元格均值。图 \@ref(fig:figure25-14) 给出了三向交互作用的最小二乘均值以及表示显著差异的字母标记。查看三向交互作用均值两两比较的另一种方式是绘制所有两两差值，从最大到最小排序，并附带差异的置信区间。图 \@ref(fig:figure25-15) 是两两差异的展示。由于差异的方差取决于比较类型，所以置信区间的宽度不全相等。相比于表 \@ref(tab:table25-4) 中展示的条区设计的不同方差数量，条-条区设计的不同方差数量要多得多。

<div class="figure" style="text-align: center">
<img src="figure/figure 25.11.png" alt="第 \@ref(sec25-8) 节中示例的三步过程的数据集" width="539" />
<p class="caption">(\#fig:figure25-11)第 \@ref(sec25-8) 节中示例的三步过程的数据集</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 25.12.png" alt="第 \@ref(sec25-8) 节中示例的三步过程的 fit model 界面" width="618" />
<p class="caption">(\#fig:figure25-12)第 \@ref(sec25-8) 节中示例的三步过程的 fit model 界面</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 25.13.png" alt="三步过程的方差分量估计和固定效应检验" width="652" />
<p class="caption">(\#fig:figure25-13)三步过程的方差分量估计和固定效应检验</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 25.14.png" alt="带字母标记的三向最小二乘均值" width="364" />
<p class="caption">(\#fig:figure25-14)带字母标记的三向最小二乘均值</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 25.15.png" alt="三步过程的三向最小二乘均值差异的图形表示" width="549" />
<p class="caption">(\#fig:figure25-15)三步过程的三向最小二乘均值差异的图形表示</p>
</div>

## 结束语 {#sec25-9}

条区设计是一种基本的设计结构，本章讨论的例子说明了可能出现的一些复杂性。包含多个例子是为了让读者能够掌握一种针对复杂设计进行恰当分析的方法。**条区设计结构被用来展示仅仅通过向实验单元分配因子水平的方式，就能产生多种尺寸的实验单元**。对基本条区进行了分析，以演示在进行各种均值比较时计算标准误估计的复杂性。显然，随着设计的复杂化，有许多不同的标准误需要评估，条区和裂区的组合设计以及条条区设计都证明了这一点。样本量和功效计算可以使用第 \@ref(chap25) 章中描述的方法，并结合所关注比较的标准误估计来进行。

## 练习 {#sec25-10}
