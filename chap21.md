
# 案例研究：随机效应模型分析 {#chap21}

前三章描述了分析随机效应模型的方法，并提供了一些示例来演示各种分析技术。本章介绍了更复杂实验情况的分析，包括估计、模型构建、假设检验和置信区间估计。

## 数据集 {#sec21-1}

在这个实验中，我们研究了多家工厂 (plants) 生产线中工人的工作效率。在公司的工厂中随机挑选了三家。从每一家工厂中，我们又随机挑选了四个装配地点 (sites) 和三个工人 (workers). 我们期望每个工人在他们工厂的每个装配地点工作五次，但由于调度问题和其他优先级，每个工人实际工作的次数各不相同。工人们到达每个地点的顺序是随机的，并且他们尽可能地遵循这个安排。响应变量是组装零件的效率，这是组装单元数量和错误次数的函数。效率得分以及工厂编号、地点编号和工人编号都列在表 \@ref(tab:table21-1) 中，其中 EFF_1, EFF_2,..., EFF_5 表示工人可能工作的五天的得分。

<table>
<caption>(\#tab:table21-1)案例研究中的数据，其中 EFF_i 表示工人在某一地点工作的第 i 次时间</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2021.1.png">
</td>
  </tr>
</tbody>
</table>

研究中的所有三个因素都是随机的，地点和工人都位于工厂内。因此，地点和工人效应以及它们的交互作用效应都嵌套在工厂内。用于描述数据的模型，对于 $i=1,2,3,j=1,2,3,4,k=1,2,3$ 和 $l=1,\ldots,n_{ijk}$

\begin{align}
y_{ijkl}=\mu+p_i+s_{j(i)}+w_{k(i)}+(sw)_{jk(i)}+\varepsilon_{ijkl}
(#eq:21-1)
\end{align}

其中 $p_i$ 表示第 i 个工厂的效应，$s_{j(i)}$ 表示工厂 i 内第 j 个地点的效应，$w_{k(i)}$ 表示表示工厂 i 内第 k 个工人的效应，$(sw)_{jk(i)}$ 表示工厂 i 内地点和工人的交互效应，$\varepsilon_{ijkl}$ 为残差项。

假设

$$\begin{aligned}&p_i\thicksim i.i.d.N(0,\sigma_p^2),s_{j(i)}\thicksim i.i.d.N(0,\sigma_\mathrm{s}^2),w_{k(i)}\thicksim i.i.d.N(0,\sigma_w^2)\\
&sw_{jk(i)}\thicksim i.i.d.N(0,\sigma_{sw}^2),\varepsilon_{ijkl}\thicksim i.i.d.N(0,\sigma_\varepsilon^2),\end{aligned}$$

且所有随机变量 $p_i,s_{j(i)},w_{k(i)},(sw)_{jk(i)},\varepsilon_{ijkl}$ 独立。

## 估计 {#sec21-2}

采用 I-III 型平方和、REML、最大似然、MIVQUE0（含和不含 NOBOUND 方法）的矩法估计程序求解方差分量，结果列于表 \@ref(tab:table21-2). 同时附上用于获得 REML 估计的 SAS^®^-Mixed代码，可通过将特定名称与 `Method = ` 选项结合来选择其他估计技术。表 \@ref(tab:table21-3) 和 \@ref(tab:table21-4) 给出了 I 型和 III 型平方和的期望均方，可用于构建求取每个方差分量的矩法估计所需的方程。$\sigma_s^2$ 的 REML, ML 和 MIVQUE0 估计均为零，而包含 `nobound` 选项的 MIVQUE0 以及 I-III 型 $\sigma_s^2$ 的解均为负值。通过从不含 `nobound` 选项的 MIVQUE0 解中将 $\sigma_s^2$ 的解设为零，得到 MIVQUE0 解。通过将负解转换为 0 来获得方差分量的估计；也就是说，$\sigma_s^2$ 的估计为 $\hat\sigma_s^2=0$. 由于 $\sigma_s^2$ 的估计为零，下一步将是对模型中的方差分量进行假设检验。可以通过逐步删除过程从模型中删除随机分量，以尝试获得一个更简单的模型来描述过程中的变异。**从模型中删除一个项本质上就是将移除的方差分量设置为零。这种模型构建可以针对与处理结构中的因素相对应的任何方差分量进行，但是不会对与设计结构中的因素相对应的方差分量进行模型构建**。

<table>
<caption>(\#tab:table21-2)使用 Proc Mixed 中可用的每种方法，为模型 \@ref(eq:21-1) 方差分量生成解，其中 MVQ 表示 MIVQUE0，MVQNB 表示 `nobound` 选项的 MIVQUE。Proc Mixed 代码适用于 Method = REML 方法</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2021.2.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table21-3)带有期望均方的 I 型方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2021.3.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table21-4)带有期望均方的 III 型方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2021.4.png">
</td>
  </tr>
</tbody>
</table>

## 模型构建 {#sec21-3}

基于 I 型平方和的方差分析表用于通过检验有关方差分量的假设来构建模型。该过程从研究与最后一个平方和相关的方差分量开始，然后逐行向上到与第一个平方和相关的方差分量。零假设是方差分量等于零。I 型平方和方法观测到的以及期望均方列于表 \@ref(tab:table21-3) 中。

第一步是检验 $H_0\colon{\sigma}_{sw}^2=0\text{ vs }H_{a}\colon{\sigma}_{sw}^2>0$. 检查表 \@ref(tab:table21-3) 中的期望均方后，残差被确定为适当的除数，因此检验统计量为

$$F_{C_{\mathrm{s}w}}=\frac{MS[worker\times site(plant)]}{MS(Residual)}=21.42\left(\text{which is in Table 2}1.5\right)$$

将 $F_{C_{\mathrm{s}w}}$ 的抽样分布与分子自由度为 18、分母自由度为 82 的 $F$ 分布的分位数进行比较。观察到的显著性水平小于 0.0001，这表明 ${\sigma}_{sw}^2$ 是过程中变异性的重要来源。

下一步是检验 $H_0\colon{\sigma}_{s}^2=0\text{ vs }H_{a}\colon{\sigma}_{s}^2>0$. $E\{MS[site(plant)]\}=\sigma_\varepsilon^2+3.4789\sigma_{sw}^2+9.1928\sigma_s^2$ 并且不存在其他期望为 $\sigma_\varepsilon^2+3.4789\sigma_{sw}^2$ 的均方可用作检验统计量的除数。因此，需要根据 $MS(Residual)$ 和 $MS[worker\times site(рlant)]$ 构建均方 $Q_s^*$，使得 $E(Q_s^*)=\sigma_\varepsilon^2+3.4789\sigma_{sw}^2$. 这样的 $Q_s^*$ 为

$$\begin{aligned}
Q_s^*& =3.4789{\left[\frac{MS[worker\times site(plant)]}{2.8569}\right]}+{\left[1-\frac{3.4989}{2.8569}\right]}MS(Residual)  \\
&=1.2177MS[worker\times site(plant)]-0.2177MS(Residual) \\
&=128.921
\end{aligned}$$

卡方分布可用于近似 $r_sQ_s^*/(\sigma_\varepsilon^2+3.4789\sigma_{sw}^2)$ 的抽样分布，使用 Satterthwaite 近似将确定相关的自由度确定为

$$\begin{aligned}r_s&=\frac{(Q_s^*)^2}{\frac{\{1.2177MS[worker\times site(plant)]\}^2}{18}+\frac{[0.2177MS(Residual)]}{82}}\\&=\frac{\left(128.921\right)^2}{\frac{[1.2177\times106.738]^2}{18}+\frac{[0.2177\times4.983]^2}{82}}=17.7\end{aligned}$$

检验统计量为

$$F_{C_s}=\frac{MS[site(plant)]}{Q_s^*}=0.65\text{ (see Table 21.5)}$$

其近似分布为自由度为 9 和 17.7 的 $F$ 分布。观察到的检验显著性水平为 0.7396，这表明 ${\sigma}_{s}^2$ 是过程中可忽略不计的变异源；也就是说，无法拒绝 $H_0\colon{\sigma}_{s}^2=0\text{ vs }H_{a}\colon{\sigma}_{s}^2>0$. 表 \@ref(tab:table21-5) 和 \@ref(tab:table21-6) 包含均方、期望均方、合适的误差项、近似分母自由度以及基于 III 型平方和的检验统计量。用于检验 $H_0\colon{\sigma}_{s}^2=0\text{ vs }H_{a}\colon{\sigma}_{s}^2>0$ 的近似 $F$ 统计量为 0.67，分母自由度估计为 18.1，以及观察到的显著性水平为 0.7217. III 型分析中关于 ${\sigma}_{s}^2$ 的结论与 I 型分析中得出的结论相同。由于 *site(plant)* 是处理结构的一部分，并且 ${\sigma}_{s}^2$ 可以忽略不计，因此一种策略是将 ${\sigma}_{s}^2$ 设置为零，从模型中消除 $s_{j(t)}$ 并根据数据拟合缩减模型。

<table>
<caption>(\#tab:table21-5)方差分量假设的 I 型检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2021.5.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table21-6)方差分量假设的 III 型检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2021.6.png">
</td>
  </tr>
</tbody>
</table>

## 缩减模型 {#sec21-4}

缩减模型为

\begin{align}
y_{ijkl}&=\mu+p_i+w_{k(i)}+(sw)_{jk(i)}+\varepsilon_{ijkl}\\i&=1,2,3,~j=1,2,3,4,~k=1,2,3,~l=1,\ldots,n_{ijk}
(#eq:21-2)
\end{align}

当从模型中移除 $s_{j(i)}$ 时，SAS-Mixed 将 site(plant) 的平方和与 worker × site(plant) 的平方和合并，如表 \@ref(tab:table21-7) 中 worker×site(plant) 的 27 个自由度所示。当 ${\sigma}_{s}^2=0$ 时，由于来自表 \@ref(tab:table21-3) 和 \@ref(tab:table21-4) 的 $MS[site(plant)]$ 和 $MS[worker \times site(plant)]$ 的期望均方相似，因此这是一个合理的过程；也就是说，在表 \@ref(tab:table21-3) 中，$E\{MS[site(plant)]\}=\sigma_\varepsilon^2+3.4789\sigma_{sw}^2$ 和 $E\{MS[worker\times site(plant)]\}=\sigma_\varepsilon^2+2.8569\sigma_{sw}^2$，在 III 型分析中，系数甚至更接近。

<table>
<caption>(\#tab:table21-7)带期望均方的缩减模型的 I 型方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2021.7.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table21-7) 列出了缩减模型的自由度、I 型均方和期望均方。worker × site(plant) 的期望均方为 $\sigma_\varepsilon^2+3.0643\sigma_{sw}^2$. 此系数通过 $3.0643 = [(9 \times 3.4789) + (18\times2.8561)/27]$ 计算，这是 I 型分析中 site(plant) 和 worker × site(plant) 的预期均方的 ${\sigma}_{s}^2=0$ 的合并系数。这种等价性之所以发生，是因为 I 型平方和是序贯性的，因此合并过程是加性的。通过比较表 \@ref(tab:table21-4), \@ref(tab:table21-6) 和  \@ref(tab:table21-8)  中的 III 型分析，可以发现这种现象不会发生在其他类型的平方和上。使用 REML、ML MIVQUE0、I-III 型方法的缩减模型的方差分量估计列于表  \@ref(tab:table21-9)  中。除了 MIVIQUE0（其值最大，为 6.61）之外，所有方法的残差方差估计大约都等于 4.98。ML 估计的 ${\sigma}_{p}^2=0$ 值最小，为 29.6，而 III 型估计最大，为 58.6. 其他方法的范围在 47.6 到 53.1 之间。由于 plant 只有三个水平，${\sigma}_{p}^2$ 是最难估计的方差分量；也就是说，与模型中的其他方差分量相比，对 ${\sigma}_{p}^2$ 的了解较少。 ${\sigma}_{w}^2$ 和 ${\sigma}_{sw}^2$ 的估计在不同的估计方法中非常一致，${\hat\sigma}_{w}^2$ 的范围为 25.5 至 28.9，${\hat\sigma}_{sw}^2$ 的范围为 25.4 至 30.7。

<table>
<caption>(\#tab:table21-8)带期望均方的缩减模型的 III 型方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2021.8.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table21-9)缩减模型的方差分量估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2021.9.png">
</td>
  </tr>
</tbody>
</table>

可以通过对残差方差分量进行假设检验来继续分析，以确定模型是否可以进一步缩减。用于检验 $H_0\colon{\sigma}_{sw}^2=0\mathrm{~vs~}H_a\colon{\sigma}_{sw}^2>0$ 的统计量为

$$\begin{aligned}F_{C_{sw}}=\frac{MS[worker\times site(plant)]}{MS(Residual)}=19.90\end{aligned}$$

其分布为自由度为 27 和 82 的 $F$ 分布。观察到的显著性水平小于 0.0001. 因此，${\sigma}_{sw}^2$ 是生成数据过程中变异的重要来源并且 $(sw)_{jk(i)}$ 应保留在模型中。同样，${\sigma}_{sw}^2$ 是适应方差分量 (adaption variance component)，其中一些工人在某些地点工作更有效，而其他工人在不同地点工作更有效。根据工人的需求定制每个地点有助于减少系统的变异性。

$Q^*_w$ 需计算用于检验 $H_0\colon{\sigma}_{w}^2=0\mathrm{~vs~}H_a\colon{\sigma}_{w}^2>0$. 均方 $Q^*_w$ 根据 $MS[worker\times site(plant)$ 和 $MS(Residual)$ 构建，使得 $E(Q^*_w)=\sigma_\varepsilon^2+3.9015\sigma_{sw}^2$. 那么

$$\begin{aligned}
Q _{w}^{*}& =3.9015{\left[\frac{MS[site\times worker(plant)]}{3.0643}\right]}+{\left[1-\frac{3.9015}{3.0643}\right]}MS(Residual)  \\
&=1.2732MS[site\times worker(plant)]-0.2732MS(Residual) \\
&=124.9085
\end{aligned}$$

检验统计量为 $F_{C_w}=456.9419/124.9085=3.66$，其分布为具有 6 和 $r_w$ 自由度的 $F$ 分布，其中

$$\begin{aligned}
r_{w}& =\frac{(Q_w^*)}{\{1.2732MS[site\times worker(plant)]\}^2/27+[0.2732\times MS(Residual)]^2/82}  \\
&=\frac{(124.9085)^2}{[1.2732\times99.1752]^2/27+[0.2732\times4.9831]^2/82} \\
&=26.42
\end{aligned}$$

该检验统计量的观测显著性水平为 0.0089，表明 ${\sigma}_{w}^2$ 是系统变异的重要组成部分。

最后，需要从 $MS[worker(plant)],MS[worker \times site(plant)]$ 和 $MS(Residual)$ 构造另一个 $Q_p^*$，使得

$$E(Q_p^*)=\sigma_\varepsilon^2+3.9277\sigma_{sw}^2+13.136\sigma_w^2$$

以便可以构建统计量来检验 $H_0\colon{\sigma}_{p}^2=0\mathrm{~vs~}H_a\colon{\sigma}_{p}^2>0$. 所需的 $Q_p^*$ 为

$$\begin{aligned}
Q_p^*=& =\frac{13.136}{13.037}\{MS[worker(plant)]\}  \\
&+\left[\frac{3.9277-\frac{13.136}{13.037}(3.0915)}{3.0643}\right]\{MS[site\times worker(plant)]\} \\
&+\left[1-\frac{13.136}{13.037}-\frac{3.9277-\frac{13.136}{13.037}(3.0915)}{3.0643}\right]MS(Residual) \\
&=1.0076MS[worker(plant)]-0.0012MS[site\times worker(plant)]-0.0065MS(Residuаl) \\
&=460.2638
\end{aligned}$$

检验统计量为

$$F_{C_p}=\frac{MS(plant)}{Q_p^*}=5.027$$

其分布为具有 2 和 $r_p$ 自由度的 $F$ 分布，其中 $r_p=(Q_p^*)^2/D_p$ 以及

$$\begin{aligned}
D_{p} =&\frac{\{1.0076MS[worker(plant)]\}^2}6+\frac{\{0.0012MS[site\times worker(plant)]\}^2}{27}  \\
&+\frac{\left[0.0065MS(Residual)\right]^2}{82}=35,330.281
\end{aligned}$$

那么 $r_p=5.9961$. 与此检验相关的观察到的显著性水平为 0.0522，这表明 ${\sigma}_{p}^2$ 是系统变异的重要贡献者，但它不如 ${\sigma}_{sw}^2$ 或 ${\sigma}_{w}^2$ 那么重要。模型 \@ref(eq:21-2) 中的所有方差分量在 $\alpha < 0.10$ 时与零显著不同，并且模型无法进一步缩减。

表 \@ref(tab:table21-8) 和 \@ref(tab:table21-10) 中的结果来自 III 型平方和。 SAS-Mixed 使用期望均方来计算适当的除数，并使用 Satterthwaite 近似来计算分母自由度。用于检验 $H_0\colon{\sigma}_{p}^2=0\mathrm{~vs~}H_a\colon{\sigma}_{p}^2>0,H_0\colon{\sigma}_{w}^2=0\mathrm{~vs~}H_a\colon{\sigma}_{w}^2>0$ 和 $H_0\colon{\sigma}_{sw}^2=0\mathrm{~vs~}H_a\colon{\sigma}_{sw}^2>0$ 的 III 型平方和分析的显著性水平分别为 0.0357, 0.0040 和小于 0.0001. III 型分析的显著性水平略低于 I 型分析（对于前两个检验），但在其他问题中，III 型分析的显著性水平可能会略高于 I 型分析。表 \@ref(tab:table21-11) 包含缩减模型的 REML 分析结果。与plant, worker(plant) 和 worker × site(plant) 方差分量 $Z$ 检验统计量相关的显著性水平分别为 的 0.2118, 0.0869 和 0.0003. 这些显著性水平远大于 I 型和 III 型平方和方法的结果。出现这种情况是因为 $Z$ 检验渐近服从正态随机变量，而在这些情况下，与每个方差分量估计相关联的自由度数量都很小。除非有很多水平与方差分量的估计相关，否则与 $Z$ 值相关的信息对于检验假设没有用处。

<table>
<caption>(\#tab:table21-10)缩减模型方差分量假设的 III 型检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2021.10.png">
</td>
  </tr>
</tbody>
</table>

plant 方差分量测量工厂总体内工厂之间的差异。由于存在很大的变异性，其中一家工厂很可能制定了适当的程序，使工人能够更有效地工作。worker(plant) 方差分量的重要性表明某些工人比其他工人更有效率，因此培训计划可以帮助提高表现不如其他工人的效率。worker × site(plant) 方差分量是一个适应方差分量，这意味着一些工人更适应在某些地点工作，而在其他地点则表现不佳，另一方面，一些其他工人在之前工人表现不佳的地点上表现出色。

## 置信区间 {#sec21-5}

[^standarddeviationunits]: 原文：Table 21.13 contains the computed degrees of freedom for the REML confi dence intervals as well as intervals that are expressed in standard deviation units.

分析的下一步是在方程 \@ref(eq:21-2) 中构造关于缩减模型方差分量的置信区间。表 \@ref(tab:table21-11) 和 \@ref(tab:table21-12) 中显示的结果分别是使用 REML 和 I 型平方和对缩减模型方差分量的估计。使用第 \@ref(sec20-2-2) 节中的方法计算 REML 解的置信区间，其中 $df=2(Z\text{-}value)^2$. $\hat{{\sigma}_p^2},\hat{{\sigma}_w^2},\hat{{\sigma}_{sw}^2}$ 的自由度分别为 1.28, 3.70 和 82.10. 表 \@ref(tab:table21-13) 包含 REML 置信区间以及以标准差单位表示的区间的计算自由度[^standarddeviationunits]。I 型分析计算的置信区间是 \@ref(sec20-2-4) 节的 Wald 区间（由于区间关于估计对称，且一些下限为负值，因此可以注意到这一点），但对于 $\sigma^2_\varepsilon$，区间是使用 \@ref(sec20-2-2) 节中的结果，基于 81.93 个自由度计算得出的。基于 I 型平方和的置信区间可以使用第 \@ref(sec20-2-2) 节的方法重新计算，其中 $df=2(Z\text{-}value)^2$. 对于 ${\hat{\sigma}_p^2}$，由此得到的自由度为 1.4677，提供的 95% 置信区间为 $[11.2668\le {\hat{\sigma}_p^2} \le 5984.61]$. 对于 ${\hat{\sigma}_w^2}$，由此得到的自由度为 3.9682，提供的 95% 置信区间为 $[9.1149\le {\hat{\sigma}_p^2} \le 212.97]$. 对于 ${\hat{\sigma}_{sw}^2}$，由此得到的自由度为 20.8392，提供的 95% 置信区间为 $[18.1625\le {\hat{\sigma}_p^2} \le 62.98]$. 表 \@ref(tab:table21-14) 总结了这些重新计算的置信区间结果。

<table>
<caption>(\#tab:table21-11)缩减模型方差分量的 REML 估计。方差分量的 Satterthwaite 型置信区间</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2021.11.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table21-12)缩减模型方差分量的 I 型平方和估计。非残差方差的 Wald 型置信区间</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2021.12.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table21-13)方差分量的 REML 估计及其计算中所用的自由度，包括区间的标准差和平方根</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2021.13.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table21-14)基于 I 型平方和，使用第 \@ref(sec20-2-2) 节中的一般 Satterthwaite 近似计算置信区间</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2021.14.png">
</td>
  </tr>
</tbody>
</table>

接下来，使用通常的 Satterthwaite 近似值，根据 I 型平方和的信息，为方差分量提供置信区间。使用 \@ref(sec20-2-1) 节中的方法构造关于 ${{\sigma}_{\varepsilon}^2}$ 的置信区间。${{\sigma}_{\varepsilon}^2}$ 的 95% 置信区间为

$$\frac{82[MS(Residual)]}{\chi_{(\alpha/2),82}^2}\leq\sigma_\varepsilon^2\leq\frac{82[MS(Residual)]}{\chi_{1-(\alpha/2),82}^2}$$

或

$$\frac{82(4.9831)}{108.937}\leq\sigma_\varepsilon^2\leq\frac{82(4.9831)}{58.8446}$$

或

$$3.781\leq\sigma_\varepsilon^2\leq6.944$$

上述关于 ${{\sigma}_{\varepsilon}^2}$ 的 95% 置信区间与表 \@ref(tab:table21-11) 中 Residual 的置信区间相同。

使用表 \@ref(tab:table21-7) 中期望均方的 ${{\sigma}_{sw}^2}$ 的 Satterthwaite 置信区间，${{\sigma}_{sw}^2}$ 的矩法解（和估计）可以表示为

$$\hat{\sigma}_{sw}^2=\frac{MS[worker\times site(plant)]}{3.064}-\frac{MS(Residual)}{3.064}=30.7388$$

那么 $r_{sw}\hat{\sigma}_{sw}^2/\sigma_{sw}^2$ 近似分布为自由度为 $r_{sw}$ 的卡方随机变量，其中

$$r_{_{sw}}=\frac{(\hat{\sigma}_{sw}^2)^2}{\left[\frac{MS[worker\times site(plant)]}{3.064}\right]\Big/27+\left[\frac{MS(Residual)}{3.064}\right]\Big/82}=24.33$$

${{\sigma}_{sw}^2}$ 的一个近似 95% 置信区间为

$$\frac{r_{sw}\hat{\sigma}_{sw}^2}{\chi_{.025,24.33}^2}\leq\sigma_{sw}^2\leq\frac{r_{sw}\hat{\sigma}_{sw}^2}{\chi_{.975,24.33}^2}$$

或

$$18.80\leq\sigma_{sw}^2\leq59.17$$

<table>
<caption>(\#tab:table21-15)缩减模型方差分量假设的 I 型检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2021.15.png">
</td>
  </tr>
</tbody>
</table>

为了构造关于 $\sigma_{sw}^2$ 的置信区间，将 $\hat\sigma_{sw}^2$ 的解表示为（从表 \@ref(tab:table21-15) 的 Error 项中获得信息）

$$\begin{aligned}\hat{\sigma}_w^2&=\{MS(worker)-1.2732MS[worker\times site(plant)]+0.2732MS(Residual)\}/13.037\\&=25.4692\end{aligned}$$

与 $\hat\sigma_{sw}^2$ 相关的近似自由度为 $r_w=3.12$，这给出近似 95% 置信区间 $[8.297\le \hat\sigma_{sw}^2 \le 327.87]$.

为了获得关于 $\sigma_{p}^2$ 的 95% 置信区间，可以使用表 \@ref(tab:table21-15) 中的结果将 $\sigma_{p}^2$ 的估计表示为

$$
\begin{aligned}
\hat\sigma_p^2=&\{MS(plant)-1.0076MS[worker(plant)]+0.0012MS[worker\times site(plant)\\&+0.0065 MS(Residual)\}/38.941\\
=&47.60
\end{aligned}
$$

与 $\hat\sigma_{p}^2$ 相关的自由度为 $r_p=1.27$. 由此产生的近似 95% 置信区间为 $[10.58\le \sigma_{p}^2 \le 12078.74]$. 使用 Satterthwaite 近似计算的置信区间集如表 \@ref(tab:table21-16) 所示。

<table>
<caption>(\#tab:table21-16)基于 I 型平方和的通常 Satterthwaite 近似置信区间</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2021.16.png">
</td>
  </tr>
</tbody>
</table>

分析的最后阶段是使用第 \@ref(sec20-2) 节中的置信域法来构建方差分量 $\sigma^2_p,\sigma^2_w,\sigma^2_{sw}$ 的置信区间。$\sigma^2_{sw}$ 95% 置信区间的端点 $c_{sw},d_{sw}$ 为

$$c_{sw}=\frac{27MS[worker\times site(plant)]/\chi_{.025,27}^2-82MS(Residual)/\chi_{.975,82}^2}{3.064}=17.96$$

和

$$d_{sw}=\frac{27MS[worker\times site(plant)]/\chi_{975,27}^2-82MS(Residual)/\chi_{.025,82}^2}{3.064}=58.74$$

将这两个表达式简化为一个得到置信区间 $17.96\le \sigma^2_{sw} \le 58.74$.

关于 $\sigma^2_{w}$ 的 95% 置信区间的两个端点是 $c_{w},d_{w}$，其中

$$c_w=\frac{6MS[worker(plant)]/\chi_{.025,6}^2-r_wQ_w^*/\chi_{.975,26.42}^2}{13.037}=-3.335$$

和

$$d_w=\frac{6MS[worker(plant)]/\chi_{.975,6}^2-r_wQ_w^*/\chi_{.025,26.42}^2}{13.037}=164.00$$

其中 $r_w=26.42$ 且 $Q^*_w=124.9085$. 简化这些表达式可以得到区间 $0\leq\sigma_w^2\leq164.00$.

关于 $\sigma^2_{p}$ 的 95% 置信区间端点是 $c_p,d_p$，其中

$$c_p=\frac{2MS(plant)]/\chi_{.025,2}^2-r_pQ_p^*/\chi_{.975,5.9961}^2}{38.941}=-41.2471$$

和

$$d_p=\frac{2MS(plant)]/\chi_{.975,2}^2-r_pQ_p^*/\chi_{.025,5.9961}^2}{38.941}=2341.94$$

其中 $r_p=5.9961$ 且 $Q^*_p=460.2638$. 简化以得到区间 $0\leq\sigma_p^2\leq 2431.94$.

表 \@ref(tab:table21-17) 展示了使用置信域法计算的置信区间。

<table>
<caption>(\#tab:table21-17)使用第 \@ref(sec20-2-5) 节置信域法计算的置信区间总结</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2021.17.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table21-18) 总结了上述四种计算置信区间的方法，其中表中的条目是各自方差分量的下限和上限。

<table>
<caption>(\#tab:table21-18)四种类型的 95% 置信区间及其上限和下限总结</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2021.18.png">
</td>
  </tr>
</tbody>
</table>

## 使用 JMP 进行计算 {#sec21-6}

图 \@ref(fig:figure21-1) 包含了示例 [21](#chap21) 中 JMP^®^ 数据表的一部分，其中 plant, site 和 worker 被指定为名义变量。图 \@ref(fig:figure21-2) 包含 fit model 界面，其中 plant, site(plant), worker(plant) 和 site×worker(plant) 被指定为随机效应，并选择 REML 方法。图 \@ref(fig:figure21-3) 显示了完全模型的方差分量估计，包括标准误估计和 95% 置信区间，这与表 \@ref(tab:table21-2) 中的 REML 估计相对应。图 \@ref(fig:figure21-4) 显示了缩减模型的 fit model 界面，其中从图 \@ref(fig:figure21-2) 的 fit model 界面中移除了 site(plant) 项。图 \@ref(fig:figure21-5) 中 JMP 的 REML 估计与表 \@ref(tab:table21-11) 显示的 SAS-Mixed 估计相同。JMP 中的 fit model 过程为随机效应模型提供了出色的结果。

<div class="figure" style="text-align: center">
<img src="figure/figure 21.1.png" alt="示例 [21](#chap21) JMP 数据表的一部分" width="730" />
<p class="caption">(\#fig:figure21-1)示例 [21](#chap21) JMP 数据表的一部分</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 21.2.png" alt="示例 [21](#chap21) 完全模型的 JMP fit model 界面" width="984" />
<p class="caption">(\#fig:figure21-2)示例 [21](#chap21) 完全模型的 JMP fit model 界面</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 21.3.png" alt="示例 [21](#chap21) 的 JMP REML 结果" width="837" />
<p class="caption">(\#fig:figure21-3)示例 [21](#chap21) 的 JMP REML 结果</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 21.4.png" alt="示例 [21](#chap21) 缩减模型的 JMP fit model 界面" width="838" />
<p class="caption">(\#fig:figure21-4)示例 [21](#chap21) 缩减模型的 JMP fit model 界面</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 21.5.png" alt="示例 [21](#chap21) 缩减模型的 JMP REML 结果" width="886" />
<p class="caption">(\#fig:figure21-5)示例 [21](#chap21) 缩减模型的 JMP REML 结果</p>
</div>

## 结束语 {#sec21-7}  

[^quitevariable]: 原文：The REML and usual Satterthwaite approximation are quite similar where the general Satterthwaite intervals using the type I sums of squares and confidence regions for the plant and worker variance components are quite variable.

通过这四种方法获得的置信区间存在很大差异，其中方差分量自由度较小的情况下的差异最大。REML 和通常的 Satterthwaite 近似值非常相似，其中对于工厂和工人方差分量，使用 I 型平方和的一般 Satterthwaite 区间和置信域是非常多变的[^quitevariable]。在本章中，使用第 \@ref(chap19) 章中描述的六种估计方法详细分析了复杂的非均衡随机效应模型，以估计方差分量。展示了第 \@ref(chap20) 章中描述的检验假设和构造置信区间的方法。发现原始模型中的一个方差分量不显著，因此采用模型构建方法提供了充分模型，其中其余方差分量的估计是有意义的。

为了确定不同方法之间可能存在多少差异，应该至少使用两种方法来估计方差分量和构造置信区间，以分析非均衡数据集。如果两种方法都得出相似的结果，那么使用哪种方法都无所谓，可以信任其中任何一种方法。**如果结果差异很大，那么应该调查设计的哪些特性导致了这些差异**。分析结果使用 SAS-Mixed 和 JMP 的 fit model 过程进行描述。

## 练习 {#sec21-8}
