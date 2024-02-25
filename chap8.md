
# 案例研究：均衡双向实验的完整分析 {#chap8}

在上一章中，我们假设当处理 T 和处理 B 在实验中发生交互时，实验者将想要比较在处理 B 的每一种可能性下处理 T 的效应，反之亦然。在许多情况下，交互作用并非在实验中随处发生——通常只有一两种处理组合会导致交互作用。在其他情况下，一种处理的可能性可能与第二种处理的可能性发生交互，而第一种处理的所有其他可能性都不与第二种处理的任何可能性发生交互。

为了对具有交互作用的数据进行更完整的分析，确定交互作用在数据中的位置是有帮助的。例如，如果已知实验中的所有交互作用都是由处理 T 的其中一种可能性（通常是对照）引起的，那么在对处理 B 的所有可能性进行平均后，仍然可以比较处理 T 的所有其他可能性。这将导致更大功效的检验，用于比较不与处理 B 交互的处理 T 的可能性。

## 主效应均值对比 {#sec8-1}

通常，处理 T 和处理 B 的可能性会提示进行主效应对比，这对实验者来说会特别有趣。这样的主效应对比会提示特殊类型的交互作用对比，这对实验者来说也应该很有趣。此外，这样的交互作用对比通常很容易解释。接下来定义了主效应对比和正交主效应对比。

::: {.definition #8-1}
♦

若 $\sum_{i=1}^tc_i=0$，则称 $\bar \mu_{i\cdot}$ 的线性组合 $\sum_{i=1}^tc_i\bar{\mu}_{i\cdot}$ 为 T 主效应均值的一个**对比** (contrast).同样，若 $\sum_{i=j}^td_j=0$，则称 $\bar \mu_{\cdot j}$ 的线性组合 $\sum_{j=1}^bd_j\bar{\mu}_{\cdot j}$ 为 B 主效应均值的一个对比。
:::

::: {.definition #8-2}
♦

若 $\sum_{i=1}^tc_ic_i^{\prime}=0$，则称两对比 $\sum_{i=1}^tc_i\bar{\mu}_{i\cdot}$ 和 $\sum_{i=1}^tc_i^\prime \bar{\mu}_{i\cdot}$ 为 T 主效应均值的**正交对比** (orthogonal contrasts). 同样，若 $\sum_{j=1}^tc_jc_j^{\prime}=0$，则称两对比 $\sum_{j=1}^tc_j\bar{\mu}_{\cdot j}$ 和 $\sum_{j=1}^tc_j^\prime \bar{\mu}_{\cdot j}$ 为 B 主效应均值的正交对比。
:::

现假设

$$\begin{aligned}S_T&=\left\{\sum_{i=1}^tc_{i1}\bar{\mu}_{i\cdot},\sum_{i=1}^tc_{i2}\bar{\mu}_{i\cdot},\ldots,\sum_{i=1}^tc_{it-1}\bar{\mu}_{i\cdot}\right\}\end{aligned}$$

为 T 主效应均值的一组 t-1 个正交对比，以及

$$\begin{aligned}S_B&=\left\{\sum_{j=1}^bc_{j1}\bar{\mu}_{\cdot j},\sum_{j=1}^bc_{j2}\bar{\mu}_{\cdot j},\ldots,\sum_{j=1}^bc_{\cdot j}\bar{\mu}_{j\cdot}\right\}\end{aligned}$$

为 B 主效应均值的一组 b-1 个正交对比。$S_T$ 和 $S_B$ 中的每一个都表明两个主效应平方和的**划分** (partitioning)。即，T 的平方和的划分定义为

\begin{equation}
S_T^*=\left\{Q_p^2=\frac{nb(\sum_ic_{ip}\bar{y}_{i\cdot\cdot})^2}{\sum_ic_{ip}^2},\quad p=1,2,\ldots,t-1\right\}
(#eq:8-1)
\end{equation}

B 的平方和的划分定义为

\begin{equation}
S_B^*=\left\{Q_q^2=\frac{nt(\sum_jd_{jq}\bar{y}_{\cdot j\cdot})^2}{\sum_jd_{jq}^2},\quad q=1,2,\ldots,b-1\right\}
(#eq:8-2)
\end{equation}

也就是说，$S^*_T$ 中的每个 $Q_p^2$ 都是一个单自由度平方和，用于检验 T 主效应均值中相应的对比是否等于 0. 此外，$S^*_T$ 中 t-1 个单自由度对比的总和等于用于检验 $H_{01}\colon\bar{\mu}_{1\cdot}=\bar{\mu}_{2\cdot}=\cdots=\bar{\mu}_{t\cdot}$ 的平方和。$S^*_B$ 的元素也存在类似的情况。

不应过分强调获得基本平方和的正交划分的可取性。从数学角度看，正交划分是不错的，但从实际的角度来看可能并不那么好。通常，一组精心选择的正交对比将使实验者能够明智、清晰、完整地解释他/她的数据。然而，实验者应该真正考虑任何以及所有有意义的对比，不应过度关注所选择的对比是否正交。

## 交互对比 {#sec8-2}

本节从双向实验中交互对比的定义开始。

若对于 $j=1,2,\cdots,b$ 有 $\sum_i\omega_{ij}=0$ 且若对于 $i=1,2,\cdots,t$ 有 $\sum_j\omega_{ij}=0$，则称 $\mu_{ij}$ 的线性组合 $\sum_i\sum_j\omega_{ij}\mu_{ij}$ 为一个**交互对比** (interaction contrasts).

双向实验主效应中的对比产生了特殊类型的交互对比。设 $\sum_{i=1}^tc_i\bar{\mu}_{i\cdot}$ 为 T 主效应均值的一个对比以及 $\sum_{j=1}^bd_j\bar{\mu}_{\cdot j}$ 为 B 主效应均值的一个对比，则 $\sum_i\sum_jc_id_j{\mu}_{i{j}}$ 为交互对比。即，对于所有的 $i,j$，取 $\omega_{ij}=c_id_j$ 即得到一个交互对比。若 $\sum_i\sum_j\omega_{ij}\omega_{ij}^{\prime}=0$，则两交互对比 $\sum_i\sum_j\omega_{ij}\mu_{ij}$ 和 $\sum_i\sum_j\omega_{ij}^{\prime}\mu_{ij}$ 为**正交交互对比** (orthogonal interaction contrasts). 两组主效应的正交对比导致交互效应的正交对比。设 $\sum_{i=1}^tc_i\bar{\mu}_{i\cdot}$ 和 $\sum_{i=1}^tc_i^{\prime}\bar{\mu}_{i\cdot}$ 为 T 主效应均值的两个对比以及 $\sum_{j=1}^bd_j\bar{\mu}_{\cdot j}$ 和 $\sum_{j=1}^bd_j^{\prime}\bar{\mu}_{\cdot j}$ 为 B 主效应均值的两个对比，若 $\sum_{i=1}^tc_ic_i^{\prime}=0$ **[或]{.underline}** $\sum_{j=1}^bd_jd_j^{\prime}=0$，则 $\sum_i\sum_jc_id_j\mu_{ij}$ 和 $\sum_i\sum_jc_i^{\prime}d^{\prime}_j\mu_{ij}$ 为正交交互对比。也就是说，如果主效应对比对中的至少一个是正交对比对，则两个交互对比将彼此正交。接下来假设 $S^*_T$ 和 $S^*_B$ 由 \@ref(sec8-1) 节定义，令 $S^*_{T\times B}$ 定义为

\begin{equation}
S_{T\times B}^*=\left\{Q_{pq}^2=\frac{n(\sum_i\sum_jc_{ip}d_{jq}\bar{y}_{ij\cdot})^2}{\sum_ic_i^2\sum_jd_j^2},\quad p=1,2,\ldots,t-1;q=1,2,\ldots,b-1\right\}
(#eq:8-3)
\end{equation}

则 $S_{T\times B}^*$ 定义了交互作用平方和的一个划分。也就是说，$\sum_p\sum_qQ_{pq}^2=T{\times}B$ 个交互作用平方和，以及式 \@ref(eq:8-3) 中的所有 (t-1)(b-1) 个单自由度平方和具有独立的概率分布。下一节将讨论一个示例，以说明上一节中描述的思想。

## 油漆铺路示例 {#sec8-3}

考虑表 \@ref(tab:table8-1) 中的实验，该实验给出了每种 paint × paving 处理组合的三次独立重复的平均值。进行这项实验是为了比较由两家不同公司生产的两种颜色的油漆 (paints) 在三种类型路面 (paving surfaces) 上的使用寿命，以周为单位。该实验的误差平方和为 455.04，具有 24 个自由度，因此 $\hat \sigma^2=18.96$. 本实验通常的方差分析表如表 \@ref(tab:table8-2) 所示。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;border-bottom: 0;">
<caption>(\#tab:table8-1)油漆铺路单元格均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%208.1.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> 颜色后跟随的 I 和 II 表示生产厂家，Asphalt I, Asphalt II, Concrete 表示三种类型的路面</td></tr></tfoot>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table8-2)油漆铺路示例方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%208.2.png">
</td>
  </tr>
</tbody>
</table>

每种油漆和铺路 (paving) 处理的可能值都会在主效应中产生两组正交对比，这可能是令人感兴趣的。这些在表 \@ref(tab:table8-3) 中给出。这两组主效应的正交对比表明，在交互效应中存在六个正交对比。这些数据如表 \@ref(tab:table8-4) 所示。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table8-3)油漆铺路实验主效应对比</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%208.3.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;border-bottom: 0;">
<caption>(\#tab:table8-4)油漆铺路数据交互作用假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%208.4.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> Type 是指 asphalt vs concrete</td></tr></tfoot>
</table>

表\@ref(tab:table8-5) 中给出了使用表 \@ref(tab:table8-3) 和 \@ref(tab:table8-4) 中表明的划分对这些数据进行更完整的分析。下面给出了从表 \@ref(tab:table8-5) 中获得三个选定对比的平方和所需的计算细节。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table8-5)油漆铺路实验方差分析表，包括单自由度检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%208.5.png">
</td>
  </tr>
</tbody>
</table>

根据式 \@ref(eq:8-1)，用于比较两种白色油漆均值的单自由度平方和为，

$$Q_1^2=\frac{3\cdot3[(1)29.0+(-1)35.0]^2}{(1)^2+(-1)^2}=162.0$$

根据式 \@ref(eq:8-2)，用于比较沥青 (asphalt) 和混凝土 (concrete) 的单自由度平方和为，

$$Q_2^2=\frac{3\cdot4[(1)(26.5)+(1)(27.5)+(-2)(29.5)]^2}{1^2+1^2+(-2)^2}=40.5$$

根据式 \@ref(eq:8-3)，而用于比较 type × type 交互作用的单自由度平方和为，

$$Q_3^2=\frac{3\left[(1)(30)+(1)(28)+(-2)(29)+(-1)(34)+(-1)(35)+(2)(36)\right]^2}{\left(1\right)^2+\left(1\right)^2+\left(-2\right)^2+\left(-1\right)^2+\left(-1\right)^2+\left(2\right)^2}=2.25$$

根据表 \@ref(tab:table8-5) 的分析，可以得出以下结论：

1.  实验中的所有交互作用都是由两种黄色油漆在两种类型的表面上的不同作用引起的，因为其交互对比是交互作用中唯一显著的单自由度平方和。
2.  因为我们现在知道了数据中交互作用的位置，我们可以进行以下观测：
    a)  由于沥青和油漆之间没有交互作用，因此可以在对所有油漆进行平均后对这两种沥青进行比较。该比较的 $F$ 统计量的值为 $F=0.32$；因此，沥青 I 和 II 之间没有显著差异。
    b)  由于白色油漆和三种路面之间没有交互作用，因此可以在对所有路面进行平均后对两种白色涂料进行比较。该比较的 $F$ 统计量值为 $F=8.54$，这表明白色油漆 II 与 I 有显著差异。从表 \@ref(tab:table8-1) 中我们可以看出，白色油漆 II 的持续时间更长。
    c)  尽管比较黄色油漆 I 和黄色油漆 II 的统计量是显著的（$F=4.45$），但在进行解释时必须小心，因为黄色油漆的品牌和路面类型之间存在显著的交互作用。
    d)  尽管沥青与混凝土的比较并不显著（$F=2.14$），但在进行解释时必须再次小心，因为油漆品牌和路面类型之间存在显著的交互作用。
3. 为了完成对这些数据的分析，我们还应该检查：
    a)  沥青上的黄色油漆 I 与黄色油漆 II，即检验 $\mu_{11}+\mu_{12}-\mu_{21}-\mu_{22}=0$.
    b)  混凝土上的黄色油漆 I 与黄色油漆 II.
    c)  对于黄色油漆 I，比较混凝土和沥青。
    d)  对于黄色油漆 II，比较混凝土和沥青。
    e)  三种路面上的白色油漆。

结果如表 \@ref(tab:table8-6) 所示。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table8-6)结论 3 中的假设检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%208.6.png">
</td>
  </tr>
</tbody>
</table>

对表 \@ref(tab:table8-6) 中的结果和表 \@ref(tab:table8-1) 中的平均值的检查表明，1）黄色油漆 II 在沥青上明显优于黄色油漆 I，但 2）黄色油漆 I 在混凝土上明显优于黄色油漆 II；3） 黄色油漆 II 在混凝土上的持续时间明显长于在沥青上；4） 黄色油漆 II 在沥青上的持续时间明显长于在混凝土上；最后，5）白色油漆在所有三种路面上持续的时间大致相同。

我们对该示例的分析所获得的所有结果都可以使用许多统计计算包获得，如 SAS^®^ 和 SPSS，使用它们的对比选项。其中一些程序可能要求将 12 种处理组合视为单向处理结构。

## 分析定量处理因素 {#sec8-4}

在本节中，假设实验的两个因素的水平都是定量的 (quantitative). 在这种情况下，可以定义对比以衡量每组处理主效应均值的曲线趋势。感兴趣的趋势通常是线性的、二次的、三次的等等。相应的正交对比将主效应平方和划分为测量线性、二次、三次等等趋势的效应，然后可以用来构建交互效应的正交对比。由此产生的对比称为 Lin T × Lin B（T 的线性效应乘以 B 的线性效应）、Lin T × Quad B （T 的线性效应乘以 B 的二次效应）等等。对于 3×4 的实验，其中两种处理具有相等的间隔水平（例如，因子 T 为 5, 10, 15，因子 B 为 2, 4, 6, 8），Lin T 对比由 $-\bar{\mu}_{1\cdot}+0\bar \mu_{2\cdot}+\bar{\mu}_{3\cdot}=0$ 定义，Lin B 对比由 $-3\bar{{\mu}}_{\cdot1}-\bar{{\mu}}_{\cdot2}+\bar{{\mu}}_{\cdot3}+3\bar{{\mu}}_{\cdot4}=0$.

注意，对于 Lin T 对比，在 $\sum_{i=1}^3c_i\bar{\mu}_{i\cdot}=0$ 中 $c_1=-1,c_2=0,c_3=1$. 对于 Lin B 对比，在 $\sum_{j=1}^4d_j\bar{\mu}_{\cdot j}=0$ 中 $d_1=-3,d_2=-1,d_3=1,d_4=3$. 因此 Lin T × Lin B 交互对比定义为 $3\mu_{11}+\mu_{12}-\mu_{13}-3\mu_{14}-3\mu_{31}-\mu_{32}+\mu_{33}+3\mu_{34}=0$. 在正交多项式系数表中可以找到在主效应均值中用于定义这些类型的正交对比的 $c_i$ 和 $d_j$ 的值。参见 Beyer (1968). 表 \@ref(tab:table8-7) 中再现了 3×4 实验的结果。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table8-7)3×4 实验的正交对比系数</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%208.7.png">
</td>
  </tr>
</tbody>
</table>

设 $x_1,x_2,\cdots,x_t$ 表示因素 T 的水平，$z_1,z_2,\cdots,z_b$ 表示因素 B 的水平。总是存在参数 $\alpha_{kh},k=0,1,2,\cdots,t-1,h=0,1,2,\cdots,b-1$ 使得单元格均值参数可以表示为 $x_i$ 和 $z_j$ 的多项式函数。也就是说，存在这样的 $\alpha_{kh}$ 使得

\begin{equation}
\mu_{ij}=\sum_{k=0}^{t-1}\sum_{h=0}^{b-1}\alpha_{kh}x_i^kz_j^h
(#eq:8-4)
\end{equation}

展开式 \@ref(eq:8-4)，可以得到

$$\begin{aligned}\mu_{ij}=&\alpha_{00}+\alpha_{10}x_i+\alpha_{20}x_i^2+\cdots+\alpha_{t-10}x_i^{t-1}+\alpha_{01}z_j+\alpha_{02}z_j^2+\cdots+\alpha_{0b-1}z_j^{b-1}\\&+\alpha_{11}x_iz_j+\alpha_{12}x_iz_j^2+\cdots+\alpha_{t-1b-1}x_i^{t-1}z_j^{b-1}\end{aligned}$$

根据方程 \@ref(eq:8-4) 中的 $\alpha_{kh}$ 值，表 \@ref(tab:table8-8) 给出了对于 3×4 实验的主效应期望值和交互对比。在构建此表时，我们假设 $x$ 的三个水平分别编码为 -1、0 和 1，$z$ 的四个水平分别编码为 -3, -1, 1 和 3. 提供表 \@ref(tab:table8-8) 的目的是指出当研究这种类型的对比时，实际上正在检验哪些假设。例如，Lin T 与 Lin B 的对比检验了假设 $40\alpha_{11}+328\alpha_{13}=0$ . 因此，如果此效应显著，那么可能是因为 $\alpha_{11} ≠ 0$ 或者 $\alpha_{13} ≠ 0$，而不仅仅是因为 $\alpha_{11} ≠ 0$（如**许多数据分析师可能认为的那样**）。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table8-8)3×4 实验中正交多项式的期望值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%208.8.png">
</td>
  </tr>
</tbody>
</table>

如果要考察正交多项式，建议先观察最高阶项的系数，然后按照降序考虑剩余项的阶数。一旦确定某个高阶项在模型中，那么所有两个分量阶数均低于该显著项的项，无论它们是否显著，也都应该被包含在模型中。例如，如果决定 $x^2z^2$ 应包括在模型中，那么模型也应该包括 $xz^2,z^2,x^2z,x^2,xz,x,z$. 这是因为正交多项式总是表示定量变量的编码值。因此，如果 $\alpha_{22}$ 非零，这意味着

$$\left(\frac{x-h_1}{c_1}\right)^2\left(\frac{z-h_2}{c_2}\right)^2$$

在模型中，其中 $\left(\frac{x-h_1}{c_1}\right)$ 和 $\left(\frac{z-h_2}{c_2}\right)$ 为 $x$ 和 $z$ 的编码值。式 $\left(\frac{x-h_1}{c_1}\right)^2\left(\frac{z-h_2}{c_2}\right)^2$ 的展开说明项 $xz^2,z^2,x^2z,x^2,xz,x,z$ 也在模型中，即使其他低阶正交多项式可能不显著。

## 多重检验 {#sec8-5}

第 \@ref(chap3) 章中讨论的任何多重比较程序都可以使用，只需进行一些非常小的调整，就可以对双因素实验的主效应进行多重比较。调整要求将 $n$ 和 $n_i$ 替换为：为估计被比较的主效应均值而进行的观测总数。在本章中，B 主效应均值的样本量为 $nt$，T 主效应均值的样本量为 $nb$. 我们对主效应均值多重比较的建议与第 \@ref(sec3-2) 节中给出的建议相同。第 \@ref(chap3) 章中给出的容易推广到交互对比的程序是 LSD procedure、Bonferroni's method、多元 $t$ 方法、模拟方法和 Scheffé's procedure. 我们发现，Scheffé's procedure 不是很令人满意，因为所需的临界点太大，程序太保守。我们对交互对比多重比较的建议如下：

1) 对交互进行 $F$ 检验。
2) 如果 $F$ 统计量显著，则使用 LSD procedure（或等价地，上一节中给出的对比程序）进行任何计划的比较。对于数据窥探和计划外比较，请使用 Johnson (1976) 给出的程序，此处不讨论。
3) 如果交互作用的 $F$ 检验不显着，实验者仍应检查她计划考虑的任何单个交互对比，但要使用多元 $t$ 方法或 Bonferroni's method. 只要所选对比是线性独立的，就使用多元 $t$ 方法；否则应使用 Bonferroni's method.

## 结束语 {#sec8-6}

在本章中，我们通过举例介绍了从实验中获得最大信息量的方法。其中包括发现实验中交互作用发生位置的方法。了解实验中交互作用发生的位置对于确定可能提出的问题的最佳答案很有价值。本章介绍的技术应该可以帮助实验者更好地分析他们的实验。还考虑了定量处理因素的分析，包括如何确定哪些趋势可能与不同水平的处理因素有关。

## 练习 {#sec8-7}
