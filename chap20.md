
# 方差分量的推断方法 {#chap20}

当研究者设计一个涉及随机效应因素的实验时，她通常希望根据模型中指定的特定方差分量做出推断。特别是，如果 $\sigma_u^2$ 是对应于因素 A 水平分布的方差分量，实验者可能希望确定是否有足够的证据得出 $\sigma_u^2>0$ 的结论。可以通过以下方式做出适当的决定：1）检验假设 $H_0{:}\sigma_u^2=0\mathrm{~vs~}{H_a}{:}\sigma_u^2>0$；2）构造关于 $\sigma_u^2$ 的置信区间；或 3）为 $\sigma_u^2$ 构造置信下限。本章讨论了随机效应模型的这些推断程序，其中假设检验方法在第 \@ref(sec20-1) 节中描述，置信区间（下界）的构造在第 \@ref(sec20-2) 节中描述。方差分量的置信区间构造一直是研究的一个活跃领域，许多作者已经为方差分量参数的特定函数开发了专门的置信区间。本章中描述的方法可在现行软件中使用，有些可用于特定问题。讨论并不全面，而是指出了已经解决的置信区间类型。在 Burdick and Graybill (1992) 以及当前统计期刊的论文中，有更为完整的讨论。

## 假设检验 {#sec20-1}

有两种基本技术来检验关于方差分量的假设。第一种技术使用方差分析表中的平方和来构造 $F$ 统计量。对于大多数均衡模型，$F$ 统计量的分布与 $F$ 分布完全相同，而对于不均衡模型，分布由 $F$ 分布近似，随着设计变得越来越不均衡，近似值变得越来越差。第二种技术基于似然比检验 (likelihood rato test)，渐近分布于卡方分布。对于均衡设计，$F$ 统计量方法可能比似然比检验更好，而对于不均衡设计，没有明确的选择。读者可能希望在决定使用哪种方法来检验感兴趣的假设之前，使用与感兴趣的数据集类似的数据结构进行模拟实验，以研究检验统计量的分布。

### 使用方差分析表 {#sec20-1-1}

[^multiplesofchi-square]: 原文："then sums of squares obtained by the usual analysis of variance are independently distributed as scalar multiples of chi-square random variables."

如果数据集是均衡的，那么通过常规方差分析获得的平方和独立分布于卡方随机变量的标量倍数[^multiplesofchi-square]。令 $Q$ 表示基于自由度 v 的平方和，其中其期望均方是四个方差分量的函数。也就是说，假设

$$E(Q/v)=\sigma_\varepsilon^2+k_1\sigma_1^2+k_2\sigma_2^2+k_3\sigma_3^2$$

然后，假设数据服从正态分布，

$$W=\frac Q{\sigma_\varepsilon^2+k_1\sigma_1^2+k_2\sigma_2^2+k_3\sigma_3^2}$$

通常分布为具有 n 个自由度的卡方随机变量。对于许多形如 $H_0{:}\sigma_1^2=0\mathrm{~vs~}{H_a}{:}\sigma_1^2>0$ 的假设，有两个独立的平方和，用 $Q_1$ 和 $Q_2$ 表示，分别基于 v~1~ 和 v~2~ 自由度，期望为

$$\begin{aligned}E(Q_1/n_1)&=\sigma_\varepsilon^2+k_1\sigma_1^2+k_2\sigma_2^2+k_3\sigma_3^2\\E(Q_2/n_2)&=\sigma_\varepsilon^2+k_2\sigma_2^2+k_3\sigma_3^2\end{aligned}$$

假设 $H_0{:}\sigma_1^2=0\mathrm{~vs~}{H_a}{:}\sigma_1^2>0$ 等价于

$$\begin{aligned}H_0\colon E\left(Q_1/v_1\right)&=E(Q_2/v_2)\text{ vs }H_a\colon E(Q_1/v_1)>E(Q_2/v_2)\end{aligned}$$

用于检验该假设的统计量为 $F = (Q_1/v_1)/(Q_2/v_2)$，在 $H_0$ 条件下，该统计量通常分布为具有 v~1~ 和 v~2~ 自由度的中心 $F$ 分布。对于较大的 $F$ 值，该假设被拒绝。此过程涉及获取平方和，然后使用其期望均方来确定每个感兴趣假设的适当除数。以下两个示例演示了此程序。

### 示例 20.1：完全随机设计结构中的双向随机效应检验统计量 {#sec20-1-2}

在完全随机设计结构中，两个因素均为随机的双向处理结构的模型为

$$y_{ijk}=\mu+a_i+b_j+c_{ij}+\varepsilon_{ijk}\quad\mathrm{~for~}i=1,2,\ldots,a,j=1,2,\ldots,b,\mathrm{~and~}k=1,2,\ldots,n$$

其中 $a_i\thicksim i.i.d.\,N(0,\sigma_a^2),b_j\thicksim i.i.d.\,N(0,\sigma_b^2),c_{ij}\thicksim i.i.d.\,N(0,\sigma_c^2),\varepsilon_{ijk}\thicksim i.i.d.\,N(0,\sigma_\varepsilon^2)$，且随机变量 $a_i,b_j,c_{ij},\varepsilon_{ijk}$ 独立分布。

表 \@ref(tab:table20-1) 显示了模型的方差分析表以及模型的平方和和期望均方。通过检查期望均方来选择适当的分子和分母来构造检验统计量。用于检验假设 $H_0{:}\sigma_a^2=0\mathrm{~vs~}{H_u}{:}\sigma_a^2>0$ 的统计量是通过在 $MSA$ 的期望均方中设置 $\sigma_a^2=0$ 来构造的。接下来，当 $H_0$ 为真时，找到与 $MSA$ 具有相同期望均方的另一个均方，并使用该均方作为除数。要检验 $H_0{:}\sigma_a^2=0\mathrm{~vs~}{H_u}{:}\sigma_a^2>0$，适当的除数是 $MSAB$；要检验 $H_0{:}\sigma_b^2=0\mathrm{~vs~}{H_u}{:}\sigma_b^2>0$，适当的除数是 $MSAB$；要检验 $H_0{:}\sigma_c^2=0\mathrm{~vs~}{H_u}{:}\sigma_c^2>0$，适当的除数是 $MSResidual$. 决策规则是拒绝 $H_0{:}\sigma_a^2=0\mathrm{~vs~}{H_u}{:}\sigma_a^2>0$ 如果 $F = MSA/MSAB>F_{\alpha,(a-1),(a-1)(b-1)}$，其中 $\alpha$ 为所选的 I 类错误率。对于 $\sigma_b^2$ 和 $\sigma_c^2$ 可以类似地确定检验统计量。表 \@ref(tab:table20-1) 包含假设列表和相应的检验统计量。很可能，当 $F$ 统计量不超过指定的分位数时，结论并不是方差分量为零，而是与系统中其他变异源相比，方差分量的大小可以忽略不计。

<table>
<caption>(\#tab:table20-1)示例 [20.1](#sec20-1-2) 双向随机效应模型方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2020.1.png">
</td>
  </tr>
</tbody>
</table>

### 示例 20.2：复杂三向随机效应检验统计量 {#sec20-1-3}

<table>
<caption>(\#tab:table20-2)示例 [20.2](#sec20-1-3) 数据</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2020.2.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table20-2) 中的数据来自一个设计，其中处理结构中的三个因素的水平是随机效应，A 的水平与 B 的水平交叉，C 的水平嵌套在 B 的水平内，所有这些都在一个完全随机设计结构中。一个描述与表 \@ref(tab:table20-2) 结构相似的大型数据集的一般模型是

$$\begin{aligned}&y_{ijkm}=\mu+a_i+b_j+(ab)_{ij}+c_{k(j)}+(ac)_{ik(j)}+\varepsilon_{ijkm}\\&\mathrm{for~}i=1,2,\ldots,a,j=1,2,\ldots,b,k=1,2,\ldots,c,\mathrm{~and~}m=1,2,\ldots,n\end{aligned}$$

参数 $\mu$ 表示总体均值，$a_i$ 表示因子 A 水平 i 的效应，$b_j$ 表示因子 B 水平 j 的效应，$(ab)_{ij}$ 表示因子 A 和因子 B 水平之间的交互作用，$c_{k(j)}$ 表示嵌套在因子 B 第 j 个水平内的因子 C 水平 k 的效应，$(ac)_{ik(j)}$ 表示因子 A 与嵌套在因子 B 内的因子 C 水平之间的交互作用，以及 $\varepsilon_{ijkm}$ 表示实验单元或抽样误差。在理想条件下，$a_i\thicksim i.i.d.~N(0,\sigma_a^2),~b_j\thicksim i.i.d.~N(0,~\sigma_b^2),~(ab)_{ij}\thicksim i.i.d.~N(0,~\sigma_{ab}^2),~(ac)_{ik(j)}\thicksim i.i.d.~N(0,~\sigma_{ac(b)}^2)$ 以及 $\varepsilon_{ijkm}\thicksim i.i.d.~N(0,~\sigma_{\varepsilon}^2)$. 此外， $a_i,b_j,(ab)_{ij},c_{k(j)},(ac)_{ik(j)}$ 和 $\varepsilon_{ijkm}$ 独立分布。一般情况下的期望均方方差分析表如表 \@ref(tab:table20-3) 所示。

<table>
<caption>(\#tab:table20-3)示例 [20.2](#sec20-1-3) 方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2020.3.png">
</td>
  </tr>
</tbody>
</table>

通过检查期望均方，可以构造 $F$ 统计量来检验每个方差分量的假设。用于检验以下假设的统计量是

1. 要检验 $H_0\colon{\sigma}_a^2=0\mathrm{~vs~}H_a\colon{\sigma}_a^2>0$，则 $F_a=MSA/MSAB$.
2. 要检验 $H_0\colon{\sigma}_{ab}^2=0\mathrm{~vs~}H_a\colon{\sigma}_{ab}^2>0$，则 $F_{ab}=MSAB/MSAC(B)$.
3. 要检验 $H_0\colon{\sigma}_{c(b)}^2=0\mathrm{~vs~}H_a\colon{\sigma}_{c(b)}^2>0$，则 $F_{c(b)}=MSC(B)/MSAC(B)$.
4. 要检验 $H_0\colon{\sigma}_{ac(b)}^2=0\mathrm{~vs~}H_a\colon{\sigma}_{ac(b)}^2>0$，则 $F_{ac(b)}=MSAC(B)/MSResidual$.

然而，没有 $F$ 统计量来检验 $H_0\colon{\sigma}_b^2=0\mathrm{~vs~}H_b\colon{\sigma}_a^2>0$，因为不涉及 $\sigma^2_b$ 的均方都不具有预期值 $\sigma_\varepsilon^2+n\sigma_{ac(b)}^2+na\sigma_{c(b)}^2+nc\sigma_{ab}^2$，这是当 $\sigma^2_b=0$ 时 $MSB$ 的期望值。但存在一个均方的线性组合（不包括 $MSB$）具有所需的期望值，即，$E[MSC(B)+MSAB-MSAC(B)]=\sigma_\varepsilon^2+n\sigma_{ac(b)}^2+na\sigma_{c(b)}^2+nc\sigma_{ab}^2$. 令 $Q=MSC(B)+MSAB-MSAC(B)$，那么用于检验 $H_0\colon{\sigma}_b^2=0\mathrm{~vs~}H_b\colon{\sigma}_a^2>0$ 的统计量为 $F_b=MSB/Q$. $F_b$ 的抽样分布可以用自由度为 b-1 和 r 的 $F$ 分布来近似。分母自由度 r 的确定是通过使用第 \@ref(chap2) 章中讨论的 Satterthwaite (1946) 近似，将 $rQ/E(Q)$ 的分布近似为卡方分布来完成的。Satterthwaite 近似用于近似 $Q = q_1MS_1 + q_2MS_2 + \cdots + q_kMS_k$ 的抽样分布，其中 $MS_i$ 表示基于自由度为 $f_i$ 的均方，均方独立分布，$q_i$ 是已知常数。那么 $rQ/E(Q)$ 近似分布为基于自由度为 r 的中心卡方随机变量，其中

$$r=\frac{(Q)^2}{\sum_{i=1}^k\frac{\left(q_iMS_i\right)^2}{f_i}}$$

假设 $U$ 是基于 f 个自由度的均方，独立分布于 $MS_1,MS_2,\cdots,MS_k$，期望为 $E(U) = E(Q) + k_0\sigma^2_0$. 检验 $H_0\colon{\sigma}_0^2=0\mathrm{~vs~}H_a\colon{\sigma}_0^2>0$ 的统计量为 $F = U/Q$，其近似分布为具有 f 和 r 自由度的 $F$ 分布。

检验 $H_0\colon{\sigma}_b^2=0\mathrm{~vs~}H_a\colon{\sigma}_b^2>0$ 的统计量为 $F_b = MSB/Q$，近似分布为具有 b-1 和 r 自由度的 $F$ 分布，其中

$$r=\frac{(Q)^2}{\frac{[MSC(B)]^2}{b(c-1)}+\frac{[MSAB]^2}{(a-1)(b-1)}+\frac{[MSAC(B)]^2}{b(a-1)(c-1)}}$$

<table>
<caption>(\#tab:table20-4)关示示例 [20.2](#sec20-1-3) 的数据，使用 I 型平方和的 Proc Mixed 代码的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2020.4.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table20-4) 包含表 \@ref(tab:table20-3) 中数据的方差分析表，其中包括除数（用误差项表示）和用于检验各个假设的 $F$ 统计量。要检验 $H_0\colon{\sigma}_b^2=0\mathrm{~vs~}H_a\colon{\sigma}_b^2>0$，令

$$\begin{aligned}
Q& =MSC(B)+MSAB-MSAC(B)  \\
&=12.0625+10.5625-0.8125=21.8125
\end{aligned}$$

$Q$ 对应的自由度为

$$\begin{aligned}
r& =\frac{(21.8125)^2}{(12.0625)^2/2+(10.5625)^2/1+(0.8125)^2/2}  \\
&=\frac{475.7852}{184.9785}=2.57
\end{aligned}$$

检验统计量 $F = 770.0625/21.8125 = 35.304$ 基于 1 和 2.57 自由度。检验的显著性水平为 0.0144，表明有证据表明 ${\sigma}_b^2 > 0$，或者由于因子 B 水平总体引起的变异是系统总变异的重要组成部分。

**为了检验均衡设计中方差分量的假设，应尽可能使用根据两个均方之比构造的 $F$ 检验。当无法使用两个均方之比时，Satterthwaite 近似是可接受的替代方案**。

**当设计不均衡时，几乎总是需要某种 Satterthwaite 近似来检验有关方差分量的假设**。此外，方差分析表中的平方和可能不具有独立分布，尽管平方和集对于某些特殊情况可能是独立的。残差或误差平方和始终独立于方差分析表中的其他平方和。因此，对于任何期望为 ${\sigma}_{{\varepsilon}}^2+{k}_0{\sigma}_0^2$ 的均方 $U$，统计量 $F_0 = U/MSResidual$ 提供了对假设 $H_0\colon{\sigma}_0^2=0\mathrm{~vs~}H_a\colon{\sigma}_0^2>0$ 的检验。在 $H_0$ 的条件下，$F$ 分布为具有 u 和 v 自由度的中心 $F$ 分布，其中 u 是与 $U$ 相关的自由度，v 是与 $MSResidual$ 相关的自由度。

期望值涉及超过两个方差分量的均方通常不能用于获得具有精确 $F$ 抽样分布的单个方差分量的检验统计量。某些均衡设计会出现精确的 $F$ 分布，如前两个示例所示。比率不精确分布于 $F$ 的一个原因是各自的均方不是独立分布的。**如果设计不太不均衡，那么使用 $F$ 分布作为近似应该是足够的**。此外，当设计不均衡时，平方和（残差除外）不会分布为卡方分布的标量倍数。

一般来说，为了检验 $H_0\colon{\sigma}_0^2=0\mathrm{~vs~}H_a\colon{\sigma}_0^2>0$，将有一个均方，记为 $U_1$，期望为

$$E(U_1)=\sigma_\varepsilon^2+k_{1a}\sigma_a^2+k_{1b}\sigma_b^2+k_{1c}\sigma_c^2$$

但不会有其他均方具有期望 $\sigma_\varepsilon^2+k_{1b}\sigma_b^2+k_{1c}\sigma_c^2$；也就是说，没有一个均方是合适的除数。方法是找到其他均方的线性组合，如 $Q=\sum_{i=1}^kq_iMS_i$ 其中 $E(Q)=\sigma_\varepsilon^2+k_{1b}\sigma_b^2+k_{1c}\sigma_c^2$. Satterthwaite 近似可用于近似 $Q$ 的抽样分布，即，求 r，使得 $rQ/E(Q)$ 近似分布为具有 r 个自由度的卡方随机变量。该近似是双重的，因为 1) 自由度是近似的，2) 组成 $Q$ 的均方不一定按照近似所要求的那样独立分布为卡方随机变量。

SAS^®^-Mixed 代码和使用 III 型平方和得出的方差分析表，用于示例 [19.2](#sec19-1-2) 中的小麦虫害数据，如表 \@ref(tab:table20-5) 所示。variaty 均方的期望值是 $\sigma_\varepsilon^2+3.1795\sigma_\mathrm{var}^2$. 为了检验假设 $H_0\colon\sigma_{\text{var}}^2=0\text{ vs }H_a\colon\sigma_{\text{var}}^2>0$，适当的除数是残差均方，它提供的 $F$ 统计量为 4.79. 将计算出的 $F$ 统计量与具有 3 和 9 个自由度的 $F$ 分布进行比较；其显著性水平为 0.0293. 由于这是一个单向实验，因此 I 型分析与 III 型分析相同。

<table>
<caption>(\#tab:table20-5)示例 [19.2](#sec19-1-2) 中数据含期望均方和 $F$ 统计量的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2020.5.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table20-6)示例 [19.3](#sec19-1-3) 的双向随机效应数据的 I 型分析，其中预期均方和误差项用于计算 $F$ 统计量</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2020.6.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table20-6) 显示了根据示例 [19.3](#sec19-1-3) 的双向随机效应数据的 I 型平方和构造的 SAS-Mixed 代码和方差分析表。为了检验假设 $H_0\colon{\sigma}_{\mathrm{row}\times\mathrm{col}}^2=0\text{ vs }H_a\colon{\sigma}_{\mathrm{row}\times\mathrm{col}}^2>0$，适当的除数是残差均方，它提供的 $F$ 统计量为 14.24. 将计算出的 $F$ 统计量与自由度为 2 和 8 的 $F$ 分布进行比较，提供的显著性水平为 0.0023. 没有精确的检验可用于检验 $H_0\colon{\sigma^2}_{\text{row}}=0\text{ vs }H_a\colon{\sigma^2}_{\text{row}}>0$ 和 $H_0\colon{\sigma^2}_{\text{col}}=0\text{ vs }H_a\colon{\sigma^2}_{\text{col}}>0$，因此需要构造近似检验。用于检验 $H_0\colon{\sigma^2}_{\text{row}}=0\text{ vs }H_a\colon{\sigma^2}_{\text{row}}>0$ 的 $MS_{\text{Row}}$ 的适当除数计算如下

$$\begin{aligned}
Q_{row} =&\,\frac{0.1429}{4.5714}MSCol+\frac1{2.2588}{\left[2.4286-\frac{0.1429}{4.5714}\times2.3126\right]}MSRow\times Col  \\
&+\left[1-\frac{0.1429}{4.5714}-\frac1{2.2588}{\left(2.4286-\frac{0.1429}{4.5714}\times2.3126\right)}\right]MSResidual \\   
=&\,0.0313\times MSCol+1.0432 \,MSRow\times Col-0.0744\,MSResidual\\
=&\,55.7806
\end{aligned}$$

与 $Q_{\mathrm{row}}$ 相关的 Satterthwaite 近似自由度计算如下

$$\begin{aligned}df_{Q_{\mathrm{row}}}&=\frac{(Q_{\mathrm{row}})^2}{\frac{(0.0313\times MSCol)^2}2+\frac{(1.0432\times MSRow\times Col)^2}2+\frac{(0.0744\times MSResidual)^2}8}\\&=1.9961\end{aligned}$$

生成的 $F$ 统计量为 $F_{\mathrm{row}}=MSRow/Q_{\mathrm{row}}=0.0113$，显著性水平为 0.9251. 
用于检验 $H_0\colon{\sigma^2}_{\text{col}}=0\text{ vs }H_a\colon{\sigma^2}_{\text{col}}>0$ 的 $MSCol$ 的适当除数计算为

$$\begin{aligned}
Q_{\text{col}}& =\frac{2.3126}{2.2588}MSRow\times Col+\left[1-\frac{2.3126}{2.2588}\right]MSResidual  \\
&=1.0238\,MSRow\times Col-0.0238\,MSResidual \\
&=56.8729
\end{aligned}$$

与 $Q_{\mathrm{col}}$ 相关的 Satterthwaite 近似自由度计算如下

$$\begin{aligned}df_{Q_{\mathrm{col}}}&=\frac{(Q_{\mathrm{col}})^2}{\frac{(1.0238\times MSRow\times Col)^2}2+\frac{(0.0238\times MSResidual)^2}8}\\&=1.9935\end{aligned}$$

生成的 $F$ 统计量为 $F_{\mathrm{col}}=MSCol/Q_{\mathrm{col}}=0.1323$，显著性水平为 0.8832.

### 似然比检验 {#sec20-1-4}

检验方差分量假设的第二种方法基于似然比程序，该程序涉及评估完全模型的似然函数值以及评估在 $H_0$ 条件下模型的似然函数值。

方程 \@ref(eq:18-3) 的一般随机模型是

$$\boldsymbol y=\boldsymbol j_n\mu+\boldsymbol Z_1 \boldsymbol u_1+\boldsymbol Z_2\boldsymbol u_2+\cdots+\boldsymbol Z_k\boldsymbol u_k+\boldsymbol\varepsilon $$

其中 $\boldsymbol{u}_1\thicksim N(0,{\sigma}_1^2\boldsymbol{I}_{t_1}),\boldsymbol{u}_2\thicksim N(0,{\sigma}_2^2\boldsymbol{I}_{t_2}),...,\boldsymbol{u}_r\thicksim N(0,{\sigma}_r^2\boldsymbol{I}_{t_r}),{\varepsilon}\thicksim N(0,{\sigma}_\varepsilon^2\boldsymbol{I}_N)$ 且这些随机变量独立分布。分布假设意味着 $\boldsymbol y$ 的边际分布为 $N(\boldsymbol j_n\mu,\boldsymbol \Sigma)$ 其中 $\boldsymbol \Sigma=\sigma_\varepsilon^2\boldsymbol I_n+\sigma_1^2\boldsymbol Z_1\boldsymbol Z_1^{\prime}+\sigma_2^2\boldsymbol Z_2\boldsymbol Z_2^{\prime}+\cdots+\sigma_k^2\boldsymbol Z_k\boldsymbol Z_k^{\prime}$. 似然方程是

$$L(\mu,\sigma_\varepsilon^2,\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2|\boldsymbol y)=(2\pi)^{-n/2}|\boldsymbol \Sigma|^{-1/2}\exp\left[-\frac12(\boldsymbol y-\boldsymbol j_n\mu)^{\prime}\boldsymbol\Sigma^{-1}(\boldsymbol y-\boldsymbol j_n\mu)\right]$$

$H_0{:{\sigma_1^2}}=0$ 条件下的似然函数为

$$L_0(\mu,\sigma_\varepsilon^2,0,\sigma_2^2,\ldots,\sigma_k^2|\boldsymbol y)=(2\pi)^{-n/2}\left|\boldsymbol \Sigma_0\right|^{-1/2}\exp\left[-\frac12(\boldsymbol y-\boldsymbol j_n\mu)^{\prime}\boldsymbol \Sigma_0^{-1}(\boldsymbol y-\boldsymbol j_n\mu)\right]$$

其中 $\boldsymbol \Sigma_0=\sigma_\varepsilon^2\boldsymbol I_n+\sigma_2^2\boldsymbol Z_2\boldsymbol Z_2^{\prime}+\sigma_3^2\boldsymbol Z_3\boldsymbol Z_3^{\prime}+\cdots+\sigma_k^2\boldsymbol Z_k\boldsymbol Z_k^{\prime}$.

该过程是获得两个似然函数的参数的最大似然估计，并根据其估计的值评估每个似然函数。似然比检验统计量为

$$LR(\sigma_1^2=0)=\frac{L_0(\hat{\mu}_0,\hat{\sigma}_{\varepsilon0}^2,0,\hat{\sigma}_{20}^2,\hat{\sigma}_{30}^2,\ldots,\hat{\sigma}_{k0}^2|\boldsymbol y)}{L(\hat{\mu},\hat{\sigma}_{\varepsilon}^2,\hat{\sigma}_1^2,\hat{\sigma}_2^2,\ldots,\hat{\sigma}_k^2|\boldsymbol y)}$$

其中 $\hat{\sigma}_{i0}^2$ 表示，在 $H_0{:{\sigma_1^2}}=0$ 条件下从似然函数中得到的 ${\sigma}_{i0}^2$ 的最大似然估计。当 $H_0{:{\sigma_1^2}}=0$ 为真，

$$\begin{aligned}-2\log[LR(\sigma_1^2=0)]=&-2\log_{\mathrm{e}}[L_0(\hat{\mu}_0,\hat{\sigma}_{\varepsilon0}^2,0,\hat{\sigma}_{20}^2,\hat{\sigma}_{30}^2,\ldots,\hat{\sigma}_{k0}^2| \boldsymbol y)]\\&+2\log_{\mathrm{e}}[L(\hat{\mu},\hat{\sigma}_{\varepsilon}^2,\hat{\sigma}_1^2,\hat{\sigma}_2^2,\ldots,\hat{\sigma}_k^2|\boldsymbol y)]\end{aligned}$$

的渐近抽样分布为自由度为 1 的中心卡方分布。存在一个自由度的原因是 $L_0(·)$ 中的参数比 $L_1(·)$ 中少一个。**决策规则是如果 $-2\log[LR(\sigma_1^2=0)]>\chi_{\alpha,1}^2$，则拒绝 $H_0$**. 似然比检验统计量可以使用 SAS-Mixed 计算，其中 METHOD = ML 用作方差分量估计程序。

### 示例 20.3：小麦品种——单向随机效应模型 {#sec20-1-5}

表 \@ref(tab:table20-7) 中给出了 SAS-Mixed 代码以及获得示例 [19.2](#sec19-1-2) 完全模型参数的最大似然估计结果。描述示例 [19.2](#sec19-1-2) 数据的模型的参数最大似然估计为 $\hat{\mu}=3.9909,\hat{\sigma}_\varepsilon^2=0.05749$ 和 $\hat{\sigma}_\mathrm{var}^2=0.04855$ 以及 $-2\log_{\mathrm{e}}(\hat{\mu},\hat{\sigma}_{\varepsilon^{}}^2,\hat{\sigma}_{\mathrm{var}}^2|\boldsymbol y)$ 的值为 4.96762. 缩减模型使用表 \@ref(tab:table20-8) 中的 SAS-Mixed 代码拟合数据，其中从表 \@ref(tab:table20-7) 的模型中删除了语句 "Random Variety;". 在$H_0{:{\sigma_1^2}}=0$ 条件下，参数的最大似然估计为 $\hat{\mu}=4.0269,\hat{\sigma}_\varepsilon^2=0.1014$ 和 $\hat{\sigma}_\mathrm{var}^2=0$ 以及 $-2\log_{\mathrm{e}}(\hat{\mu},\hat{\sigma}_{\varepsilon^{}}^2,\hat{\sigma}_{\mathrm{var}}^2|\boldsymbol y)$ 的值为 7.13832. 检验 $H_0$ 的似然比检验的 -2log~e~ 值为 $7.13832 - 4.96762 = 2.1707$. 值 2.171 与具有一个自由度的中心卡方分布的分位数进行比较。此检验的显著性水平为 0.1407。使用表 \@ref(tab:table20-5) 中的期望均方构造的 $F$ 检验提供了 0.0293 的显著性水平。**似然比检验的抽样分布是渐近分布，对于大样本量来说是可以接受的**，而在本例中，样本量很小。前面给出的 $F$ 检验的抽样分布是准确的。**对于涉及超过两个方差分量的模型，其他方差分量的 $F$ 检验是小样本量的近似检验，通常相当充分，并且对于小样本量的情况，可能比基于渐近分布的检验更好**。

<table>
<caption>(\#tab:table20-7)示例 [19.2](#sec19-1-2) 使用 Method = ML 评估似然函数的 Proc Mixed 代码和拟合完全模型的结果</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2020.7.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table20-8)示例 [19.2](#sec19-1-2) 拟合缩减模型的 Proc Mixed 代码和结果，使用 Method = ML 评估当 $\sigma^2_{\text{var}} = 0$ 时的似然函数</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2020.8.png">
</td>
  </tr>
</tbody>
</table>

### 示例 20.4：不均衡双向 {#sec20-1-6}

为检验示例 [19.3](#sec19-1-3) 中的双向随机效应数据的假设 $H_0$，似然比统计量是通过拟合一个在 random 语句中包含了所有项的模型，然后拟合一个在 random 语句中不包含 row × col 项的模型来获得的，如表 \@ref(tab:table20-9) 所示。用于拟合这两个模型的 SAS-Mixed 代码包含在表 \@ref(tab:table20-9) 中。不包含行效应和列效应的模型拟合结果也包含在内，但没有提供相应的 SAS-Mixed 代码。给出了方差分量的最大似然估计、截距估计和 -2log(likelihood) 值。用于假设 $H_0:\sigma_{\mathrm{~row\times col}}^2=0\mathrm{~vs~}H_a{:}\sigma_{\mathrm{~row\times col}}^2>0$ 的似然比统计量为 $73.4-68.7=4.7$，其在 $H_0$ 条件下为具有单自由度的卡方抽样分布。显著性水平为 0.030，表明有足够的信息去相信 $\sigma_{\mathrm{~row\times col}}^2>0$. 对于检验行和列方差分量是否等于零的假设，似然比统计量的值均为零。这是因为两个方差分量的最大似然估计均为零。因此，在检验行方差分量是否为零时，无论行是否包含在 random 语句中，-2log(likelihood) 的值都保持不变。

<table>
<caption>(\#tab:table20-9)各种模型的方差分量最大似然估计，以便计算似然比检验统计量来检验每个单独方差分量为零的假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2020.9.png">
</td>
  </tr>
</tbody>
</table>

## 构造置信区间 {#sec20-2}

有些程序可以为某些模型中的某些方差分量提供精确的置信区间，但大多数情况下，所获得的置信区间都是近似的，并依赖于某种类型的近似。

### 残差方差 $\sigma^2_\varepsilon$ {#sec20-2-1}

对于一般随机模型，$\sigma^2_\varepsilon$ 的 $(1 - \alpha)100\%$ 置信区间为

$$\frac{v\hat{\sigma}_\varepsilon^2}{\chi_{(\alpha/2),v}^2}\leq\sigma_\varepsilon^2\leq\frac{v\hat{\sigma}_\varepsilon^2}{\chi_{1-(\alpha/2),v}^2}$$

其中 $\hat{\sigma}_\varepsilon^2=\text{SSRESIDUAL}/v$，v 表示与 $\hat{\sigma}_\varepsilon^2$ 相关的自由度。$\chi_{1-(\alpha/2),v}^2,\chi_{\alpha/2,v}^2$ 表示自由度为 v 的卡方分布的下 $\alpha/2$ 和上 $\alpha/2$ 分位点。

### 一般 Satterthwaite 近似 {#sec20-2-2}

方差分量抽样分布的**一般 Satterthwaite 近似** (general Satterthwaite approximation) 是通过将 $r\hat{{\sigma}}_i^2/E(\hat{{\sigma}}_i^2)$ 的前两阶矩与基于 r 个自由度的卡方分布的前两个阶相等，然后解出 r 来获得的。卡方分布的一阶矩等于其自由度 r. 因此，令 $r\hat{{\sigma}}_i^2/E(\hat{{\sigma}}_i^2)$ 的一阶矩等于自由度为 r 的卡方分布的一阶矩，并不能提供有关 r 的信息。自由度为 r 的卡方分布的方差等于 2r. 令 $r\hat{{\sigma}}_i^2/E(\hat{{\sigma}}_i^2)$ 的方差等于自由度为 r 的卡方分布的方差，得到

$$\mathrm{Var}\left(\frac{r\hat{\sigma}_i^2}{E(\hat{\sigma}_i^2)}\right)=2r\quad\mathrm{or}\quad\frac{r^2}{\left[E(\hat{\sigma}_i^2)\right]^2}\mathrm{Var}(\hat{\sigma}_i^2)=2r$$

这意味着

$$r=\frac{2[\operatorname{E}(\hat{\sigma}_i^2)]^2}{\operatorname{Var}(\hat{\sigma}_i^2)}$$

r 的值是通过将 $E(\hat{\sigma}_i^2)$ 替换为 $\hat{\sigma}_i^2$ ，并将 ${\operatorname{Var}(\hat{\sigma}_i^2)}$ 替换为其方差的估计值来估计的。如果获得了 REML 或 ML 解，则可以使用信息矩阵的逆来提供方差分量估计的方差的估计 (estimates of the variances of the estimated variance  components). SAS-Mixed 使用在方差分量估计值处评估的信息矩阵的逆来提供估计方差分量的方差和协方差的渐近估计值。

表 \@ref(tab:table20-10) 包含使用 `covtest`, `cl` 和 `asycov` 选项的 SAS-Mixed 代码，以生成方差分量估计的标准误估计，将每个估计值与其标准误估计的比值计算为 Z 值，将 Satterthwaite 类型的置信区间标记为 Lower 和 Upper，以及方差分量估计值的渐近协方差阵。表 \@ref(tab:table20-10) 中的 data step 计算了与每个方差分量估计相关的近似自由度，然后使用这些自由度来计算 Satterthwaite 置信区间。$\hat{\sigma}_{\text{var}}^2$ 和 $\hat{\sigma}_{\varepsilon}^2$ 的近似自由度分别为 1.758 和 8.828。重新计算的置信区间与 SAS-Mixed 提供的置信区间相同。

<table>
<caption>(\#tab:table20-10)SAS-Mixed 代码，采用 Method = REML 和 Data Step 来计算示例 [19.2](#sec19-1-2) 中小麦昆虫损害的方差分量的 Satterthwaite 类型置信区间</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2020.10.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table20-11) 中的 SAS-Mixed 代码为示例 [20.1](#sec20-1-2) 中的数据提供了方差分量的 REML 估计。方差分量的置信区间是使用 $df = 2(\text{Z-value})^2$ 计算的，分别提供了 0.50, 0.95, 0.84, 1.73, 0 和 8 个自由度。应该注意到，当与随机效应相关的水平数非常少，因此自由度非常小，所得到的置信区间将非常宽。

<table>
<caption>(\#tab:table20-11)用于提供方差分量的 REML 估计和 Satterwaite 置信区间的 SAS-Mixed 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2020.11.png">
</td>
  </tr>
</tbody>
</table>

### 方差分量函数的近似置信区间 {#sec20-2-3}

研究人员通常希望构造关于模型中方差分量的某些函数的置信区间。如下一节所述，涉及均衡模型的某些情况下存在精确的置信区间。

在本节中，使用一阶泰勒级数 (Kendall and Stuart, 1973) 来估计方差函数的方差，然后使用它来构造 Satterthwaite 类型的置信区间。设 $\sigma$ 表示方差分量的向量，$\text{Var}(\sigma)$ 表示方差分量的方差估计和协方差阵估计。假设感兴趣的方差分量的函数为 ${\varphi}({\sigma})$，且 ${\varphi}({\sigma})$ 具有连续的一阶导数。${\varphi}({\sigma})$ 的最大似然估计是 ${\varphi}(\hat{\sigma})$，其中 $\hat{\sigma}$ 表示 $\sigma$ 的最大似然估计。与方差分量对应的信息矩阵的逆的估计是 $\hat{{V}}(\hat{{\sigma}})$，它提供了 $\text{Var}(\sigma)$ 的估计。计算 ${\varphi}({\sigma})$ 关于每个方差分量的导数，并在 $\hat{\sigma}$ 处评估它们。然后令

${\varphi}(\hat{\sigma})$ 的近似方差为 ${\sigma}_{{\varphi}(\hat{\sigma})}^2=\hat{\boldsymbol{f}}^{\prime}\hat{\boldsymbol{V}}(\hat{{\sigma}})\hat{\boldsymbol{f}}$. 与该估计对应的 Satterthwaite 近似自由度计算如下

$$r=\frac{2[\varphi(\hat{\sigma})]^2}{\hat{\sigma}_{\varphi(\hat{\sigma})}^2}$$

假设研究人员有兴趣估计示例 [19.2](#sec19-1-3) 的昆虫损害研究数据的总变异性；也就是说，研究人员想要估计 $\varphi(\sigma)=\sigma^2_{\text{var}}+\sigma^2_{\varepsilon}$. $\varphi(\sigma)$ 的一阶导数向量为

$$\boldsymbol f=\begin{bmatrix}\frac{\partial\varphi(\sigma)}{\partial\sigma_\mathrm{var}^2}\\\frac{\partial\varphi(\sigma)}{\partial\sigma_\varepsilon^2}\end{bmatrix}=\begin{bmatrix}1\\1\end{bmatrix}$$

方差分量估计的协方差阵为

$$\hat{\boldsymbol V}(\hat{\sigma})=\begin{bmatrix}\hat{\sigma}_{\hat{\sigma}_{\mathrm{var}}^2}^2&\hat{\sigma}_{\hat{\sigma}_{\mathrm{var}}^2,\hat{\sigma}_{\hat{\varepsilon}}^2}^2\\\hat{\sigma}_{\hat{\sigma}_{\mathrm{var}}^2,\hat{\sigma}_{\hat{\varepsilon}}^2}^2&\hat{\sigma}_{\hat{\sigma}_{\hat{\varepsilon}}^2}^2\end{bmatrix}$$

因此 $\varphi(\hat{\sigma})$ 的近似方差为 $$\hat{\sigma}_{\varphi(\hat{\sigma})}^2=\hat{\sigma}_{\hat{\sigma}_\mathrm{var}^2}^2+\hat{\sigma}_{\hat{\sigma}_\varepsilon^2}^2+2\hat{\sigma}_{\hat{\sigma}_\mathrm{var}^2\hat{\sigma}_\varepsilon^2}$$. 使用表 \@ref(tab:table20-10) 中的信息可以得到

$$\varphi(\hat{\sigma})=0.07316+0.05700=0.13010$$

以及

$$\hat{\sigma}_{\varphi(\hat{\sigma})}^2=0.006087+0.000736-2(0.00031)=0.006203$$

方差估计对应的近似自由度为

$$r=\frac{2(0.13010)^2}{0.006203}=5.46$$

最后，$\varphi(\sigma)=\sigma_{\mathrm{var}}^2+\sigma_{\varepsilon}^2$ 的 95% 置信区间为

$$0.052287<\sigma_\mathrm{var}^2+\sigma_\varepsilon^2<0.70253$$

作为第二个例子，组内相关系数定义为

$$\rho=\frac{\sigma_\mathrm{var}^2}{\sigma_\varepsilon^2+\sigma_\mathrm{var}^2}$$

可以使用泰勒级数方法构造近似置信区间，其中

$$\varphi(\rho)=\frac{\sigma_\mathrm{var}^2}{\sigma_\varepsilon^2+\sigma_\mathrm{var}^2}$$

$\rho=\frac{\sigma_\mathrm{var}^2}{\sigma_\varepsilon^2+\sigma_\mathrm{var}^2}$ 导数的估计向量的转置为

$$\hat{\boldsymbol f}^{\prime}=\left[\frac{(\hat{\sigma}_\varepsilon^2+\hat{\sigma}_\mathrm{var}^2)-(\hat{\sigma}_\mathrm{var}^2)}{(\hat{\sigma}_\varepsilon^2+\hat{\sigma}_\mathrm{var}^2)^2},\frac{-(\hat{\sigma}_\mathrm{var}^2)}{(\hat{\sigma}_\varepsilon^2+\hat{\sigma}_\mathrm{var}^2)^2}\right]$$

使用表 \@ref(tab:table20-10) 中的信息，组内相关性的估计值为 $\hat\rho = 0.5621$，其方差估计为 0.2156。相应的近似自由度为 $r = 2.93$. 得到的组内相关系数的近似置信区间为 (0.1785, 8.214)，可以简化为 (0.1785, 1)，因为类内相关系数不能大于 1. **对于方差分量的线性函数，一阶泰勒级数方法非常有效，而对于如组内相关性这样的非线性函数，其效果较差**。组内相关性的置信区间更好，将在第 \@ref(sec20-2-3) 节中描述。

当使用平方和的方法获得方差分量的估计值时，可以使用**通常的 Satterthwaite 近似** (usual Satterthwaite approximation) 来构造置信区间。假设每个方差分量的估计值可以表示为方差分析表中给出的均方值的线性组合，如 $\hat{{\sigma}}_s^2=\sum_{i=1}^kc_{is}MS_i+c_{\varepsilon s}MSResidual$. 然后自由度的近似数计算为

$$r_s=\frac{(\hat{\sigma}_s^2)^2}{\sum_{i=1}^k\frac{\left(c_{is}MS_i\right)^2}{df_{_{MS_i}}}+\frac{\left(c_{_{\mathbf{\varepsilon}s}}MSResidual\right)^2}{df_{_{MSResidual}}}}$$

$\sigma^2_s$ 的近似 $(1 - \alpha)100\%$ 置信区间为

$$\frac{r_s\hat{\sigma}_s^2}{\chi_{\alpha/2,r_s}^2}\leq\sigma_s^2\leq\frac{r_s\hat{\sigma}_s^2}{\chi_{1-(\alpha/2),r_s}^2}$$

当使用平方和方法估计方差分量时，可以很容易地从期望均方和 SAS-Mixed 提供的误差项获得均方系数。例如，示例 [20.2](#sec20-1-3) 中 $\sigma^2_a$ 的估计为

$$\hat{\sigma}_a^2=\frac18{\left[MSA-MS(A\times B)\right]}=\frac18(39.0625-10.5625)=3.5625$$

自由度的近似数为

$$\begin{aligned}r&=\frac{(3.5625)^2}{\frac{[(1/8)39.0625]^2}1+\frac{[(1/8)10.5625]^2}1}=0.305\end{aligned}$$

$\sigma^2_a$ 的最终置信区间为 (0.4226, 27556199203.09)，这是由于自由度数非常低 (0.305) 而导致的非常宽的区间。另外，示例 [20.2](#sec20-1-3) 中 $\sigma^2_b$ 的估计为

$$\begin{aligned}
\hat{\sigma}_b^2& =\frac18\{MSB-[MS(A\times B)+MS(C(B))-MS(A\times C(B))]\}  \\
&=\frac18(770.0625-12.0625-10.5625+0.8125) \\
&=93.5313
\end{aligned}$$

自由度的近似数为

$$\begin{aligned}r&=\frac{(93.5313)^2}{\frac{[(1/8)770.0625]^2}1+\frac{[(1/8)12.0625]^2}2\frac{[(1/8)10.5625]^2}{1}\frac{[(1/8)0.8125]^2}2}\\ &=0.944 \end{aligned}$$

$\sigma^2_b$ 的最终置信区间为 (18.1347, 141493.19)，这也是由于自由度数较少而导致的较宽区间。**基本上，当只有少量自由度可用于估计方差分量时，对该方差分量将会知之甚少**。

### 方差分量的 Wald 型置信区间 {#sec20-2-4}

Wald 型置信区间可以使用最大似然估计的渐近正态性来计算。方差分量 $\sigma^2_s$ 的 $(1 - \alpha)100\%$ Wald 置信区间为 $\hat{{\sigma}}_s^2-Z_{\alpha/2}\sqrt{\hat{\sigma}_{\hat{{\sigma}}_s^2}^2}\leq{\sigma}_s^2\leq\hat{{\sigma}}_s^2+Z_{\alpha/2}\sqrt{\hat{\sigma}_{\hat{{\sigma}}_s^2}^2}$. Wald 置信区间关于 $\hat{{\sigma}}_s^2$ 对称。在数据正态性的假设下，与方差相关的抽样分布是不对称的卡方分布。但随着自由度数量的增加，卡方分布的形状变得更加对称。因此，**当与方差分量相关的自由度较大时，Wald 置信区间应该足够，但当自由度数量较小时，Satterthwaite 型置信区间将更接近于现实**。当在 SAS-Mixed 中使用 METHOD = TYPEx 并请求置信区间时，除残差方差分量计算的是 Satterthwaite 置信区间外，所有方差分量都提供 Wald 置信区间。

### 一些精确的置信区间 {#sec20-2-5}

由均衡模型的方差分析表中的所有平方和（不包括 $SSRESIDUAL$）与残差平方和无关，因此可以使用下一个结果来构造方差分量的置信区间，如 $\sigma^2_1$，当模型中存在期望为 $\sigma_\varepsilon^2+a\sigma_1^2$ 的均方时，其置信系数至少为 $1-\alpha$. 令 $Q_1=MSRESIDUAL$ 并假设其基于 u~1~ 个自由度，并令 $Q_2$ 为基于 u~2~ 个自由度，期望为 $\sigma_\varepsilon^2+a\sigma_1^2$. 关于 $\sigma_\varepsilon^2$ 和 $\sigma_\varepsilon^2+a\sigma_1^2$ 的一组精确同时 $(1 - \alpha)100\%$ 置信区间由下式给出

$$\begin{aligned}&\frac{u_1Q_1}{\chi_{\rho/2,u_1}^2}\leq\sigma_\varepsilon^2\leq\frac{u_1Q_1}{\chi_{1-(\rho/2),u_1}^2}\\&\frac{u_2Q_2}{\chi_{\rho/2,u_2}^2}\leq\sigma_\varepsilon^2+a\sigma_1^2\leq\frac{u_2Q_2}{\chi_{1-(\rho/2),u_2}^2}\end{aligned}$$

其中 $\rho=1-\sqrt{1-\alpha}$. 这两个区域在 $(\sigma_{\varepsilon^{}}^2,\sigma_1^2)$ 平面上的交集为 $(\sigma_{\varepsilon^{}}^2,\sigma_1^2)$ 提供了 $(1 - \alpha)100\%$ 同时置信域 (simultaneous confidence region). 置信域图如图 \@ref(fig:figure20-1) 所示。通过确定置信域内 $\sigma_1^2$ 的最大值和最小值，可以获得关于 $\sigma_1^2$ 的 $(1 - \alpha)100\%$ 置信区间。图 \@ref(fig:figure20-1) 中交集区域中 $\sigma_1^2$ 的最小值由 c 给出，最大值由 d 给出。 c 和 d 的值可以通过求解相应的线相交处的 $\sigma_1^2$ 的值来确定。两组线相交于

$$\begin{aligned}c&=\frac{u_2Q_2/\chi_{\rho/2,u_2}^2-u_2Q_2/\chi_{1-(\rho/2),u_2}^2}a\\d&=\frac{u_2Q_2/\chi_{1-(\rho/2),u_2}^2-u_1Q_1/\chi_{\rho/2,u_2}^2}a\end{aligned}$$

从而提供 $c\leq\sigma_1^2\leq d$ 作为关于 $\sigma_1^2$ 的 $(1 - \alpha)100\%$ 置信区间。c 的值对应于如下两条线相交时 $\sigma_1^2$ 的值

$$\sigma_\varepsilon^2=u_1Q_1/\chi_{1-(\rho/2),u_1}^2\quad\mathrm{and}\quad\sigma_\varepsilon^2=-a\sigma_1^2+u_2Q_2/\chi_{\rho/2,u_2}^2$$

d 的值对应于如下两条线相交时 $\sigma_1^2$ 的值

$$\sigma_\varepsilon^2=u_1Q_1/\chi_{\rho/2,u_1}^2\quad\mathrm{~and~}\quad\sigma_\varepsilon^2=-a\sigma_1^2+u_2Q_2/\chi_{1-(\rho/2),u_2}^2$$

<div class="figure" style="text-align: center">
<img src="figure/figure 20.1.png" alt="$\sigma_{\varepsilon}^2,\sigma_1^2$ 的置信域图" width="1065" />
<p class="caption">(\#fig:figure20-1)$\sigma_{\varepsilon}^2,\sigma_1^2$ 的置信域图</p>
</div>

上述关于 $(\sigma_{\varepsilon^{}}^2,\sigma_1^2)$ 的联合置信域可用于构造任意关于 $\sigma_{\varepsilon^{}}^2$ 和 $\sigma_1^2$ 的连续函数的置信区间。令 $\varphi(\sigma_{\varepsilon^{}}^2,\sigma_1^2)$ 表示 $\sigma_{\varepsilon^{}}^2$ 和 $\sigma_1^2$ 在 $(\sigma_{\varepsilon^{}}^2,\sigma_1^2)$ 空间上的连续函数，并令 $\mathfrak{R}_\alpha(\sigma_{\varepsilon^{}}^2,\sigma_1^2)$ 表示 $\sigma_{\varepsilon^{}}^2$ 和 $\sigma_1^2$ 的联合置信域 (joint confidence region). $\varphi(\sigma_{\varepsilon^{}}^2,\sigma_1^2)$ 的置信下限为

$$L=\min_{(\sigma_\varepsilon^2,\sigma_1^2)\in\Re_\alpha(\sigma_\varepsilon^2,\sigma_1^2)}[\varphi(\sigma_\varepsilon^2,\sigma_1^2)]$$

$\varphi(\sigma_{\varepsilon^{}}^2,\sigma_1^2)$ 的置信上限为

$$U=\max_{(\sigma_\varepsilon^2,\sigma_1^2)\in\Re_\alpha(\sigma_\varepsilon^2,\sigma_1^2)}[\varphi(\sigma_\varepsilon^2,\sigma_1^2)]$$

例如，如果

$$\varphi(\sigma_\varepsilon^2,\sigma_1^2)=\frac{\sigma_1^2}{\sigma_\varepsilon^2+\sigma_1^2}$$

最大值出现在 $(u_1Q_1/\chi_{1-(\alpha/2),{u_1}}^2C)$，最小值出现在 $(u_1Q_1/\chi_{\alpha/2,u_1}^2,d)$. 因此，使用示例 [20.1](#sec20-1-2) 中的数据，得到的 $\frac{\sigma_1^2}{\sigma_\varepsilon^2+\sigma_1^2}$ 的置信区间约为

$$-0.048511<\frac{\sigma_1^2}{\sigma_\varepsilon^2+\sigma_1^2}<1.89748$$

或

$$0<\frac{\sigma_1^2}{\sigma_\varepsilon^2+\sigma_1^2}<1$$

因为众所周知 $\frac{\sigma_1^2}{\sigma_\varepsilon^2+\sigma_1^2}$ 必须大于 0 且小于 1. 当自由度数很小时，最小化和最大化方法可能会导致疯狂的结果，特别是当 $\varphi(\sigma_\varepsilon^2,\sigma_1^2)$ 是非线性函数时。

Williams (1962) 提出的 $\sigma_1^2$ 的另一个置信区间具有由下式给出的置信下限和上限：

$$\begin{aligned}c_1&=\frac{u_2[Q_2-Q_1F_{\alpha/2,u_2,u_1}]}{a\chi_{\alpha/2,u_2}^2}\quad\text{and}\quad d_1=\frac{u_2[Q_2-Q_1F_{1-(\alpha/2),u_2,u_1}]}{a\chi_{1-(\alpha/2),u_2}^2}\end{aligned}$$

其置信系数至少为 $1 - 2\alpha$. 模拟研究（Boardman, 1974；以及作者未发表的研究）表明，如果取 $\rho=\alpha$ 为而不是 $\rho=1-√(1-\alpha)$，如上所述，第一种程序为 $\sigma_1^2s$ 提供的置信区间的置信系数不小于 $1-\alpha$. 类似的模拟表明，Williams 提出的区间的置信系数不小于 $1-\alpha$. Williams 的区间比 $c\leq\sigma_1^2\leq d$ 刻画的区间稍短。因此 (c, d) 区间比 Williams 区间稍微保守一些。

### 示例 20.5：均衡单向随机效应处理结构 {#sec20-2-6}

<table>
<caption>(\#tab:table20-12)示例 [20.5](#sec20-2-6) 数据</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2020.12.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table20-12) 中的数据来自单向处理结构，其中从工厂工作的工人总体中随机选择了五名工人。工人的水平是随机效应，研究是在完全随机设计结构下进行的，其中五名工人，每名工人三个观测值，响应是在固定时间段内装配的单元数。可用于刻画这些数据的模型是

$$y_{ij}=\mu+iv_i+\varepsilon_{ij},\quad i=1,2,\ldots,5,j=1,2,3,w_i\thicksim i.i.d.N(0,\sigma_{w}^2)\mathrm{~and~}\varepsilon_{ij}\thicksim i.i.d.N(0,\sigma_{\varepsilon}^2)$$
  
<table>
<caption>(\#tab:table20-13)示例 [20.5](#sec20-2-6) 工人数据的 SAS-Mixed 代码和方差分析表结果</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2020.13.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table20-13) 提供了包括期望均方的方差分析表。$\sigma_\varepsilon^2$ 的 95% 置信区间为

$$\frac{16}{20.483}\leq\sigma_\varepsilon^2\leq\frac{16}{3.247}\quad\mathrm{or}\quad0.781\leq\sigma_\varepsilon^2\leq4.928$$

$\sigma_\varepsilon^2+3\sigma_w^2$ 的 95% 置信区间为

$$\frac{128}{11.143}\leq\sigma_\varepsilon^2+3\sigma_w^2\leq\frac{128}{0.484}\quad\mathrm{or}\quad11.487\leq\sigma_\varepsilon^2+3\sigma_w^2\leq264.463$$

通过使用这两个区间相交产生的联合置信域，关于 $\sigma_w^2$ 的近似 95% 置信区间为 $c\leq\sigma_w^2\leq d$，其中

$$\begin{aligned}c&=\frac{11.487-4.928}3=2.186\quad\mathrm{and}\quad d=\frac{264.463-0.781}3=87.818\end{aligned}$$

因此，可以得到 $2.186\leq\sigma_w^2\leq 87.818$ 作为 $\sigma_w^2$ 的近似 95% 置信区间。

这种程序可用于构造关于任何方差分量（例如 $\sigma_1^2$）的置信区间，当存在两个使得 $E(Q_2) = E(Q_1) + \sigma_1^2$ 的独立均方 $Q_1$ 和 $Q_2$，分别基于 u~1~ 和 u~2~ 自由度。令 ${\sigma}_0^2=E(Q_1)$ ，并在之前的开发中用 $\sigma^2_0$ 替换 $\sigma^2_\varepsilon$；那么 $\sigma_1^2$ 的 $(1 - \alpha)100\%$ 近似置信区间为 $c\leq\sigma_1^2\leq d$. Williams (1962) 的结果也适用于这种情况 (Graybill, 1976, 定理 15.3.5).

### 示例 20.6 {#sec20-2-7}

为了演示上述更一般的程序，对示例 [20.2](#sec20-1-3) 中的数据构造了关于 $\sigma_a^2$ 的 90% 置信区间，其中数据在表 \@ref(tab:table20-2) 中，期望均方在表 \@ref(tab:table20-3) 中，并且方差分析表数据见表 \@ref(tab:table20-4). 令 $Q_2=MSA$ 且 $Q_1=MSAB$，两者都基于 1 个自由度。它们的期望分别是

$$E(Q_1)=\sigma_\varepsilon^2+2\sigma_{ac(b)}^2+4\sigma_{ab}^2\quad\mathrm{~and~}\quad E(Q_2)=\sigma_\varepsilon^2+2\sigma_{ac(b)}^2+4\sigma_{ab}^2+8\sigma_a^2,$$

令 $\sigma_0^2=\sigma_\varepsilon^2+2\sigma_{ac(b)}^2+4\sigma_{ab}^2,\sigma_1^2=\sigma_a^2$ 以及 $a=8$. 那么 $\sigma_0^2$ 的置信限为

$$\begin{aligned}c&=\frac{(1\times39.0625/3.8415)-(1\times10.5625/0.00393)}8=-334.69\\d&=\frac{(1\times39.0625/0.00393)-(1\times10.5625/3.8415)}8=1242.10\end{aligned}$$

或

$$-334.69\leq\sigma_a^2\leq1242.10$$

在此示例中，图 \@ref(fig:figure20-1) 中确定 c 值的两条线在参数空间之外相交，因此下限被截断为零，因此置信区间仅包括参数空间中 $\sigma_a^2$ 的值。因此，$\sigma_a^2$ 的置信区间为 $0 \le \sigma_a^2 \le 1242.10$. 类似的技术可用于构造关于 $\sigma_{ab}^2,\sigma_{ac(b)}^2$ 和 $\sigma_{c(b)}^2$ 的置信区间，但不能用于获取关于 $\sigma_b^2$ 的置信区间。

Burdick and Graybill (1992) 提出了几种可以获得方差分量精确置信区间的情况。他们还提出了方差分量几个函数的近似置信区间。对于每个单元格中观测值数量相等且没有缺失单元格的均衡设计，通常可以获得精确的置信区间。这些区间基于分析中的一组平方和，这些平方和作为独立的卡方随机变量进行分布。Burdick and Graybill (1992) 介绍了 Graybill and Wang (1980), Lu et al. (1989) 以及 Ting et al. (1990) 为方差分量的各种函数开发的置信区间。以下讨论是关于与在完全随机设计结构中进行的均衡单向随机效应处理结构相关的置信区间，Burdick and Graybill 称之为单层嵌套设计 (one-fold nested design). 该模型是

$$\begin{aligned}&y_{ij}=\mu+u_i+\varepsilon_{ij},\quad i=1,2,\ldots,t,\quad j=1,2,\ldots,n\\&u_i\thicksim i.i.d.N(0,\sigma_u^2),\varepsilon_{ij}\thicksim i.i.d.N(0,\sigma_\varepsilon^2)\end{aligned}$$

其中 $u_i,\varepsilon_{ij}$ 为独立随机变量。

该模型的方差分析表基于两个平方和

$$\begin{aligned}\mathrm{SSBetween}&=n\sum_{i=1}^t{(\bar{y}_{i\cdot}-\bar{y}_{\cdot\cdot})^2}=(t-1)Q_1\\\text{SSWithin}&=\sum_{i=1}^t\sum_{j=1}^n{(y_{ij}-\bar{y}_{i\cdot})^2}=t(n-1)Q_2\end{aligned}$$

其中 $Q_1$ 和 $Q_2$ 为各自的均方。与平方和对应的期望均方为

$$\begin{aligned}E(\text{MSBetween})&=E(Q_1)=\sigma_\varepsilon^2+n\sigma_u^2\quad\mathrm{~and~}\quad E(\text{MSWITHIN})=E(Q_2)=\sigma_\varepsilon^2\end{aligned}$$

两个平方和是具有如下抽样分布的独立随机变量

$$\frac{(t-1)Q_1}{\sigma_\varepsilon^2+n\sigma_u^2}\thicksim\chi_{t-1}^2\quad\mathrm{and}\quad\frac{t(n-1)Q_2}{\sigma_\varepsilon^2}\thicksim\chi_{t(n-1)}^2$$

$\sigma^2_\varepsilon$ 的 95% 精确置信区间为

$$\left[\frac{Q_2}{F_{\alpha/2,t(n-1),\infty}}\leq\sigma_\varepsilon^2\leq\frac{Q_2}{F_{1-(\alpha/2),t(n-1),\infty}}\right]$$

其中 $F_{\alpha/2,t(n-1),\infty}$ 和 $F_{1-(\alpha/2),t(n-1),\infty}$ 为分子分母自由度分别为 $t(n-1)$ 和 $\infty$ 的 $F$ 分布的上分位数和下分位数。可以注意到 $F_{\alpha/2,t(n-1),\infty}=\chi_{\alpha/2,t(n-1)}^2/[t(n-1)]$ 以及 $F_{1-(\alpha/2),t(n-1),\infty}=\chi_{1-(\alpha/2),t(n-1)}^2/[t(n-1)]$. 使用这些分位数表示法，上述关于 $\sigma^2_\varepsilon$ 的置信区间与第 \@ref(sec20-2-1) 节中的区间相同。

由 Burdick and Graybill (1992) 给出的关于 $\sigma^2_u$ 的 $(1-\alpha)100\%$ 近似置信区间为

$$\left[\frac{Q_1-Q_2-\sqrt{V_L}}n\leq\sigma_u^2\leq\frac{Q_1-Q_2+\sqrt{V_U}}n\right]$$

其中 

$$V_L=G_1^2Q_1^2+H_2^2Q_2^2+G_{12}Q_1Q_2\quad V_U=H_1^2Q_1^2+G_2^2Q_2^2+H_{12}Q_1Q_2$$

其中

$$\begin{aligned}
G_{1} &=1-\frac1{F_{\alpha/2,t-1,\infty}},\quad G_2=1-\frac1{F_{\alpha/2,t(n-1),\infty}} \\
H_{1} &=\frac1{F_{1-(\alpha/2),t-1,\infty}}-1,~H_2=\frac1{F_{1-(\alpha/2),t(n-1),\infty}}-1 \\
G_{12} &=\frac{(F_{\alpha/2,t-1,t(n-1)}-1)^2-G_1^2F_{\alpha/2,t-1,t(n-1)}^2-H_2^2}{F_{\alpha/2,t-1,t(n-1)}} \\
H_{12} &=\frac{(1-F_{1-(\alpha/2),t-1,t(n-1)})^2-H_1^2F_{1-(\alpha/2),t-1,t(n-1)}^2-G_2^2}{F_{1-(\alpha/2),t-1,t(n-1)}} 
\end{aligned}$$

当下界为负时，下限设置为 0. 当 $Q_1/Q_2<F_{\alpha/2,t-1,t(n-1)}$ 时下界为负。 使用 $Q_1/Q_2>F_{\alpha,t-1,t(n-1)}$ 作为决策规则来拒绝 $H_0{:}\sigma_u^2=0\mathrm{~vs~}{H_a}{:}\sigma_u^2>0$，这提供了一个精确检验，是一个**一致最有效无偏检验** (uniformly most powerful unbiased test) (Lehmann, 1986).

由 Burdick and Graybill (1992) 给出的关于 $\sigma_\varepsilon^2+\sigma_u^2$ 的近似 $(1 - \alpha)100\%$ 置信区间为

$$\left[\hat{\sigma}_\varepsilon^2+\hat{\sigma}_u^2-\frac{\sqrt{G_1^2Q_1^2+G_2^2(n-1)Q_2^2}}n\leq\sigma_\varepsilon^2+\sigma_u^2\leq\hat{\sigma}_\varepsilon^2+\hat{\sigma}_u^2+\frac{\sqrt{H_1^2Q_1^2+H_2^2(n-1)Q_2^2}}n\right]$$

其中

$$\begin{aligned}G_1&=1-\frac1{F_{\alpha/2,t-1,\infty}},\quad G_2=1-\frac1{F_{\alpha/2,t(n-1),\infty}}\\H_1&=\frac1{F_{1-(\alpha/2),t-1,\infty}}-1,\quad H_2=\frac1{F_{1-(\alpha/2),t(n-1),\infty}}-1\end{aligned}$$

组内相关性 $\rho=\frac{\sigma_u^2}{\sigma_e^2+\sigma_u^2}$ 的近似 $(1 - \alpha)100\%$ 置信区间为

$$\begin{aligned}\left[\frac{L^*-1}{L^*-1+n}\le\rho\le\frac{U^*-1}{U^*-1+n}\right]\end{aligned}$$

其中

$$L^*=\frac{Q_1}{Q_2F_{\alpha/2,t-1,t(n-1)}}\quad\mathrm{and}\quad U^*=\frac{Q_1}{Q_2F_{1-(\alpha/2),t-1,t(n-1)}}$$

比率 $\sigma_{u}^2/\sigma_{\varepsilon}^2$ 的 $(1 - \alpha)100\%$ 置信区间由下式给出

$$\left[\frac{L^*-1}n\leq\frac{\sigma_u^2}{\sigma_\varepsilon^2}\leq\frac{U^*-1}n\right]$$

Burdick and Graybill (1992) 还描述了当模型中存在两个或多个方差分量时，为均衡和不均衡设计构造关于方差分量和方差分量函数的置信区间的方法。本节的讨论仅限于具有两个方差分量的模型。以下示例用于演示这些置信区间的计算。

### 示例 20.6 （续） {#sec20-2-8}

示例 [20.5](#sec20-2-6) 的数据用于演示 $\sigma_\varepsilon^2,\sigma_w^2,\sigma_w^2/\sigma_\varepsilon^2$ 和 $\rho$ 置信区间的计算。两个平方和为 $Q_1=32.0,Q_2=1.6$，且 $n=3, t-1=4,t(n-1)=10,F_{0.025,3,10}=4.46834,F_{0.975,3,10}=0.11307,F_{0.025,3,\infty}=2.78582,F_{0.025,10,\infty}=2.04832, -0.11208,H_{12}=-8.15882,V_\text{L}{ = 4 2 6 . 1 2 9 , V _\text{U}{ = 5 3 5 1 5 7 1 , L }^{*}=}4.47593$ 以及 $U^{*}=176.878$. $\sigma_\varepsilon^2,\sigma_w^2,\sigma_w^2/\sigma_\varepsilon^2$ 和 $\rho$ 的 95% 置信区间分别为 $[0.78113\leq\sigma_{\varepsilon}^2\leq4.92767],[3.25237\leq\sigma_{w^{}}^2\leq87.2449],[4.78382\leq\sigma_{\varepsilon}^2+\sigma_{u}^2\leq89.1765],[1.15864\leq(\sigma_{w^{}}^2/\sigma_{\varepsilon}^2)\leq58.6259]$ 和 $[0.53675\leq\rho\leq0.98323]$.

REML 解的 Satterthwaite 类型置信区间显示在表 \@ref(tab:table20-14) 中。使用三种方法构造关于 $\sigma^2_w$ 的置信区间。三个 95% 置信区间是 Burdick and Graybill $(3.25237\leq\sigma_w^2\leq87.2449)$, REML $(3.4964\leq\sigma_w^2\leq99.2870)$ 以及联合置信域 $(2.186\leq\sigma_w^2\leq 87.818)$. 三个区间之间没有太大差异，并且随着方差分量对应的随机效应水平的增加，差异会减小。

<table>
<caption>(\#tab:table20-14)示例 [20.5](#sec20-2-6) 工人数据的 SAS-Mixed 代码和方差分量的 REML 估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2020.14.png">
</td>
  </tr>
</tbody>
</table>

## 模拟研究 {#sec20-3}

由于计算置信区间的选项众多，可以使用给定数据的结构进行模拟研究，以提供有关置信区间程序性能的信息。例如，表 \@ref(tab:table20-15) 中的信息包含 SAS data step 代码，用于对示例 [20.5](#sec20-2-6) 中数据集的结构进行模拟，其中每个数据集有五名工人，每个工人有三个观测值。模拟中使用的方差分量值对应于表 \@ref(tab:table20-15) 中示例 [20.5](#sec20-2-6) 的 REML 估计值（请注意，在这种情况下，矩估计与 REML 估计值相同）。模拟结果表明，96.63% 的置信区间覆盖 ${\sigma_w^2}$ 的真实值，95.26% 的置信区间覆盖 ${\sigma_\varepsilon^2}$ 的真实值。这两个经验置信率都非常接近人们期望的 95% 覆盖率。当数据服从正态分布时，对于这种尺寸的设计，每个 Satterthwaite 置信区间的覆盖率都相当充足。

<table>
<caption>(\#tab:table20-15)用于单向随机效应 Satterthwaite 置信区间覆盖范围检查的模拟代码和结果</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2020.15.png">
</td>
  </tr>
</tbody>
</table>

对于不均衡模型，方差分析表中的平方和并不一定具有独立的卡方分布。因此，用于均衡模型的技术不能直接应用于不均衡模型，否则会违反假设。如果模型不太不均衡，平方和将几乎是独立的卡方随机变量，均衡模型技术应该提供有效的置信区间。对于不均衡模型和大样本量，可以使用最大似然估计及其渐近性质来构造关于方差分量的置信区间。上述示例说明了对于任何给定问题，不同的置信区间可能会有何不同。当设计不均衡且样本量不大时，一种建议是使用最宽的置信区间。评估特定方法有效性的另一种技术是使用模拟。

## 结束语 {#sec20-4}

本章介绍关于方差分量和方差分量函数的假设检验和构造置信区间的方法。假设检验是使用平方和方法完成的，其中使用期望均方或使用似然比检验构造检验统计量。基于平方和方法的检验似乎比基于似然比的检验对于小样本量实验具有更好的覆盖范围。SAS-Mixed 提供的置信区间以及其他特殊情况下的置信区间也得到了描述。如第 \@ref(sec20-3) 节所述，**可以通过模拟研究来评估特定置信区间对于给定数据结构的实用性**。针对不均衡模型和一些均衡模型的情况，描述了围绕 Satterthwaite 近似值的近似方法。还包括了几个示例来演示这些技术。SAS-Mixed 提供的置信区间对于各个方差分量似乎是足够的，但对于方差分量的非线性函数（如组内相关系数或两个方差分量的比率）的置信区间，则需要专门的方法。

## 练习 {#sec20-5}
