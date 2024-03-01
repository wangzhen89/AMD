
# 方差分量的推断方法

当研究者设计一个涉及随机效应因素的实验时，她通常希望根据模型中指定的特定方差分量做出推断。特别是，如果 $\sigma_u^2$ 是对应于因素 A 水平分布的方差分量，实验者可能希望确定是否有足够的证据得出 $\sigma_u^2>0$ 的结论。可以通过以下方式做出适当的决定：1）检验假设 $H_0{:}\sigma_u^2=0\mathrm{~vs~}{H_a}{:}\sigma_u^2>0$；2）构建关于 $\sigma_u^2$ 的置信区间；或 3）为 $\sigma_u^2$ 构建置信下限。本章讨论了随机效应模型的这些推断程序，其中假设检验方法在第 \@ref(sec20-1) 节中描述，置信区间（下界）的构造在第 \@ref(sec20-2) 节中描述。方差分量的置信区间构建一直是研究的一个活跃领域，许多作者已经为方差分量参数的特定函数开发了专门的置信区间。本章中描述的方法可在当前软件中使用，有些可用于特定问题。讨论并不全面，而是指出了已经解决的置信区间类型。在 Burdick and Graybill (1992) 以及当前统计期刊的论文中，有更为完整的讨论。

## 假设检验 {#sec20-1}

有两种基本技术来检验关于方差分量的假设。第一种技术使用方差分析表中的平方和来构建 $F$ 统计量。对于大多数均衡模型，$F$ 统计量的分布与 $F$ 分布完全相同，而对于不均衡模型，分布由 $F$ 分布近似，随着设计变得越来越不均衡，近似值变得越来越差。第二种技术基于似然比检验 (likelihood rato test)，渐近分布于卡方分布。对于均衡设计，$F$ 统计量方法可能比似然比检验更好，而对于不均衡设计，没有明确的选择。读者可能希望在决定使用哪种方法来检验感兴趣的假设之前，使用与感兴趣的数据集类似的数据结构进行模拟实验，以研究检验统计量的分布。

### 使用方差分析表 {#sec20-1-1}

[^multiplesofchi-square]: 原文："then sums of squares obtained by the usual analysis of variance are independently distributed as scalar multiples of chi-square random variables."

如果数据集是均衡的，那么通过常规方差分析获得的平方和独立分布于卡方随机变量的标量倍数[^multiplesofchi-square]。令 $Q$ 表示基于自由度为 v 的平方和，其中其预期均方是四个方差分量的函数。也就是说，假设

$$E(Q/\nu)=\sigma_\varepsilon^2+k_1\sigma_1^2+k_2\sigma_2^2+k_3\sigma_3^2$$

然后，假设数据服从正态分布，

$$W=\frac Q{\sigma_\varepsilon^2+k_1\sigma_1^2+k_2\sigma_2^2+k_3\sigma_3^2}$$

通常分布为具有 n 个自由度的卡方随机变量。对于许多形如 $H_0{:}\sigma_1^2=0\mathrm{~vs~}{H_a}{:}\sigma_1^2>0$ 的假设，有两个独立的平方和，用 $Q_1$ 和 $Q_2$ 表示，分别基于 v~1~ 和 v~2~ 自由度，期望为

$$\begin{aligned}E(Q_1/n_1)&=\sigma_\varepsilon^2+k_1\sigma_1^2+k_2\sigma_2^2+k_3\sigma_3^2\\E(Q_2/n_2)&=\sigma_\varepsilon^2+k_2\sigma_2^2+k_3\sigma_3^2\end{aligned}$$

假设 $H_0{:}\sigma_1^2=0\mathrm{~vs~}{H_a}{:}\sigma_1^2>0$ 等价于

$$\begin{aligned}H_0\colon E\left(Q_1/v_1\right)&=E(Q_2/v_2)\text{ vs }H_a\colon E(Q_1/v_1)>E(Q_2/v_2)\end{aligned}$$

用于检验该假设的统计量为 $F = (Q_1/v_1)/(Q_2/v_2)$，在 $H_0$ 条件下，该统计量通常分布为具有 v~1~ 和 v~2~ 自由度的中心 $F$ 分布。对于较大的 $F$ 值，该假设被拒绝。此过程涉及获取平方和，然后使用其预期均方来确定每个感兴趣假设的适当除数。以下两个示例演示了此程序。

### 示例 20.1：完全随机设计结构中的双向随机效应检验统计量 {#sec20-1-2}

在完全随机设计结构中，两个因素均为随机的双向处理结构的模型为

$$y_{ijk}=\mu+a_i+b_j+c_{ij}+\varepsilon_{ijk}\quad\mathrm{~for~}i=1,2,\ldots,a,j=1,2,\ldots,b,\mathrm{~and~}k=1,2,\ldots,n$$

其中 $a_i\thicksim i.i.d.\,N(0,\sigma_a^2),b_j\thicksim i.i.d.\,N(0,\sigma_b^2),c_{ij}\thicksim i.i.d.\,N(0,\sigma_c^2),\varepsilon_{ijk}\thicksim i.i.d.\,N(0,\sigma_\varepsilon^2)$，且随机变量 $a_i,b_j,c_{ij},\varepsilon_{ijk}$ 独立分布。

表 \@ref(tab:table20-1) 显示了模型的方差分析表以及模型的平方和和预期均方。通过检查预期均方来选择适当的分子和分母来构建检验统计量。用于检验假设 $H_0{:}\sigma_a^2=0\mathrm{~vs~}{H_u}{:}\sigma_a^2>0$ 的统计量是通过在 $MSA$ 的预期均方中设置 $\sigma_a^2=0$ 来构建的。接下来，当 $H_0$ 为真时，找到与 $MSA$ 具有相同预期均方的另一个均方，并使用该均方作为除数。要检验 $H_0{:}\sigma_a^2=0\mathrm{~vs~}{H_u}{:}\sigma_a^2>0$，适当的除数是 $MSAB$；要检验 $H_0{:}\sigma_b^2=0\mathrm{~vs~}{H_u}{:}\sigma_b^2>0$，适当的除数是 $MSAB$；要检验 $H_0{:}\sigma_c^2=0\mathrm{~vs~}{H_u}{:}\sigma_c^2>0$，适当的除数是 $MSResidual$. 决策规则是拒绝 $H_0{:}\sigma_a^2=0\mathrm{~vs~}{H_u}{:}\sigma_a^2>0$ 如果 $F = MSA/MSAB>F_{\alpha,(a-1),(a-1)(b-1)}$，其中 $\alpha$ 为所选的 I 类错误率。对于 $\sigma_b^2$ 和 $\sigma_c^2$ 可以类似地确定检验统计量。表 \@ref(tab:table20-1) 包含假设列表和相应的检验统计量。很可能，当 $F$ 统计量不超过指定的分位数时，结论并不是方差分量为零，而是与系统中其他变异源相比，方差分量的大小可以忽略不计。

### 示例 20.2：复杂三向随机效应检验统计量 {#sec20-1-3}

表 \@ref(tab:table20-2) 中的数据来自一个设计，其中处理结构中的三个因素的水平是随机效应，A 的水平与 B 的水平交叉，C 的水平嵌套在 B 的水平内，所有这些都在一个完全随机设计结构中。一个描述与表 \@ref(tab:table20-2) 结构相似的大型数据集的一般模型是

$$\begin{aligned}&y_{ijkm}=\mu+a_i+b_j+(ab)_{ij}+c_{k(j)}+(ac)_{ik(j)}+\varepsilon_{ijkm}\\&\mathrm{for~}i=1,2,\ldots,a,j=1,2,\ldots,b,k=1,2,\ldots,c,\mathrm{~and~}m=1,2,\ldots,n\end{aligned}$$

参数 $\mu$ 表示总体均值，$a_i$ 表示因子 A 水平 i 的效应，$b_j$ 表示因子 B 水平 j 的效应，$(ab)_{ij}$ 表示因子 A 和因子 B 水平之间的交互作用，$c_{k(j)}$ 表示嵌套在因子 B 第 j 个水平内的因子 C 水平 k 的效应，$(ac)_{ik(j)}$ 表示嵌套在 B 水平内的因子 A 和因子 C 水平之间的交互作用，以及 $\varepsilon_{ijkm}$ 表示实验单元或抽样误差。在理想条件下，$a_i\thicksim i.i.d.~N(0,\sigma_a^2),~b_j\thicksim i.i.d.~N(0,~\sigma_b^2),~(ab)_{ij}\thicksim i.i.d.~N(0,~\sigma_{ab}^2),~(ac)_{ik(j)}\thicksim i.i.d.~N(0,~\sigma_{ac(b)}^2)$ 以及 $\varepsilon_{ijkm}\thicksim i.i.d.~N(0,~\sigma_{\varepsilon}^2)$. 此外， $a_i,b_j,(ab)_{ij},c_{k(j)},(ac)_{ik(j)}$ 和 $\varepsilon_{ijkm}$ 独立分布。一般情况下的期望均方方差分析表如表 \@ref(tab:table20-3) 所示。通过检查期望均方，可以构建 $F$ 统计量来检验每个方差分量的假设。用于检验以下假设的统计量是

1. 要检验 $H_0\colon{\sigma}_a^2=0\mathrm{~vs~}H_a\colon{\sigma}_a^2>0$，则 $F_a=MSA/MSAB$.
2. 要检验 $H_0\colon{\sigma}_{ab}^2=0\mathrm{~vs~}H_a\colon{\sigma}_{ab}^2>0$，则 $F_{ab}=MSAB/MSAC(B)$.
3. 要检验 $H_0\colon{\sigma}_{c(b)}^2=0\mathrm{~vs~}H_a\colon{\sigma}_{c(b)}^2>0$，则 $F_{c(b)}=MSC(B)/MSAC(B)$.
4. 要检验 $H_0\colon{\sigma}_{ac(b)}^2=0\mathrm{~vs~}H_a\colon{\sigma}_{ac(b)}^2>0$，则 $F_{ac(b)}=MSAC(B)/MSResidual$.

然而，没有 $F$ 统计量来检验 $H_0\colon{\sigma}_b^2=0\mathrm{~vs~}H_b\colon{\sigma}_a^2>0$，因为不涉及 $\sigma^2_b$ 的均方都不具有预期值 $\sigma_\varepsilon^2+n\sigma_{ac(b)}^2+na\sigma_{c(b)}^2+nc\sigma_{ab}^2$，这是当 $\sigma^2_b=0$ 时 $MSB$ 的期望值。但存在一个均方的线性组合（不包括 $MSB$）具有所需的期望值，即，$E[MSC(B)+MSAB-MSAC(B)]=\sigma_\varepsilon^2+n\sigma_{ac(b)}^2+na\sigma_{c(b)}^2+nc\sigma_{ab}^2$. 令 $Q=MSC(B)+MSAB-MSAC(B)$，那么用于检验 $H_0\colon{\sigma}_b^2=0\mathrm{~vs~}H_b\colon{\sigma}_a^2>0$ 的统计量为 $F_b=MSB/Q$. $F_b$ 的抽样分布可以用自由度为 b-1 和 r 的 $F$ 分布来近似。分母自由度 r 的确定是通过使用第 \@ref(chap2) 章中讨论的 Satterthwaite (1946) 近似，将 $rQ/E(Q)$ 的分布近似为卡方分布来完成的。Satterthwaite 近似用于近似 $Q = q_1MS_1 + q_2MS_2 + \cdots + q_kMS_k$ 的抽样分布，其中 $MS_i$ 表示基于自由度为 $f_i$ 的均方，均方独立分布，$q_i$ 是已知常数。那么 $rQ/E(Q)$ 近似分布为基于自由度为 r 的中心卡方随机变量，其中

$$r=\frac{(Q)^2}{\sum_{i=1}^k\frac{\left(q_iMS_i\right)^2}{f_i}}$$

假设 $U$ 是基于 f 个自由度的均方，独立分布于 $MS_1,MS_2,\cdots,MS_k$，期望为 $E(U) = E(Q) + k_0\sigma^2_0$. 检验 $H_0\colon{\sigma}_0^2=0\mathrm{~vs~}H_a\colon{\sigma}_0^2>0$ 的统计量为 $F = U/Q$，其近似分布为具有 f 和 r 自由度的 $F$ 分布。

检验 $H_0\colon{\sigma}_b^2=0\mathrm{~vs~}H_a\colon{\sigma}_b^2>0$ 的统计量为 $F_b = MSB/Q$，近似分布为具有 b-1 和 r 自由度的 $F$ 分布，其中

$$r=\frac{(Q)^2}{\frac{[MSC(B)]^2}{b(c-1)}+\frac{[MSAB]^2}{(a-1)(b-1)}+\frac{[MSAC(B)]^2}{b(a-1)(c-1)}}$$

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

为了检验均衡设计中方差分量的假设，应尽可能使用根据两个均方之比构建的 $F$ 检验。当无法使用两个均方之比时，Satterthwaite 近似是可接受的替代方案。

当设计不均衡时，几乎总是需要某种 Satterthwaite 近似来检验有关方差分量的假设。此外，方差分析表中的平方和可能不具有独立分布，尽管平方和集对于某些特殊情况可能是独立的。残差或误差平方和始终独立于方差分析表中的其他平方和。因此，对于任何期望为 ${\sigma}_{{\varepsilon}}^2+{k}_0{\sigma}_0^2$ 的均方 $U$，统计量 $F_0 = U/MSResidual$ 提供了对假设 $H_0\colon{\sigma}_0^2=0\mathrm{~vs~}H_a\colon{\sigma}_0^2>0$ 的检验。在 $H_0$ 的条件下，$F$ 分布为具有 u 和 v 自由度的中心 $F$ 分布，其中 u 是与 $U$ 相关的自由度，v 是与 $MSResidual$ 相关的自由度。

期望值涉及超过两个方差分量的均方通常不能用于获得具有精确 $F$ 抽样分布的单个方差分量的检验统计量。某些均衡设计会出现精确的 $F$ 分布，如前两个示例所示。比率不精确分布于 $F$ 的一个原因是各自的均方不是独立分布的。如果设计不太不均衡，那么使用 $F$ 分布作为近似应该是足够的。此外，当设计不均衡时，平方和（残差除外）不会分布为卡方分布的标量倍数。

一般来说，为了检验 $H_0\colon{\sigma}_0^2=0\mathrm{~vs~}H_a\colon{\sigma}_0^2>0$，将有一个均方，记为 $U_1$，期望为

$$E(U_1)=\sigma_\varepsilon^2+k_{1a}\sigma_a^2+k_{1b}\sigma_b^2+k_{1c}\sigma_c^2$$

但不会有其他均方具有期望 $\sigma_\varepsilon^2+k_{1b}\sigma_b^2+k_{1c}\sigma_c^2$；也就是说，没有一个均方是合适的除数。方法是找到其他均方的线性组合，如 $Q=\sum_{i=1}^kq_iMS_i$ 其中 $E(Q)=\sigma_\varepsilon^2+k_{1b}\sigma_b^2+k_{1c}\sigma_c^2$. Satterthwaite 近似可用于近似 $Q$ 的抽样分布，即，求 r，使得 $rQ/E(Q)$ 近似分布为具有 r 个自由度的卡方随机变量。该近似是双重的，因为 1) 自由度是近似的，2) 组成 $Q$ 的均方不一定按照近似所要求的那样独立分布为卡方随机变量。

SAS^®^-Mixed 代码和使用 III 型平方和得出的方差分析表，用于示例 [19.2](#sec19-1-2) 中的小麦虫害数据，如表 \@ref(tab:table20-5) 所示。variaty 均方的期望值是 $\sigma_\varepsilon^2+3.1795\sigma_\mathrm{var}^2$》 为了检验假设 $H_0\colon\sigma_{\text{var}}^2=0\text{ vs }H_a\colon\sigma_{\text{var}}^2>0$，适当的除数是残差均方，它提供的 $F$ 统计量为 4.79. 将计算出的 $F$ 统计量与具有 3 和 9 个自由度的 $F$ 分布进行比较；其显著性水平为 0.0293. 由于这是一个单向实验，因此 I 型分析与 III 型分析相同。

表 \@ref(tab:table20-6) 显示了根据示例 [19.3](#sec19-1-3) 的双向随机效应数据的 I 型平方和构建的 SAS-Mixed 代码和方差分析表。为了检验假设 $H_0\colon{\sigma}_{\mathrm{row}\times\mathrm{col}}^2=0\text{ vs }H_a\colon{\sigma}_{\mathrm{row}\times\mathrm{col}}^2>0$，适当的除数是残差均方，它提供的 $F$ 统计量为 14.24. 将计算出的 $F$ 统计量与自由度为 2 和 8 的 $F$ 分布进行比较，提供的显著性水平为 0.0023. 没有精确的检验可用于检验 $H_0\colon{\sigma^2}_{\text{row}}=0\text{ vs }H_a\colon{\sigma^2}_{\text{row}}>0$ 和 $H_0\colon{\sigma^2}_{\text{col}}=0\text{ vs }H_a\colon{\sigma^2}_{\text{col}}>0$，因此需要构建近似测试。用于检验 $H_0\colon{\sigma^2}_{\text{row}}=0\text{ vs }H_a\colon{\sigma^2}_{\text{row}}>0$ 的 $MS_{\text{Row}}$ 的适当除数计算如下

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



