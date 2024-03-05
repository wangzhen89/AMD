
# 分析三向和高阶处理结构 {#chap16}

> "Statistics are the triumph of the quantitative method, and the quantitative method is the victory of sterility and death" - Hilaire Belloc

在第 \@ref(chap7) - \@ref(chap15) 章中，我们讨论了双向处理结构的分析。这 \@ref(chap9) 章中给出的方法和结果可以推广到更复杂的处理结构；随着处理结构复杂性的增加，这种分析只会变得稍微复杂一些。我们通过具体分析三向处理结构来说明泛化方法。

在第 \@ref(sec16-1) 节中，我们给出了分析高阶处理结构时应遵循的一般策略。第 \@ref(sec16-2) 节讨论了均衡和不均衡处理结构的分析。对不均衡实验的讨论包括每个处理组合至少观测一次的情况和一些处理组合缺失的情况。

## 一般策略 {#sec16-1}

假设处理 $T_i,B_j,C_k$ 同时应用于同一实验单元。设 $\mu_{ijk}$ 表示对处理组合 $(T_i,B_j,C_k)$ 的预期响应，$i=1,2,\cdots,t;j=1,2,\cdots,b;k=1,2,\cdots,c$. 这些处理组合之间不存在三向相互作用的前提是

$$(\mu_{ijk}-\mu_{i^{\prime}jk}-\mu_{ij^{\prime}k}+\mu_{i^{\prime}j^{\prime}k})-(\mu_{ijk^{\prime}}-\mu_{i^{\prime}jk^{\prime}}-\mu_{ij^{\prime}k^{\prime}}+\mu_{i^{\prime}j^{\prime}k^{\prime}})=0\quad\mathrm{for~all~}i,i^{\prime},j,j^{\prime},k,\mathrm{and~}k^{\prime}$$

这意味着对于所有 k 和 k' 值，因素 C 在 k 水平上的 T × B 交互作用与因素 C 在 k' 水平上的 T × B 交互作用相同。类似地，T × C 交互作用在因素 B 的所有水平上都是相同的，以及 B × C 交互作用在因素 T 的所有水平上都是相同的。无交互作用陈述的等价表达为

1. $\mu_{ijk}-\bar{\mu}_{ij\cdot}-\bar{\mu}_{i\cdot k}-\bar{\mu}_{\cdot jk}+\bar{\mu}_{i\cdot\cdot}+\bar{\mu}_{\cdot j\cdot}+\bar{\mu}_{\cdot\cdot k}-\bar{\mu}_{\cdot\cdot\cdot}=0\quad\text{ for all }i,j,\mathrm{~and~}k$
2. 存在参数 $\mu,\tau_1,\tau_2,\cdots,\tau_t,\beta_1,\beta_2,\cdots,\beta_b,\xi_1,\xi_2,\ldots,\xi_c,\gamma_{11},\gamma_{12},\cdots,\gamma_{tb},\eta_{11},\eta_{12},\cdots,\eta_{tc},\theta_{11},\theta_{12},\cdots,\theta_{bc}$ 使得
$$\mu_{ijk}=\mu+\tau_i+\beta_j+\xi_k+\gamma_{ij}+\eta_{ik}+\theta_{jk}\quad\text{for all }i,j,\mathrm{~and~}k$$

也就是说，$\mu_{ijk}$ 可以用主效应和双因素交互效应来描述。**在分析三向处理结构时，第一步也是最重要的一步是确定是否存在三因素交互作用，尽管实验者可能对此不感兴趣**。如果没有三因素交互作用，那么第二步是确定是否存在任何双因素交互作用。如果也没有双因素交互作用，那么可以分析每个主效应。如果存在三因素交互作用，实验者应该在所选的第三个处理因素（通常是最不感兴趣的因素）的每个水平上分析两个处理因素的双向处理结构。显然，这些双向分析可以通过让每个处理成为所选的第三个处理来进行。使用统计计算包可以获得的分析类型与可用于双向处理结构的分析类型相似。

图 \@ref(fig:figure16-1) 显示了分析三向处理结构的一般策略。这种策略也可应用于四向和高阶处理结构。如果对所有处理组合进行相同次数的观测，则可以使用许多不同类型的统计软件来分析结果数据。

<div class="figure" style="text-align: center">
<img src="figure/figure 16.1.png" alt="三因素实验的分析策略" width="648" />
<p class="caption">(\#fig:figure16-1)三因素实验的分析策略</p>
</div>

## 均衡和不均衡实验 {#sec16-2}

如果观测到所有处理组合，但观测次数不相等，则可以使用 III 型分析。如果每种处理组合都至少观测一次，那么所有的主效应和交互作用假设仍然可以得到检验，并且回答的问题与在处处完全均衡时回答的问题相同。

如果某些处理组合缺失，那么就像第 \@ref(chap13) 章中的情况一样，无法检验涉及缺失处理组合的假设。实验者应该在感兴趣的处理之间指定她自己的 IV 型假设。这种假设可以通过使用第 \@ref(chap1) 章中描述的矩阵程序或通过许多统计计算包中可用的 Contrast 语句进行检验。

第 \@ref(chap17) 章中的示例演示了可能需要的一些步骤，以便获得对缺失处理组合数据的完整分析。

## I 型和 II 型分析 {#sec16-3}

第 \@ref(chap10) 章描述了双向处理结构实验的 I 型和 II 型分析。这两种分析都使用模型比较程序为每个效应产生平方和。I 型分析序贯式地 (sequentially) 拟合模型，每个效应都会针对模型中之前的所有其他效应进行调整。在第 II 型分析中，针对相同或更低水平的所有其他效应调整每个效应的平方和。考虑由下式给出的三向模型

$$\begin{aligned}y_{ijk\ell}&=\mu+T_i+B_j+(TB)_{ij}+C_k+(TC)_{ik}+(BC)_{jk}+(TBC)_{ijk}+\varepsilon_{ijk\ell}\\&\mathrm{for~}i=1,2,\ldots,t,~j=1,2,\ldots,b;~k=1,2,\ldots,c;~\mathrm{and}~\ell=1,2,\ldots,n_{ijk}\end{aligned}$$

假设对于所有 $i,j,k$ 有 $n_{ijk}>0$. 也就是说，每个三向单元格至少观测一次。表 \@ref(tab:table16-1) 显示了使用第 \@ref(chap10) 章中描述的缩减符号的 I 型和 II 型平方和及其相应的自由度。

<table>
<caption>(\#tab:table16-1)在每个单元格中至少有一次观测的三向实验的 I 型和 II 型平方和</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2016.1.png">
</td>
  </tr>
</tbody>
</table>

## 结束语 {#sec16-4}

本章讨论了三向及高阶处理结构实验的分析。给出了一个流程图，提供了分析此类实验的一般策略。

首先研究最高阶交互作用效应是很重要的。许多实验者避免考虑高阶交互作用，因为他们通常不太确定如何处理这些交互作用。应该避免这种诱惑。尽管经验告诉我们，非常高阶的交互作用很少是有意义的，但无论何时，都必须处理它们。第 \@ref(chap8) 章中讨论的技术可以推广到三向和更高阶的处理结构，并有助于确定哪些处理组合导致了交互作用。交点的原因可能是研究中确定的最重要的信息。

## 练习 {#sec16-5}
