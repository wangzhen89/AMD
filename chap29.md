
# 交叉设计分析 {#chap29}

交叉设计 (crossover design) 用于比较按序列给予实验单元（如动物或人）的处理。也就是说，每个实验单元按照预定的序列进行每种处理。交叉设计的目的是通过观察应用于同一实验单元的处理，消除比较处理时实验单元之间的差异。

尽管交叉设计消除了处理比较中实验单元之间的差异，但其他问题可能会以**延滞** (carryover) 或**残余** (residual) 效应的形式出现。例如，当首先给予处理 $A$ 并且在应用处理 $B$ 时其效应尚未消失，就会出现延滞效应。如果 $A$  的残留效应干扰了个体对处理 $B$ 的响应（或正向或负向），那么处理 $A$ 对处理 $B$ 有残余效应。

交叉设计模型必须包含序列效应、时期效应、处理效应、延滞效应、实验单元误差项和时间间隔误差项。第一节讨论了一般模型及其假设，最后两节分别讨论了两个时期和两个以上时期两个处理的交叉设计。

在一般的交叉设计中，比较 $t$ 个处理，其中在每个实验单元上观察每个处理；也就是说，这些处理以特定的序列应用于实验单元。实验者构建 $t$ 个处理的 $s$ 个序列，并将实验单元随机分配给第 $i$ 个序列。表 \@ref(tab:table29-1) 包含一组可应用于实验单元的三种处理的序列。将序列分配给个体（可能的实验单元）意味着个体是序列的实验单元；将处理分配给时间间隔意味着该时间间隔是处理的实验单元。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-1)对实验单元应用三种处理（A、B 和 C）的可能序列集</caption>
 <thead>
  <tr>
   <th style="text-align:left;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;">  <img src="table/table%2029.1.png">
</td>
  </tr>
</tbody>
</table>

## 定义，假设和模型 {#sec29-1}

描述在第 $i$ 个时期内分配给第 $j$ 个序列的第 $\ell$ 只动物的观测结果的响应的模型为

\begin{align}
y_{ijk\ell}&=\mu+S_i+\delta_{i\ell}+P_j+T_k+\varepsilon_{ijk\ell}\\i&=1,2,\ldots,s,j=1,2,\ldots,p,k=1,2,\ldots,t,\mathrm{~and~}\ell=1,2,\ldots,n_j
(#eq:29-1)
\end{align}

在上述模型中，$S_i$ 是第 $i$ 个序列的效应，$P_j$ 是第 $j$ 个周期的效应，$T_k$ 是第 $k$ 个处理的效应，其中 $k$ 的值由第 $i$ 个序列和第 $j$ 个时期的 $i,j$ 组合确定。在理想条件下，将有 $\delta_{i\ell}\sim i.i.d.N(0,\sigma_{\delta}^{2}),\varepsilon_{ijk\ell}\sim i.i.d.N(0,\sigma_{\varepsilon}^{2})$，并且所有 $\delta_{i\ell}$ 和 $\varepsilon_{ijk\ell}$ 彼此独立。

误差项 $\delta_{i\ell}$ 和 $\varepsilon_{ijk\ell}$ 的理想条件对于进行适当的分析非常重要。由于实验单位是随机分配给序列的，因此适当的假设是 $\delta_{i\ell}\sim i.i.d.N(0,\sigma_{\delta}^{2})$. 由于 $\varepsilon_{ijk\ell}$ 是时间间隔的误差，并且在某种意义上是重复测量误差，因此这些误差可能不是独立分布的。上面给出的理想条件相当于假设时间周期误差的协方差阵满足复合对称条件（参见 \@ref(sec27-1) 节）。在这种情况下，可以使用通常的独立的 $\delta_{i\ell}$ 和 $\varepsilon_{ijk\ell}$ 方差分析方法来分析观测数据。复合对称性假设对于两周期设计是必须满足的，对于许多三周期设计也可能是满足的。时间序列误差结构可能更适合具有三个或更多周期的交叉设计，并且应该使用第 \@ref(chap27) 章中讨论的方法。

两周期/两处理交叉设计 (two-period/two-treatment crossover design) 在 \@ref(sec29-2) 节中描述，具有两个以上处理或两个以上时期的设计在 \@ref(sec29-3) 和 \@ref(sec29-4) 节中描述。

## 两周期/两处理交叉设计 {#sec29-2}

考虑具有两个处理序列（$AB$ 和 $BA$）的两周期/两处理交叉设计，其中 $AB$ 意味着处理 $A$ 在第 $1$ 时期分配到一个实验单元，处理 $B$ 在第 $2$ 时期分配到同一实验单元。同样，$BA$ 序列意味着处理 $B$ 在时期 $1$ 被分配到实验单元，处理 $A$ 在时期 $2$ 被分配给相同的实验单元。还假设有 $n_1$ 个实验单元分配给 $AB$ 序列，且有 $n_2$ 个实验单元分配给 $BA$ 序列。令观测到的响应的模型为

\begin{align}
y_{ijk\ell}&=\mu+S_i+\delta_{i\ell}+P_j+T_k+\varepsilon_{ijk\ell}\\i&=AB,BA,~j=1,2,~k=A,B,~\mathrm{and}~\ell=1,2,\ldots,n_i
(#eq:29-2)
\end{align}

其中下标 $k$ 的值由序列 $i$ 和周期 $j$ 的组合确定。令

$$\mu_{ij}=\mu+S_i+P_j+T_k$$

是对应于第 $i$ 个序列和第 $j$ 个时期的单元均值的模型。令 $\hat \mu_{ij}$ 为第 $(i,j)$ 单元格中所有观测值的观测均值。注意到

\begin{align}
\hat \mu_{11}=\bar y_{11\cdot\cdot}\,\mathrm{estimates}\,\mu+S_{AB}+P_1+T_A\\
\hat \mu_{12}=\bar y_{12\cdot\cdot}\,\mathrm{estimates}\,\mu+S_{AB}+P_2+T_B\\
\hat \mu_{21}=\bar y_{21\cdot\cdot}\,\mathrm{estimates}\,\mu+S_{BA}+P_1+T_B\\
\hat \mu_{22}=\bar y_{22\cdot\cdot}\,\mathrm{estimates}\,\mu+S_{BA}+P_2+T_A
(#eq:29-3)
\end{align}

俄
