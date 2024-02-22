
# 交叉设计分析 {#chap29}

交叉设计 (crossover design) 用于比较按序列给予实验单元（如动物或人）的处理。也就是说，每个实验单元按照预定的序列进行每种处理。交叉设计的目的是通过观察应用于同一实验单元的处理，消除比较处理时实验单元之间的差异。

尽管交叉设计消除了处理比较中实验单元之间的差异，但其他问题可能会以**延滞** (carryover) 或**残余** (residual) 效应的形式出现。例如，当首先给予处理 $A$ 并且在应用处理 $B$ 时其效应尚未消失，就会出现延滞效应。如果 $A$ 的滞留效应 (lingering effect) 干扰了个体对处理 $B$ 的响应（或正向或负向），那么处理 $A$ 对处理 $B$ 有延滞效应。

交叉设计模型必须包含序列效应、时期效应、处理效应、延滞效应、实验单元误差项和时间间隔误差项。第一节讨论了一般模型及其假设，最后两节分别讨论了两个时期和两个以上时期两个处理的交叉设计。

在一般的交叉设计中，比较 $t$ 个处理，其中在每个实验单元上观察每个处理；也就是说，这些处理以特定的序列应用于实验单元。实验者构建 $t$ 个处理的 $s$ 个序列，并将实验单元随机分配给第 $i$ 个序列。表 \@ref(tab:table29-1) 包含一组可应用于实验单元的三种处理的序列。将序列分配给个体（可能的实验单元）意味着个体是序列的实验单元；将处理分配给时间间隔意味着该时间间隔是处理的实验单元。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-1)对实验单元应用三种处理（A、B 和 C）的可能序列集</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.1.png">
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

在上述模型中，$S_i$ 是第 $i$ 个序列的效应，$P_j$ 是第 $j$ 个时期的效应，$T_k$ 是第 $k$ 个处理的效应，其中 $k$ 的值由第 $i$ 个序列和第 $j$ 个时期的 $i,j$ 组合确定。在理想条件下，将有 $\delta_{i\ell}\sim i.i.d.N(0,\sigma_{\delta}^{2}),\varepsilon_{ijk\ell}\sim i.i.d.N(0,\sigma_{\varepsilon}^{2})$，并且所有 $\delta_{i\ell}$ 和 $\varepsilon_{ijk\ell}$ 彼此独立。

误差项 $\delta_{i\ell}$ 和 $\varepsilon_{ijk\ell}$ 的理想条件对于进行适当的分析非常重要。由于实验单元是随机分配给序列的，因此适当的假设是 $\delta_{i\ell}\sim i.i.d.N(0,\sigma_{\delta}^{2})$. 由于 $\varepsilon_{ijk\ell}$ 是时间间隔的误差，并且在某种意义上是重复测量误差，因此这些误差可能不是独立分布的。上面给出的理想条件相当于假设时间时期误差的协方差阵满足复合对称条件（参见 \@ref(sec27-1) 节）。在这种情况下，可以使用通常的独立的 $\delta_{i\ell}$ 和 $\varepsilon_{ijk\ell}$ 方差分析方法来分析观测数据。复合对称性假设对于两时期设计是必须满足的，对于许多三时期设计也可能是满足的。时间序列误差结构可能更适合具有三个或更多时期的交叉设计，并且应该使用第 \@ref(chap27) 章中讨论的方法。

两时期/两处理交叉设计 (two-period/two-treatment crossover design) 在 \@ref(sec29-2) 节中描述，具有两个以上处理或两个以上时期的设计在 \@ref(sec29-3) 和 \@ref(sec29-4) 节中描述。

## 两时期/两处理交叉设计 {#sec29-2}

考虑具有两个处理序列（$AB$ 和 $BA$）的两时期/两处理交叉设计，其中 $AB$ 意味着处理 $A$ 在第 $1$ 时期分配到一个实验单元，处理 $B$ 在第 $2$ 时期分配到同一实验单元。同样，$BA$ 序列意味着处理 $B$ 在时期 $1$ 被分配到实验单元，处理 $A$ 在时期 $2$ 被分配给相同的实验单元。还假设有 $n_1$ 个实验单元分配给 $AB$ 序列，且有 $n_2$ 个实验单元分配给 $BA$ 序列。令观测到的响应的模型为

\begin{align}
y_{ijk\ell}&=\mu+S_i+\delta_{i\ell}+P_j+T_k+\varepsilon_{ijk\ell}\\i&=AB,BA,~j=1,2,~k=A,B,~\mathrm{and}~\ell=1,2,\ldots,n_i
(#eq:29-2)
\end{align}

其中下标 $k$ 的值由序列 $i$ 和时期 $j$ 的组合确定。令

$$\mu_{ij}=\mu+S_i+P_j+T_k$$

是对应于第 $i$ 个序列和第 $j$ 个时期的单元均值的模型。令 $\hat \mu_{ij}$ 为第 $(i,j)$ 单元格中所有观测值的观测均值。注意到

\begin{align}
\hat \mu_{11}=\bar y_{11\cdot\cdot}\,\mathrm{estimates}\,\mu+S_{AB}+P_1+T_A\\
\hat \mu_{12}=\bar y_{12\cdot\cdot}\,\mathrm{estimates}\,\mu+S_{AB}+P_2+T_B\\
\hat \mu_{21}=\bar y_{21\cdot\cdot}\,\mathrm{estimates}\,\mu+S_{BA}+P_1+T_B\\
\hat \mu_{22}=\bar y_{22\cdot\cdot}\,\mathrm{estimates}\,\mu+S_{BA}+P_2+T_A
(#eq:29-3)
\end{align}

从式 \@ref(eq:29-3) 可看出

\begin{align}
\hat {\bar\mu}_{1\cdot}&=\frac{\bar y_{11\cdot\cdot}+\bar y_{12\cdot\cdot}}{2}\,\mathrm{estimates}\,\mu+S_{AB}+\bar P_\cdot+\frac{T_A+T_B}{2}\\
\hat {\bar\mu}_{2\cdot}&=\frac{\bar y_{21\cdot\cdot}+\bar y_{22\cdot\cdot}}{2}\,\mathrm{estimates}\,\mu+S_{BA}+\bar P_\cdot+\frac{T_A+T_B}{2}\\
\hat {\bar\mu}_{\cdot 1}&=\frac{\bar y_{11\cdot\cdot}+\bar y_{21\cdot\cdot}}{2}\,\mathrm{estimates}\,\mu+\bar S_{\cdot}+\bar P_1+\frac{T_A+T_B}{2}\\
\hat {\bar\mu}_{\cdot 2}&=\frac{\bar y_{12\cdot\cdot}+\bar y_{22\cdot\cdot}}{2}\,\mathrm{estimates}\,\mu+\bar S_{\cdot}+\bar P_2+\frac{T_A+T_B}{2}
(#eq:29-4)
\end{align}

另请注意，两个序列的均值差 $\hat {\bar\mu}_{1\cdot}-\hat {\bar\mu}_{2\cdot}$ 估计 $S_{AB} - S_{BA}$，即两个序列参数之差。时期 $1$ 和时期 $2$ 的均值差 $\hat {\bar\mu}_{\cdot1}-\hat {\bar\mu}_{\cdot2}$ 估计 $P_1 - P_2$. 最后，为了估计两种处理之差 $T_A - T_B$，可以使用 $(\hat {\bar\mu}_{11}-\hat {\bar\mu}_{12}-\hat {\bar\mu}_{21}+\hat {\bar\mu}_{22})/2$ 。在理想条件下，可以证明

\begin{equation}
\begin{bmatrix}\hat{\mu}_{11}\\\hat{\mu}_{12}\\\hat{\mu}_{21}\\\hat{\mu}_{22}\end{bmatrix}\thicksim N\left(\begin{bmatrix}\mu_{11}\\\mu_{12}\\\mu_{21}\\\mu_{22}\end{bmatrix},\,\begin{bmatrix}\frac{\sigma_\varepsilon^2+\sigma_\delta^2}{n_1}&\frac{\sigma_\delta^2}{n_1}&0&0\\\\\frac{\sigma_\delta^2}{n_1}&\frac{\sigma_\epsilon^2+\sigma_\delta^2}{n_1}&0&0\\\\0&0&\frac{\sigma_\epsilon^2+\sigma_\delta^2}{n_2}&\frac{\sigma_\delta^2}{n_2}\\\\0&0&\frac{\sigma_\delta^2}{n_2}&\frac{\sigma_\epsilon^2+\sigma_\delta^2}{n_2}\end{bmatrix}\right)
(#eq:29-5)
\end{equation}

从式 \@ref(eq:29-5) 可证明

$$\begin{aligned}&\hat{\bar{\mu}}_{1\cdot}-\hat{\bar{\mu}}_{2\cdot}\thicksim N\bigg[S_{AB}-S_{BA},\bigg(\frac{\sigma_{\varepsilon}^{2}+2\sigma_{\delta}^{2}}{2}\bigg)\bigg(\frac{1}{n_{1}}+\frac{1}{n_{2}}\bigg)\bigg]\\&\hat{\bar{\mu}}_{\cdot1}-\hat{\bar{\mu}}_{\cdot2}\thicksim N\bigg[P_1-P_2,\frac{\sigma_{\varepsilon}^{2}}{2}\bigg(\frac{1}{n_1}+\frac{1}{n_2}\bigg)\bigg]\\&\frac{\hat\mu_{11}-\hat\mu_{12}-\hat\mu_{21}+\hat\mu_{22}}{2}\thicksim N\bigg[T_A-T_B,\frac{\sigma_{\varepsilon}^{2}}{2}\bigg(\frac{1}{n_1}+\frac{1}{n_2}\bigg)\bigg]\end{aligned}$$

与本节讨论的两期/两处理交叉设计相对应的方差分析表的形式如表 \@ref(tab:table29-2) 所示。用于检验 $H_0: T_A = T_B$ 的 $F$ 统计量由 $F = TMS/WSEMS$ 给出，其中 $TMS$ 是表 \@ref(tab:table29-2) 中的处理均方 (treatment mean square)，$WSEMS$ 是个体内误差均方 (within subject error mea square). 若 $F>F_{\alpha,1,n_1+n_2-2}$ 则拒绝 $H_0$. $T_A-T_B$ 的 $(1-\alpha)100\%$ 置信区间由

$$\frac{\hat{\mu}_{11}-\hat{\mu}_{12}-\hat{\mu}_{21}+\hat{\mu}_{22}}2\pm t_{\alpha,v}\sqrt{\frac{\hat\sigma_{\varepsilon}^2}2{\left(\frac1{n_1}+\frac1{n_2}\right)}}$$

给出。其中 $v=n_1+n_2-n$ 且 $\hat\sigma_\varepsilon^2=WSEMS$. 作为一个例子，考虑 Grizzle (1965) 使用的数据。数据如表 \@ref(tab:table29-3) 所示。使用表 \@ref(tab:table29-4) 中所示的命令通过 SAS^®^-GLM 程序分析数据。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-2)两期/两处理交叉实验的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.2.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-3)Grizzle (1965) 数据</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> Seq </th>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> Period </th>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> Trt </th>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> Person </th>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> Y </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> AB </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 0.2 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> AB </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> AB </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> AB </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> -0.7 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> AB </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> -0.8 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> AB </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 0.2 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> AB </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 14 </td>
   <td style="text-align:center;"> 0.6 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> AB </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 14 </td>
   <td style="text-align:center;"> 1.1 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> AB </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 15 </td>
   <td style="text-align:center;"> 0.3 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> AB </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 15 </td>
   <td style="text-align:center;"> 0.4 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> AB </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 16 </td>
   <td style="text-align:center;"> 1.5 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> AB </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 16 </td>
   <td style="text-align:center;"> 1.2 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BA </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 21 </td>
   <td style="text-align:center;"> 1.3 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BA </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 21 </td>
   <td style="text-align:center;"> 0.9 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BA </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 22 </td>
   <td style="text-align:center;"> -2.3 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BA </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 22 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BA </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 23 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BA </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 23 </td>
   <td style="text-align:center;"> 0.6 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BA </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 24 </td>
   <td style="text-align:center;"> -0.8 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BA </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 24 </td>
   <td style="text-align:center;"> -0.3 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BA </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 25 </td>
   <td style="text-align:center;"> -0.4 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BA </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 25 </td>
   <td style="text-align:center;"> -1 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BA </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 26 </td>
   <td style="text-align:center;"> -2.9 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BA </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 26 </td>
   <td style="text-align:center;"> 1.7 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BA </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 27 </td>
   <td style="text-align:center;"> -1.9 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BA </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 27 </td>
   <td style="text-align:center;"> -0.3 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BA </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 28 </td>
   <td style="text-align:center;"> -2.9 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BA </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 28 </td>
   <td style="text-align:center;"> 0.9 </td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-4)用于分析表 \@ref(tab:table29-3) 中数据的 SAS-GLM 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> TITLE 'CRSOVER EXAMPLE 29.1 – A TWO PERIOD/TWO TREATMENT DESIGN';\
DATA GRIZ;\
INPUT SEQ $ PERIOD TRT $ PERSON Y;\
LINES;\
AB 1 A 11 0.2\
AB 2 B 11 1.0\
AB 1 A 12 0.0\
AB 2 B 12 -0.7\
AB 1 A 13 -0.8\
AB 2 B 13 0.2\
AB 1 A 14 0.6\
AB 2 B 14 1.1\
AB 1 A 15 0.3\
AB 2 B 15 0.4\
AB 1 A 16 1.5\
AB 2 B 16 1.2\
BA 1 B 21 1.3\
BA 2 A 21 0.9\
BA 1 B 22 -2.3\
BA 2 A 22 1.0\
BA 1 B 23 0.0\
BA 2 A 23 0.6\
BA 1 B 24 -0.8\
BA 2 A 24 -0.3\
BA 1 B 25 -0.4\
BA 2 A 25 -1.0\
BA 1 B 26 -2.9\
BA 2 A 26 1.7\
BA 1 B 27 -1.9\
BA 2 A 27 -0.3\
BA 1 B 28 -2.9\
BA 2 A 28 0.9\
;\
PROC GLM;\
 TITLE2 'STATISTICAL ANALYSIS USING SAS-GLM';\
 CLASSES SEQ TRT PERIOD PERSON;\
 MODEL Y=SEQ PERSON(SEQ) TRT PERIOD;\
 LSMEANS TRT/PDIFF;\
 RANDOM PERSON(SEQ);\
 RUN; </td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table29-5) 显示了模型的方差分析，并给出了 $WSEMS$ 的值 $\hat \sigma_\varepsilon^2 = 1.245$. 表 \@ref(tab:table29-6) 给出了 $H_0:T_A = T_B$ 的检验： 以及个体间误差均方 $BSEMS = \hat \sigma_{\varepsilon} ^2 + 2\hat \sigma_{\delta}^2 = 1.001$. $H_0$ 观察到的显著性水平为 $\hat \alpha = 0.1165$. 表 \@ref(tab:table29-7) 中显示了处理 $A$ 和 $B$ 的最小二乘均值，以及比较处理 $A$ 和处理 $B$ 时观察到的显著性水平 $\hat \alpha = 0.1165$。方差分析表中各行的预期均方如表 \@ref(tab:table29-8) 所示。如果想求出双向序列×时期的估计值，可以使用表 \@ref(tab:table29-9) 中的命令。双向均值如表 \@ref(tab:table29-10) 所示。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-5)具有个体内误差均方的模型方差分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.5.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-6)表 \@ref(tab:table29-3) 中数据的 III 型方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.6.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-7)处理主效应均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.7.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-8)表 \@ref(tab:table29-3) 中数据的预期均方表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.8.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-9)用于分析表 \@ref(tab:table29-3) 中数据的 SAS-GLM 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> PROC GLM;\
 TITLE2 'STATISTICAL ANALYSIS USING SAS-GLM';\
 CLASSES SEQ TRT PERIOD PERSON;\
 MODEL Y=SEQ PERSON(SEQ) PERIOD SEQ*PERIOD;\
 LSMEANS SEQ*PERIOD;\
 RANDOM PERSON(SEQ);\
RUN; </td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-10)序列×时期均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.8.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table29-3) 中的数据也可以使用表 \@ref(tab:table29-11) 中所示的命令通过 SAS Mixed 程序进行分析。由于 GLM 分析中受试者间误差均方小于受试者内误差均方，因此结果与 GLM 分析获得的结果略有不同。因此，Mixed 程序将与个体相对应的方差分量估计为 $\hat \sigma _\delta ^2 = 0$，而 GLM 将使用矩准则方法估计相同的方差分量：

$$\hat\sigma_\delta^2=\frac{BSEMS-WSEMS}2=\frac{1.005-1.245}2=\frac{-0.240}2=-0.120$$

目前未给出使用表 \@ref(tab:table29-11) 中的代码的结果，感兴趣的读者需要运行表 \@ref(tab:table29-11) 中的代码。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-11)用于分析表 \@ref(tab:table29-3) 中数据的 SAS-Mixed 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> PROC MIXED;\
 TITLE2 ‘STATISTICAL ANALYSIS USING SAS-MIXED’;\
 CLASSES SEQ TRT PERIOD PERSON;\
 MODEL Y=SEQ PERIOD TRT;\
 LSMEANS TRT/PDIFF;\
 RANDOM PERSON(SEQ);\
 RUN;\
 PROC MIXED;\
 TITLE2 ‘STATISTICAL ANALYSIS USING SAS-MIXED’;\
 CLASSES SEQ TRT PERIOD PERSON;\
 MODEL Y=SEQ PERIOD SEQ*PERIOD;\
 LSMEANS SEQ*PERIOD;\
 RANDOM PERSON(SEQ);\
 RUN; </td>
  </tr>
</tbody>
</table>

上述分析假设第一时期的处理不会影响第二时期的观测结果。接下来考虑存在延滞的情况。当存在延滞时，可以假设

\begin{align}
\hat \mu_{11}=\bar y_{11\cdot\cdot}\,\mathrm{estimates}\,\mu+S_{AB}+P_1+T_A\\
\hat \mu_{12}=\bar y_{12\cdot\cdot}\,\mathrm{estimates}\,\mu+S_{AB}+P_2+T_B+\lambda_A\\
\hat \mu_{21}=\bar y_{21\cdot\cdot}\,\mathrm{estimates}\,\mu+S_{BA}+P_1+T_B\\
\hat \mu_{22}=\bar y_{22\cdot\cdot}\,\mathrm{estimates}\,\mu+S_{BA}+P_2+T_A+\lambda_B
(#eq:29-6)
\end{align}

其中 $\lambda_A$ 是序列 $AB$ 的参数，对应于从时期 $1$ 施加的处理 $A$ 到时期 $2$ 的延滞，而 $\lambda_B$ 是序列 $BA$ 的参数，对应于从时期 $1$ 施加的处理 $B$ 到时期 $2$ 的延滞。

从式 \@ref(eq:29-6) 可看出

\begin{align}
\hat {\bar\mu}_{1\cdot}&=\frac{\bar y_{11\cdot\cdot}+\bar y_{12\cdot\cdot}}{2}\,\mathrm{estimates}\,\mu+S_{AB}+\bar P_\cdot+\frac{T_A+T_B}{2}+\frac{\lambda_A}{2}\\
\hat {\bar\mu}_{2\cdot}&=\frac{\bar y_{21\cdot\cdot}+\bar y_{22\cdot\cdot}}{2}\,\mathrm{estimates}\,\mu+S_{BA}+\bar P_\cdot+\frac{T_A+T_B}{2}+\frac{\lambda_B}{2}\\
\hat {\bar\mu}_{\cdot 1}&=\frac{\bar y_{11\cdot\cdot}+\bar y_{21\cdot\cdot}}{2}\,\mathrm{estimates}\,\mu+\bar S_{\cdot}+\bar P_1+\frac{T_A+T_B}{2}\\
\hat {\bar\mu}_{\cdot 2}&=\frac{\bar y_{12\cdot\cdot}+\bar y_{22\cdot\cdot}}{2}\,\mathrm{estimates}\,\mu+\bar S_{\cdot}+\bar P_2+\frac{T_A+T_B}{2}+\frac{\lambda_A+\lambda_B}{2}
(#eq:29-7)
\end{align}

当 $\lambda_A\ne \lambda_B$ 时，延滞就存在。也就是说，在无延滞的情况下，$\lambda_A$ 和 $\lambda_B$ 实际上不必等于零，它们只需要彼此相等。在 $\lambda_A= \lambda_B$ 的情况下，$\lambda_A$ 和 $\lambda_B$ 与时期效应混淆 (confounded)，不需要将参数包括在模型中 \@ref(eq:29-2) 中。在存在延滞的情况下，两个序列差的均值 $\hat {\bar\mu}_{1\cdot}-\hat {\bar\mu}_{2\cdot}$ 估计 $S_{AB}-S_{BA}+[(\lambda_A-\lambda_B)/2]$。时期 $1$ 和时期 $2$ 差的均值 $\hat {\bar\mu}_{\cdot1}-\hat {\bar\mu}_{\cdot2}$ 估计 $P_1-P_2+[(\lambda_A+\lambda_B)/2]$ 。最后，$(\hat{\mu}_{11}-\hat{\mu}_{12}-\hat{\mu}_{21}+\hat{\mu}_{22})/2$ 估计 $T_A-T_B+[(\lambda_A-\lambda_B)/2]$ 。在存在延滞的情况下，时期 $2$ 的数据无助于估计两种处理的差异。当存在延滞时，必须通过 $\hat \mu_{11}-\hat \mu_{21}$ 来估计 $T_A-T_B$，并且该估计仅取决于时期 $1$ 的数据。还可以注意到，如果仅使用时期 $1$ 的数据来估计处理效应，则时期 $1$ 的模型可以简化为

\begin{equation}
y_{k\ell}=\mu+T_k+\varepsilon_{k\ell}^*,~k=A,B,~\mathrm{and~}~\ell=1,2,\ldots,n_i
(#eq:29-8)
\end{equation}

其中 $\varepsilon_{k\ell}^{*}=\delta_{k\ell}+\varepsilon_{k\ell}$. 而且，

$$\mathrm{Var}(\hat{\mu}_{11}-\hat{\mu}_{21})=(\sigma_\varepsilon^2+\sigma_\delta^2){\left(\frac1{n_1}+\frac1{n_2}\right)}$$

表 \@ref(tab:table29-12) 显示了当存在延滞时，与两时期/两处理交叉设计相对应的方差分析表的形式。该表与表 \@ref(tab:table29-2) 具有相同的形式，除了预期均方列。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-12)两时期/两处理交叉实验的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.12.png">
</td>
  </tr>
</tbody>
</table>

如果不存在延滞，那么个体被分配到哪个序列应该无关紧要。也就是说，从哲学上讲，$S_{AB}-S_{BA}$ 应该等于零。因此，$F=SeqMS/BSEM$ 给出的 $F$ 统计量提供了是否存在延滞的测试。如果 $F>F_{\alpha,1,n1+n2-2}$，则可以得出结论，存在显著的延滞效应。需要注意的是，延滞检验是个体间的比较，因此该检验的功效不如基于个体内的比较。如果认为可能存在延滞效应，建议不要使用两时期/两处理的交叉设计。在某些情况下，实验者可能能够在交叉设计的两个阶段之间包括所谓的“洗脱” (wash-out) 阶段。洗脱期是指一段足够长的时间，以便在对实验单元进行第二种处理之前，第一时期中给予的处理的任何延滞效应将被消除或清除。同样极其重要的是，要注意在两时期/两处理交叉设计中的延滞问题在一些具有两个以上处理和/或两个以上时期的交叉设计中减少和/或消除。此类设计将在下一节中进行考虑。

作为可能存在延滞的两时期/两处理交叉实验的示例，再次考虑表 \@ref(tab:table29-3) 中显示的 Grizzle (1965) 数据。将使用 SAS-Mixed 程序（使用表 \@ref(tab:table29-13) 中显示的命令）分析数据。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-13)用于分析表 \@ref(tab:table29-3) 中数据的 SAS-Mixed 代码，当存在延滞时</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> PROC MIXED;\
 TITLE2 'STATISTICAL ANALYSIS USING SAS-MIXED';\
 CLASSES SEQ TRT PERIOD PERSON;\
 MODEL Y=SEQ PERIOD TRT;\
 LSMEANS TRT/PDIFF;\
 RANDOM PERSON(SEQ);\
RUN;\
PROC MIXED;\
TITLE2 ‘STATISTICAL ANALYSIS USING SAS-MIXED’;\
 CLASSES SEQ TRT PERIOD PERSON;\
 MODEL Y=SEQ PERIOD SEQ\*PERIOD;\
 LSMEANS SEQ\*PERIOD;\
 ESTIMATE ‘A-B FROM PERIOD 1’ SEQ 1 -1 SEQ*PERIOD 1 0 -1 0;\
 RANDOM PERSON(SEQ);\
 RUN; </td>
  </tr>
</tbody>
</table>

第一组 Mixed 命令用于获得延滞检验和处理效应检验（如果没有任何显著的延滞）。这些命令还给出了处理主效应均值，如果没有延滞，则可以使用该均值。第二组命令使用不同的 random 语句，以获得序列×时期双向均值，包括 `ESTIMATE` 选项，以便从时期 $1$ 数据中获得处理效应估计值。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-14)固定效应检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.14.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table29-14) 给出了对表 \@ref(tab:table29-13) 中第一组 Mixed 命令的固定效应的检验。请注意，序列效应在 $0.0665$ 水平下显著，表明存在一些延滞的证据。还请注意，处理效应仅在 $0.1002$ 水平下显著。
处理主效应平均值见表 \@ref(tab:table29-15)，该表还比较了两种处理平均值。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-15)处理均值及两种处理均值差的检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.15-1.png">
</td>
  </tr>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.15-2.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table29-16) 给出了表 \@ref(tab:table29-13) 中第二组 Mixed 命令的固定效应检验。请注意，序列的检验与表 \@ref(tab:table29-14) 中给出的检验相同。还请注意，序列×时期的检验与表 \@ref(tab:table29-14) 中处理的检验相同。表 \@ref(tab:table29-17) 给出了序列×时期均值，表 \@ref(tab:table29-18) 给出了仅来自时期 $1$ 数据比较两种处理的检验。表 \@ref(tab:table29-18) 的检查显示，在表 \@ref(tab:table29-7) 中给出的分析中，延滞效应明显掩盖了部分处理效应，其中假设无延滞，因为处理效应不显著，而仅使用时期 $1$ 数据的比较发现处理存在显著差异。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-16)固定效应的额外检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.16.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-17)序列×时期双向均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.17.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-18)处理 $A$ 和处理 $B$，来自时期 $1$ 的估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.18.png">
</td>
  </tr>
</tbody>
</table>

## 具有两个以上时期的交叉设计 {#sec29-3}

当没有延滞时，具有两个以上时期和/或两个以上处理的交叉设计的统计分析是简单的。当满足个体间误差和个体内误差的理想条件时，式 \@ref(eq:29-1) 中的模型是合适的。表 \@ref(tab:table29-4)，\@ref(tab:table29-9)，和 \@ref(tab:table29-11) 仍然可以用来获得有用和有趣的统计量。由于没有延滞时的分析是直接的，本节集中讨论存在延滞的情况。

考虑的第一种情况是仍然有两种处理的情况，但涉及具有三个时期的序列。一种可能性是使用两个序列 $ABA$ 和 $BAB$. 对于 $ABA$ 序列，实验单元将在时期 $1$ 接受处理 $A$，在时期 $2$ 接受处理 $B$，并在时期 $3$再次接受处理 $A$. $BAB$ 序列将进行类似处理。另一种可能性是使用两个序列 $ABB$ 和 $BAA$. 第三种可能性是使用所有这四个序列。也就是说，$ABA$、$BAB$、$ABB$ 和 $BAA$ 都将被使用。本节重点介绍使用两个序列 $ABA$ 和 $BAB$ 的交叉设计。其他的可能性也可以类似地分析。

表 \@ref(tab:table29-19) 确定了与 $ABA$ 和 $BAB$ 交叉设计的序列与时期组合相关的模型参数。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-19)三时期/两处理交叉设计的单元个均值参数</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.19.png">
</td>
  </tr>
</tbody>
</table>

假设 $n_1$ 个个体已被分配到 $ABA$ 序列，$n2$ 个个体已被分配到 $BAB$ 序列。令 $y_{ij\ell}$ 为在从序列 $i$ 时期 $j$ 的个体 $\ell$ 中观察到的响应。令 $\hat\mu _{ij} = \bar y_{ij\cdot} , i = 1, 2; j = 1, 2, 3$. 在式 \@ref(eq:29-1) 中给出的模型误差项的理想条件下，可以证明： 

\begin{equation}
\begin{bmatrix}\hat{\mu}_{11}\\\hat{\mu}_{12}\\\hat{\mu}_{13}\\\hat{\mu}_{21}\\\hat{\mu}_{22}\\\hat{\mu}_{23}\end{bmatrix}\sim N\left(\begin{bmatrix}\mu_{11}\\\mu_{12}\\\mu_{13}\\\mu_{21}\\\mu_{21}\\\mu_{23}\end{bmatrix},\begin{bmatrix}\frac{\sigma_\varepsilon^2+\sigma_\delta^2}{n_1}&\frac{\sigma_\delta^2}{n_1}&\frac{\sigma_\delta^2}{n_1}&0&0&0\\\frac{\sigma_\delta^2}{n_1}&\frac{\sigma_\varepsilon^2+\sigma_\delta^2}{n_1}&\frac{\sigma_\delta^2}{n_1}&0&0&0\\\frac{\sigma_\delta^2}{n_1}&\frac{\sigma_\delta^2}{n_1}&\frac{\sigma_\epsilon^2+\sigma_\delta^2}{n_1}&0&0&0\\0&0&0&\frac{\sigma_\varepsilon^2+\sigma_\delta^2}{n_2}&\frac{\sigma_\delta^2}{n_2}&\frac{\sigma_\delta^2}{n_2}\\0&0&0&\frac{\sigma_\delta^2}{n_2}&\frac{\sigma_\varepsilon^2+\sigma_\delta^2}{n_2}&\frac{\sigma_\delta^2}{n_2}\\0&0&0&\frac{\sigma_\delta^2}{n_2}&\frac{\sigma_\delta^2}{n_2}&\frac{\sigma_\varepsilon^2+\sigma_\delta^2}{n_2}\end{bmatrix}\right)
(#eq:29-9)
\end{equation}

注意到对比

$$\begin{aligned}

\begin{aligned}\mu_{11}-\frac12\mu_{12}-\frac12\mu_{13}-\mu_{21}+\frac12\mu_{22}+\frac12\mu_{23}\end{aligned}& =(\mu+S_{ABA}+P_1-T_A)-\frac12(\mu+S_{ABA}+P_2+T_B+\lambda_A)  \\
&-\frac12(\mu+S_{ABA}+P_3+T_A+\lambda_B)-(\mu+S_{BAB}+P_1+T_B) \\
&+\frac12(\mu+S_{BAB}+P_2+T_A+\lambda_B)+\frac12(\mu+S_{BAB}+P_3+T_B+\lambda_A) \\
&=T_{A}-T_{B}
\end{aligned}$$

因此，即使 $ABA$、$BAB$ 交叉设计中存在延滞，两种处理的直接差异也可以通过下式来估计

\begin{equation}
\hat{\mu}_{11}-\frac12\hat{\mu}_{12}-\frac12\hat{\mu}_{13}-\hat{\mu}_{21}+\frac12\hat{\mu}_{22}+\frac12\hat{\mu}_{23}
(#eq:29-10)
\end{equation}

还注意到

\begin{equation}
\mathrm{Var}\left(\hat{\mu}_{11}-\frac12\hat{\mu}_{12}-\frac12\hat{\mu}_{13}-\hat{\mu}_{21}+\frac12\hat{\mu}_{22}+\frac12\hat{\mu}_{23}\right)=\frac32\sigma_{\varepsilon}^2\left(\frac1{n_1}+\frac1{n_2}\right)
(#eq:29-11)
\end{equation}

因此，式 \@ref(eq:29-11) 中的对比是个体内对比，其方差仅取决于 $\sigma^2_\varepsilon$.

一个可能感兴趣的问题是如何定义处理主效应均值？一种可能性是通过对表 \@ref(tab:table29-19) 中接受处理 $A$ 的所有单元进行平均来估计 $A$ 的主效应均值。也就是说，$A$ 平均值为 $(\mu_{11} + \mu_{22} + \mu_{13})/3$，$B$ 平均值为 $(\mu_{21} + \mu_{12} + \mu_{23})/3$。根据表 \@ref(tab:table29-19) 中定义的效应模型参数，这两个函数分别等于

$$\begin{aligned}\mu+\frac23S_{ABA}+\frac13S_{BAB}+\bar{P}_\cdot+T_A+\frac13(\lambda_B+\lambda_A)\\\mu+\frac13S_{ABA}+\frac23S_{BAB}+\bar{P}_\cdot+T_B+\frac13(\lambda_A+\lambda_B)\end{aligned}$$

这样的定义没有意义，因为这样的 $A$ 均值和 $B$ 均值的差等于 $\frac{1}{3}S_{ABA}-\frac{1}{3}S_{BAB}+T_A-T_B$，而不等于 $T_A-T_B$. 也就是说，这样定义的 $A$ 均值和 $B$ 均值与序列效应混叠 (aliased). $A$ 主效应均值的第二个定义是

$$\mu+\bar{S}_\cdot+\bar{P}_\cdot+T_A+\frac13(\lambda_B+\lambda_A)$$

$B$ 主效应均值类似地定义为

$$\mu+\bar{S}_\cdot+\bar{P}_\cdot+T_B+\frac13(\lambda_B+\lambda_A)$$

这样的定义是合理的，因为第一个等于表 \@ref(tab:table29-19) 中单元格均值的下列线性函数

$$\frac23\mu_{11}-\frac1{12}\mu_{12}-\frac1{12}\mu_{13}-\frac13\mu_{21}+\frac5{12}\mu_{22}+\frac5{12}\mu_{23}$$

第二个等于

$$-\frac13\mu_{11}+\frac5{12}\mu_{12}+\frac5{12}\mu_{13}+\frac23\mu_{21}-\frac1{12}\mu_{22}-\frac1{12}\mu_{23}$$

此外，这两个均值之差为

$$\mu_{11}-\frac12\mu_{12}-\frac12\mu_{13}-\mu_{21}+\frac12\mu_{22}+\frac12\mu_{23}=T_A-T_B$$

请注意，这与公式 \@ref(eq:29-10) 中给出的单元格均值对比相同。单元格均值内的对比意味着延滞效应是

\begin{equation}
\mu_{11}-\mu_{13}-\mu_{21}+\mu_{23}=\lambda_A-\lambda_B
(#eq:29-12)
\end{equation}

该对比估计的方差是

\begin{equation}
\mathrm{Var}(\hat{\mu}_{11}-\hat{\mu}_{13}-\hat{\mu}_{21}+\hat{\mu}_{23})=2\sigma_{\varepsilon}^2\left(\frac1{n_1}+\frac1{n_2}\right)
(#eq:29-13)
\end{equation}

且此对比也是个体内对比。

为了举例说明，请考虑表 \@ref(tab:table29-20) 中的数据，其中三名个体分配给 $ABA$ 序列，另外三名个体分配给 $BAB$ 序列。因此 $n_1 = n_2 = 3$. 读取数据并计算序列×时期单元格均值的 SAS 代码如表 \@ref(tab:table29-21) 所示。单元平均值如表 \@ref(tab:table29-22) 所示。

因此，$A$ 主效应平均值的估计由下式给出

$$\begin{aligned}
\frac23\hat{\mu}_{11}& -\frac1{12}\hat{\mu}_{12}-\frac1{12}\hat{\mu}_{13}-\frac13\hat{\mu}_{21}+\frac5{12}\hat{\mu}_{22}+\frac5{12}\hat{\mu}_{23}  \\
&=\frac23(24.133)-\frac1{12}(26.533)-\frac1{12}(23.933)-\frac13(26.367)+\frac5{12}(26.233)+\frac5{12}(24.833)=24.372
\end{aligned}$$

$B$ 主效应平均值的估计由下式给出

$$\begin{aligned}-\frac13\hat{\mu}_{11}&+\frac5{12}\hat{\mu}_{12}+\frac5{12}\hat{\mu}_{15}+\frac23\hat{\mu}_{21}-\frac1{12}\hat{\mu}_{22}-\frac1{12}\hat{\mu}_{23}\\&=-\frac13(24.133)+\frac5{12}(26.533)+\frac5{12}(23.933)+\frac23(26.367)-\frac1{12}(26.233)-\frac1{12}(24.833)=26.306\end{aligned}$$

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-20) $ABA/BAB$ 交叉实验的数据</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> Seq </th>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> Per </th>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> Trt </th>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> Person </th>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> Y </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> ABA </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 25.1 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> ABA </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 27.6 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> ABA </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 24.5 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> ABA </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 22 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> ABA </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 24.3 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> ABA </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 21.6 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> ABA </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 25.3 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> ABA </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 27.7 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> ABA </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 25.7 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BAB </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;"> 25.5 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BAB </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;"> 23.7 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BAB </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;"> 24.9 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BAB </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> 27.4 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BAB </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> 27.9 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BAB </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> 24.6 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BAB </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 26.2 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BAB </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> A </td>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 27.1 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> BAB </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> B </td>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 25 </td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-21)用于获取表 \@ref(tab:table29-20) 中数据的序列×时期单元格均值的 SAS 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> TITLE 'CRSOVR EXAMPLE #29.2 AN ABA/BAB DESIGN';\
DATA CRS;\
 INPUT SEQ $ PER TRT $ PERSON Y;\
 LINES;\
ABA 1 A 1 25.1\
ABA 2 B 1 27.6\
ABA 3 A 1 24.5\
ABA 1 A 2 22.0\
ABA 2 B 2 24.3\
ABA 3 A 2 21.6\
ABA 1 A 3 25.3\
ABA 2 B 3 27.7\
ABA 3 A 3 25.7\
BAB 1 B 4 25.5\
BAB 2 A 4 23.7\
BAB 3 B 4 24.9\
BAB 1 B 5 27.4\
BAB 2 A 5 27.9\
BAB 3 B 5 24.6\
BAB 1 B 6 26.2\
BAB 2 A 6 27.1\
BAB 3 B 6 25.0\
;\
PROC MEANS;\
CLASS SEQ PER;\
VAR Y;\
RUN; </td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-22)序列×时期单元格均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.22.png">
</td>
  </tr>
</tbody>
</table>

$A$ 均值与 $B$ 均值之差为

$$\hat{\tau}_A-\hat{\tau}_B=24.372-26.306=-1.934$$

人们可以使用 SAS-Mixed 程序来估计处理均值、检验处理差异、检验延滞以及估计两处理/三时期交叉设计的方差分量。例如，为了分析表 \@ref(tab:table29-20) 中的数据，必须定义一个新变量来标识前一时期给予的处理。为了不消除第一时期数据，还必须为第一时期数据定义这个新变量。对于 $ABA$ 序列，可以将这个新变量定义为在时期 $1$、$2$ 和 $3$ 中分别取值为 $O$、$A$ 和 $B$. 对于 $BAB$ 序列，新变量将分别取值 $O$、$B$ 和 $A$. 表 \@ref(tab:table29-23) 中的 SAS 代码使用了一个名为 `PRIORTRT` 的新变量，该变量具有前几句中给出的属性。定义了这个新变量的数据如表 \@ref(tab:table29-24) 所示。可用于分析表 \@ref(tab:table29-20) 中数据的 SAS-Mixed 命令如表 \@ref(tab:table29-25) 所示。结果输出如表 \@ref(tab:table29-20) 至 \@ref(tab:table29-29) 所示。表 \@ref(tab:table29-26) 给出了方差分量的估计。从该表中可以看出，$\sigma^2_\delta=2.1378$ 和 $\sigma^2_\varepsilon=0.7872$. 表 \@ref(tab:table29-27) 给出了处理 $A$ 和处理 $B$ 之差的检验。表 \@ref(tab:table29-27) 中标记为 Trt 的行是检验 $H_{01}:\tau_A=\tau_B$，标记为 PRIORTRT 的行是检验 $H_{02}:\lambda_A=\lambda_B$. 请注意，表 \@ref(tab:table29-27) 中标记为 PER 的行只有一个自由度。这可能看起来很奇怪，因为有三个时期。然而，在时期 $1$ 和 $\lambda_0$ 之间存在完全混淆 (completely confounding)，因此标记为 PER 的行仅比较时期 $2$ 和时期 $3$. 也就是说，标记为 PER 的行检验 $H_{03}:P_2=P_3$. 表 \@ref(tab:table29-25) 中的 `CONTRAST` 选项给出了一个检验统计量，将所有三个时期相互比较；即该 `CONTRAST` 选项检验 $H_{03}:P_1+\lambda_0=P_2=P_3$ 的结果。该选项的结果附在表 \@ref(tab:table29-27) 的底部。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-23)为表 \@ref(tab:table29-20) 中的数据创建 `PRIORTRT` 变量的 SAS 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> DATA CRS2; SET CRS;\
IF SEQ='ABA' AND PER=2 THEN PRIORTRT='A';\
ELSE IF SEQ='ABA' AND PER=3 THEN PRIORTRT='B';\
ELSE IF SEQ='BAB' AND PER=2 THEN PRIORTRT='B';\
ELSE IF SEQ='BAB' AND PER=3 THEN PRIORTRT='A';\
ELSE PRIORTRT='O';\
RUN;\
PROC PRINT;\
TITLE2 'A PRINT OF THE DATA WITH PRIORTRT VARIABLE INCLUDED';\
RUN; </td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-24)打印定义了 `PRIORTRT` 变量的新数据</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.24.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-25)分析表 \@ref(tab:table29-24) 中数据的 SAS 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> PROC MIXED;\
 TITLE3 'AN ANALYSIS USING MIXED';\
 CLASSES SEQ PERSON TRT PER PRIORTRT;\
 MODEL Y=SEQ TRT PER PRIORTRT/DDFM=KR;\
 LSMEANS TRT/PDIFF;\
 ESTIMATE 'TRT DIFF' TRT 1 -1;\
 ESTIMATE 'CARRYOVER' PRIORTRT 1 -1;\
 CONTRAST 'PERIOD DIFF' PER 0 1 -1, PER 1 -.5 -.5 PRIORTRT -.5 -.5 1;\
 RANDOM PERSON(SEQ);\
RUN; </td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-26)方差分量估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.26.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-27)固定效应检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.27-1.png">
</td>
  </tr>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.27-2.png">
</td>
  </tr>
</tbody>
</table>

处理主效应均值的估计值如表 \@ref(tab:table29-28) 所示。请注意，SAS-Mixed 程序给出的值与表 \@ref(tab:table29-22) 中使用主效应平均值的第二定义计算的估计值相同。表 \@ref(tab:table29-25)  中两个 ESTIMATE 语句的结果如表 \@ref(tab:table29-29) 所示。这些估计的第一个 $\tau_A-\tau_B$ 给出了估计值的标准误估计，并给出了用于检验 $H_{01}:\tau_A=\tau_B$ 的 $t$ 统计量以及检验观察到的显著性水平。第二个 ESTIMATE 语句估计 $\lambda_A-\lambda_B$ 给出了估计值的标准误估计，并给出了用于检验 $H_{02}:\lambda_A=\lambda_B$ 的 $t$ 统计量以及检验观察到的显著性水平。读者应注意，表 \@ref(tab:table29-29) 中观察到的处理差异和延滞差异的显著性水平与表 \@ref(tab:table29-27) 中相应的观察到的显著性水平相同。处理差异估计的标准误估计可以通过将 \@ref(eq:29-11) 中的 $\sigma_\varepsilon^2$ 替换为 $\hat\sigma_\varepsilon^2$ 来计算，即

$$\begin{aligned}
\widehat{\mathrm{Var}}\left(\hat{\mu}_{11}-\frac12\hat{\mu}_{12}-\frac12\hat{\mu}_{13}-\hat{\mu}_{21}+\frac12\hat{\mu}_{22}+\frac12\hat{\mu}_{23}\right)& =\frac32\hat\sigma_\varepsilon^2\left(\frac1{n_1}+\frac1{n_2}\right)  \\
&=\left(\frac32\right)(0.7872)\left(\frac13+\frac13\right)=0.7872
\end{aligned}$$

因此 $\hat \tau_A-\hat \tau_B$ 标准误估计为 $√0.7872=0.8872$.

类似地，可以通过将方程 \@ref(eq:29-13) 中的 $\sigma_\varepsilon^2$ 替换为 $\hat\sigma_\varepsilon^2$ 并取结果的平方根，来获得 $\hat \lambda_A-\hat \lambda_B$ 的标准误估计。与这些的标准误估计中的每一个相关的自由度为 $2(n_1+n_2-2)=8$.

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-28)处理主效应均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.28.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-29)从 Estimate 选项中检验处理和延滞效应</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.29.png">
</td>
  </tr>
</tbody>
</table>

## 具有两种以上处理的交叉设计 {#sec29-4}

接下来，考虑六个序列中的三处理/三时期交叉设计。设计序列为 $ABC,ACB,BAC,BCA,CAB,CBA$. 所有这六个序列都是必需的，为了获得不会与延滞混叠的处理均值的直接比较。Williams (1949) 开发的设计在延滞效应方面是均衡的。Williams 交叉设计在延滞效应方面是均衡的，其特点是每次处理后都会进行完全相同次数的其他处理。

上面给出的六序列设计是 Williams 设计。在这六个序列的设计中，每个处理（字母）跟随另一个处理两次。还要注意，虽然具有序列 $ABC,BCA,CAB$ 的三序列设计是拉丁方设计，但它不是 Williams 设计。在这种设计中，$A$ 后面跟着 $B$ 两次，但 $A$ 后面永远不会跟着 $C$. 同样，$B$ 后面跟着 $C$ 两次，但 $B$ 后面从来没有跟着 $A$；$C$ 后面跟着 $A$ 两次，但 $C$ 后面从来没有跟着 $B$。三个时期内三种处理的 Williams 设计必须使用 $ABC,ACB,BAC,BCA,CAB,CBA$ 的所有六个序列。

假设有四处理和四时期，显然，24 序列设计具有处理 $A,B,C,D$ 的所有 $24$ 种可能的排序，这将是一个 Williams 设计。另一种四处理四时期的 Williams 设计的序列为 $ABCD,BDAC,CADB,DCBA$. 设 $t$ 为交叉设计中要使用的处理数和时期数。如果 $t$ 是偶数，则 Williams 设计可以用一个特殊的 $t×t$ 拉丁方设计来构建，如果 $t$ 是奇数，Williams 设计可以通过使用两个特殊的拉丁方设计来构建。

表 \@ref(tab:table29-30) 给出了三时期/三处理 Williams 设计的单元格均值参数。为了表明当存在延滞时可以估计处理差异，可以证明对于表 \@ref(tab:table29-30) 中的单元格均值

\begin{align}
&(5\mu_{11}-2\mu_{12}-3\mu_{13}+4\mu_{21}+2\mu_{22}-6\mu_{23}-5\mu_{31}+2\mu_{32}+3\mu_{33}-4\mu_{41}-2\mu_{41}\\&+6\mu_{43}-\mu_{51}+4\mu_{52}-3\mu_{53}+\mu_{61}-4\mu_{62}+3\mu_{63})/24=\tau_A-\tau_B
(#eq:29-14)
\end{align}

还注意到，也可以比较延滞效应，

\begin{align}
&(\mu_{11}+2\mu_{12}-3\mu_{13}+0\mu_{21}+2\mu_{22}-2\mu_{23}-\mu_{31}-2\mu_{32}+3\mu_{33}+0\mu_{41}-2\mu_{41}\\&+2\mu_{43}-\mu_{51}+0\mu_{52}+\mu_{53}+\mu_{61}+0\mu_{62}-\mu_{63})/8=\lambda_A-\lambda_B
(#eq:29-15)
\end{align}

还可以获得简化为 $\tau_A-\tau_C$ 和 $\tau_B-\tau_C$ 以及 $\lambda_A-\lambda_C$ 和 $\lambda_B-\lambda_C$ 的单元格均值参数的类似函数。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-30)六个序列的三时期/三处理交叉设计单元格均值参数</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.30.png">
</td>
  </tr>
</tbody>
</table>

例如，考虑表 \@ref(tab:table29-31) 中的数据。这些数据来自六个序列的三时期/三处理交叉设计。有少数个体缺失数据。使用表 \@ref(tab:table29-31) 中的命令创建 SAS 数据集。这些命令还创建了一个延滞参数，该参数给出了上一时期的处理值。对于第一时期，该参数的值为“O”。数据最初使用 SAS-Mixed 软件进行分析，使用表 \@ref(tab:table29-33) 中的命令。表 \@ref(tab:table29-34) 给出了个体间和个体内方差分量的估计，即，$\hat\sigma^2_\delta = 3.2278$ 和 $\hat\sigma^2_\varepsilon = 0.8934$.

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-31)三时期/三处理交叉设计的数据</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.31.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-32)分析表 \@ref(tab:table29-31) 中数据的 SAS 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.32.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-33)分析表 \@ref(tab:table29-31) 中数据的 SAS 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.33.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-34)方差分量的估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.34.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table29-35) 给出了固定效应的检验。对表 \@ref(tab:table29-35) 的检验表明，有一个很显著处理效应（$p < 0.0001$）和显著的延滞效应（$p = 0.0306$）。表 \@ref(tab:table29-36) 和 \@ref(tab:table29-37) 给予 `LSMEANS` 选项的结果。值得注意的是，无论是处理还是先前处理的最小二乘均值都不可估 (not estimable)。原因是 SAS-Mixed 程序生成的模型参数函数是不可估的。

例如，SAS-Mixed 程序将处理 $A$ 最小二乘均值定义为

$$\mu+\frac{S_{ABC}+S_{ACB}+S_{BAC}+S_{BCA}+S_{CAB}+S_{CBA}}6+\frac{P_1+P_2+P_3}3+\tau_A+\frac{\lambda_A+\lambda_B+\lambda_C+\lambda_O}4$$

该函数不可估，因为 $P_1$ 和 $\lambda_O$ 在模型中彼此完全混叠，因此，它们必须具有相同的乘数。在上式中，它们不具有相同的乘数，因为 $P_1$ 的乘数等于 $1/3$，而 $\lambda_O$ 的乘数等于 $1/4$. 处理 $A$ 最小二乘均值的更好定义为

\begin{align}
\mu&+\frac{S_{ABC}+S_{ACB}+S_{BAC}+S_{BCA}+S_{CAB}+S_{CBA}}6+\frac{P_1+P_2+P_3}3\\&+\tau_A+\frac{2\lambda_A+2\lambda_B+2\lambda_C+3\lambda_O}9
(#eq:29-16)
\end{align}

在该定义中，$P_1$ 和 $\lambda_O$ 的乘数均等于 $1/3$. 此外，由于四个延滞参数的乘数之和必须为 $1$，因此 $\lambda_A,\lambda_B,\lambda_C$ 的乘数必须等于 $2/9$. 式 \@ref(eq:29-16) 中的函数是可估计。以类似的方式，可以定义处理 $B$ 和处理 $C$ 的最小二乘均值，并且这些平均值之间的成对差异将不包括延滞参数。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-35)固定效应的检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.35.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-36)处理最小二乘均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.36.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-37)处理和延滞参数之间的成对比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.37.png">
</td>
  </tr>
</tbody>
</table>

尽管在表 \@ref(tab:table29-36) 中处理最小二乘均值不可估，但成对差异是可估计的，处理和延滞参数之间的成对比较见表 \@ref(tab:table29-37) 。从表 \@ref(tab:table29-37) 可以看出，

$$\hat{\tau}_{A}-\hat{\tau}_{B}=0.8195,\hat{\tau}_{A}-\hat{\tau}_{C}=-3.0861,\mathrm{and}\,\hat{\tau}_{B}-\hat{\tau}_{C}=-3.9056$$

相应的观察到的显著性水平为 $p = 0.0041$、$p < 0.0001$ 和 $p < 0.0001$.
因此，所有处理均值彼此之间存在显著差异。

从表 \@ref(tab:table29-37) 中还可以看出，

$$\hat{\lambda}_A-\hat{\lambda}_B=-0.7706,\hat{\lambda}_A-\hat{\lambda}_C=0.1842,\mathrm{~and~}\hat{\lambda}_B-\hat{\lambda}_C=0.9548$$

相应的观察到的显著性水平为 $p = 0.0557$、$p = 0.6380$ 和 $p = 0.0118$.
因此，处理 $A$ 和 $C$ 的延滞效应彼此相似，且均与处理 $B$ 的延滞显著不同。

表 \@ref(tab:table29-38) 给出了额外的 SAS 代码，该代码使用与 \@ref(eq:29-16) 中给出的处理 $A$ 相似的定义计算处理最小二乘均值。表 \@ref(tab:table29-38) 中 SAS 命令的结果见表 \@ref(tab:table29-39). 对表 \@ref(tab:table29-39) 的检查表明，处理 $C$ 的平均值最大，其次是处理 $A$ 和 $B$. 还注意到，表 \@ref(tab:table29-39) 中给出的处理平均值和延滞效应之间的成对比较与表 \@ref(tab:table29-37) 中的相应比较完全相同。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-38)计算处理可估最小二乘均值的额外 SAS 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.38.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table29-39)表 \@ref(tab:table29-38) 中给出的命令的结果</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2029.39.png">
</td>
  </tr>
</tbody>
</table>

## 小结 {#sec29-5}

本章考虑了交叉设计，即每个个体随着时间的推移接受一系列处理。特定个体接受的一系列特殊处理被称为序列。当使用交叉设计时，需要考虑延滞问题。如果怀疑可能存在延滞，可以在进行下一次处理之前纳入洗脱期，以尽可能减少延滞的影响。当两时期/两处理交叉设计中存在延滞时，不能使用第二个时期的数据来估计处理差异。然而，当能够使用两个以上时期和/或两种以上处理时，无论是否存在延滞，都可以获得处理差异的估计。
