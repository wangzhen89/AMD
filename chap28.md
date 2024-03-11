
# 案例研究：重复测量的复杂例子 {#chap28}

本章考虑了几个例子及其使用混合模型程序的统计分析。考虑的示例包括具有重复测量的裂区设计、重复测量中嵌套的重复测量以及多地点研究。

##  复杂舒适度实验 {#sec28-1}

工程师在三种环境 (environments) 中测试两种类型的服装 (clothing). 由于男性和女性对环境的反应也不同，因此人 (person) 的性别被列为一个因素。四个人（两男两女）被放入一个环境室 (chambers)（分配了三种环境之一）。一名男性和一名女性着类型 1 的服装，另一名男性和女性着类型 2 的服装。在 1、2 和 3 小时结束时记录每个人的舒适度得分。本实验的数据如表 \@ref(tab:table28-1) 所示。实验单位分为三种尺寸。最大的实验单位是一个室，或者相当于一个四人小组。室实验单元的实验设计是完全随机设计结构中的单向处理结构（环境），在环境的每一水平具有三次重复。中型实验单位是一个人。一个人的实验设计是在九个区组（每个区组包含四个实验单元，人）中的随机完全区组设计结构中的双向处理结构（性别 × 服装）。最小的实验单位是 1 小时的时间间隔 (Hour)，这是一个重复的测量。小时的实验设计是在 36 个区组（每个区组包含 3 个实验单元，1 小时的时间间隔）的随机完全区组设计结构中的单向处理结构 (Time).

<table>
<caption>(\#tab:table28-1)复杂舒适度实验数据</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.1.png">
</td>
  </tr>
</tbody>
</table>

该实验的模型为（其中模型被分为三个部分，对应于三种尺寸的实验单元，重复测量被视为<u>时间的裂区</u>因子）

\begin{align}
y_{ijkm\ell}=&\,[\mu+E_i+\eta_{i\ell\ell}]+[S_j+C_k+(SC)_{jk}+(ES)_{ij}+(EC)_{ik}+(ESC)_{ijk}+\delta_{ijk\ell}]+[T_m+(ET)_{im}\\&+(ST)_{jm}+(CT)_{km}+(SCT)_{jkm}+(EST)_{ijm}+(ECT)_{ikm}+(ESCT)_{ijkm}+\varepsilon_{ijkm\ell}]
(#eq:28-1)
\end{align}

其中 $E_i$ 表示第 i 个环境，$S_j$ 表示第 j 个性别，$C_k$ 表示第 k 种服装，$T_m$ 表示第 m 个时间点，$\eta_{i\ell}$ 表示随机室效应且假定 $\eta_{i\ell}\thicksim i.i.d.N(0,\sigma_{\eta}^{2})$，$\delta_{ijk\ell}$ 表示随机人效应且假定 ${\delta_{ijk\ell}}\sim N(0,{\sigma_\delta^2})$，$\varepsilon_{ijkm\ell}$ 表示给定小时的随机测量误差且假设 ${\varepsilon_{ijkm\ell}}\sim i.i.d.N(0,{\sigma_{\varepsilon}^2})$. 还假设所有随机效应独立分布。还要注意，式 \@ref(eq:28-1) 的第一个括号是模型的室部分，第二个括号是模型的人部分，而第三个括号是这个模型的小时部分。如果时间<u>时间的裂区</u>假设不合适，那么模型将是

\begin{align}
y_{ijkm\ell}=&\,\mu+E_i+\eta_{i\ell}+S_j+C_k+(SC)_{jk}+(ES)_{ij}+(EC)_{ik}+(ESC)_{ijk}+T_m+(ET)_{im}+(ST)_{jm}\\&+(CT)_{km}+(SCT)_{jkm}+(EST)_{ijm}+(ECT)_{ikm}+(ESCT)_{ijkm}+\varepsilon_{ijkm\ell}
(#eq:28-2)
\end{align}

其中假设 $\eta_{i\ell}\thicksim i.i.d.N(0,\sigma_\eta^2)$ 且

$$\left.{\boldsymbol\varepsilon}_{ijk\ell}=\begin{bmatrix}{\varepsilon}_{ijk1\ell}\\{\varepsilon}_{ijk2\ell}\\{\varepsilon}_{ijk3\ell}\end{bmatrix}\thicksim N_3\left(\begin{bmatrix}0\\0\\0\end{bmatrix}\right.,\begin{bmatrix}{\sigma}_{11}&{\sigma}_{12}&{\sigma}_{13}\\{\sigma}_{21}&{\sigma}_{22}&{\sigma}_{23}\\{\sigma}_{31}&{\sigma}_{32}&{\sigma}_{33}\end{bmatrix}\right)\thicksim N_3(0,{\boldsymbol\Sigma})\text{ (say)}$$

并且分析将取决于为 $\boldsymbol\Sigma$ 假设的结构。注意，式 \@ref(eq:28-2) 中的模型与 \@ref(eq:28-1) 中的模型相同，只是从模型中删除了对应于人的随机效应 $\delta_{ijk\ell}$。对于某些 $\boldsymbol\Sigma$ 有时可能希望包括人随机效应，例如当 $\boldsymbol\Sigma$ 具有 AR(1) 结构时。

如果希望执行<u>时间的裂区</u>分析，可使用的基本 SAS^®^-Mixed 代码如表 \@ref(tab:table28-2) 所示。这样的分析结果就留给读者去追寻了。这样的分析假设实验设计是裂-裂区设计，以室为整区实验单位，以人为子区实验单位，以时间测量为次子区实验单位。

<table>
<caption>(\#tab:table28-2)用于分析表 \@ref(tab:table28-1) 中数据的 SAS Mixed 代码，将重复测量因子视为<u>时间的裂区</u>因子</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.2.png">
</td>
  </tr>
</tbody>
</table>

这里将采取的方法是将时间因子视为重复测量，然后考虑重复测量的各种协方差结构。然后，该设计将被描述为对子区因子（人）进行重复测量的裂区设计。表 \@ref(tab:table28-3) 显示了在假设重复测量的协方差结构具有复合对称性时使用的基本 SAS-Mixed 代码。请注意，Sex × Clo × Rep × Env 的组合对于每个人来说都是唯一的，因此此效应可用于识别出现重复测量的人。

<table>
<caption>(\#tab:table28-3)使用具有复合对称协方差结构的 `repeated`` 选项分析表 \@ref(tab:table28-1) 中数据的 SAS Mixed 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.3.png">
</td>
  </tr>
</tbody>
</table>

对于这些数据，考虑了重复测量的协方差阵的四种结构。所考虑的协方差结构是：复合对称，AR(1)，一个非结构化协方差阵以及人也被包括为随机效应的 AR(1) 重复测量协方差结构。对于其中的前三个，只需要更改表 \@ref(tab:table28-3) 中 `repeated`选项中的 "Type = "语句。为了将人也包括为随机因素，可以使用下面给出的 random 和 `repeated` 选项


```sas
RANDOM REP(ENV) SEX*CLO*REP*ENV; 
REPEATED TIME/SUBJECT=SEX*CLO*REP*ENV TYPE=AR(1);
```

表 \@ref(tab:table28-4) 中给出了所考虑的四个协方差结构中每一个的最终拟合统计量。

<table>
<caption>(\#tab:table28-4)四种协方差结构的拟合统计量</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.4.png">
</td>
  </tr>
</tbody>
</table>

对表 \@ref(tab:table28-4) 中的 AIC, AICC 和 BIC 拟合统计量的检查表明，对于非结构化协方差矩阵假设，每一个都有其最小值。将复合对称与非结构化进行比较的卡方统计量为 $\chi^2=117.7-107.5=10.2$，自由度为 6-2=4。相应的观察到的显著性水平为 $\hat\alpha=0.0372$. 因此，非结构化明显优于复合对称。该示例的其余分析将假设重复测量的非结构化协方差阵。表\@ref(tab:table28-5) 中给出了产生固定效应分析表的 SAS-Mixed 代码。表 \@ref(tab:table28-5) 中所述分析的协方差参数估计如表 \@ref(tab:table28-6) 所示。从表 \@ref(tab:table28-6) 中可以看出 $\hat{{\sigma}}_{{\eta}}^{2}=2.4047$，并且

$$\boldsymbol{\Sigma}=\begin{bmatrix}0.2157&0.2256&0.1248\\0.2256&0.3525&0.1894\\0.1248&0.1894&0.1462\end{bmatrix}$$

<table>
<caption>(\#tab:table28-5)使用非结构化协方差结构的 `repeated`` 选项分析表 \@ref(tab:table28-1) 中数据的 SAS Mixed 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.5.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table28-6)表 \@ref(tab:table28-5) 中分析的协方差参数估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.6.png">
</td>
  </tr>
</tbody>
</table>

固定效应分析表见表 \@ref(tab:table28-7)。对表 \@ref(tab:table28-7) 的检查显示，显著的效应（$\hat\alpha<0.05$）为：Sex, Clo, Env × Clo, Sex × Clo, Env × Sex × Clo, Time 和 Env × Time. 由于三向交互 Env × Sex × Clo 是显著的，我们应该研究这些三向均值之间有趣的成对比较。唯一涉及时间的显著性交互是环境与时间的交互，因此我们也应该研究这些双向均值之间的成对比较。除了上面描述的成对比较，分析还考虑了每种环境下时间均值之间的线性和二次比较。这里这样做是为了说明，读者需要自己确定这种比较是否对类似于本研究的研究感兴趣。需要附加到表 \@ref(tab:table28-5) 中的额外 SAS-Mixed 语句如表 \@ref(tab:table28-8) 所示。这些附加语句的结果如表 \@ref(tab:table28-9) - \@ref(tab:table28-16) 所示。

<table>
<caption>(\#tab:table28-7)固定效应分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.7.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table28-8)用于获得感兴趣比较的额外 SAS-Mixed 语句</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.8.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table28-9) 包含了每种环境下 Sex × Clothing 的四种组合和 Reps 的时间均值的线性和二次对比。对表 \@ref(tab:table28-9) 的检查表明，关于环境 1 和 2，线性和二次对比都具有统计学意义，环境 3 都不具有统计学意义。

<table>
<caption>(\#tab:table28-9)每种环境下时间的线性和二次对比</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.9.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table28-10) 包含 Environment × Sex × Clothing 的最小二乘均值。这些均值是在 Reps 和 Time 上取均值。表 \@ref(tab:table28-11) 包含 Environment × Time 最小二乘均值，这些均值是在 Reps, Clothing 和 Sex 上取均值。有 12 个 Environment × Sex × Clothing 均值。如果要在这 12 个均值之间进行所有成对比较，则需要 66 个比较。然而，许多这样的成对比较对实验者来说并不有趣。通常情况下，我们感兴趣的那些比较往往是针对每两种因素的组合，对比另一个因素的各种可能性。也就是说，人们通常想要比较每个 Environment × Sex 组合的 Clothing 类型（六个比较）、每个 Environment × Clothing 组合的 Sex（六个对比）以及每个 Sex × Clothing 组合所处的 Environment（12个比较）。因此，在 66 个可能的比较中，总共有 24 个可能是实验者感兴趣的。这些比较的结果分别显示在表 \@ref(tab:table28-12) - \@ref(tab:table28-14) 中。

<table>
<caption>(\#tab:table28-10)Environment × Sex × Clothing 三向最小二乘均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.10.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table28-11)Environment 和 Time 的双向最小二乘均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.11.png">
</td>
  </tr>
</tbody>
</table>

通过查看表 \@ref(tab:table28-12) 中的比较可以发现，在不同 Environment × Sex 的每一种组合下，Clothing 类型都有显著差异；进一步分析表 \@ref(tab:table28-13) 表明，在每一种 Environment × Clothing 类型的组合中，男性和女性之间存在显著差异。最后，对表 \@ref(tab:table28-14) 的考察揭示了只有在 Sex = 1 且 Clo = 2 的情况下，Environment 1 和 2 才有显著差异；Sex = 2 且 Clo = 1 的情况下也是如此。

<table>
<caption>(\#tab:table28-12)不同 Environment × Sex 组合下 Clothing 类型的比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.12.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table28-13)不同 Environment × Clothing 类型组合下 Sex 的比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.13.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table28-14)不同 Sex × Clothing 类型组合下 Environment 的比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.14.png">
</td>
  </tr>
</tbody>
</table>

如果要对 Environment × Time 这九种组合进行所有可能的成对比较，总共需要进行 36 次比较。但其中可能只有 18 次比较对实验者有意义。表 \@ref(tab:table28-15) 包含了在每个时间段内（这些均值是对 reps, Sex 和 clothing 进行平均后的结果）各 Environment 间的均值比较，而表 \@ref(tab:table28-16) 则包含了在每个环境下（这些均值同样是对 reps, Sex 和 clothing 进行平均后的结果）各个时间段的均值比较。通过观察表 \@ref(tab:table28-15) 可知，在 time 1 时，Environment 1 与 2, 3 都有显著差异，而 Environment 2 与 3 之间没有显著差异。在 time 2 和 3 时，这三个 Environment 彼此间均无显著差异。进一步分析表 \@ref(tab:table28-16) 可知，对于 Environment 1和 2 来说，所有 time 之间的差异都是显著的；而对于 Environment 3 而言，其在任何 time 之间都没有显著差异。

<table>
<caption>(\#tab:table28-15)每个时间点环境间的比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.15.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table28-16)每环境下时间之间的比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.16.png">
</td>
  </tr>
</tbody>
</table>

显示每次比较的 time 均值的图如图 \@ref(fig:figure28-1) 所示。该图是由表 \@ref(tab:table28-17) 中所示的 SAS 命令创建的。对图的检查说明了有关表 \@ref(tab:table28-16) 的 time 比较的陈述。

<div class="figure" style="text-align: center">
<img src="figure/figure 28.1.png" alt="每种环境下响应随时间的变化" width="622" />
<p class="caption">(\#fig:figure28-1)每种环境下响应随时间的变化</p>
</div>

<table>
<caption>(\#tab:table28-17)获得每种环境下时间比较的图形的额外 SAS 语句</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.17.png">
</td>
  </tr>
</tbody>
</table>

## 家庭态度实验 {#sec28-2}

农村和城市环境中家庭的态度在三个时间段内每六个月测量一次。这些数据来自七个农村家庭和十个城市家庭，每个家庭都由一个儿子、一个父亲和一个母亲组成。数据在表 \@ref(tab:table28-18) 中给出。注意，在这个例子中，有两组重复测量数据，即三个时间段和三个家庭成员。家庭成员是一个重复测量，因为不能随机地将父亲、母亲和儿子分配给三个家庭成员。这个实验可以描述为一个完全随机设计的单向实验，包含两组重复测量数据，即家庭成员和嵌套在家庭成员中的时间。

<table>
<caption>(\#tab:table28-18)家庭态度研究数据</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.18.png">
</td>
  </tr>
</tbody>
</table>

本实验的处理结构是三向处理结构，因素为地区（A：农村与城市）、家庭成员（M：儿子、父亲、母亲）和时间（T：1, 2, 3）。设计结构是完全随机设计，有两组重复测量。可代表该实验的一个模型由下式给出

\begin{align}
y_{ijk\ell}&=\mu+A_i+\eta_{i\ell}+M_j+(AM)_{ij}+\delta_{ij\ell}+T_k+(AT)_{ik}+(MT)_{jk}+(AMT)_{ijk}+\varepsilon_{ijk\ell}\\&i=1,2;~j=1,2,3;~k=1,2,3
(#eq:28-3)
\end{align}

在两组重复测量的<u>时间的裂区</u>假设下，令 $\eta_{i\ell}$ 表示随机家庭效应且假设 $\eta_{i\ell}\thicksim i.i.d.~N(0,~\sigma_\eta^2)$，${\delta_{ijk\ell}}$ 表示随机人效应且假设 ${\delta_{ijk\ell}}\sim{N}(0,{\sigma_\delta^2})$，${\varepsilon}_{ijk\ell}$ 表示给定时间段的测量误差且假设 ${\varepsilon_{ijk\ell}}\sim i.i.d.N(0,{\sigma_{\varepsilon}^{2}})$. 在上述假设下执行分析的 SAS-Mixed 命令如表 \@ref(tab:table28-19) 所示。然而，如果读者愿意，则由他们来执行此分析。使用表 \@ref(tab:table28-19) 中的命令进行的分析与裂区实验的分析类似。参见示例 [24.5](#sec24-6-3)。

<table>
<caption>(\#tab:table28-19)用于分析表 \@ref(tab:table28-18) 中日期的 SAS 命令，假设两组重复测量都满足<u>时间的裂区</u>假设</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.19.png">
</td>
  </tr>
</tbody>
</table>

本节中考虑的分析将考虑模型

\begin{align}
y_{ijk\ell}&=\mu+A_i+M_j+(AM)_{ij}+T_k+(AT)_{ik}+(MT)_{jk}+(AMT)_{ijk}+\varepsilon_{ijk\ell}\\&i=1,2;~j=1,2,3;~k=1,2,3
(#eq:28-4)
\end{align}

令

$$
\boldsymbol \varepsilon_{i\ell}=\begin{bmatrix}\varepsilon_{i11\ell}\\
\varepsilon_{i12\ell}\\
\varepsilon_{i13\ell}\\
\varepsilon_{i21\ell}\\
\varepsilon_{i22\ell}\\
\varepsilon_{i23\ell}\\
\varepsilon_{i31\ell}\\
\varepsilon_{i32\ell}\\
\varepsilon_{i33\ell}\end{bmatrix}
$$

为第 i 个地区第 $\ell$ 个家庭的误差向量。

假设 $\boldsymbol \varepsilon_{i\ell}\sim N(\boldsymbol 0,\boldsymbol \Sigma)$，统计分析将依赖于 $\boldsymbol \Sigma$ 的结构。协方差阵 $\boldsymbol \Sigma$ 为 9 × 9 的矩阵。这种设置忽略了重复测量的嵌套结构，并简单地假设有九个重复测量，并且第 \@ref(chap27) 章中考虑的所有结构在这里都是可能的。考虑到重复测量的嵌套的结构作了如下假设

$$
\boldsymbol \varepsilon_{i\ell}=\operatorname{Var}\left(\begin{bmatrix}\varepsilon_{i11\ell}\\
\varepsilon_{i12\ell}\\
\varepsilon_{i13\ell}\\
\varepsilon_{i21\ell}\\
\varepsilon_{i22\ell}\\
\varepsilon_{i23\ell}\\
\varepsilon_{i31\ell}\\
\varepsilon_{i32\ell}\\
\varepsilon_{i33\ell}\end{bmatrix}\right)=\begin{bmatrix}{\theta}_{11}\boldsymbol{V}&{\theta}_{12}\boldsymbol{V}&{\theta}_{13}\boldsymbol{V}\\{\theta}_{21}\boldsymbol{V}&{\theta}_{22}\boldsymbol{V}&{\theta}_{23}\boldsymbol V\\{\theta}_{31}\boldsymbol{V}&{\theta}_{32}\boldsymbol{V}&{\theta}_{33}\boldsymbol{V}\end{bmatrix}=\boldsymbol{\Theta}\otimes\boldsymbol{V}\mathrm{~(say)}
$$

其中 $\boldsymbol{\Theta}\otimes\boldsymbol{V}$ 表示如下两个方差协方差阵的直积 (direct product)：

$$\boldsymbol{\Theta}=\begin{bmatrix}\theta_{11}&\theta_{12}&\theta_{13}\\\theta_{21}&\theta_{22}&\theta_{23}\\\theta_{31}&\theta_{32}&\theta_{33}\end{bmatrix}\quad\mathrm{and}\quad\boldsymbol{V}=\begin{bmatrix}v_{11}&v_{12}&v_{13}\\v_{21}&v_{22}&v_{23}\\v_{31}&v_{32}&v_{33}\end{bmatrix}$$

当 $\boldsymbol{\Theta}$ 是非结构化且 $\boldsymbol{V}$ 具有复合对称性、AR(1)  或非结构化时，SAS-Mixed 程序将允许用户选择此类协方差结构。将数据如同裂区实验一样分析相当于假设 $\boldsymbol{\Theta}$ 和 $\boldsymbol{V}$ 都满足复合对称结构。

表 \@ref(tab:table28-18) 中的数据在四种不同的协方差结构下进行了分析，前三个是通过取对应于家庭成员的协方差阵和对应于时间的协方差阵的直积获得的。家庭成员考虑的结构是非结构化的，时间考虑的结构是复合对称、AR(1) 和非结构化。考虑的第四个协方差结构是非结构化 9 × 9 协方差阵，它将家庭成员和时间的所有组合视为九次重复测量。检查所得的拟合统计量，以确定在比较地区、家庭成员和时间的固定效应以及这些效应的所有可能的交互作用时要考虑的最佳协方差结构。表 \@ref(tab:table28-20) 中给出了所使用的 SAS-Mixed 命令，表 \@ref(tab:table28-21) 为续表。可能会注意到，所有分析中都使用了选项 `DDFM = BETWITHIN`，而不是 `DDFM = KR`。对于获得的拟合统计量，使用这两个选项中的哪一个并不重要。然而，如果要考虑对直积协方差结构下的固定效应进行检验，则使用 `DDFM = KR` 无法给出合理的结果，因为 9 × 9 协方差阵估计在直积结构下是奇异的。表 \@ref(tab:table28-22) 给出了拟合统计结果。在假设 $\boldsymbol{V}$ 具有复合对称性或 AR(1) 的两种情况下，Mixed 程序不会产生任何结果，因为获得了无限的似然值。发生这种情况的原因似乎是 Mixed 程序选择的协方差参数的默认起始值。包含 `PARMS` 选项以强制 Mixed 程序以不同的起始值开始。语句 `PARMS 10, 1, 10, 1, 1, 10, .5;` 为下列矩阵

$${\boldsymbol\Theta}=\begin{bmatrix}{\theta}_{11}&{\theta}_{12}&{\theta}_{13}\\{\theta}_{21}&{\theta}_{22}&{\theta}_{23}\\{\theta}_{31}&{\theta}_{32}&{\theta}_{33}\end{bmatrix}\quad\mathrm{and}\quad{\boldsymbol V}=\begin{bmatrix}1&{\rho}&{\rho}\\{\rho}&1&{\rho}\\{\rho}&{\rho}&1\end{bmatrix}$$

选择参数的起始值：$\theta_{11}=10,\theta_{21}=1,\theta_{22}=10,\theta_{31}=1,\theta_{32}=1,\theta_{33}=10,\rho=0.5$. 不失一般性，$\boldsymbol{V}$ 中的对角线元可以视为等于 1。相同的起始值也用于 $\boldsymbol{V}$ 具有 AR(1) 结构时的情况。在这种情况下

$$\boldsymbol{V}=\begin{bmatrix}1&{\rho}&{\rho}^2\\{\rho}&1&{\rho}\\{\rho}^2&{\rho}&1\end{bmatrix}$$

同样不失一般性，$\boldsymbol{V}$ 中的对角线元素可以视为等于 1.

<table>
<caption>(\#tab:table28-20)用于比较协方差结构的 SAS 命令</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.20.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table28-21)用于比较协方差结构的 SAS 命令（续）</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.21.png">
</td>
  </tr>
</tbody>
</table>

检查表 \@ref(tab:table28-22) 中的 AIC, AICC 和 BIC 值发现，当 $\boldsymbol{\Theta}$ 是非结构化的，$\boldsymbol{V}$ 具有复合对称性时，协方差结构 $\boldsymbol{\Theta}\otimes\boldsymbol{V}$ 这三个准则在下均达到了各自的最小值。这就是用于考虑固定效应参数推断的假设。

<table>
<caption>(\#tab:table28-22)家庭态度数据的拟合统计量</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.22.png">
</td>
  </tr>
</tbody>
</table>

为了获得固定效应的检验，表 \@ref(tab:table28-23) 中给出了 SAS-Mixed 命令。由于 SAS-Mixed 程序按字母顺序对 Classes 语句中的变量进行排序，因此在此分析之前对数据进行了排序。对数据进行排序有助于将协方差参数估计与适当的家庭成员相匹配。协方差参数的估计值如表 \@ref(tab:table28-24) 所示。注意到

$$\hat{{\boldsymbol\Theta}}=\begin{bmatrix}\hat{{\theta}}_{_{FF}}&\hat{{\theta}}_{_{FM}}&\hat{{\theta}}_{_{FS}}\\\hat{{\theta}}_{_{MF}}&\hat{{\theta}}_{_{MM}}&\hat{{\theta}}_{_{MS}}\\\hat{{\theta}}_{_{SF}}&\hat{{\theta}}_{_{SM}}&\hat{{\theta}}_{_{SS}}\end{bmatrix}=\begin{bmatrix}8.1570&2.7308&1.2415\\2.7508&9.7410&2.6035\\1.2415&2.6035&14.7547\end{bmatrix}\quad\mathrm{and}\quad\hat{{\rho}}=0.9637$$

$\hat{{\boldsymbol\Theta}}$ 中的下标已更改为对应于父亲 (F)、母亲 (M) 和儿子 (S)。表 \@ref(tab:table28-25) 中给出了 9 × 9 协方差矩阵 $\boldsymbol{\Theta}\otimes\boldsymbol{V}$ 的估计结果，表 \@ref(tab:table28-26) 中展示了对固定效应的检验。

<table>
<caption>(\#tab:table28-23)用于获得固定效应检验表和协方差参数估计的 SAS 命令</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.23.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table28-24)协方差参数估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.24.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table28-25)每个家庭下的协方差阵估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.25.png">
</td>
  </tr>
</tbody>
</table>

对表 \@ref(tab:table28-26) 的检查表明，A, M, T 和 A × T 具有统计显著性。为了探究这些效应，我们应该看看家庭成员 (M) 的主效应均值和地区与时间的双向均值 (A × T).

<table>
<caption>(\#tab:table28-26)固定效应检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.26.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table28-27) 所展示了用于计算相应的最小二乘均值并在它们之间进行成对比较的 SAS-Mixed 命令。这些命令应附加到表 \@ref(tab:table28-23) 中的命令中。 M 和 A × T 的最小二乘均值分别显示在表 \@ref(tab:table28-28) 和 \@ref(tab:table28-29) 中。表 \@ref(tab:table28-30) - \@ref(tab:table28-32) 给出了这些均值之间有趣的成对比较。对表 \@ref(tab:table28-28) 和 \@ref(tab:table28-30) 的检查表明，儿子的均值显著小于父母的均值，并且父亲和母亲的均值彼此没有显著差异。

<table>
<caption>(\#tab:table28-27)用于获得最小二乘均值及其之间的成对比较的额外 SAS 命令</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.27.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table28-28)家庭成员最小二乘均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.28.png">
</td>
  </tr>
</tbody>
</table>

对表 \@ref(tab:table28-29) 的检查表明，农村家庭和城市家庭的均值随着时间的推移而增加。然而，对表 \@ref(tab:table28-31) 的检查显示，农村家庭的 0 个月和 6 个月均值没有显著差异，但均显著小于 12 个月均值。对于城市家庭来说，所有时间均值彼此均有显著差异。对表 \@ref(tab:table28-32) 的检验表明，城市家庭在 0 个月和 12 个月时的均值显著小于农村家庭。农村家庭与城市家庭 6 个月均值的差异不足以达到 0.05 水平的显著性。

<table>
<caption>(\#tab:table28-29)地区和时间组合的双向最小二乘均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.29.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table28-30)家庭成员比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.30.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table28-31)每个地区内的时间比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.31.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table28-32)每个时间下的地区比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.32.png">
</td>
  </tr>
</tbody>
</table>

## 多地点研究 {#sec28-3}

本节考虑涉及三种药物的实验，其中每个个体在三个不同的时间点重复测量。数据由三名不同的研究人员（或在三个不同的中心）收集。数据如表 \@ref(tab:table28-33) 所示。请注意，有很多缺失的数据点，如空单元格所示。

<table>
<caption>(\#tab:table28-33)多地点药物实验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.33.png">
</td>
  </tr>
</tbody>
</table>

该分析假设中心 (center) 和个体是随机效应，药物和时间是固定效应。如果假设重复测量满足<u>时间的裂区</u>假设，则表 \@ref(tab:table28-33) 中数据的模型为

\begin{equation}
y_{ijk\ell}=\mu+\eta_i+D_j+\delta_{ij\ell}+T_k+(DT)_{jk}+\varepsilon_{ijk\ell}
(#eq:28-5)
\end{equation}

其中 $D_j$ 表示第 j 种药物，$T_k$ 表示地 k 个时间效应，$(DT)_{ij}$ 表示时间与药物的交互作用。此外，$\eta_i$ 表示中心随机效应且假设 $\eta_i\thicksim i.i.d.N(0,\sigma_\eta^2)$，$\delta_{ij\ell}$ 表示第 i 个中心的第 j 种药物的第 $\ell$ 个个体的随机个体效应，且假设 ${\delta}_{ij\ell}\thicksim N(0,{\sigma}_\delta^2)$. 最后，$\varepsilon_{ijk\ell}$ 表示与第 i 个中心的第 j 种药物的第 $\ell$ 个个体第 i 个中心的第 j 种药物的第 $\ell$ 个个体的第 $\ell$ 个时间相关的随机误差，且假设 $\varepsilon_{ijk\ell}\sim i.i.d.N(0,\sigma_{\varepsilon}^2)$. 还假设所有随机效应都是独立分布的。

如果<u>时间的裂区</u>假设不合适，则模型将为

\begin{equation}
y_{ijk\ell}=\mu+\eta_i+D_j+T_k+(DT)_{jk}+\varepsilon_{ijk\ell}
(#eq:28-6)
\end{equation}

其中假设 $\eta_i\thicksim i.i.d.N(0,\sigma_\eta^2)$ 以及

$$\left.{\boldsymbol\varepsilon}_{ij\ell}=\begin{bmatrix}{\varepsilon}_{ij1\ell}\\{\varepsilon}_{ij2\ell}\\{\varepsilon}_{ij3\ell}\end{bmatrix}\thicksim N_3\left(\begin{bmatrix}0\\0\\0\end{bmatrix}\right.,\begin{bmatrix}{\sigma}_{11}&{\sigma}_{12}&{\sigma}_{13}\\{\sigma}_{21}&{\sigma}_{22}&{\sigma}_{23}\\{\sigma}_{31}&{\sigma}_{32}&{\sigma}_{33}\end{bmatrix}\right)\thicksim N_3(\boldsymbol{0},{\boldsymbol\Sigma})\text{ (say)}$$

分析将取决于 $\boldsymbol\Sigma}$ 的假设结构。请注意，式 \@ref(eq:28-6) 中的模型与 \@ref(eq:28-5) 中的模型相同，只是与个体对应的随机效应 $\delta_{ij\ell}$ 已从模型中删除。人们可能还希望将人作为某些 $\boldsymbol\Sigma}$ 的随机效应，例如当 $\boldsymbol\Sigma}$ 具有 AR(1) 结构时。

考虑了重复测量的四种协方差结构：非结构化、复合对称、AR(1) 和 AR(1) 加随机个体效应。对于非结构化假设，BIC 准则最小，而对于 AR(1) 加随机个体效应结构，AIC 和 AICC 都是最小的。对于后一种结构，$\rho$ 的估计为 –0.3141，从哲学角度来看这似乎没有意义，因为这意味着第一个和第二个时间点之间的相关性为负，但第一和第三个时间点间的相关性是正的，这似乎不合理。因此，通过假设重复测量具有非结构化协方差结构，可以获得对固定效应的推断。

表 \@ref(tab:table28-34) 中给出了用于获得固定效应推断的 SAS-Mixed 命令。协方差参数的估计如表 \@ref(tab:table28-35) 所示，混合效应的假设检验如表 \@ref(tab:table28-36) 所示。对表 \@ref(tab:table28-36) 的检查表明，药物和时间的主效应都是显著的，而药物与时间交互效应不显著。因此，可以在药物和时间的主效应均值之间进行比较。表 \@ref(tab:table28-37) 给出了药物主效应均值，表 \@ref(tab:table28-38) 给出了药物主效应均值之间的成对比较。对这两个表的检查表明，药物 1 的均值显著小于药物 2 的均值和药物 3 的均值。此外，药物 2 和药物 3 的均值之间没有显著差异。表 \@ref(tab:table28-39) 给出了时间主效应均值，表 \@ref(tab:table28-40) 给出了时间均值之间的成对比较。时间均值着随着时间的推移而增加，所有成对比较都具有统计学意义。

<table>
<caption>(\#tab:table28-34)用于分析表 \@ref(tab:table28-33) 中数据的 SAS 命令</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.34.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table28-35)协方差参数估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.35.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table28-36)固定效应假设检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.36.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table28-37)药物最小二乘均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.37.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table28-38)药物主效应均值的成对比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.38.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table28-39)时间最小二乘均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.39.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table28-40)时间主效应均值的成对比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2028.40.png">
</td>
  </tr>
</tbody>
</table>

## 练习 {#sec28-4}

