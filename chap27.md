
# 不满足理想条件时重复测量实验的分析 {#chap27}

> "He uses statistics as a drunken man uses lamp-posts--for support rather than illumination." - Andrew Lang

[^subject]: 本书将 "subject" 译作个体，原因见介绍页的[术语规范](#术语规范)。

重复测量设计涉及一个或多个步骤，其中研究人员不能将一个或更多因素的水平随机分配给实验单元。使用时间作为因素是无法使用随机化的最常见实验情况。例如，当在同一实验单元的多个时间点收集数据时，无法随机化时间点的顺序。时间 1 必须是第一个，时间 2 必须是第二个，依此类推。重复测量因子的这种非随机分配会影响实验单元之间的方差和协方差，第 \@ref(chap26) 章中描述的理想条件可能无效。本章介绍了当第 \@ref(chap26) 章中给出的理想条件不成立时，分析重复测量实验数据的策略。此外，还描述了允许检查第 \@ref(chap26) 章中描述的理想条件是否满足的程序。

## 介绍 {#sec27-1}

考虑类似于表 \@ref(tab:table26-1) 中描述的实验情况。设 $y_{ijk}$ 代表处理组 i 中个体[^subject] k 在时间 j 观测到的响应，并设

$$\boldsymbol{y}_{ik}=\begin{bmatrix}{y}_{i1k}\\{y}_{i2k}\\\vdots\\{y}_{ipk}\end{bmatrix}$$

表示处理组 i 中个体 k 的响应向量。

可以用来描述这些数据的模型是

$$y_{ijk}=\mu+\alpha_i+\tau_j+\gamma_{ij}+\varepsilon_{ijk}^*,i=1,2,\ldots,t;~j=1,2,\ldots,p;~k=1,2,\ldots,n_i$$

令

$$\boldsymbol{\varepsilon}_{ik}^*=\begin{bmatrix}{\varepsilon}_{i1k}^*\\{\varepsilon}_{i2k}^*\\\vdots\\{\varepsilon}_{ipk}^*\end{bmatrix}$$

表示处理组 i 中个体 k 的误差向量。假设 $\boldsymbol{\varepsilon}_{ik}^*$ 独立同分布于均值为 $\boldsymbol 0$ 和协方差矩阵为 $\boldsymbol \Sigma$ 的 p 元正态分布。即 $\boldsymbol \varepsilon_{ik}^*\sim i.i.d.N(\boldsymbol{0},\boldsymbol{\Sigma}),i=1,2,\ldots,t;k=1,2,\ldots,n_i$

令

$$\boldsymbol{\Sigma}=\begin{bmatrix}{\sigma}_{11}&{\sigma}_{12}&\cdots&{\sigma}_{1p}\\{\sigma}_{21}&{\sigma}_{22}&\cdots&{\sigma}_{2p}\\\vdots&\vdots&\ddots&\vdots\\{\sigma}_{p1}&{\sigma}_{p2}&\cdots&{\sigma}_{pp}\end{bmatrix}$$

表示重复测量向量的协方差阵。

::: {.definition #27-1}
♦

若 $\boldsymbol{\Sigma}={\lambda}\boldsymbol{I}_t+\boldsymbol{\eta}\boldsymbol{j}'+\boldsymbol{j}\boldsymbol{\eta}'$，其中 $\boldsymbol j$ 是全 1 的 p × 1 向量，$\boldsymbol \eta$ 是由常数组成的 p × 1 向量，则重复测量被认为满足 Huynh–Feldt (H–F) 条件（参见 Huynh and Feldt, 1970）。
:::

当定义 27.1 成立时，$\boldsymbol{\Sigma}$ 具有如下形式

$$\boldsymbol{\Sigma}=\begin{bmatrix}\lambda+2{\eta}_1&{\eta}_1+{\eta}_2&\cdots&{\eta}_1+{\eta}_p\\{\eta}_2+{\eta}_1&\lambda+2{\eta}_2&\cdots&{\eta}_2+{\eta}_p\\\vdots&\vdots&\ddots&\vdots\\{\eta}_p+{\eta}_1&{\eta}_p+{\eta}_2&\cdots&\lambda+2{\eta}_p\end{bmatrix}$$

H-F 条件的一个特例是当重复测量具有**复合对称协方差结构** (compound 
symmetry covariance structure, **CS**)，则对于 $\sigma^2,\rho$ 的某种取值，

$$\boldsymbol{\Sigma}={\sigma}^2{\begin{bmatrix}1&{\rho}&\cdots&{\rho}\\{\rho}&1&\cdots&{\rho}\\\vdots&\vdots&\ddots&\vdots\\{\rho}&{\rho}&\cdots&1\end{bmatrix}}$$

请注意，第 \@ref(sec26-1) 节中讨论的<u>时间的裂区分析</u> (split-plot-in-time analysis) 描述的理想条件是具有复合对称性的重复测量的协方差阵的特例，其中 $\sigma^2=\sigma_\delta^2+\sigma_\varepsilon^2$ 以及 $\rho=\sigma_\delta^2/(\sigma_\delta^2+\sigma_\varepsilon^2)$. **复合对称结构比时间的裂区结构更通用，因为复合对称结构中的 $\rho$ 可以为负**。

:::: {.rmdnote data-latex="{规则}"}
如果满足 H-F 条件，那么许多涉及时间比较的重要问题可以通过分析重复测量实验来回答，就像分析满足第 \@ref(sec26-1) 节给出的理想条件的重复测量实验一样，并且如果重复测量满足复合对称性，则可以使用第 \@ref(chap26) 章中给出的<u>时间的裂区分析</u>方法。特别是，**当且仅当重复测量满足定义 27.1 中给出的 H-F 条件时，时间主效应和时间 × 处理交互效应的<u>时间的裂区检验</u>才可以在统计上证明有效**。此外，使用<u>时间的裂区分析</u>时，比较时间主效应的对比以及比较处理变量特定值内时间效应的对比都具有统计有效性。如果关注双向均值或边际均值，<u>时间的裂区分析</u>可以提供均值的正确估计，但标准误的估计是不正确的。
::::

:::: {.rmdnote data-latex="{规则}"}
如果重复测量具有复合对称性且 $\rho>0$，则<u>时间的裂区分析</u>给出的所有结果都是正确的。在这种情况下，可以说 $\varepsilon_{ik}^*=\delta_{ik}+\varepsilon_{ijk}$，其中 $\delta_{ik}\sim i.i.d.N(0,\sigma_\delta^2),\varepsilon_{ijk}\sim i.i.d.N(0,\sigma_\varepsilon^2)$，以及 $\delta_{ik}$ 和 $\varepsilon_{ijk}$ 是独立的。
::::

**提问：如果 H-F 条件不满足呢？**

[^adjustpvalue]: 原文：but adjust the p-values by adjusting the degrees of freedom corresponding to relevant effect mean squares.

在不满足 H–F 条件的情况下，可以考虑几种分析方法。一种总是合适的方法是将重复测量的向量视为多元响应向量，并使用**多元方差分析** (multivariate analysis of variance, MANOVA) 方法。第二种是使用<u>时间的裂区分析</u>，但通过调整与相关效应均方相对应的自由度来调整 $p$ 值[^adjustpvalue]。第三种方法是使用 SAS^®^-Mixed 和 SAS^®^-Gilimix 程序中可用的混合模型方法，并对协方差结构进行建模。

多元方差分析方法在第 \@ref(sec27-2) 节中描述，调整自由度方法在第 \@ref(sec27-3) 节中讨论，混合模型方法将在第 \@ref(sec27-4) 节中讨论。

## MANOVA 法 {#sec27-2}

本节考虑使用多元方差分析法来分析重复测量实验。这些方法**总是**适用于以下情况

$$\boldsymbol\varepsilon_{ik}^*\sim i.i.d.\mathrm{~}N(\boldsymbol0,\boldsymbol{\Sigma}),\quad i=1,2,\ldots,t;\mathrm{~}k=1,2,\ldots,n_i$$

**MANOVA 方法只能应用于给定实验单元上的所有重复测量都存在或全部缺失的实验**，因为所使用的统计软件将自动从分析中删除一个或多个重复测量值缺失的任何实验单元。然而，对于其他两种方法，这并不存在问题。（对于它们，）分配给每个处理的实验单元数量无需均衡。要充分利用这一章节的内容，读者应需理解第 \@chap(6) 章讨论的模型的矩阵形式。

具有类似于表 \@ref(tab:table26-1) 中给出的结构的重复测量实验的一般多元模型是式 \@ref(eq:26-1) 中定义的模型的矩阵形式的推广。多元模型可以表示为

\begin{equation}
\boldsymbol Y=\boldsymbol {XB}+\boldsymbol E
(#eq:27-1)
\end{equation}

其中 $\boldsymbol Y$ 表示在实验中测量的所有数据。数据矩阵 $\boldsymbol Y$ 的每一行对应于特定的实验单元，每一列对应于一个重复测量中。因此 $\boldsymbol Y$ 是 N × p 矩阵，其中 $N=\sum_{i=1}^t n_i$. 矩阵 $\boldsymbol X$ 是假设秩为 t 的 N × r 设计矩阵。$\boldsymbol B$ 的每一列都是未知参数的 r × 1 向量，每一列对应于特定的重复测量。矩阵 $\boldsymbol E$ 是不可观测随机误差的 N × p 矩阵。假设 $\boldsymbol E$ 的行独立分布于 $N(\boldsymbol{0},\boldsymbol{\Sigma})$. 因此，当 $\boldsymbol E$ 中的行是独立的时，一行中的元素可以彼此相关 (correlated) 且可以具有不同的方差。

对于多元模型 \@ref(eq:27-1)，可以检验以下形式的一般假设

\begin{equation}
H_0{:\boldsymbol{CBM}}=0\mathrm{~vs~}H_a{:\boldsymbol{CBM}}\neq0
(#eq:27-2)
\end{equation}

其中 $\boldsymbol C$ 是秩为 g 的 g × r 矩阵，$\boldsymbol M$ 是秩为 q 的 p × q 矩阵。

[^crossproduct]: 原文：one first needs the least squares estimates of the parameters in B and an observed residual sum-of-squares and cross-products matrix. 译者在此看得云里雾里的，这句话其实不重要，只需理解接下来的式表达即可。

为了检验式 \@ref(eq:27-2) 中的假设，首先需要求得矩阵 $\boldsymbol B$ 中参数的最小二乘估计以及观测到的残差平方和与叉乘矩阵 (cross-products matrix)[^crossproduct]。分别记作 $\hat{\boldsymbol B}$ 和 $\hat{\boldsymbol E}$ ，并由以下式给出：

\begin{equation}
\hat{\boldsymbol B}=(\boldsymbol X^{\prime}\boldsymbol X)^-\boldsymbol X^{\prime}\boldsymbol Y\quad\mathrm{and}\quad\hat{\boldsymbol E}=\boldsymbol Y^{\prime}[\boldsymbol I-\boldsymbol X(\boldsymbol X^{\prime}\boldsymbol X)^-\boldsymbol X^{\prime}]\boldsymbol Y
(#eq:27-3)
\end{equation}

用于检验等式 \@ref(eq:27-2) 中的假设，似然比检验统计量 (likelihood ratio test statistic) 由下式给出

\begin{equation}
\Lambda=\frac{|\boldsymbol R|}{\left|\boldsymbol H+\boldsymbol R\right|}
(#eq:27-4)
\end{equation}

其中

$$\boldsymbol R=\boldsymbol M'\hat{\boldsymbol E}\boldsymbol M,\quad \boldsymbol H=\boldsymbol M^{\prime}\hat{\boldsymbol B}\boldsymbol C^{\prime}\left[\boldsymbol C(\boldsymbol X^{\prime}\boldsymbol X)^{-}\boldsymbol C^{\prime}\right]^{-1}\boldsymbol C \hat{\boldsymbol B} \boldsymbol M$$

以及 $|\boldsymbol W|$ 表示矩阵 $\boldsymbol W$ 的行列式。

这种统计量被称为 Wilks' 似然比准则 (likelihood ratio criterion) (Morrison, 1976). $\Lambda$ 的抽样分布相当复杂，但对于大多数实际目的，可以获得一个近似 $\alpha$ 水平的检验，拒绝 $H_0$ 当

$$-\left(N-t-\frac{|q-g|+1}2\right)\mathrm{log}_{\mathrm{e}}(\Lambda)>\chi_{\alpha,qg}^2$$

只有当 $q$ 和 $g$ 都大于 $2$ 才可使用的一个更好的近似是，拒绝 $H_0$ 当

$$F>F_{\alpha,qg,ab-c}$$

其中

$$F=\frac{(1-\Lambda^{1/b})(ab-c)}{qg\Lambda^{1/b}}$$

以及

\begin{align}
&a=N-t-\frac{|q-s|+1}2 \\
&b=\left(\frac{q^{2}s^{2}-4}{q^{2}+s^{2}-5}\right)^{1/2} \\
& c=\frac{qs-2}2  \\
&s=\min(q,g)
(#eq:27-5)
\end{align}

而只要 $q = 1, 2$ 或 $g = 1, 2$，式 \@ref(eq:27-2) 就存在精确的 $F$ 检验。这些检验如下

1. 对于 $g = 1$ 以及任何 $q$，拒绝 $H_0$ 如果
\begin{equation}
F=\left(\frac{1-\Lambda}\Lambda\right)\left(\frac{N-t-q+1}q\right)>F_{\alpha,q,N-t-q+1}
(#eq:27-6)
\end{equation}

2. 对于 $q = 1$ 以及任何 $g$，拒绝 $H_0$ 如果
\begin{equation}
F=\left(\frac{1-\Lambda}\Lambda\right)\left(\frac{N-t}g\right)>F_{\alpha,g,N-t}
(#eq:27-7)
\end{equation}

3. 对于 $g = 2$ 以及任何 $q > 1$，拒绝 $H_0$ 如果
\begin{equation}
F=\left(\frac{1-\sqrt{\Lambda}}{\sqrt{\Lambda}}\right){\left(\frac{N-t-q+1}q\right)}{>F_{\alpha,2q,2(N-t-q+1)}}
(#eq:27-8)
\end{equation}

4. 对于 $q = 2$ 以及任何 $g > 1$，拒绝 $H_0$ 如果
\begin{equation}
F=\left(\frac{1-\sqrt{\Lambda}}{\sqrt{\Lambda}}\right)\left(\frac{N-t-1}g\right)>F_{\alpha,2g,2(N-t-1)}
(#eq:27-9)
\end{equation}

多元方法的一个缺点是必须满足 $p < N – t$. 当 $p \ge N – t$ 时，通常可以将相邻的重复测量组合成 $p^*$ 个新变量，或者仅分析重复测量尺寸为 $p^*$ 的子集，其中 $p^* < N – t$.

[^leafareaindex]: 译者注，叶面积指数 = 总绿叶面积 / 土地面积。
[^addictivemodel]: 译者注，加性模型见第 \@ref(chap7) 章。

为了说明本节中描述的分析，考虑进行一项实验，研究四种高粱 (sorghum) 品种和五种肥料 (fertilizer) 水平在叶面积指数[^leafareaindex] (leaf area index) 上的差异，其中四种高粱品种用 V1, V2, V3 和 V4 表示，五个肥料水平用 1, 2, 3, 4 和 5 表示。还假设这 20 个品种 × 肥料的组合被随机分配到田间 20 个区 (plots). 对于此示例，假设肥料水平和品种之间不存在交互作用，并且可以使用基本的双向加性模型[^addictivemodel]来分析数据。最后，假设从植物出苗后两周开始，在五周内的每周对每个品种 × 肥料的区中进行叶面积指数的测量。获得的数据列于表 \@ref(tab:table27-1) 中。

<table>
<caption>(\#tab:table27-1)四个高粱品种的叶面积指数</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.1.png">
</td>
  </tr>
</tbody>
</table>

对于表 \@ref(tab:table27-1) 中的数据，数据矩阵为

$$\boldsymbol Y=\begin{bmatrix}5.00&4.84&4.02&3.75&3.13\\4.42&4.30&3.67&3.29&2.83\\4.42&4.10&3.46&3.09&2.82\\4.01&3.89&3.21&2.89&2.56\\3.36&3.10&2.67&2.47&2.16\\5.82&5.60&5.05&4.72&4.46\\5.73&5.59&5.00&4.65&4.42\\5.31&5.19&4.86&4.44&4.22\\4.92&4.66&4.56&4.16&3.99\\3.96&3.86&3.50&3.13&2.95\\5.65&5.97&5.27&5.07&4.52\\5.39&5.49&5.08&4.87&4.32\\5.15&5.28&4.93&4.67&4.15\\4.50&4.89&4.74&4.49&4.10\\3.75&3.74&3.55&3.28&3.00\\5.86&5.60&5.37&5.00&4.37\\5.82&5.55&5.29&4.95&4.07\\5.26&5.06&4.76&4.48&3.94\\4.87&4.75&4.55&4.33&3.83\\3.96&3.76&3.56&3.18&2.96\end{bmatrix}$$

参数矩阵为

$$\boldsymbol{B}=\begin{bmatrix}
\mu^{(1)}&\mu^{(2)}&\mu^{(3)}&\mu^{(4)}&\mu^{(5)}\\
\tau_1^{(1)}&\tau_1^{(2)}&\tau_1^{(3)}&\tau_1^{(4)}&\tau_1^{(5)}\\
\tau_2^{(1)}&\tau_2^{(2)}&\tau_2^{(3)}&\tau_2^{(4)}&\tau_2^{(5)}\\
\tau_3^{(1)}&\tau_3^{(2)}&\tau_3^{(3)}&\tau_3^{(4)}&\tau_3^{(5)}\\
\tau_4^{(1)}&\tau_4^{(2)}&\tau_4^{(3)}&\tau_4^{(4)}&\tau_4^{(5)}\\
\beta_1^{(1)}&\beta_1^{(2)}&\beta_1^{(3)}&\beta_1^{(4)}&\beta_1^{(5)}\\
\beta_2^{(1)}&\beta_2^{(2)}&\beta_2^{(3)}&\beta_2^{(4)}&\beta_2^{(5)}\\
\beta_3^{(1)}&\beta_3^{(2)}&\beta_3^{(3)}&\beta_3^{(4)}&\beta_3^{(5)}\\
\beta_4^{(1)}&\beta_4^{(2)}&\beta_4^{(3)}&\beta_4^{(4)}&\beta_4^{(5)}\\
\beta_5^{(1)}&\beta_5^{(2)}&\beta_5^{(3)}&\beta_5^{(4)}&\beta_5^{(5)}\\
\end{bmatrix}$$

其中 $\tau$ 对应不同品种，$\beta$ 对应不同肥料水平，设计矩阵为

$$\boldsymbol X=\begin{bmatrix}
1 & 1 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
1 & 1 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
1 & 1 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 \\
1 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 \\
1 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 \\
1 & 0 & 1 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
1 & 0 & 1 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
1 & 0 & 1 & 0 & 0 & 0 & 0 & 1 & 0 & 0 \\
1 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 1 & 0 \\
1 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 1 \\
1 & 0 & 0 & 1 & 0 & 1 & 0 & 0 & 0 & 0 \\
1 & 0 & 0 & 1 & 0 & 0 & 1 & 0 & 0 & 0 \\
1 & 0 & 0 & 1 & 0 & 0 & 0 & 1 & 0 & 0 \\
1 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 1 & 0 \\
1 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 1 \\
1 & 0 & 0 & 0 & 1 & 1 & 0 & 0 & 0 & 0 \\
1 & 0 & 0 & 0 & 1 & 0 & 1 & 0 & 0 & 0 \\
1 & 0 & 0 & 0 & 1 & 0 & 0 & 1 & 0 & 0 \\
1 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 1 & 0 \\
1 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 1 
\end{bmatrix} $$

注意，$\boldsymbol B$ 的每一列代表数据矩阵 $\boldsymbol Y$ 的第 $j(j=1,2,3,4,5)$ 个响应列的双向加性模型所需的参数，设计矩阵 $\boldsymbol X$ 是一个 20 × 10 矩阵且秩等于 8；因此 $N = 20,p = 10$ 且 $t = 8$。式 \@ref(eq:27-3) 给出的 $\hat{\boldsymbol B}$ 的值为

$$\hat{\boldsymbol{B}}=\begin{bmatrix}3.350&3.283&3.003&2.788&2.150\\0.222&0.106&-0.198&-0.260&-0.312\\1.128&1.040&0.990&0.874&0.996\\0.868&1.134&1.110&1.130&1.006\\1.134&1.004&1.102&1.042&0.822\\1.395&1.398&1.173&1.115&0.982\\1.152&1.128&1.006&0.940&0.772\\0.847&0.803&0.748&0.685&0.645\\0.387&0.443&0.511&0.483&0.482\\-0.430&-0.489&-0.434&-0.470&-0.370\end{bmatrix}$$

以及式 \@ref(eq:27-3) 中的 $\boldsymbol E$ 的值为

$$\hat{\boldsymbol{E}}=\begin{bmatrix}0.237&0.171&0.162&0.228&0.129\\0.171&0.247&0.163&0.231&0.135\\0.162&0.163&0.268&0.303&0.184\\0.228&0.231&0.303&0.392&0.241\\0.129&0.135&0.184&0.241&0.247\end{bmatrix}$$

等品种主效应均值检验由式 \@ref(eq:27-4) 得出，通过取

$$\boldsymbol{C}=\begin{bmatrix}0&1&-1&0&0&0&0&0&0&0\\0&1&0&-1&0&0&0&0&0&0\\0&1&0&0&-1&0&0&0&0&0\end{bmatrix}\quad\mathrm{and}\quad\boldsymbol{M}=\begin{bmatrix}1\\1\\1\\1\\1\\1\end{bmatrix}$$

那么 $\Lambda= 0.04345,g = 3,q = 1$. 由于 $q = 1$，可以使用式 \@ref(eq:27-7)，得到

$$F=\frac{1-0.04345}{0.04345}\times\frac{12}3=88.06$$

具有 3 个和 12 个自由度。观察到的显著性水平 $\hat\alpha$ 小于 0.0001.

等时间主效应均值的检验由式 \@ref(eq:27-4) 得出，通过取

$$\boldsymbol{C}=\begin{bmatrix}1&\frac14&\frac14&\frac14&\frac14&\frac15&\frac15&\frac15&\frac15&\frac15\\\end{bmatrix}\quad\mathrm{and}\quad\boldsymbol{M}=\begin{bmatrix}1&1&1&1\\-1&0&0&0\\0&-1&0&0\\0&0&-1&0\\0&0&0&-1\end{bmatrix}$$

可得到 $\Lambda= 0.04345,g=1,q=4$. 由于 $g = 1$，因此可以使用式 \@ref(eq:27-6)，并得到

$$\begin{aligned}F=&\frac{1-0.00502}{0.00502}\times\frac94=445.96\end{aligned}$$

具有 4 个和 9 个自由度。观察到的显著性水平 $\hat\alpha$ 小于 0.0001.

品种 × 时间交互作用的检验由式 \@ref(eq:27-4) 得出，通过取

$$\boldsymbol{C}=\begin{bmatrix}0&1&-1&0&0&0&0&0&0&0\\0&1&0&-1&0&0&0&0&0&0\\0&1&0&0&-1&0&0&0&0&0\end{bmatrix}\quad\mathrm{and}\quad\boldsymbol{M}=\begin{bmatrix}1&1&1&1\\-1&0&0&0\\0&-1&0&0\\0&0&-1&0\\0&0&0&-1\end{bmatrix}$$

可得到 $\Lambda= 0.01426,g=1,q=4$. 由于 g 和 q 都大于 1，可以使用式 \@ref(eq:27-5)，其中 $s = 3,a = 11,b = 2.646$ 以及 $c = 5$. 那么

$$F=\frac{(1-0.01426^{1/2.646})(11\times2.646-5)}{4\times3\times0.01426^{1/2.646}}=\frac{(1-0.2006)(24.1060)}{12(0.2006)}=8.00$$

具有 12 个和 24.1 个自由度。观察到的显著性水平 $\hat\alpha$ 小于 0.0001.

请注意，等品种均值的检验与通过<u>时间的裂区分析</u>获得的检验相同。这是三个检验中唯一与通过<u>时间的裂区分析</u>获得的相应检验相同的检验。以类似的方式可以得到肥料主效应（$F=94.36$）和肥料 × 时间交互效应（$F=1.91$）的检验统计量。

$\boldsymbol{B}$ 中参数的许多特殊函数是可估的，并且可以对这些可估函数进行推断。$\boldsymbol{B}$ 中参数的大多数有趣的线性函数都可以写成 $\boldsymbol{c'Bm}$ 的形式，其中 $\boldsymbol{c}$ 是 r × 1 向量，$\boldsymbol{m}$ 是 p × 1 向量。正如第 \@ref(chap6) 章中那样，当且仅当存在一个向量 $\boldsymbol{u}$ 使得 $\boldsymbol{X^{\prime}}\boldsymbol{X}\boldsymbol{u}=\boldsymbol{c}$ 时，$\boldsymbol{c'Bm}$ 才是可估的。对于 $\boldsymbol m$ 没有限制。

$\boldsymbol{c'Bm}$ 的最佳估计是 $\boldsymbol{c'\hat Bm}$，其中 $\boldsymbol{\hat B}$ 由式 \@ref(eq:27-3) 中给出。 $\boldsymbol{c'\hat Bm}$ 的标准误估计由下式给出

$$\widehat{s.e.}(\boldsymbol c^{\prime}\hat{\boldsymbol B}\boldsymbol m)=\sqrt{\boldsymbol c^{\prime}(\boldsymbol X^{\prime}\boldsymbol X)^-\boldsymbol c\cdot\frac{{\boldsymbol m^{\prime}\hat{\boldsymbol E}\boldsymbol m}}{N-t}}$$

相应的自由度为 N-t. 因此，$\boldsymbol{c'Bm}$ 的 $(1 - \alpha)100\%$ 置信区间由下式给出

\begin{equation}
\boldsymbol{c^{\prime}\hat{B}m}\pm t_{\alpha,N-t}\left[\widehat{s.e.}(\boldsymbol{c^{\prime}\hat{B}m})\right]
(#eq:27-10)
\end{equation}

用于检验 $H_0:\boldsymbol{c'Bm}=a_0$ 的 $t$ 统计量由下式给出

\begin{equation}
t=\frac{\boldsymbol{c'\hat Bm}-a_0}{\widehat{s.e.}(\boldsymbol{c'\hat Bm})}
(#eq:27-11)
\end{equation}

如果 $t>t_{\alpha/2,N-t}$ 则拒绝 $H_0$.

例如，考虑估计 V1 的边际均值。对于此边际均值，

$$\boldsymbol c'=\begin{bmatrix}1&1&0&0&0&0.2&0.2&0.2&0.2&0.2\end{bmatrix}\quad\mathrm{and}\quad\boldsymbol{m}=\begin{bmatrix}0.2\\0.2\\0.2\\0.2\\0.2\end{bmatrix}$$

$\boldsymbol{c'\hat Bm}$ 的值为 3.496，其标准误估计为 0.0594. 利用式 \@ref(eq:27-10) 可以得到 V1 边际均值的 95% 置信区间为 3.496 ± (2.179) (0.0594).

作为第二个例子，考虑估计表 \@ref(tab:table27-1) 中数据的时间 1 的 边际均值。对于此边际均值，

$$\boldsymbol c^{\prime}=[1\quad0.25\quad0.25\quad0.25\quad0.2\quad0.2\quad0.2\quad0.2\quad0.2]\quad\mathrm{and}\quad\boldsymbol{m}=\begin{bmatrix}1\\0\\0\\0\\0\end{bmatrix}$$

$\boldsymbol{c'\hat Bm}$ 的值为 4.858，其标准误估计为 0.0315.

作为第三个示例，考虑估计 V1 和 V2 边际均值之差。此时

$$\boldsymbol c^{\prime}=\begin{bmatrix}0&1&-1&0&0&0&0&0&0&0\end{bmatrix}\quad\mathrm{and}\quad\boldsymbol{m}=\begin{bmatrix}0.2\\0.2\\0.2\\0.2\\0.2\end{bmatrix}$$

$\boldsymbol{c'\hat Bm}$ 的值为 -1.094，其标准误估计为 0.0840. 用于比较这两个边际平均值的 $t$ 统计量为 $t = -1.094/0.0840 = 13.02$，其观察到的显著性水平为 $\hat\alpha < 0.0001$.

当使用 `MANOVA` 选项及其 `M = ` 选项时，可以从 SAS^®^-GLM 程序获得上面给出的许多推断结果。为了说明这一点，使用表 \@ref(tab:table27-2) 中给出的 SAS 命令重新分析表 \@ref(tab:table27-1) 中的数据。第一个 `MANOVA` 选项用于获取 E^，第二个 `MANOVA` 选项用于获取比较品种主效应均值和肥料主效应均值的检验【注意，M = (0.2 0.2 0.2 0.2 0.2) 告诉 `MANOVA` 选项计算五次重复的时间测量平均值】。肥料 × 时间交互作用和品种 × 时间交互作用的检验是从第三个 `MANOVA` 选项获得的，其中 $\boldsymbol M$ 是时间对比的 4 × 5 矩阵。

<table>
<caption>(\#tab:table27-2)使用 MANOVA 分析表 \@ref(tab:table27-1) 中的数据的 SAS-GLM 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.2.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table27-3) 给出了从表 \@ref(tab:table27-2) 的多元方差分析中获得的 $\hat{\boldsymbol E}$ 值。表 \@ref(tab:table27-4) 给出了肥料和品种主效应的 MANOVA 检验，表 \@ref(tab:table27-5) 给出了 MANOVA 交互检验。

<table>
<caption>(\#tab:table27-3)误差平方和和叉乘矩阵</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.3.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table27-4)肥料和品种主效应检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.4.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table27-5)MANOVA 交互检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.5.png">
</td>
  </tr>
</tbody>
</table>

`MANOVA` 和 `CONTRAST` 选项都可以用来检验形如 $H_0{:\boldsymbol{c'Bm}}=0$ 的假设。例如，通过将以下三个语句添加到表 \@ref(tab:table27-2) 中的 SAS 命令中


```sas
CONTRAST ‘V1-V2’ VARIETY 1 -1 0 0;
MANOVA M=(.2 .2 .2 .2 .2) / PRINTE;
MANOVA M=(1 0 0 0 0) / PRINTE;
```

我们可以从第一个 MANOVA 语句中获得对品种 1 和 2 主效应均值比较的 Wilks 统计量，以及从第二个 MANOVA 语句中对时间 1 处的品种 1 和品种 2 进行比较的 Wilks 统计量。通过在上述每个 MANOVA 命令中包含 `PRINTE` 选项，还可以获得 $\boldsymbol {m'\hat E m}$ 的值。表 \@ref(tab:table27-6) 给出了比较品种 1 主效应均值与品种 2 主效应均值的检验统计量 ($F = 169.4$)。请注意，169.4 是本节前面给出的 $t$ 统计量的平方。另请注意，当 ${\boldsymbol m}' = [0.2 0.2 0.2 0.2 0.2]$ 时，$\boldsymbol {m'\hat E m} = 0.2115$.

<table>
<caption>(\#tab:table27-6)比较品种 1 和品种 2 主效应均值的检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.6.png">
</td>
  </tr>
</tbody>
</table>

## $p$ 值调整法 {#sec27-3}

回想一下，当重复测量满足复合对称性假设时，重复测量实验的分析可以通过第 \@ref(chap26) 章中描述的<u>时间的裂区</u>方法获得。如果只考虑涉及时间因素差异的检验，那么当重复测量满足 H-F 条件时，<u>时间的裂区</u>检验也有效。第二种分析重复测量实验的方法已被证明比第 \@ref(sec27-2) 节中描述的基于 MANOVA 方法的分析更强大。该方法是在不满足 H-F 条件的情况下执行<u>时间的裂区分析</u>，并调整检验统计量的分子和分母自由度，从而产生调整的 $p$ 值。

设

$$\boldsymbol{\Sigma}=\begin{bmatrix}{\sigma}_{11}&{\sigma}_{12}&\cdots&{\sigma}_{1p}\\{\sigma}_{21}&{\sigma}_{22}&\cdots&{\sigma}_{2p}\\\vdots&\vdots&\ddots&\vdots\\{\sigma}_{p1}&{\sigma}_{p2}&\cdots&{\sigma}_{pp}\end{bmatrix}$$

为重复测量的协方差阵。令

\begin{equation}
\theta=\frac{p^2(\bar{\sigma}_{ii}-\bar{\sigma}_{\cdot\cdot})^2}{(p-1)\left[\sum_{i=1}^p\sum_{j=1}^p\sigma_{ij}^2-2p\sum_{i=1}^p\bar{\sigma}_{i\cdot}^2+p^2\bar{\sigma}_{\cdot\cdot}^2\right]}
(#eq:27-12)
\end{equation}

其中

$$\bar{\sigma}_{ii}=\frac1p{\sum_{j=1}^p\sigma_{ij}},\quad\bar{\sigma}_{i\cdot}=\frac1p{\sum_{j=1}^p\sigma_{ij}},\quad\mathrm{and}\quad\bar{\sigma}_{\cdot\cdot}=\frac1{p^2}{\sum_{i=1}^p{\sum_{j=1}^p\sigma_{ij}}}$$

Box (1954) 提出用 $\theta$ 来衡量 $\boldsymbol \Sigma$ 偏离复合对称性的程度，并证明了 $1/(p - 1) < \theta < 1$. $\theta$的值越小，$\boldsymbol \Sigma$ 偏离复合对称性就越远。假设 $F_{TIME}$ 和 $F_{TIME × TRT}$ 分别是时间主效应和时间 × 处理交互作用的<u>时间的裂区</u>检验统计量，适用于表 \@ref(tab:table26-1) 所描述的重复测量场景。Box 表明，当 $\boldsymbol \Sigma$ 关于复合对称性偏离了 $\theta$ 时，则，当没有时间效应时

\begin{equation}
F_{TIME}\text{ is approximately distributed as }F[\theta(p-1),\theta(N-t)(p-1)]
(#eq:27-13)
\end{equation}

以及当无交互效应时

\begin{equation}
F_{TIME\times TRT}\text{ is approximately distributed as F}[\theta(t-1)(p-1),\theta(N-t)(p-1)]
(#eq:27-14)
\end{equation}

不幸的是，由于 $\sigma_{ij}$ 未知，所以 $\theta$ 未知。如果可以估计 $\sigma_{ij}$，则可以估计 $\theta$.

已有 $\theta$ 的三种估计。第一个被称为 Box 保守估计。Box (1954) 建议保守的方法是取 $\theta$ 为其可能值的最小值。即，取 $\theta = 1/(p - 1)$. 这可能太保守，不予推荐。

Greenhouse 和 Geisser (1959) 提出了第二种可能性。令 $\hat{\boldsymbol{Q}}=\boldsymbol{C\hat{\boldsymbol{\Sigma}}}\boldsymbol{C'}$ 其中 $\boldsymbol{C}$ 是满足 $\boldsymbol{Cj} = \boldsymbol 0$ 且 $\boldsymbol{CC'} = \boldsymbol I_{p-1}$ 的任何 $(p - 1) × p$ 矩阵，其中 $\boldsymbol j$ 是全 1 的 p × 1 向量，而 $\boldsymbol I_{p-1}$ 是 (p - 1) × ( p - 1) 的单位阵。然后

$$\hat{{\theta}}~=~\frac{\left(\sum_{i=1}^{p-1}\hat{{q}}_{ii}\right)^2}{(p-1)\sum_{i=1}^{p-1}\sum_{j=1}^{p-1}\hat{{q}}_{ij}^2}=\frac{\left[\operatorname{tr}(\hat{\boldsymbol{Q}})\right]^2}{(p-1)\operatorname{tr}(\hat{\boldsymbol{Q}}\hat{\boldsymbol{Q}}^{\prime})}$$

是 $\theta$ 的 Greenhouse and Geisser (G–G) 估计，其中 $\text{tr}(\boldsymbol B)$ 是矩阵 $\boldsymbol B$ 的迹。

Huynh and Feldt (1976) 提出了第三种估计 $\theta$ 的方法。他们的方法将 $\theta$ 估计为

$$\tilde{{\theta}}=\frac{N(p-1)\hat{{\theta}}-2}{(p-1)(N-r-(p-1)\hat{{\theta}})}$$

其中 N 为总样本量，N - r 为<u>时间的裂区分析</u>中误差项的自由度。

可以注意到，Box 的修正是三种调整方法中最保守的，而 H–F 修正是三个调整方法中最不保守的。还应注意的是，如果 $\theta$ 的任何估计值大于 1，则在计算调整的 $p$ 值时，在式 \@ref(eq:27-13) 和 \@ref(eq:27-14) 中用 1 代替 $\theta$，因为不希望增加任何效应的自由度。

可以注意到，如果 $\boldsymbol{\Sigma}$ 满足 H-F 条件，则对于某 $\lambda$，$\boldsymbol{C\Sigma C'}=\lambda \boldsymbol{I}$. 当这为真时，我们说 $\boldsymbol{C\Sigma C'}$ 满足球形条件。此时 $H_0{:\boldsymbol{C\Sigma C'}}=\lambda\boldsymbol{I}$ 的似然比检验是可用的。似然比检验统计量由下式给出：

\begin{equation}
\boldsymbol{\Lambda}=\frac{|\hat{\boldsymbol{Q}}|}{\left[\frac1{p-1}\mathrm{tr}(\hat{\boldsymbol{Q}})\right]^{p-1}}
(#eq:27-15)
\end{equation}

若 $-2\log_{\mathrm{e}}(\Lambda)>\chi_{\alpha,p(p-1)/2-1}^2$ 则拒绝 $H_0$. 存在许多矩阵 $\boldsymbol C$ 满足 $\boldsymbol{Cj} = \boldsymbol 0 $ 和 $\boldsymbol{CC'} = \boldsymbol I_{p-1}$. 应当注意，$\Lambda$ 的值不取决于可能选择哪个可行的 $\boldsymbol C$ 矩阵。

例如，考虑表 \@ref(tab:table27-1) 中的数据。对于这些数据

$$\hat{\boldsymbol E}=\begin{bmatrix}0.237&0.171&0.162&0.228&0.129\\0.171&0.247&0.163&0.231&0.135\\0.162&0.163&0.268&0.303&0.184\\0.228&0.231&0.303&0.392&0.241\\0.129&0.135&0.184&0.241&0.247\end{bmatrix}$$

基于 12 个自由度。

因此，

$$\begin{aligned}
\boldsymbol{\hat{\Sigma}}=\frac1{12}\boldsymbol{\hat{E}}& =\frac{1}{12}\begin{bmatrix}0.237&0.171&0.162&0.228&0.129\\0.171&0.247&0.163&0.231&0.135\\0.162&0.163&0.268&0.303&0.184\\0.228&0.231&0.303&0.392&0.241\\0.129&0.135&0.184&0.241&0.247\end{bmatrix}  \\\\
&=\begin{bmatrix}0.0198&0.0143&0.0135&0.0190&0.0108\\0.0143&0.0206&0.0136&0.0193&0.0113\\0.0135&0.0136&0.0223&0.0253&0.0153\\0.0190&0.0193&0.0253&0.0327&0.0201\\0.0108&0.0113&0.0153&0.0201&0.0206\end{bmatrix}
\end{aligned}$$

取

$$\boldsymbol{C}=\begin{bmatrix}\frac1{\sqrt{2}}&-\frac1{\sqrt{2}}&0&0&0\\\frac1{\sqrt{6}}&\frac1{\sqrt{6}}&-\frac2{\sqrt{6}}&0&0\\\frac1{\sqrt{12}}&\frac1{\sqrt{12}}&\frac1{\sqrt{12}}&-\frac3{\sqrt{12}}&0\\\frac1{\sqrt{20}}&\frac1{\sqrt{20}}&\frac1{\sqrt{20}}&\frac1{\sqrt{20}}&-\frac4{\sqrt{20}}\end{bmatrix}$$

可得

$$\boldsymbol{Q}=\boldsymbol{C}\boldsymbol{\hat{\Sigma}}\boldsymbol{C}^{\prime}=\begin{bmatrix}0.00592&-0.00019&-0.00003&0.00013\\-0.00019&0.00830&0.00399&0.00178\\-0.00003&0.00399&0.00486&0.00078\\0.00013&0.00178&0.00078&0.00875\end{bmatrix}$$

然后得到 $\hat{{\theta}}=0.795$ 和 $\tilde{{\theta}}=1.747$. 用于检验 $H_0{:\boldsymbol{C\Sigma C'}}=\lambda\boldsymbol{I}$ 的似然比检验统计量的值

$$\Lambda=\frac{|\hat{\boldsymbol{Q}}|}{\left[\frac1{p-1}\mathrm{tr}(\hat{\boldsymbol{Q}})\right]^{p-1}}=\frac{1.20895(10)^{-9}}{\left(\frac{0.027831}4\right)^4}=0.5159$$

因此 $-2\log_e(\Lambda)=1.32$，将其与 $\chi_{0.05,9}^2=16.919$ 比较，因此根据这些数据无法拒绝 $H_0{:\boldsymbol{C\Sigma C'}}=\lambda\boldsymbol{I}$.

假设表 \@ref(tab:table27-1) 中的数据是使用 SAS-GLM 程序进行<u>时间的裂区分析</u>。表 \@ref(tab:table27-7) 中给出了所使用的 SAS 命令。

<table>
<caption>(\#tab:table27-7)使用时间的裂区分析来分析表 \@ref(tab:table27-1) 中的数据的 SAS-GLM 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.7.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table27-8) 中给出了肥料和品种的主效应检验，表 \@ref(tab:table27-9) 中给出了时间主效应以及时间 × 肥料和时间 × 品种交互效应的检验以及 G–G 调整的自由度。调整的自由度是通过将原始的<u>时间的裂区分析</u>中的自由度分别与$\hat\theta$ 相乘来计算的。表 \@ref(tab:table27-8) 中的检验无需调整。请注意，表 \@ref(tab:table27-8) 中的 $F$ 统计量与前面给出的品种和肥料主效应的 MANOVA 检验统计量一致。即使使用了调整的自由度，表 \@ref(tab:table27-9) 中的所有显著概率仍小于 0.0001.

<table>
<caption>(\#tab:table27-8)时间的裂区主效应均值检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.8.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table27-9)具有 G-G 调整自由度的时间的裂区检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.9.png">
</td>
  </tr>
</tbody>
</table>

SAS-GLM 程序及其 `repeated` 选项也可用于产生本节所述的分析以及 MANOVA 分析的一些结果。表 \@ref(tab:table27-10) 给出了可用于表 \@ref(tab:table27-1) 中数据分析的 SAS 命令。

<table>
<caption>(\#tab:table27-10)使用 `repeated` 选项分析表 \@ref(tab:table27-1) 中的数据的 SAS-GLM 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.10.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table27-11) 中标记为 "Orthogonal components" 的一行给出了 $H_0{:\boldsymbol{C\Sigma C'}}=\lambda\boldsymbol{I}$ 的检验。也就是说，这是对 H–F 条件是否满足的检验。表 \@ref(tab:table27-12) 给出了比较时间主效应的 MANOVA 检验，表  \@ref(tab:table27-13)  给出了时间与品种交互作用和时间与肥料交互作用的 MANOVA 检验，表  \@ref(tab:table27-14)  给出了品种与肥料主效应的检验。最后，表  \@ref(tab:table27-15)  给出了时间主效应的<u>时间的裂区分析</u>，以及时间与品种和时间与肥料交互作用。给出了三组不同的 $p$ 值。它们是未调整的<u>时间的裂区分析</u> $p$ 值, Greenhouse–Geisser (G-G) 调整的 $p$ 值和 Huyhn–Feldt (H–F) 调整的 $p$ 值。还给出了 $\hat \theta=0.7954$ 和 $\tilde\theta=1.7473$.

<table>
<caption>(\#tab:table27-11)检验 $H_0{:\boldsymbol{C\Sigma C'}}=\lambda\boldsymbol{I}$</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.11.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table27-12)时间主效应均值的 MONOVA 检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.12.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table27-13)时间与肥料与品种的交互作用的 MONOVA 检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.13.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table27-14)肥料与品种主效应检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.14.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table27-15)调整 $p$ 值的<u>时间的裂区分析</u>
</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.15.png">
</td>
  </tr>
</tbody>
</table>

## 混合模型法 {#sec27-4}

再次考虑第 \@ref(sec27-1) 节中描述的符号，其中

$$\boldsymbol y_{ik}=\begin{bmatrix}y_{i1k}\\y_{i2k}\\\vdots\\y_{ipk}\end{bmatrix}$$

是处理组 i 中个体 k 的响应向量。假设误差向量

\begin{equation}
\boldsymbol \varepsilon_{ik}^*\thicksim\mathrm{independent~}N(\boldsymbol 0,\boldsymbol{\Sigma}_i),\quad i=1,2,\ldots,t;\quad k=1,2,\ldots,n_i
(#eq:27-16)
\end{equation}

设 $\boldsymbol y$ 是一个包含所有数据向量的向量。即

$$\boldsymbol y^{\prime}=[\boldsymbol y_{11}^{\prime}\quad \boldsymbol y_{12}^{\prime}\quad\ldots\quad \boldsymbol y_{1n_{1}}^{\prime}\quad \boldsymbol y_{21}^{\prime}\quad \boldsymbol y_{22}^{\prime}\quad\ldots\quad \boldsymbol y_{2n_{2}}^{\prime}\quad\ldots\quad \boldsymbol y_{t1}^{\prime}\quad \boldsymbol y_{t2}^{\prime}\quad \boldsymbol y_{t2}^{\prime}\quad\ldots\quad \boldsymbol y_{tn_{t}}^{\prime}]$$

在式 \@ref(eq:27-16) 的条件下，

$$\boldsymbol V=\mathrm{Cov}(\boldsymbol y)=\begin{bmatrix}\boldsymbol{\Sigma}_1&\boldsymbol{0}&\cdots&\boldsymbol{0}&\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{0}&\cdots&\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{0}\\
\boldsymbol{0}&\boldsymbol{\Sigma}_1&\cdots&\boldsymbol{0}&\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{0}&\cdots&\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{0}\\
\vdots&\vdots&\ddots&\vdots&\vdots&\vdots&\cdots&\vdots&\cdots&\vdots&\vdots&\cdots&\vdots\\
\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{\Sigma}_1&\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{0}&\cdots&\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{0}\\
\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{0}&\boldsymbol{\Sigma}_2&\boldsymbol{0}&\cdots&\boldsymbol{0}&\cdots&\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{0}\\\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{0}&\boldsymbol{0}&\boldsymbol{\Sigma}_2&\cdots&\boldsymbol{0}&\cdots&\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{0}\\\vdots&\vdots&\cdots&\boldsymbol{0}&\boldsymbol{0}&\boldsymbol{0}&\ddots&\vdots&\vdots&\vdots&\vdots&\vdots\\\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{0}&\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{\Sigma}_2&\cdots&\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{0}\\\vdots&\vdots&\cdots&\vdots&\vdots&\vdots&\cdots&\vdots&\ddots&\vdots&\vdots&\vdots&\vdots\\\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{0}&\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{0}&\cdots&\boldsymbol{\Sigma}_t&\boldsymbol{0}&\cdots&\boldsymbol{0}\\\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{0}&\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{0}&\cdots&\boldsymbol{0}&\boldsymbol{\Sigma}_t&\cdots&\boldsymbol{0}\\\vdots&\vdots&\vdots&\vdots&\vdots&\vdots&\cdots&\vdots&\vdots&\vdots&\vdots&\ddots&\vdots\\\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{0}&\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{0}&\boldsymbol{0}&\boldsymbol{0}&\boldsymbol{0}&\cdots&\boldsymbol{\Sigma}_t\end{bmatrix}$$

设 $\boldsymbol{y}\sim N(\boldsymbol{X}\boldsymbol{\beta},V)$，或等价地 $\boldsymbol{y}=\boldsymbol{X\beta}+\boldsymbol{\varepsilon}$ 其中 $\boldsymbol{\varepsilon}\sim N(\boldsymbol{0},V)$. 现假设 $\boldsymbol V$ 已知且 $\boldsymbol\ell'\boldsymbol \beta$ 可估。在这种情况下，可以证明 $\boldsymbol\ell'\boldsymbol \beta$ 的最佳估计是 $\boldsymbol\ell'\hat{\boldsymbol \beta}_\boldsymbol V$ 其中 $\boldsymbol{\hat{\beta}}_V=(\boldsymbol{X}^{\prime}\boldsymbol{V}^{-1}\boldsymbol{X})^{-}\boldsymbol{X'}\boldsymbol{V}^{-1}\boldsymbol{y}$. 此外，$\boldsymbol\ell'\hat{\boldsymbol \beta}_\boldsymbol V$ 的标准误可证明等于 $\boldsymbol\ell'(\boldsymbol{X'V}^{-1}\boldsymbol{X})^-\boldsymbol\ell$. 估计 $\boldsymbol\ell'\hat{\boldsymbol \beta}_\boldsymbol V$ 称为 $\boldsymbol\ell'\boldsymbol \beta$ 的**广义最小二乘估计** (generalized least squares estimator). 不幸的是，$\boldsymbol V$ 很少是已知的。但假设 $\boldsymbol V$ 可被估计为 $\hat {\boldsymbol V}$，那么 $\boldsymbol\ell'\boldsymbol \beta$ 可被估计为 $\boldsymbol\ell'\hat{\boldsymbol \beta}_\boldsymbol V$，其中 $\hat{\boldsymbol \beta}_\boldsymbol V=\left(\boldsymbol{X}^{\prime}\hat {\boldsymbol{V}}^{-1}\boldsymbol{X}\right)^{-}\boldsymbol{X'}\hat {\boldsymbol{V}}^{-1}\boldsymbol{y}$，那么 $\boldsymbol\ell'\hat{\boldsymbol \beta}_\boldsymbol V$ 的一个标准误估计为 $\widehat{s.e.}(\boldsymbol\ell'\hat{\boldsymbol{\beta}}_{\hat{\boldsymbol V}})=\boldsymbol\ell'\left(\boldsymbol{X'}\hat{\boldsymbol{V}}^{-1}\boldsymbol{X}\right)^{-}\boldsymbol\ell$. 估计 $\boldsymbol\ell'\hat{\boldsymbol \beta}_\boldsymbol V$ 称为 $\boldsymbol\ell'\boldsymbol \beta$ 的**广义最小二乘估计的估计** (estimated generalized least squares estimator). 用于检验 $\boldsymbol{\ell'\beta}=0$ 的近似 $t$ 统计量（其自由度也必须近似）为

\begin{equation}
t=\frac{\boldsymbol\ell^{\prime}\hat{\boldsymbol\beta}_{\hat{\boldsymbol V}}}{\boldsymbol\ell^{\prime}(\boldsymbol X^{\prime}\hat{\boldsymbol V}^{-1}\boldsymbol X)\boldsymbol\ell}
(#eq:27-17)
\end{equation}

假设要检验 $H_0{:\boldsymbol{H\beta}}=\boldsymbol0\mathrm{~vs~}H_a{:\boldsymbol{H\beta}}\ne\boldsymbol0$，其中 $\boldsymbol{H}$ 是秩为 q 的 q × p 矩阵。用于检验 $H_0$ 的一个近似 $F$ 统计量为

\begin{equation}
F=\frac{(\boldsymbol H\boldsymbol{\hat{\beta}}_{\hat{\boldsymbol V}})^{\prime}\left[\boldsymbol H(\boldsymbol X\hat{\boldsymbol V}^{-1}\boldsymbol{\boldsymbol X})^{-}\boldsymbol{H}^{\prime}\right]^{-1}(\boldsymbol H{\hat{\boldsymbol \beta}}_{\hat{\boldsymbol V}})}q
(#eq:27-18)
\end{equation}

这个 $F$ 统计量的分子自由度由 q 给出，分母自由度必须近似。近似方法超出了本书的范围，但可以推荐两种方法。如果没有缺失数据值，可以使用一种称为 Satterthwaite's method 的方法。这是第二章中讨论的 Satterthwaite's method 的推广。有关更多详细信息，请参见 Giesbrecht and Burns (1985), McLean and Sanders (1988) 和 Fai and Cornelius (1996). 当有缺失值时，应使用 Kenward–Roger's method. 详见 Kenward and Roger (1997).

混合模型程序可用于获得 $\boldsymbol\Sigma_i$ 的估计，由此可估计 $\boldsymbol V$. 所使用的估计方法是基于式 \@ref(eq:27-16) 中的假设的似然方法，即 $\boldsymbol y$ 具有多元正态分布。这些方法超出了本书的范围，只是后面会有一些一般性的评论。$\boldsymbol\Sigma_i$ 的估计取决于 $\boldsymbol\Sigma_i$ 是否具有任何结构。可以考虑各种结构。相对流行的结构包括复合对称 (CS) 结构 、Hyuhn–Feldt 结构、AR(1) 结构、异质复合对称 (heterogeneous compound symmetry, HCS) 和异质 AR(1). 混合模型程序也可以处理 $\boldsymbol\Sigma_i$ 毫无结构的**非结构化** (unstructured, **UN**) 情况。一旦估计了 $\boldsymbol V$，混合模型程序就可以计算 $\boldsymbol\ell'\hat{\boldsymbol \beta}_\boldsymbol V$ 及其标准误估计。表 \@ref(tab:table27-16) 给出了一些相对流行的用于重复测量的协方差结构。表 \@ref(tab:table27-16) 中的结构假设 t 个处理具有相同的参数。即，它假设所有的 $\boldsymbol\Sigma_i$ 是相等的。也可以分别估计每个处理的协方差参数。

<table>
<caption>(\#tab:table27-16)用于重复测量实验的一些流行的协方差结构</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.16.png">
</td>
  </tr>
</tbody>
</table>

综上所述，重复测量模型为

\begin{equation}
y_{ijk}=\mu+\alpha_i+\tau_j+\gamma_{ij}+\varepsilon_{ijk}^*,\quad i=1,2,\ldots,t;~j=1,2,\ldots,p;~k=1,2,\ldots,n_i
(#eq:27-19)
\end{equation}

其中

$$\boldsymbol{\varepsilon}_{ik}^*=\begin{bmatrix}{\varepsilon}_{i1k}^*\\{\varepsilon}_{i2k}^*\\\vdots\\{\varepsilon}_{ipk}^*\end{bmatrix}$$

是处理组 i 中个体 k 的误差向量，以及对于 $i=1,2,\ldots,t;k=1,2,\ldots,n_i$，$\boldsymbol{\varepsilon}_{ik}^*$ 独立分布于 $N(\boldsymbol 0,\boldsymbol \Sigma_i)$. 通过向式 \@ref(eq:27-19) 中的模型添加随机个体分量 (random subject component)，可以获得其他的协方差结构。这样的模型将由下式给出

\begin{equation}
y_{ijk}=\mu+\alpha_i+\delta_{ik}+\tau_j+\gamma_{ij}+\varepsilon_{ijk}^*,i=1,2,\ldots,t;~j=1,2,\ldots,p;~k=1,2,\ldots,n_i
(#eq:27-20)
\end{equation}

其中 $\delta_{ik}\thicksim i.i.d.N(0,\sigma_\delta^2),i=1,2,\ldots,t;k=1,2,\ldots,n_i$.

在模型中添加随机个体分量会将重复测量的协方差结构改为

$$\left.\mathrm{Cov}(y_{ik})=\mathrm{Cov}\left(\begin{bmatrix}y_{i1k}\\y_{i2k}\\\vdots\\y_{ipk}\end{bmatrix}\right.\right)=\sigma_\delta^2\boldsymbol{J}_p+\boldsymbol{\Sigma}_i$$

其中 $\boldsymbol{J}_p$ 是全 1 的 p × p 矩阵。当 $\Sigma_i$ 满足复合对称性或 $\Sigma_i$ 具有非结构化形式时，添加随机个体分量没有任何影响，但对于许多其他协方差结构添加随机个体分量可能非常有用。

接下来讨论了用于拟合模型的两种方法，如式 \@ref(eq:27-16) 和 \@ref(eq:27-17) 中的方法。一种是确定模型参数的**最大似然估计** (maximum likelihood, **ML**)，另一种是定义模型参数的**受限最大似然估计** (restricted maximum likelihood estimates, **REML**). 这两种方法都是数值密集型的，通过这些方法拟合模型超出了本书的范围。尽管如此，下面将讨论这两种方法所涉及的一些问题。

### 最大似然法 {#sec27-4-1}

可以基于 $\boldsymbol y$ 的分布形成似然函数，并将似然视为协方差矩阵 $\boldsymbol V$ 中的参数以及 $\boldsymbol \beta$ 中固定效应的参数的函数。这种方法是数值密集型的，得到的估计可能不会产生感兴趣参数函数（如 $\boldsymbol P'\boldsymbol\beta$）的无偏估计，参数函数估计的分布性质可能是未知的，除非是渐近的。求解似然方程需要一个迭代过程，该迭代过程可能收敛，也可能不收敛，它可以收敛于局部最大值而不是全局最大值。ML 法的最大缺点是它倾向于低估方差-协方差参数，这导致固定效应估计的标准误估计太小。这导致 I 类错误率，并且置信区间未达到预期置信水平。有关双向混合模型的讨论，请参阅第 \@ref(chap22) 章。

### 受限最大似然法 {#sec27-4-2}

受限最大似然法也是数值密集型的，并且除了渐近之外，估计的分布特性是未知的。然而，它优于 ML 方法，因为所得的感兴趣参数函数的估计具有较小的偏差，并且 REML 方法不会像 ML 方法那样低估方差-协方差参数。这导致固定效应估计的标准误估计更合适，从而得到更理想的 I 类错误率，以及往往更接近所需置信水平的置信区间。

考虑模型的矩阵形式，$\boldsymbol{y}=\boldsymbol{X\beta}+\boldsymbol{\varepsilon}$，其中 $\boldsymbol{\varepsilon}\sim N(\boldsymbol{0},\boldsymbol{V})$. 令 $\boldsymbol L$ 为满足 $\boldsymbol{LX} = \boldsymbol 0$ 的行满秩矩阵，并且有 $\text{rank}(\boldsymbol L) = n-\text{rank}(\boldsymbol X)$，其中 n 为 $\boldsymbol y$ 的维度。

设 $\boldsymbol y^* = \boldsymbol {Ly}$. 那么 $\boldsymbol y^*\sim N(\boldsymbol 0,\boldsymbol{LVL}^{\prime})$. 由 $\boldsymbol y^*$ 形成的似然函数仅取决于方差-协方差参数。方差-协方差参数的 REML 估计是使受限似然函数达到最大化的参数值，该函数基于 $\boldsymbol y^*$ 的分布。一旦找到方差-协方差参数的估计，就可以估计 $\boldsymbol V$，然后可以根据式 \@ref(eq:27-17) 和 \@ref(eq:27-18) 计算检验统计量。

作为使用混合模型方法的一种策略，通常建议首先利用 REML 确定重复测量值的方差-协方差矩阵的合适结构。一旦获得协方差阵结构后，就可以考虑关于固定效应参数可估函数的推断。有关双向混合模型的 REML 讨论，请参阅第 \@ref(chap22) 章。接下来将考虑选择重复测量协方差阵结构的一些可能性。

诸如 SAS-Mixed 等统计软件会计算几种统计量，它们提供了有用信息，有助于选择重复测量的协方差结构。第一个这样的统计量由 Akaike (1974) 提出，即著名 Akaiki's Information Criterion (AIC). 第二个统计量由 Schwarz (1978) 提出，被称为 Schwarz' Bayesian Criterion (BIC). 第三个由 Hurvich and Tsai (1989) 给出，通常表示为 AICC. 这三个准则都是 $-2\log_{\mathrm{e}}(\hat{{L}})$ 的函数，其中 $\hat{\boldsymbol{L}}$ 是受限最大似然函数的最大值。这些准则的不同之处仅在于随着协方差结构中参数数量的增加，对 $-2\log_{\mathrm{e}}(\hat{{L}})$ 的惩罚方式不同。对于这三种准则，一般认为数值越小越好。令 $d$ 表示所考虑协方差结构中的参数数量，$N$ 表示个体数量，$N^*$ 表示总观测次数减去设计矩阵 $\boldsymbol X$ 的秩。那么，这三种准则可以分别定义为：

\begin{align}
\mathrm{AIC}&=-2\log_\mathrm{e}(\hat{\boldsymbol{L}})+2d\\\mathrm{BIC}&=-2\log_\mathrm{e}(\hat{\boldsymbol{L}})+d[\log_\mathrm{e}(N)]\\\mathrm{AICC}&=-2\log_\mathrm{e}(\hat{\boldsymbol{L}})+2d[N^*/(N^*-d-1)]
(#eq:27-21)
\end{align}

[^lrtcovariance]: 译者注：这里 “特例” 的意思是，若以 “x<y” 表示 x 是 y 的特例，则 “复合对称 < H-F”，“复合对称 < 非结构化”，“H-F < 非结构化”，“AR(1) < 非结构化”，但 AR(1) 与复合对称结构无特例之关系，AR(1) 也与 H-F 结构也无特例之关系，就无法利用 LRT 进行比较。特例之关系可结合表 \@ref(tab:table27-16) 看出。

当样本量相对较大时，还可以计算似然比检验 (LRT) 统计量，该统计量可以用于比较重复测量的两种协方差结构，只要零假设下的结构是备择假设下的结构的特例。例如，可以获得将复合对称结构与 H–F 结构或非结构化结构进行比较的 LRT 统计量，也可以获得将非结构化结构与 H–F 或 AR (1) 结构相比较的 LRT 统计量，但由于 AR(1) 结构与复合对称结构或 H-F 结构彼此之间并无特例之关系，因此无法直接通过 LRT 统计量将 AR(1) 结构与复合对称结构或 H-F 结构进行比较[^lrtcovariance]。为了进行这样的比较，需要依赖于式 \@ref(eq:27-21) 中的准则。

例如，再次考虑表 \@ref(tab:table27-4) 中的药物数据。将使用 SAS-Mixed 程序对这些数据进行各种分析。可以使用的基本 SAS 命令如表 \@ref(tab:table27-17) 所示。基础模型是一个包含药物效应、时间效应和药物-时间交互作用的双向模型。Repeated 语句表明时间是一个重复因子，`Type = CS` 选项表示拟合模型假设重复测量的协方差矩阵具有复合对称性，而 `Subject = Person` 选项则表明 Person 是数据集中的一项变量，用于识别构成特定个体的重复测量构成。请注意，对于这些数据，Person变量具有 1-24 的唯一值，每个个体一个值。注意：如果让 Person 变量在三个药物组中都取值 1-8，那么为了正确识别每个个体，就需要使用 `Subject = Person(Drug)`。最后，`R` 选项被用来在输出中提供 Person = 1 的协方差阵。请注意，该分析假设所有 24 个个体具有相同的协方差阵 $\boldsymbol \Sigma$，所以只需要打印第一个。如果希望获得第 1, 9, 17 号个体的协方差阵，则需要使用 `R = 1, 9, 17` 而不是仅使用 `R`.

<table>
<caption>(\#tab:table27-17)使用 `repeated` 选项分析表 \@ref(tab:table27-4) 中数据的 SAS-Mixed 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.17.png">
</td>
  </tr>
</tbody>
</table>

本例中使用的方法是重新拟合假设了各种协方差结构的模型，以便在考虑对药物和时间因素的推断时选择将使用的协方差结构。因此，使用表 \@ref(tab:table27-18) 中给出的每个 repeated 语句对本例中的数据进行重新分析。在每次分析中，语句 "ODS RTF SELECT R;" 用于获取 RTF 输出文件中的每个 $\boldsymbol R$ 矩阵，而 "ODS output FITSTATISTICS=FIT1;" 用于创建一个数据集，该数据集包含每次分析的各种拟合统计量。所获得的拟合统计量是 AIC, AICC 和 BIC，以及每个分析的 $-2\log_{\mathrm{e}}(\hat{L})$ 值。

<table>
<caption>(\#tab:table27-18)考虑不同 `repeated` 选项的协方差结构</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.18.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table27-18) 中的前六个 repeated 语句对应于表 \@ref(tab:table27-16) 中给出的协方差结构。最后四个 repeated 语句包括 `Group = Drug` 选项，这会导致为每种药物拟合不同的协方差阵。每个药物组中的第一个人是编号为 1, 9 和 17 的人，`R = 1,9,17` 选项将这些的协方差阵估计置于输出中。

表 \@ref(tab:table27-19) 给出了假设复合对称时重复测量的协方差阵估计。请注意，$\hat{{\sigma}}^2=33.4182$ 以及

$$\hat{\rho}=\frac{25.9702}{33.4182}=0.777$$

表 \@ref(tab:table27-20) 给出了假设 H–F 条件下的协方差阵估计。此时

$$\begin{aligned}
\hat{\sigma}^2+2\hat{\eta}_1 &=28.5296 \\
\hat{\sigma}^2+2\hat{\eta}_2 &=40.2442 \\
\hat{\eta}_1+\hat{\eta}_2 &=26.9390
\end{aligned}$$

同时求解 $\hat{\sigma}^2,\hat{\eta}_1,\hat{\eta}_2$ 得到

$$\hat{\sigma}^2=7.4479,\quad\hat{\eta}_1=10.5409,\quad\mathrm{and}\quad\hat{\eta}_2=16.3982$$

<table>
<caption>(\#tab:table27-19)复合对称 (CS) 的 $\hat{\boldsymbol \Sigma}$</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.19.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table27-20)H-F 条件的 $\hat{\boldsymbol \Sigma}$</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.20.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table27-21)非结构化 (UN) 情况的 $\hat{\boldsymbol \Sigma}$</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.21.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table27-22)AR(1) 情况的 $\hat{\boldsymbol \Sigma}$</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.22.png">
</td>
  </tr>
</tbody>
</table>

根据这些，以及 $\boldsymbol \Sigma$ 的其他两个对角元，可以得到

$$\hat{\eta}_3=15.2671,\quad\mathrm{~and~}\quad\hat{\eta}_4=11.6326$$

表 \@ref(tab:table27-21) 给出了假设为非结构化协方差阵的协方差阵估计。这与在 SAS-GLM 分析中包括 `MANOVA/PRINTE` 选项时获得的估计相同。

表 \@ref(tab:table27-22) 给出了假设 AR(1) 结构的协方差阵估计。此时 $\hat{\sigma}^2=32.4945$ 以及 $\hat{\rho}=26.7309/32.4945=0.8226$.

表 \@ref(tab:table27-23) 给出了假设异质复合对称结构的协方差阵估计。在这种情况下，可以得到

$$\begin{aligned}
&\hat{\sigma}_1^2 =31.1038_{\prime}  \\
&\hat{\sigma}_2^2 =38.6218,  \\
&\hat{\sigma}_3^2=29.3752, \\
&\hat{\sigma}_4^2=34.7840
\end{aligned}$$

以及

$$\hat{\rho}=\frac{27.0585}{\sqrt{31.1038\cdot38.6218}}=0.7807$$

<table>
<caption>(\#tab:table27-23)具有异质方差的复合对称结构的 $\hat{\boldsymbol \Sigma}$</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.23.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table27-24) 显示了异质 AR(1) 结构的结果。此时

$$\begin{aligned}
\hat{\sigma}_1^2 &=30.7872 \\
\hat{\sigma}_2^2 &=39.3259 \\
\hat{\sigma}_3^2 &=31.0182 \\
\hat{\sigma}_4^2 &=32.4256 
\end{aligned}$$

以及

$$\hat{\rho}=\frac{29.0118}{\sqrt{30.7872\cdot39.3256}}=0.8338$$

<table>
<caption>(\#tab:table27-24)具有异质方差的 AR(1) 结构的 $\hat{\boldsymbol \Sigma}$</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.24.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table27-25) - \@ref(tab:table27-28) 给出了假设每种药物的协方差阵不同时的协方差阵估计。表 \@ref(tab:table27-25) 为复合对称结构，表 \@ref(tab:table27-26) 为 H–F 结构，表 \@ref(tab:table27-27) 为非结构化结构，表 \@ref(tab:table27-28) 为 AR(1) 结构。每个单独的协方差参数估计留给读者去寻找。

<table>
<caption>(\#tab:table27-25)复合对称下每种药物的三个协方差阵估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.25.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table27-26)H-F 条件下每种药物的三个协方差阵估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.26.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table27-27)非结构化条件下每种药物的三个协方差阵估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.27.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table27-28)AR (1) 条件下每种药物的三个协方差阵估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.28.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table27-29) 显示了上面考虑的每个协方差结构的 AIC 值。可以看出，AIC 的最小值为 488.603，并且该最小值出现在 AR(1) 结构中。因此，在 AIC 下，AR(1) 结构将是最佳的协方差结构选择。

<table>
<caption>(\#tab:table27-29)每种协方差结构的 AIC</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.29.png">
</td>
  </tr>
</tbody>
</table>


表 \@ref(tab:table27-30) 显示了上面考虑的每个协方差结构的 AICC 值。可以看到，AICC 的最小值为 488.751，并且 AR(1) 结构也出现该最小值。因此，在 AICC 下，AR(1) 结构将是最佳的协方差结构选择。

<table>
<caption>(\#tab:table27-30)每种协方差结构的 AICC</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.30.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table27-31) 显示了上面考虑的每个协方差结构的 BIC 值。可以看到，BIC 的最小值为 490.959，并且该最小值出现在 AR(1) 结构中。因此，在 BIC 下，AR(1) 结构将是最佳的协方差结构选择。因此，在 AIC, AICC 和 BIC 这三个准则中，重复测量的首选协方差结构是 AR(1).

<table>
<caption>(\#tab:table27-31)每种协方差结构的 BIC</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.31.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table27-32) 给出了本例中考虑的每个协方差结构的 $-2\log_{\mathrm{e}}(\hat{L})$ 值。假设想要使用似然比检验来比较复合对称结构与 H-F 结构。这是可以做到的，因为复合对称性是 H-F 结构的特例。也就是说，假设我们要检验 $H_0:\boldsymbol \Sigma$ 具有 CS 结构与 $H_a:\boldsymbol \Sigma$具有 H–F 结构。LRT 统计量的值是它们各自 $-2\log_{\mathrm{e}}(\hat{L})$ 值的差。在 $H_0$ 下，该检验统计量近似卡方分布，其自由度为协方差结构参数数量之差。当 p = 4 时，CS 结构有两个参数，而 H–F 结构有 5 个参数。因此，卡方检验统计量为 $\chi^2=488.797-487.034=1.763$，自由度为 5 - 2 = 3。结果观察到的显著性水平为 $\hat\alpha = 0.6230$，并且 $H_0$ 不能被拒绝。

<table>
<caption>(\#tab:table27-32)每种协方差结构的 $-2\log_{\mathrm{e}}(\hat{L})$</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.32.png">
</td>
  </tr>
</tbody>
</table>

作为使用 LRT 的第二个示例，请考虑将 AR(1) 结构与 ARH(1) 结构进行比较。那么当 p = 4 时，这两个结构中的参数数量分别为 2 和 5。这里 $\chi^2=484.603-482.767=1.836$，具有 3 个自由度。卡方临界点是 ${\chi}_{0.05,3}^2=7.815$，并且由于 1.836 < 7.815，因此无法拒绝 AR(1) 结构而选择 ARH(1) 结构。

作为第三个示例，考虑将 AR(1) 结构与 UN 结构进行比较。AR(1) 结构中的参数数量为 2，当 p = 4 时，UN 结构中的参数数量为 p(p + 1)/2 = 10。这里 $\chi^2=484.603-477.372=7.231$ 其中 10 - 2 = 8 个自由度。由此产生的观察到的显著性水平为 $\hat \alpha = 0.5119$，并且无法拒绝 AR(1) 结构而选择 UN 结构。

上述三个检验是可能的，因为其中一个协方差结构是另一个协方差结构的特例。无法给出将 AR(1) 与 CS 进行比较的 LRT，因为两者都不是另一者的特例。同样，人们无法获得将 ARH(1) 与 CSH 进行比较的 LRT.

在考虑固定效应检验之前，假设拟合一个协方差结构为 AR(1) 的模型，但将与个体对应的随机分量添加到模型中。也就是说，假设重复测量满足 AR(1) 结构，则对式 \@ref(eq:27-20) 中的模型进行拟合。这相当于为重复测量选择协方差结构，该结构等于

$$\boldsymbol{\Sigma}={\sigma}_\delta^2\begin{bmatrix}1&1&1&1\\1&1&1&1\\1&1&1&1\\1&1&1&1\\1&1&1&1\end{bmatrix}+{\sigma}^2\begin{bmatrix}1&{\rho}&{\rho}^2&{\rho}^3\\{\rho}&1&{\rho}&{\rho}^2\\{\rho}^2&{\rho}&1&{\rho}\\{\rho}^3&{\rho}^2&{\rho}&1\end{bmatrix}$$

表 \@ref(tab:table27-33) 给出了用于在重复测量和随机个体效应上拟合具有 AR(1) 结构的模型的 SAS-Mixed 命令，表 \@ref(tab:table27-34) 给出了使用表 \@ref(tab:table27-33) 中的命令进行分析的拟合统计量的值。请注意，AIC, AICC 和 BIC 的值均大于没有随机个体效应的 AR(1) 结构模型的值。因此，对于这些药物数据，不建议添加随机个体项。此外，将具有随机个体项的模型与不具有随机个体项的模型进行比较的 LRT 统计量为 $\chi^2=484.6-483.7=0.9$. 向模型添加随机个体项会使参数数量增加 1。因此，$\chi^2=484.6-483.7=0.9$ 将与自由度为 1 的卡方临界点进行比较。这个临界点是 $\chi_{0.05,1}^2=3.84$，并且向模型添加随机个体项并不会给出明显更好的拟合。

<table>
<caption>(\#tab:table27-33)用于向重复测量实验添加随机个体效应的 SAS-Mixed 命令</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.33.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table27-34)具有随机个体效应的重复测量模型的拟合统计量</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.34.png">
</td>
  </tr>
</tbody>
</table>

接下来分析数据以考虑有关药物和时间效应的推断。 SAS-Mixed 命令如表 \@ref(tab:table27-35) 所示。表 \@ref(tab:table27-36) 显示了式 \@ref(eq:27-18) 给出的药物和时间主效应以及药物 × 时间交互作用的检验统计量。值得注意的是，药物与时间的交互作用非常显著。因此，双向药物 × 时间平均值应在每个药物组内和每个时间水平内相互比较。最小二乘均值列于表 \@ref(tab:table27-37) 中，表 \@ref(tab:table27-38) 中给出了每个时间水平下药物之间的成对比较，表 \@ref(tab:table27-39) 中给出了每种药物中时间水平之间的成对比较。

<table>
<caption>(\#tab:table27-35)用于分析固定效应的最终 SAS-Mixed 命令</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.35.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table27-36)III 型固定效应检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.36.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table27-37)模型中每种效应的最小二乘均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.37.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table27-38)每个时间水平下药物均值之间的比较，带有 Tukey-Kramer 调整显著性水平</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.38.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table27-39)每种药物下时间均值之间的比较，带有 Tukey-Kramer 调整显著性水平</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2027.39.png">
</td>
  </tr>
</tbody>
</table>

## 总结 {#sec27-5}

本章考虑了三种替代方法，当<u>时间的裂区分析</u>不合适时，可使用这些方法来分析重复测量实验。考虑的方法包括多元方差分析法、调整 $p$ 值法和混合模型法。人们选择使用哪种方法可能取决于可以使用的统计软件。当拥有可用的混合模型软件时，这种方法可能是最容易使用的方法。如果使用混合模型程序难以收敛，那么如果拥有允许使用多元方差分析法的统计软件，则可以考虑多元方差分析法。如果软件不适用于这两种方法，则可以考虑使用调整 $p$ 值法。

## 练习 {#sec27-6}

