
# 模型的矩阵形式 {#chap6}

> "The sciences do not try to explain, they hardly even try to interpret, they mainly make models. By a model is meant a mathematical construct which, with the addition of certain verbal interpretations, describes observed phenomena. The justification of such a mathematical construct is solely and precisely that it is expected to work" - John Von Neumann

当使用不均衡 (unbalanced) 的固定效应模型、随机效应模型或混合效应模型时，求和符号变得非常费力，有时甚至几乎无法使用。这个问题可以通过使用模型的矩阵形式表达来解决。本章讨论了模型矩阵形式的构建，并描述了如何使用矩阵获得最小二乘估计、检验假设、计算最小二乘均值或总体边际均值以及构建置信区间。在 \@ref(sec6-3) 节中讨论了可估性 (estimability) 的概念。

## 基本符号 {#sec6-1}

模型的矩阵形式可以表示为

\begin{equation}
\underset{n\times1}{\boldsymbol y}=\underset{n × p}{\vphantom{\boldsymbol y}\boldsymbol X} \,\,\underset{p × 1}{\vphantom{\boldsymbol y}\boldsymbol \beta}+\underset{n × 1}{\vphantom{\boldsymbol y}\boldsymbol\varepsilon}
(#eq:6-1)
\end{equation}

其中 $\boldsymbol y$ 表示 $n \times 1$ 观测向量，$\boldsymbol X$ 表示 $n \times p$ 已知常数矩阵，称为**设计矩阵** (designed matrix)，$\boldsymbol \beta$ 表示 $p \times 1$ 未知参数向量，$\boldsymbol \varepsilon$ 表示 $n \times 1$ 未观测的误差向量。第 $i$ 个观测值（$\boldsymbol y$ 的第 $i$ 个元素）的模型具有以下形式

\begin{equation}
y_i=\beta_0+\beta_1x_{i1}+\beta_2x_{i2}+\cdots+\beta_{p-1}x_{ip-1}+\varepsilon_i,\quad i=1,2,\ldots,n
(#eq:6-2)
\end{equation}

用于将模型 \@ref(eq:6-2) 表示为矩阵模型 \@ref(eq:6-1) 的向量和矩阵是

\begin{equation}
\boldsymbol{y}=\begin{bmatrix}y_1\\y_2\\\vdots\\y_n\end{bmatrix},\quad\boldsymbol{X}=\begin{bmatrix}1&x_{11}&x_{12}&\cdots&x_{1p-1}\\1&x_{21}&x_{22}&\cdots&x_{2p-1}\\\vdots&\vdots&\vdots&\ddots&\vdots\\1&x_{n1}&x_{n2}&\cdots&x_{np-1}\end{bmatrix},\quad\boldsymbol{\beta}=\begin{bmatrix}\beta_0\\\beta_1\\\beta_2\\\vdots\\\beta_{p-1}\end{bmatrix},\quad\mathrm{and}\quad\boldsymbol{\varepsilon}=\begin{bmatrix}\varepsilon_1\\\varepsilon_2\\\vdots\\\varepsilon_n\end{bmatrix}
(#eq:6-3)
\end{equation}

式 \@ref(eq:6-3) 中的矩阵可用于表示多种模型类型，包括单向模型、双向模型、析因模型和部分析因模型等设计模型以及回归模型、协方差分析模型、随机模型等。通过指定 $\boldsymbol X$ 的适当元素以及 $\boldsymbol \beta$ 和 $\boldsymbol \varepsilon$ 的适当假设，可以构建效应模型、混合效应模型、裂区模型、重复测量模型和随机系数回归模型。以下各节介绍了各种实验情况下的一些矩阵模型。

### 简单线性回归模型 {#sec6-1-1}

简单线性回归模型可以表示为 $y_i=\beta_0+\beta_1x_i+\varepsilon_i,i=1,2,\ldots,n$ 并且可以用矩阵形式表示为

$$\begin{bmatrix}y_1\\y_2\\\vdots\\y_n\end{bmatrix}=\begin{bmatrix}1&x_1\\1&x_2\\\vdots&\vdots\\1&x_n\end{bmatrix}\begin{bmatrix}{\beta}_0\\{\beta}_1\end{bmatrix}+\begin{bmatrix}{\varepsilon}_1\\{\varepsilon}_2\\\vdots\\{\varepsilon}_n\end{bmatrix}$$

$\boldsymbol X$ 中的全 1 列对应于回归模型的截距 $\boldsymbol \beta_0$，$x_i$ 列对应于回归模型的斜率。

### 单向处理结构模型 {#sec6-1-2}

为了在完全随机设计结构中表示具有 $t$ 个处理的单向处理结构的模型，其中第 $i$ 个处理有 $n_i$ 个观测结果，设自变量 $x_{ij}$ 定义为

$$x_{kij}=\begin{cases}0&\text{ if the }ij\text{th observation is not from the }k\text{th treatment}\\1&\text{ if the }ij\text{th observation is from the }k\text{th treatment}&\end{cases}$$

其中 $i=1,2,\cdots,t$，$j=1,2,\cdots,n_i$. 变量 $x_{kij}$ 被称为**指示变量** (indicator variable)，因为当它的值为 1 时，它表明观测结果来自处理 k. 当它的值为 0 时，它表明观测结果不来自处理 k. 

**均值模型** (means model) 可表示为 $y_{ij}={\mu}_1x_{1ij}+{\mu}_2x_{2ij}+\cdots+{\mu}_tx_{tij}+{\varepsilon}_{ij}$，其中 $i=1,2,\cdots,t$，$j=1,2,\cdots,n_i$. 或者以矩阵表示法

$$\begin{bmatrix}{y}_{11}\\{y}_{12}\\\vdots\\{y}_{1n_1}\\{y}_{22}\\\vdots\\{y}_{2{n}_2}\\\vdots\\{y}_{t1}\\\vdots\\{y}_{t{n}_t}\end{bmatrix}=\begin{bmatrix}1&0&\cdots&0\\1&0&\cdots&0\\\vdots&\vdots&\vdots&\vdots\\1&0&\cdots&0\\0&1&\cdots&0\\0&1&\cdots&0\\\vdots&\vdots&\vdots&\vdots\\0&1&\cdots&0\\\vdots&\vdots&\vdots&\vdots\\0&0&\cdots&1\\\vdots&\vdots&\ddots&\vdots\\0&0&\cdots&1\end{bmatrix}\begin{bmatrix}\mu_1\\\mu_2\\\vdots\\\mu_t\end{bmatrix}+\boldsymbol \varepsilon$$

均值模型通常表示为 $y_{ij}=\mu_i+\varepsilon_{ij}$，其中 $i=1,2,\cdots,t$，$j=1,2,\cdots,n_i$.

**效应模型** (effect model) 可表示为 $y_{ij}=\mu+\tau_1x_{1ij}+\tau_2x_{2ij}+\cdots+\tau_tx_{tij}+\varepsilon_{ij}$，其中 $i=1,2,\cdots,t$，$j=1,2,\cdots,n_i$. 或者以矩阵表示法


$$\begin{bmatrix}{y}_{11}\\{y}_{12}\\\vdots\\{y}_{1{n}_1}\\{y}_{21}\\{y}_{22}\\\vdots\\{y}_{2{n}_2}\\\vdots\\{y}_{t1}\\\vdots\\{y}_{t{n}_t}\end{bmatrix}=\begin{bmatrix}1&1&0&\cdots&0\\1&1&0&\cdots&0\\\vdots&\vdots&\vdots&\vdots&\vdots\\1&1&0&\cdots&0\\1&0&1&\cdots&0\\1&0&1&\cdots&0\\\vdots&\vdots&\vdots&\vdots&\vdots\\1&0&1&\cdots&0\\\vdots&\vdots&\vdots&\vdots&\vdots\\1&0&0&\cdots&1\\\vdots&\vdots&\vdots&\ddots&\vdots\\1&0&0&\cdots&1\end{bmatrix}\begin{bmatrix}\mu\\\tau_1\\\tau_2\\\vdots\\\tau_t\end{bmatrix}+\boldsymbol \varepsilon $$

效应模型通常表示为 $y_{ij}=\mu+\tau_i+\varepsilon_{ij}$，其中 $i=1,2,\cdots,t$，$j=1,2,\cdots,n_i$. 

均值模型和效应模型之间的区别在于，效应模型的设计矩阵包含模型截距 $\mu$ 的全 1 列，而均值模型不包含全 1 列。

### 双向处理结构模型 {#sec6-1-3}

在完全随机设计结构中，对于具有 $t$ 行处理和 $b$ 列处理的双向处理结构，模型的一种形式为

\begin{equation}
y_{ijk}=\mu_{ij}+\varepsilon_{ijk}\quad i=1,2,\ldots,t,j=1,2,\ldots,b,\mathrm{~and~}k=1,2,\ldots,n_{ij}
(#eq:6-4)
\end{equation}

式子 \@ref(eq:6-4) 中使用的模型称为均值模型，可以用矩阵形式表示为

$$\begin{bmatrix}y_{111}\\y_{112}\\\vdots\\y_{11n_1}\\y_{121}\\\vdots\\y_{12n_2}\\\vdots\\y_{1b1}\\\vdots\\y_{1bn_{1b}}
\\y_{211}\\\vdots\\y_{21n_{21}}\\\vdots\\y_{tb1}\\\vdots\\y_{tbn_{tb}}
\end{bmatrix}=\begin{bmatrix}
1&0&\cdots&0&0&\cdots&0\\
1&0&\cdots&0&0&\cdots&0\\
\vdots&\vdots&\cdots&\vdots&\vdots&\cdots&\vdots\\
1&0&\cdots&0&0&\cdots&0\\
0&1&\cdots&0&0&\cdots&0\\
\vdots&\vdots&\cdots&\vdots&\vdots&\cdots&\vdots\\
0&1&\cdots&0&0&\cdots&0\\
\vdots&\vdots&\cdots&\vdots&\vdots&\cdots&\vdots\\
0&0&\cdots&1&0&\cdots&0\\
\vdots&\vdots&\cdots&\vdots&\vdots&\cdots&\vdots\\
0&0&\cdots&1&0&\cdots&0\\
0&0&\cdots&0&1&\cdots&0\\
\vdots&\vdots&\cdots&\vdots&\vdots&\cdots&\vdots\\
0&0&\cdots&0&1&\cdots&0\\
\vdots&\vdots&\cdots&\vdots&\vdots&\cdots&\vdots\\
0&0&\cdots&0&0&\cdots&1\\
\vdots&\vdots&\cdots&\vdots&\vdots&\cdots&\vdots\\
0&0&\cdots&0&0&\cdots&1\end{bmatrix}
\begin{bmatrix}\mu_{11}\\\mu_{12}\\\vdots\\\mu_{1b}\\\mu_{21}\\\vdots\\\mu_{tb}\end{bmatrix} +
\boldsymbol \varepsilon $$

双向效应模型可表示为

\begin{equation}
y_{ijk}=\mu+\tau_i+\beta_j+\gamma_{ij}+\varepsilon_{ijk}\quad i=1,2,\ldots,t,j=1,2,\ldots,b,k=1,2,\ldots,n_{ij}
(#eq:6-5)
\end{equation}

双向效应模型的矩阵形式为

$$
\begin{bmatrix}y_{111}\\y_{112}\\\vdots\\y_{11n_1}\\y_{121}\\\vdots\\y_{12n_2}\\\vdots\\y_{1b1}\\\vdots\\y_{1bn_{1b}}
\\y_{211}\\\vdots\\y_{21n_{21}}\\\vdots\\y_{tb1}\\\vdots\\y_{tbn_{tb}}
\end{bmatrix}=\begin{bmatrix}
1 & 0 & \cdots & 0 & 0 & \cdots & 0 \\
1 & 0 & \cdots & 0 & 0 & \cdots & 0 \\
\vdots & \vdots & \cdots & \vdots & \vdots & \cdots & \vdots \\
1 & 0 & \cdots & 0 & 0 & \cdots & 0 \\
0 & 1 & \cdots & 0 & 0 & \cdots & 0 \\
\vdots & \vdots & \cdots & \vdots & \vdots & \cdots & \vdots \\
0 & 1 & \cdots & 0 & 0 & \cdots & 0 \\
\vdots & \vdots & \cdots & \vdots & \vdots & \cdots & \vdots \\
0 & 0 & \cdots & 1 & 0 & \cdots & 0 \\
\vdots & \vdots & \cdots & \vdots & \vdots & \cdots & \vdots \\
0 & 0 & \cdots & 1 & 0 & \cdots & 0 \\
0 & 0 & \cdots & 0 & 1 & \cdots & 0 \\
\vdots & \vdots & \cdots & \vdots & \vdots & \cdots & \vdots \\
0 & 0 & \cdots & 0 & 1 & \cdots & 0 \\
\vdots & \vdots & \cdots & \vdots & \vdots & \cdots & \vdots \\
0 & 0 & \cdots & 0 & 0 & \cdots & 1 \\
\vdots & \vdots & \cdots & \vdots & \vdots & \cdots & \vdots \\
0 & 0 & \cdots & 0 & 0 & \cdots & 1 \\
\end{bmatrix}  \begin{bmatrix}\mu_{11}\\\mu_{12}\\\vdots\\\mu_{1b}\\\mu_{21}\\\vdots\\\mu_{tb}\end{bmatrix}+
\boldsymbol \varepsilon
$$

本书强调与实验设计情境相对应的模型，而不是纯粹的回归情境。接下来的两个例子将展示如何从数据结构构建这样的模型。

### 示例 6.1：双向处理结构的均值模型 {#sec6-1-4}

表 \@ref(tab:table6-1) 中的信息代表了一个完全随机设计结构中的双向处理结构的数据，其中有三行处理和三列处理，每个单元格有一个或两个观测值。表 \@ref(tab:table6-1) 中数据的均值模型的矩阵形式为：

$$\begin{bmatrix}3\\6\\9\\10\\2\\5\\3\\8\\4\\2\\6\end{bmatrix}=\begin{bmatrix}
1&0&0&0&0&0&0&0&0\\
1&0&0&0&0&0&0&0&0\\
0&1&0&0&0&0&0&0&0\\
0&0&1&0&0&0&0&0&0\\
0&0&0&1&0&0&0&0&0\\
0&0&0&0&1&0&0&0&0\\
0&0&0&0&1&0&0&0&0\\
0&0&0&0&0&1&0&0&0\\
0&0&0&0&0&0&1&0&0\\
0&0&0&0&0&0&0&1&0\\
0&0&0&0&0&0&0&0&1\end{bmatrix}\begin{bmatrix}\mu_{11}\\\mu_{12}\\\mu_{13}\\\mu_{21}\\\mu_{22}\\\mu_{23}\\\mu_{31}\\\mu_{32}\\\mu_{33}\end{bmatrix}+\boldsymbol \varepsilon$$

表 \@ref(tab:table6-1) 中数据的效应模型的矩阵形式为

$$
\begin{bmatrix}3\\6\\9\\10\\2\\5\\3\\8\\4\\2\\6\end{bmatrix}=\begin{bmatrix}
1 & 1 & 0 & 0 & 1 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
1 & 1 & 0 & 0 & 1 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
1 & 1 & 0 & 0 & 0 & 1 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
1 & 1 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
1 & 0 & 1 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
1 & 0 & 1 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
1 & 0 & 1 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
1 & 0 & 1 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
1 & 0 & 0 & 1 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 \\
1 & 0 & 0 & 1 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 \\
1 & 0 & 0 & 1 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 
\end{bmatrix} 
\begin{bmatrix}
\mu_1 \\
\tau_1 \\
\tau_2 \\
\tau_3 \\
\beta_1 \\
\beta_2 \\
\beta_3 \\
\gamma_{11} \\
\gamma_{12} \\
\gamma_{13} \\
\gamma_{21} \\
\gamma_{22} \\
\gamma_{23} \\
\gamma_{31} \\
\gamma_{32} \\
\gamma_{33} 
\end{bmatrix}
+\boldsymbol \varepsilon
$$

或

$$\boldsymbol y=\boldsymbol j\boldsymbol\mu+\boldsymbol X_1\boldsymbol \tau+\boldsymbol X_2\boldsymbol\beta+\boldsymbol X_3\boldsymbol \gamma+\boldsymbol\varepsilon $$

其中，$\boldsymbol j$ 是与上述设计矩阵的第一列对应的 11×1 向量，$\boldsymbol X_1$ 是与第 2-4 列对应的 11×3 矩阵，$\boldsymbol X_2$ 是与第 5-7 列对应的 11×3 矩阵，$\boldsymbol X_3$ 是与以上设计矩阵的最后九列对应的 11×9 矩阵。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table6-1)CRD 结构中双向处理结构的数据</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%206.1.png">
</td>
  </tr>
</tbody>
</table>

其他处理和设计结构的设计矩阵也是以类似的方式构建的。幸运的是，大多数用于拟合非均衡数据结构的软件都会使用上述类型的表示形式，并在模型中指定分类效应 (categorical effects) 时自动生成设计矩阵的必要列。

## 最小二乘估计 {#sec6-2}

一旦模型以矩阵形式指定，分析的下一步是获得参数向量 $\boldsymbol \beta$ 的**最小二乘估计** (least squares estimator)。最小二乘法可用于估计模型的参数。要使用此方法，假设模型可以表示为

\begin{equation}
y_i=f(x_i;\boldsymbol{\beta})+\varepsilon_i\quad\mathrm{~for~}i=1,2,\ldots,n
(#eq:6-6)
\end{equation}

其中 $y_i=f(x_i;\boldsymbol{\beta})$ 是由 $x_i$ 表示的设计变量向量的函数，且依赖于参数向量 $\boldsymbol \beta$. $\boldsymbol \beta$ 的最小二乘估计通常用 $\hat {\boldsymbol \beta}$ 表示，它能最小化平方和

\begin{equation}
SS(\boldsymbol{\beta})=\sum_{i=1}^n[y_i-f(x_i;\boldsymbol{\beta})]^2
(#eq:6-7)
\end{equation}

如果除了假设模型的形式为 \@ref(eq:6-6) 之外，还额外假设 $\varepsilon_i \sim i.i.d\,N(0,\sigma^2),i＝1,2,\cdots,n$，则 $\boldsymbol \beta$ 的最小二乘估计也是**最大似然估计** (maximum likelihood estimator).

例如，完全随机设计中单向处理结构的均值模型的模型函数为

$$f(x_{ij};\boldsymbol{\beta})=\mu_i,\quad i=1,2,\ldots,t;\quad j=1,2,\ldots,n_i$$

$\mu_i$ 的最小二乘估计记作 $\hat \mu_1,\cdots,\hat \mu_t$，它们是这样的值，即能够最小化

$$\begin{aligned}SS({\mu})=\sum_{i=1}^t\sum_{j=1}^{n_j}(y_{ij}-\mu_i)^2\end{aligned}$$

完全随机设计中双向处理结构的均值模型的模型函数为

$$f(x_{ijk};\boldsymbol{\beta})=\mu_{ij}\quad i=1,2,\ldots,t;j=1,2,\ldots,b;k=1,2,\ldots,n_{ij}$$

$\mu_{ij}$ 的最小二乘估计记作 $\hat \mu_{11},\cdots,\hat \mu_{tb}$，它们是这样的值，即能够最小化

$$\begin{aligned}SS({\mu})=\sum_{i=1}^{t}\sum_{j=1}^{b}\sum_{k=1}^{n_{ij}}(y_{ijk}-\mu_{ij})^2\end{aligned}$$

完全随机设计中双向处理结构的效应模型的模型函数为

$$f(x_{ijk};\boldsymbol{\beta})=\mu+\tau_i+{\beta}_j+\gamma_{ij},\quad i=1,2,\ldots,t;j=1,2,\ldots,b;k=1,2,\ldots,n_{ij}$$

$\mu,\tau_i,\beta_j,\gamma_{ij}$ 的最小二乘估计是通过最小化下式来获得的

$$SS(\mu,\tau_i,\beta_j,\gamma_{ij})=\sum_{i=1}^t\sum_{j=1}^b\sum_{k=1}^{n_{ij}}(y_{ijk}-\mu-\tau_i-\beta_j-\gamma_{ij})^2$$

通常，模型可以写成矩阵形式，如式 \@ref(eq:6-1) 中的矩阵形式，并且 $\boldsymbol \beta$ 的最小二乘估计是这样的 $\hat {\boldsymbol \beta}$ ，即能最小化平方和

\begin{equation}
SS(\boldsymbol\beta)=(\boldsymbol y-\boldsymbol X \boldsymbol\beta)^{\prime}(\boldsymbol y-\boldsymbol X \boldsymbol\beta)
(#eq:6-8)
\end{equation}

### 最小二乘方程组 {#sec6-2-1}

矩阵表示和计算可用于确定能使残差平方和最小化的 $\boldsymbol\beta$ 值。当进行最小化时，会得到一组 $\hat {\boldsymbol \beta}$ 必须满足的方程，这些方程被称为模型的**最小二乘方程组** (least squares 
equations) 或**正规方程组** (normal equations). 对于模型 \@ref(eq:6-1)，正规方程组由下式给出：

\begin{equation}
\boldsymbol X^{\prime} \boldsymbol X\hat{\boldsymbol{\beta}}=\boldsymbol X^{\prime}\boldsymbol{y}
(#eq:6-9)
\end{equation}

任何满足正规方程组的向量 $\hat {\boldsymbol \beta}$ 都是 $\boldsymbol\beta$ 的最小二乘估计。对于一些模型，最小二乘估计不必是唯一的。为了帮助读者更加熟悉正规方程组，下面提供了第 \@ref(sec6-1) 节中讨论的模型的正规方程组。

单向均值模型的正规方程组为

$$\begin{bmatrix}n_1&0&\cdots&0\\0&n_2&\cdots&0\\\vdots&\vdots&\ddots&\vdots\\0&0&\cdots&n_t\end{bmatrix}\begin{bmatrix}\hat{\mu}_1\\\hat{\mu}_2\\\vdots\\\hat{\mu}_t\end{bmatrix}=\begin{bmatrix}y_{1\cdot}\\y_{2\cdot}\\\vdots\\y_{t\cdot}\end{bmatrix}\mathrm{~where~}y_{i\cdot}=\sum_{i=1}^{n_i}y_{ij}$$

单向效应模型的正规方程组为

$$\begin{bmatrix}n_{\cdot}&n_1&n_2&\cdots&n_t\\n_1&n_1&0&\cdots&0\\n_2&0&n_2&\cdots&0\\\vdots&\vdots&\vdots&\ddots&\vdots\\n_t&0&0&\cdots&n_t\end{bmatrix}\begin{bmatrix}\hat{\mu}\\\hat{\tau}_1\\\hat \tau_2\\\vdots\\\hat{\tau}_t\end{bmatrix}=\begin{bmatrix}y_{\cdot\cdot}\\y_{1\cdot}\\y_{2\cdot}\\\vdots\\y_{t\cdot}\end{bmatrix}\quad\mathrm{where~}y_{\cdot\cdot}=\sum_{i=1}^t\sum_{j=1}^{n_i}y_{ij}\quad\mathrm{and}\quad n_{\cdot}=\sum_{i=1}^tn_t$$

使用表 \@ref(tab:table6-1) 中数据的双向均值模型的正规方程组为

$$\begin{bmatrix}2&0&0&0&0&0&0&0&0\\0&1&0&0&0&0&0&0&0\\0&0&1&0&0&0&0&0&0\\0&0&0&1&0&0&0&0&0\\0&0&0&0&2&0&0&0&0\\0&0&0&0&0&1&0&0&0\\0&0&0&0&0&0&1&0&0\\0&0&0&0&0&0&0&1&0\\0&0&0&0&0&0&0&0&1\end{bmatrix}\begin{bmatrix}\hat{\mu}_{11}\\\hat{\mu}_{12}\\\hat{\mu}_{13}\\\hat{\mu}_{21}\\\hat{\mu}_{22}\\\hat{\mu}_{23}\\\hat{\mu}_{31}\\\hat{\mu}_{32}\\\hat{\mu}_{33}\end{bmatrix}=\begin{bmatrix}y_{11\cdot}\\y_{12\cdot}\\y_{13\cdot}\\y_{21\cdot}\\y_{22\cdot}\\y_{23\cdot}\\y_{31\cdot}\\y_{32\cdot}\\y_{33\cdot}\end{bmatrix}$$

其中

$$y_{ij\cdot}=\sum_{k=1}^{n_{ij}}y_{ijk}$$

与表 \@ref(tab:table6-1) 中的数据相对应的双向效应模型的正规方程组为

$$$$

$$\begin{bmatrix}
11&4&4&3&4&4&3&2&1&1&1&2&1&1&1&1\\
4&4&0&0&2&1&1&2&1&1&0&0&0&0&0&0\\
4&0&4&0&1&2&1&0&0&0&1&2&1&0&0&0\\
3&0&0&3&1&1&1&0&0&0&0&0&0&1&1&1\\
4&2&1&1&4&0&0&2&0&0&1&0&0&1&0&0\\
4&1&2&1&0&4&0&0&1&0&0&2&0&0&1&0\\
3&1&1&1&0&0&3&0&0&1&0&0&1&0&0&1\\
2&2&0&0&2&0&0&2&0&1&0&0&0&0&0&0\\
1&1&0&0&0&1&0&0&1&0&0&0&0&0&0&0\\
1&1&0&0&0&0&1&0&0&1&0&0&0&0&0&0\\
1&0&1&0&1&0&0&0&0&0&1&0&0&0&0&0\\
2&0&2&0&0&2&0&0&0&0&0&2&0&0&0&0\\
1&0&1&0&0&0&1&0&0&0&0&0&1&0&0&0\\
1&0&0&1&1&0&0&0&0&0&0&0&0&1&0&0\\
1&0&0&1&0&1&0&0&0&0&0&0&0&0&1&0\\
1&0&0&1&0&0&1&0&0&0&0&0&0&0&0&1\end{bmatrix}
\begin{bmatrix}\hat\mu\\\hat\tau_{1}\\\hat{\tau}_2\\\hat\tau_{3}\\\hat\beta_{1}\\\hat\beta_{2}\\\hat\beta_{3}\\\hat\gamma_{11}\\\hat\gamma_{12}\\\hat\gamma_{13}\\\hat\gamma_{21}\\\hat\gamma_{22}\\\hat\gamma_{23}\\\hat\gamma_{31}\\\hat\gamma_{32}\\\hat\gamma_{33}\end{bmatrix}=\begin{bmatrix}y_{\cdot\cdot\cdot}\\y_{1\cdot\cdot}\\y_{2\cdot\cdot}\\y_{3\cdot\cdot}\\y_{\cdot1\cdot}\\y_{\cdot2\cdot}\\y_{\cdot3\cdot}\\y_{11\cdot}\\y_{12\cdot}\\y_{13\cdot}\\y_{21\cdot}\\y_{22\cdot}\\y_{23\cdot}\\y_{31\cdot}\\y_{32\cdot}\\y_{33\cdot}\end{bmatrix}$$

其中

$$y_{\cdot\cdot\cdot}=\sum_{i=1}^t\sum_{j=1}^b\sum_{k=1}^{n_{ij}}y_{ijk},\quad y_{i\cdot\cdot}=\sum_{j=1}^b\sum_{k=1}^{n_{ij}}y_{ijk},\quad y_{\cdot j\cdot}=\sum_{i=1}^t\sum_{k=1}^{n_{ij}}y_{ijk},\quad\mathrm{and}\quad y_{ij\cdot}=\sum_{k=1}^{n_{ij}}y_{ijk}$$

当 $\boldsymbol X^\prime \boldsymbol X$ 满秩时 (Graybill, 1976)，即 $\boldsymbol X^\prime \boldsymbol X$ 是非奇异的，则 $\boldsymbol X^\prime \boldsymbol X$ 的逆存在，并且 ${\boldsymbol \beta}$ 的最小二乘估计（方程 \@ref(eq:6-9) 中 $\hat {\boldsymbol \beta}$ 的解）为

\begin{equation}
\hat {\boldsymbol \beta}=(\boldsymbol X^\prime \boldsymbol X)^{-1}\boldsymbol{X^{\prime}y}
(#eq:6-10)
\end{equation}

当 $\boldsymbol X^\prime \boldsymbol X$ 是满秩时，最小二乘估计是唯一的。计算 $\boldsymbol X^\prime \boldsymbol X$ 的逆通常不是一项容易的任务。计算软件开发的最重要方面之一是，现在统计学家可以求超大型矩阵的逆，这在计算机出现之前是不可能的。然而，当 $\boldsymbol X^\prime \boldsymbol X$ 中存在某些模式时，可以利用这些模式更容易地求逆。对于单向均值模型和双向均值模型的正规方程组，$\boldsymbol X^\prime \boldsymbol X$ 是对角矩阵（所有对角元非零，非对角元为零），并且 $\boldsymbol X^\prime \boldsymbol X$ 的逆是通过简单地将每个对角元替换为其倒数来获得的。因此，单向均值模型的 $\mu_i$ 的最小二乘估计为

$$\begin{bmatrix}\hat{{\mu}}_1\\\hat{{\mu}}_2\\\vdots\\\hat{{\mu}}_t\end{bmatrix}=\begin{bmatrix}\frac1{n_1}&0&\cdots&0\\0&\frac1{n_2}&\cdots&0\\\vdots&\vdots&\ddots&\vdots\\0&0&\cdots&\frac1{n_t}\end{bmatrix}\begin{bmatrix}y_{1\cdot}\\y_{2\cdot}\\\vdots\\y_{t\cdot}\end{bmatrix}$$

或等价地

$$\hat{\mu}_i=\frac{y_{i\cdot}}{n_i}=\bar{y}_{i\cdot}\quad i=1,2,\ldots,t$$

类似地，双向均值模型 $\mu_{ij}$ 的最小二乘估计为

$$\hat{\mu}_{ij}=\frac{y_{ij\cdot}}{n_{ij}}=\bar{y}_{ij\cdot}\quad i=1,2,\ldots,t,\quad j=1,2,\ldots,b$$

与均值模型的正规方程组不同，效应模型的 $\boldsymbol X^\prime \boldsymbol X$ 是奇异的，并且 $\boldsymbol X^\prime \boldsymbol X$ 的逆不存在。在这种情况下，正规方程组有许多解（实际上是无穷多个最小二乘解）。效应模型被称为**过度指定** (overspecified models) 或**奇异模型** (sigular models)，因为与从收集的数据中可以唯一估计的参数相比，这些模型具有更多的参数。通常使用过度指定的模型，并且有几种方法可以求解其相应的正规方程组。以下讨论涉及完全随机设计结构中的双向处理结构，但类似的技术也可用于其他因素效应模型。

理论上，广义逆可用于求解关于 ${\boldsymbol \beta}$ 的正规方程组 (Graybill, 1976)，但求解过度指定模型的正规方程组的常用方法是对模型中的参数进行限制（这实际上产生了 g-inverse 解）。对模型参数施加限制可以通过多种方式实现，这里考虑其中两种方式。

### 零和限制 {#sec6-2-2}

一种常见的技术是要求某些参数的和等于零。这种方法自方差分析之初就被用来求解正规方程组。对于使用表 \@ref(tab:table6-1) 中数据的双向效应模型，**零和限制** (Sum-to-Zero Restrictions) 为：

$$\begin{aligned}
&\begin{aligned}\sum_{i=1}^3\tau_i=0,\sum_{j=1}^3\beta_j=0,\sum_{i=1}^3\gamma_{i1}=0,\sum_{i=1}^3\gamma_{i2}=0,\sum_{i=1}^3\gamma_{i3}=0\end{aligned} \\
&\sum_{j=1}^3\gamma_{1j}=0,\sum_{j=1}^3\gamma_{2j}=0,\,\text{and }\sum_{j=1}^3\gamma_{3j}=0
\end{aligned}$$

接下来，通过在考虑限制的情况下根据其他参数求解一些参数，然后将表达式替换回模型，将这些限制合并到模型中。例如，可以替换的参数有

$$\begin{aligned}\tau_3&=-\tau_1-\tau_2,\quad\beta_3=-\beta_1-\beta_2,\quad\gamma_{13}=-\gamma_{11}-\gamma_{12}\\\gamma_{23}&=-\gamma_{21}-\gamma_{22},\quad\gamma_{31}=-\gamma_{11}-\gamma_{21},\quad\gamma_{32}=-\gamma_{12}-\gamma_{22}\\\gamma_{33}&=-\gamma_{13}-\gamma_{23}=-\gamma_{31}-\gamma_{32}=\gamma_{11}+\gamma_{12}+\gamma_{21}+\gamma_{22}\end{aligned}$$

因此，替换模型中的 $\tau_3,~\beta_3,~\gamma_{13},~\gamma_{23},~\gamma_{33},~\gamma_{31}$ 和 $\gamma_{32}$，以获得重新参数化的模型

$$
\begin{bmatrix}y_{111}\\y_{112}\\y_{121}\\y_{131}\\y_{211}\\y_{221}\\y_{222}\\y_{231}\\y_{311}\\y_{321}\\y_{331}\end{bmatrix}=
\begin{bmatrix}1&1&0&1&0&1&0&0&0\\1&1&0&1&0&1&0&0&0\\1&1&0&0&1&0&1&0&0\\1&1&0&0&1&0&1&0&0\\1&1&0&-1&-1&-1&-1&0&0\\1&0&1&1&0&0&0&1&0\\1&0&1&0&1&0&0&0&1\\1&0&1&0&1&0&0&0&1\\1&0&1&-1&-1&0&0&-1&-1\\1&-1&-1&0&0&-1&0&-1&0\\1&-1&-1&0&1&0&-1&0&-1\\1&-1&-1&-1&-1&1&1&1&1\end{bmatrix}\begin{bmatrix}\mu^*\\\tau_1^*\\\tau_2^*\\\beta_1^*\\\beta_2^*\\\gamma_{11}^*\\\gamma_{12}^*\\\gamma_{21}^*\\\gamma_{22}^*\end{bmatrix}+\boldsymbol{\varepsilon}
$$

其表示为

$$\begin{aligned}\boldsymbol{y}=\boldsymbol{X}^*\boldsymbol{\beta}^*+\boldsymbol{\varepsilon}\end{aligned}$$

表 \@ref(tab:table6-1) 中数据对应的零和限制正态方程的解为 $\hat{\boldsymbol{\beta^*}}=(\boldsymbol X^{*\prime} \boldsymbol X^*)^{-1}\boldsymbol X^{*\prime}\boldsymbol{y}$. 可得到

$$\begin{aligned}
\hat{\boldsymbol \beta}^{*\prime}&=[\hat{\mu}^*,\hat{\tau}_1^*,\hat{\tau}_2^*,\hat{\beta}_1^*,\hat{\beta}_2^*,\hat{\gamma}_{11}^*,\hat{\gamma}_{12}^*,\hat{\gamma}_{21}^*,\hat{\gamma}_{22}^*] \\
&=[5.500,2.333,-0.833,-2.000,-0.5000,-1.333,1.667,-0.667,-0.167] 
\end{aligned}$$

$\boldsymbol \beta$ 剩余元素的估计从如下限制中获得：

$$\begin{aligned}
\hat{{\tau}}_3^*&=-\hat{{\tau}}_1^*-\hat{{\tau}}_2^*=-1.500,\quad\hat{{\beta}}_3^*=-\hat{{\beta}}_1^*-\hat{{\beta}}_2^*=2.500 \\
\hat{\gamma}_{13}^*&=-\hat{\gamma}_{11}^*-\hat{\gamma}_{12}^*=-0.333,\quad\hat{\gamma}_{23}^*=-\hat{\gamma}_{21}^*-\hat{\gamma}_{22}^*=0.833 \\
\hat{\gamma}_{31}^*&=-\hat{\gamma}_{11}^*-\hat{\gamma}_{21}^*=2.000,\quad\hat{\gamma}_{32}^*=-\hat{\gamma}_{12}^*-\hat{\gamma}_{22}^*=-1.500 \\
\hat{\gamma}_{33}^*&=\hat{\gamma}_{11}^*+\hat{\gamma}_{12}^*+\hat{\gamma}_{21}^*+\hat{\gamma}_{22}^*=-0.500
\end{aligned}$$

关于均值模型参数 $\mu_{ij}$，参数 $\mu^*,t_i^*,\beta_{j}^*,\gamma^*_{ij}$ 可以通过如下定义以满足零和限制

$$\mu^*=\bar{\mu}_{\cdot\cdot},\quad\tau_i^*=\bar{\mu}_{i\cdot}-\bar{\mu}_{\cdot\cdot\cdot},\quad\beta_j^*=\bar{\mu}_{\cdot j}-\bar{\mu}_{\cdot\cdot\cdot},\quad\mathrm{and}\quad\gamma_{ij}^*=\mu_{ij}-\bar{\mu}_{i\cdot}-\bar{\mu}_{\cdot j}+\bar{\mu}_{\cdot\cdot}$$

### 置零限制 {#sec6-2-3}

另一种常用于求解双向效应模型正规方程组的重新参数化技术为，将每组中的最后一个参数**设置为零** (Set-to-Zero)（选择最后一个参数是为了方便；也可以选择第一个、第二个或任何其他参数）。对于表 \@ref(tab:table6-1) 中数据的双向效应模型，限制为

$$\tau_3=0,\beta_3=0,\gamma_{13}=0,\gamma_{23}=0,\gamma_{33}=0,\gamma_{31}=0,\mathrm{~and~}\gamma_{23}=0$$

通过将上述限制条件纳入双向效应模型，得到重新参数化后的模型为

$$\begin{bmatrix}y_{111}\\y_{112}\\y_{121}\\y_{131}\\y_{211}\\y_{221}\\y_{222}\\y_{231}\\y_{311}\\y_{321}\\y_{331}\end{bmatrix}=\begin{bmatrix}1&1&0&1&0&1&0&0&0\\1&1&0&1&0&1&0&0&0\\1&1&0&0&1&0&1&0&0\\1&1&0&0&0&0&0&0&0\\1&0&1&1&0&0&0&1&0\\1&0&1&0&1&0&0&0&1\\1&0&1&0&1&0&0&0&1\\1&0&1&0&0&0&0&0&0\\1&0&0&1&0&0&0&0&0\\1&0&0&0&1&0&0&0&0\\1&0&0&0&0&0&0&0&0\end{bmatrix}\begin{bmatrix}\mu^+\\\tau_1^+\\\tau_2^+\\\beta_1^+\\\beta_2^+\\\gamma_{11}^+\\\gamma_{12}^+\\\gamma_{21}^+\\\gamma_{22}^+\end{bmatrix}+\boldsymbol\varepsilon  $$

可以表示为 $\boldsymbol{y}=\boldsymbol{X}^+\boldsymbol{\beta}^++\boldsymbol{\varepsilon}$。矩阵 $\boldsymbol{X}^+$ 是通过从双向效应模型完整设计矩阵 $\boldsymbol X$ 中删除与 $\tau_3,\beta_3,\gamma_{31},\gamma_{32},\gamma_{33},\gamma_{13},\gamma_{23}$ 对应的列而获得的。关于重新参数化模型设计矩阵的求解过程，置零限制比零和限制简单得多。置零限制对应的正规方程组的解为

$$\hat{\boldsymbol{\beta}}^+=(\boldsymbol X^{+\prime}\boldsymbol X^{+})^{-1}\boldsymbol X^{+\prime}\boldsymbol y$$

可以得到

$$\begin{aligned}
\hat{\boldsymbol{\beta}}^{+\prime}& =[\hat{{\mu}}^+,\hat{{\tau}}_1^+,\hat{{\tau}}_2^+,\hat{{\beta}}_1^+,\hat{{\beta}}_2^+,\hat{{\gamma}}_{11}^+,\hat{{\gamma}}_{12}^+,\hat{{\gamma}}_{21}^+,\hat{{\gamma}}_{22}^+]  \\
&=[6.0,4.0,2.0,-2.0,-4.0,-3.5,3.0,-4.0,0.0]
\end{aligned}$$

其余参数的估计是零，因为它们是由置零限制指定的，即，

$$\hat{\tau}_3^+=\hat{\beta}_3^+=\hat{\gamma}_{31}^+=\hat{\gamma}_{32}^+=\hat{\gamma}_{33}^+=\hat{\gamma}_{13}^+=\hat{\gamma}_{23}^+=0$$

为了将参数置零限制条件与均值模型参数 $\mu_{ij}$ 关联起来，定义 $\mu^+,t_i^+,\beta_j^+,\gamma_{ij}^+$

$$\begin{aligned}\mu^+=\mu_{tb},&\,\,\tau_i^+=\mu_{ib}-\mu_{tb},\,\,\beta_j^+=\mu_{tj}-\mu_{tb}\,\mathrm{~and~}\,\gamma_{ij}^+=\mu_{ij}-\mu_{tj}-\mu_{ib}+\mu_{tb}\end{aligned}$$

当 $\boldsymbol X^\prime \boldsymbol X$ 的不满秩时（即 $\boldsymbol X^\prime \boldsymbol X$ 是奇异的），正规方程组有多个可能的解。这是因为模型参数过多，即模型中的参数（在双向效应模型的情况下为 16 个）比可用数据（有 9 个数据单元，因此最多可以估计 9 个参数）能够唯一确定的参数要多。可以唯一估计的参数数量可被称为**基本** (essential) 参数数量。为了应对过度参数化模型以及非唯一的最小二乘解，必须考虑可估性的概念，这是 \@ref(sec6-3) 节的主题。下一个例子将给出单向处理结构效应模型的两种可能的解。

### 示例 6.2：单向处理结构 {#sec6-2-4}

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table6-2)第 \@ref(sec6-2-4) 节中均值和效应模型的单向处理结构数据</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%206.2.png">
</td>
  </tr>
</tbody>
</table>

这是一个单向处理结构的示例，在完全随机设计结构中具有四种处理。数据如表 \@ref(tab:table6-2) 所示。 $\boldsymbol X^*$ 矩阵是通过使用零和限制重新参数化模型来构造的；即，假设 $\tau_1^*+\tau_2^*+\tau_3^*+\tau_4^*=0$.

由此产生的设计矩阵为

$$\boldsymbol X^*=\begin{bmatrix}
1 & 1 & 0 & 0 \\
1 & 1 & 0 & 0 \\
1 & 1 & 0 & 0 \\
1 & 1 & 0 & 0 \\
1 & 1 & 0 & 0 \\
1 & 1 & 0 & 0 \\
1 & 1 & 0 & 0 \\
1 & 0 & 1 & 0 \\
1 & 0 & 1 & 0 \\
1 & 0 & 1 & 0 \\
1 & 0 & 1 & 0 \\
1 & 0 & 1 & 0 \\
1 & 0 & 0 & 1 \\
1 & 0 & 0 & 1 \\
1 & 0 & 0 & 1 \\
1 & 0 & 0 & 1 \\
1 & 0 & 0 & 1 \\
1 & 0 & 0 & 1 \\
1 & -1 & -1 & -1 \\
1 & -1 & -1 & -1 \\
1 & -1 & -1 & -1 \\
1 & -1 & -1 & -1 \\
1 & -1 & -1 & -1 \\
1 & -1 & -1 & -1 \\
1 & -1 & -1 & -1 \\
1 & -1 & -1 & -1 \\
\end{bmatrix}$$

零和限制模型的正规方程组为

$$\begin{bmatrix}26&-1&-3&-2\\-1&15&8&8\\-3&8&13&8\\-2&8&8&14\end{bmatrix}\begin{bmatrix}\hat{\mu}^*\\\hat{\tau}_1^*\\\hat{\tau}_2^*\\\hat{\tau}_3^*\end{bmatrix}=\begin{bmatrix}55.1\\-0.9\\-2.5\\-5.9\end{bmatrix}$$

零和限制模型的正规方程组相应最小二乘解为

$$\hat{\boldsymbol{\beta}}^*=\begin{bmatrix}\hat{{\mu}}^*\\\hat{{\tau}}_1^*\\\hat{{\tau}}_2^*\\\hat{{\tau}}_3^*\end{bmatrix}=\begin{bmatrix}2.1510\\0.0204\\0.5690\\-0.4510\end{bmatrix}$$

以及 $\hat{\tau}_4^*=-\hat{\tau}_1^*-\hat{\tau}_2^*-\hat{\tau}_3^*=-0.1384$. 就均值模型参数而言，零和限制参数可表示为

$$\mu^*=\bar{\mu}_{\cdot},\quad\tau_1^*=\mu_1-\bar{\mu}_{\cdot},\quad\tau_2^*=\mu_2-\bar{\mu}_{\cdot},\quad\tau_3^*=\mu_3-\bar{\mu}_{\cdot},\quad\tau_4^*=\mu_4-\bar{\mu}_{\cdot}$$

具有 $\tau_ 4 ^+ = 0$ 的置零限制设计矩阵 $\boldsymbol X^+$ 是通过在 $\boldsymbol X^*$ 中用“0”替换“-1”来构建的。置零限制模型的正规方程组为

$$\begin{bmatrix}26&7&5&6\\7&7&0&0\\5&0&5&0\\6&0&0&6\end{bmatrix}\begin{bmatrix}\hat{{\mu}}^+\\\hat{{\tau}}_1^+\\\hat{{\tau}}_2^+\\\hat{{\tau}}_3^+\end{bmatrix}=\begin{bmatrix}55.1\\15.2\\13.6\\10.2\end{bmatrix}$$

置零限制模型的最小二乘解为

$$\boldsymbol{\hat{\beta}}^+=\begin{bmatrix}{\hat{\mu}}^+\\{\hat{\tau}}_1^+\\{\hat{\tau}}_2^+\\{\hat{\tau}}_3^+\end{bmatrix}=\begin{bmatrix}2.0125\\0.1589\\0.7075\\-0.3125\end{bmatrix}$$

且 $\tau_ 4 ^+ = 0$. 就均值模型参数而言，置零限制参数可表示为

$$\mu^+=\mu_4,\quad\tau_1^+=\mu_1-\mu_4,\quad\tau_2^+=\mu_2-\mu_4,\quad\tau_3^+=\mu_3-\mu_4,\quad\tau_4^+=\mu_4-\mu_4=0$$

需要估计的最后一个参数是总体方差 $\sigma^2$. 基于 $\boldsymbol \beta$ 最小二乘解的 $\sigma^2$ 的估计为

\begin{align}
\hat{\sigma}^2& =\frac{1}{n-r}(\boldsymbol y-X\hat{\boldsymbol \beta})'(\boldsymbol y-\boldsymbol X\hat{\boldsymbol \beta})  \\
&=\frac1{n-r}\sum_{i=1}^n{(y_i-\hat{\beta}_0-\hat{\beta}_1x_{i1}-\hat{\beta}_2x_{i2}-\cdots-\hat{\beta}_{p-1}x_{ip-1})^2}
(#eq:6-11)
\end{align}

其中 $r=\operatorname{rank}(\boldsymbol X)$.

如果假设误差独立分布，且前四阶矩等于正态分布的前四阶矩，那么 $\hat \sigma^2$ 是 $\sigma^2$ 的最佳二次无偏估计 (best quadratic unbiased estimate). 如果误差也服从正态分布，那么 $\hat \sigma^2$ 是 $\sigma^2$ 的最佳无偏估计 (best unbiased estimate)，且 $(n - r)\hat \sigma^2 /\sigma^2$ 的抽样分布是自由度为 $n - r$ 的中心卡方分布。

## 可估性和连通的设计 {#sec6-3}

当过度指定模型或非满秩模型用于实验情况时，存在许多不同的最小二乘解（实际上有无限多个解）。如果两名研究人员分析上述两个数据集，一个使用零和限制，另一个使用置零限制，他们可能会得出不同的结论。对于上一节中的双向示例，$\tau_2^* = -0.833, \tau_2^+ = 2.000$；因此，一位研究人员可能会说 $\tau_2$ 很可能为负，而另一位研究人员可能会说 $\tau_2$ 很可能为正——**而这两种说法都是错误的**。

### 可估函数 {#sec6-3-1}

由于两位研究人员都在分析相同的数据集，因此他们似乎应该只考虑对两个重新参数化模型具有相同估计的参数或参数的函数。这种参数函数称为**参数的可估函数** (estimable functions of the parameters).

::: {.definition #6-1}
♦

参数或参数的函数 $f(\boldsymbol \beta)$ 是可估的，当且仅当参数或参数的函数的估计值对于最小二乘解的选择是不变的；也就是说，无论使用正规方程组的哪个解，估计值都是相同的。
:::

如果两位研究人员对参数的可估函数进行了估计，即使他们采用了两种不同的最小二乘解，他们也会得到相同的估计值。因此，他们对于参数的可估函数会得出相同的结论。对于矩阵模型，$\boldsymbol \beta$ 的线性可估函数表现为参数向量的线性组合，如 $\boldsymbol a^\prime \boldsymbol \beta$，其中 $\boldsymbol a^\prime$ 是一个 p × 1 的常数向量。线性函数 $\boldsymbol a^\prime \boldsymbol \beta$ 是可估的，当且仅当存在向量 $r$，使得 $\boldsymbol a =\boldsymbol X ^\prime \boldsymbol X r$. 每个形如 $\boldsymbol x^\prime_i \boldsymbol \beta$ 的函数都是可估的，其中 $\boldsymbol x^\prime_i$ 是 $\boldsymbol X$ 的第 $i$ 行。此外，$\boldsymbol x^\prime_i \boldsymbol \beta$ 的任何线性组合都是可估的函数。考虑 \@ref(sec6-2) 节中单向示例的两个解。由于每个参数 $\mu,\tau_1,\tau_2,\tau_3,\tau_4$ 都有两个不同的解，因此这些参数被认为是**不可估的** (nonestimable). 但是，通过计算每种方法关于 $\mu + \tau_i$ 的估计值，可以发现

$$\hat{\mu}^*+\hat{\tau}_i^*=\hat{\mu}^++\hat{\tau}_i^+,\quad i=1,2,3,4$$

证明 $\mu + \tau_i$ 是参数的可估函数。

$\tau_i$ 的所有对比，如差异 ${\tau}_1{-\tau}_2{,\tau}_1{-\tau}_3{,\tau}_2{-\tau}_3$ 或者 $\sum_{i=1}^tc_i\tau_i$ 其中 $\sum_{i=1}^tc_i=0$，可证明是单向模型的可估函数。

对于双向效应模型，一些可估函数为

$\mu+\tau_i+\beta_j+\gamma_{ij},\quad\gamma_{ij}-\gamma_{ij^\prime}-\gamma_{i^\prime j}+\gamma_{i^\prime j^\prime},\quad\beta_j-\beta_j+\bar{\gamma}_{\cdot j}-\bar{\gamma}_{\cdot j^\prime},\quad\tau_i-\tau_i+\bar{\gamma}_{i\cdot}-\bar{\gamma}_{i^\prime\cdot}$

可估函数在第 \@ref(chap10) 章进行了更详细的讨论。这里需记住的重点是可估函数的定义。从数据集中进行推断时，必须只考虑可估的参数函数，因为它们是参数的函数，其估计值不依赖于选择的最小二乘解。SAS^®^-GLM 和 SAS^®^-Mixed 都会检查所请求估计的参数或参数函数是否确实是可估的。如果参数或参数函数不可估，则不提供估计值。

### 连通性 {#sec6-3-2}

与可估函数相关的另一个概念是双向处理结构的**连通性** (connectedness). 如果可以假设行处理和列处理的水平之间没有交互作用，那么可以将处理组合平均值建模为

\begin{equation}
\mu_{ij}=\mu+\tau_i+\beta_j,\quad i=1,2,\ldots,b,j=1,2,\ldots,t
(#eq:6-12)
\end{equation}

如果双向处理结构中的数据以这样的方式出现在双向单元格中，即对于模型 \@ref(eq:6-12) 中的所有 $j\ne j^\prime$ 以及 $i\ne i^\prime$，$\beta_j-\beta_{j^\prime}$ 和 $\tau_i-\tau_{i^\prime}$ 都是可估的，则该双向处理结构被称为**连通的** (connected)。图 \@ref(fig:figure6-1) 中的安排 I 是一个连通的实验，而安排 II 不是。例如，使用安排 I，

$$\beta_1-\beta_2=(\mu+\tau_1+\beta_1)-(\mu+\tau_1+\beta_5)+(\mu+\tau_2+\beta_5)-(\mu+\tau_2+\beta_2)$$

我们可以得到单元格均值的线性组合，因此 $\beta_1-\beta_2$ 是可估的。对于安排 II，没有这样的单元格均值的线性组合可以得到 $\beta_1-\beta_2$，因此 $\beta_1-\beta_2$ 是不可估的。下一节将讨论关于参数的可估函数的假设检验。

<div class="figure" style="text-align: center">
<img src="figure/figure 6.1.png" alt="连通和非连通的双向处理结构" width="972" />
<p class="caption">(\#fig:figure6-1)连通和非连通的双向处理结构</p>
</div>

## 关于线性模型参数的检验假设 {#sec6-4}

对于线性模型中参数线性函数的假设检验，有多种方法可以开发适当的统计量。这里使用的方法用矩阵表示法，相当于条件误差原理（第 \@ref(chap1) 章）和似然比统计量。讨论仅限于检验关于参数的可估函数的假设。特别是，考虑检验假设

\begin{equation}
H_0{:\boldsymbol{H\beta}}=\boldsymbol{h}\mathrm{~vs~}H_a{:\boldsymbol{H\beta}}\neq\boldsymbol{h}
(#eq:6-13)
\end{equation}

其中线性组合 $\boldsymbol{H\beta}$ 是 $\boldsymbol \beta$ 的可估函数，$\boldsymbol H$ 是秩为 $q$ 的 q × p矩阵（即，$\boldsymbol H$ 的所有行是线性无关的）。相应的检验统计量为

\begin{equation}
F_c=\frac{SSH_0/q}{\hat{\sigma}^2}
(#eq:6-14)
\end{equation}

其中 $\hat \sigma^2$ 由式 \@ref(eq:6-11) 给出，并且

\begin{equation}
SSH_0=(\boldsymbol H\hat{\boldsymbol{\beta}}-\boldsymbol h)^{\prime}\left[\boldsymbol H(\boldsymbol X^{\prime}\boldsymbol X)^{-}\boldsymbol H^{\prime}\right]^{-1}(\boldsymbol H\hat{\boldsymbol{\beta}}-\boldsymbol h)
(#eq:6-15)
\end{equation}

这被称为**因偏离零假设而产生的平方和** (sum of squares due to deviations from the null hypothesis)（记号 "$(\boldsymbol X^{\prime}\boldsymbol X)^{-}$" 表示矩阵 $\boldsymbol X^{\prime}\boldsymbol X$ 的广义逆，Graybill, 1976）。在误差向量元素为 $i.i.d.\,N(0,\sigma^2)$ 的假设下，$F_c$ 服从自由度为 $q,n-r$ 的 $F$ 分布。

式 \@ref(eq:6-13) 中的假设总是可以用重新参数模型 $\boldsymbol{y}=\boldsymbol{X}^*\boldsymbol{\beta}^*+\boldsymbol \varepsilon$ 来等价地表述为下式，其中 $\boldsymbol X^{*\prime}\boldsymbol X^*$ 是非奇异的

$$H_0{:}H^*\boldsymbol{\beta}^*=\boldsymbol{h}^*\mathrm{~vs~}H_a{:}H^*\boldsymbol{\beta}^*\neq\boldsymbol{h}^*$$

则式 \@ref(eq:6-15) 的 $SSH_0$ 可计算为

$$SSH_0=(\boldsymbol H^*\hat{\boldsymbol{\beta}}^*-\boldsymbol{h}^*)^{\prime}[\boldsymbol H^*(\boldsymbol X^{*\prime}\boldsymbol X^*)^{-1}(\boldsymbol{H}^*{}^{\prime}\hat{\boldsymbol{\beta}}^*-\boldsymbol{h}^*)]$$

对于单向模型，在原始效应模型中检验假设

$$\begin{aligned}H_0\colon\tau_1=\tau_2=\cdots=\tau_t\text{ vs }H_a\colon\tau_i\neq\tau_{i'}\quad\text{for some }i\neq i'\end{aligned}$$

相当于在置零重新参数化模型中检验

$$\begin{aligned}H_0\colon\tau_1^*=\tau_2^*=&\cdots=\tau_{t-1}^*=0\text{ vs }H_a{:}\tau_i^*\neq0\quad\text{for some }i\leq t\text{ -1}\end{aligned}$$

对于示例 [6.2](#sec6-2-4) 关于 $\boldsymbol \beta^*$ 的零假设是

$$H_0{:}\begin{bmatrix}0&1&0&0\\0&0&1&0\\0&0&0&1\end{bmatrix}\begin{bmatrix}\mu^*\\\tau_1^*\\\tau_2^*\\\tau_3^*\end{bmatrix}=0\quad\mathrm{or}\quad \boldsymbol H^*\boldsymbol \beta^*=0$$

偏离 $H_0$ 而产生的平方和为

$$SSH_0=(\hat{\tau}_1^*\quad\hat{\tau}_1^*\quad\hat{\tau}_1^*)\,\boldsymbol{Z}^{-1}\begin{bmatrix}\hat{\tau}_1^*\\\hat{\tau}_2^*\\\hat{\tau}_3^*\end{bmatrix}$$

其中 $\boldsymbol{Z}=\boldsymbol{H}^*(\boldsymbol{X}^*\boldsymbol{'}\boldsymbol{X}^*)^{-1}\boldsymbol{H}^*\boldsymbol{'}$. 这是矩阵 $(\boldsymbol{X}^*\boldsymbol{'}\boldsymbol{X}^*)^{-1}$ 中与 $\tau_1,\tau_2,\tau_3)$ 相关的行和列的部分。

可估函数 $\boldsymbol a^\prime \boldsymbol \beta$ 的一个 $(1-\alpha)100\%$ 置信区间为

$$\boldsymbol a^{\prime}\hat{\boldsymbol \beta}-[t_{\alpha/2,n-p}]S_{\boldsymbol a^{\prime}\hat{\boldsymbol \beta}}\leq a^{\prime}\boldsymbol{\beta}\leq a^{\prime}\hat{\beta}+[t_{\alpha/2,n-p}]S_{\boldsymbol a^{\prime}\hat{\boldsymbol \beta}}$$

其中 $S_{\boldsymbol a'\boldsymbol{\beta}}^2=\hat{\sigma}^2\boldsymbol{a'}(\boldsymbol{X'}\boldsymbol{X})^{-1}\boldsymbol{a}$. 通过使用第 \@ref(chap3) 章中讨论的多重比较程序，可以构造关于几个可估函数的同时置信区间 (simultaneous confidence intervals).

## 总体边际均值 {#sec6-5}

在通过方差分析对交叉分类数据集进行分析后，实验者通常对估计特定效应或单元格均值感兴趣。**总体边际均值** (population marginal means) 被定义为在特定类别上平均化的参数线性组合，就好像每个单元格中都有一个观测值一样 (Searle et al., 1980). 如果每个单元格中至少有一个观测值，那么所有的总体边际均值都是可估的，而如果某些单元格为空，则它们不一定可估。这个定义不依赖于单元格中的样本量。如果数据代表单元格的比例抽样，那么实验者可能希望考虑单元格均值的加权平均值，其中权重由样本量给出（参见第 \@ref(chap10) 章）。或者，可能需要使用其他加权方案来对各个单元格进行平均。本节末尾提供了一个示例，以演示不同的可能性。

对于完全随机设计的单向处理结构，第 $i$ 个处理的总体边际平均值为 $\mu+\tau_i=\mu_i$，并由下式估计

$$\widehat{\mu+\tau_i}=\hat{\mu}+\hat{\tau_i}$$

其中 $\hat \mu,\hat \tau_i$ 是从正规方程组的任何解得到的。这些估计值被称为**总体边际均值估计** (estimated population marginal means).

对于完全随机设计中的双向处理结构，第 $(i,j)$ 个单元格的总体边际平均值为 $\begin{aligned}\mu_{ij}=\mu+\tau_i+\beta_j+\gamma_{ij}\end{aligned}$. 第 $i$ 行的总体边际均值是该行中 $\mu_{ij}$ 的平均值，或者

$$\bar{\mu}_{i\cdot}=\sum_{j=1}^b\frac{\mu_{ij}}b=\mu+\tau_i+\bar{\beta}_{\cdot}+\bar{\gamma}_{i\cdot}$$

第 $j$ 列的总体边际均值是该列中 $\mu_{ij}$ 的平均值，或者

$$\bar{{\mu}}_{\cdot j}=\sum_{i=1}^t\frac{{\mu}_{ij}}t={\mu}+\bar{{\tau}}_{\cdot}+{\beta}_j+\bar{{\gamma}}_{\cdot j}$$

总体边际均值估计分别为

$$\begin{aligned}
\hat{{\mu}}_{ij} &=\hat{\mu}+\hat{\tau}_i+\hat{\beta}_j+\hat{\gamma}_{ij},  \\
\hat{\bar{{\mu}}}_{i\cdot} &=\hat{\mu}+\hat{\tau}_i+\sum_{j=1}^b\frac{\hat{\beta}_j}b+\sum_{j=1}^b\frac{\hat{\gamma}_{ij}}b,\\
\hat{\bar{{\mu}}}_{\cdot j}&=\hat{{\mu}}+\sum_{i=1}^t\frac{\hat{{\tau}}_i}t+\hat{{\beta}}_j+\sum_{i=1}^t\frac{\hat{{\gamma}}_{ij}}t
\end{aligned}$$

这些估计对于任何可能的最小二乘解都是唯一的

$$\hat{\mu},\hat{\tau}_1,\hat{\tau}_2,\ldots,\hat{\tau}_t,\hat{\beta}_1,\hat{\beta}_2,\ldots,\hat{\beta}_b,\hat{\gamma}_{11},\hat{\gamma}_{12},\ldots,\hat{\gamma}_{tb}$$

根据零和限制模型，表 \@ref(tab:table6-1) 中双向示例的总体边际平均值 $\hat{\bar\mu}_{1\cdot}$ 的估计为

$$\begin{aligned}
\hat{\bar{\mu}}_{1\cdot}& =\hat{\mu}^*+\hat{\tau}_1^*+\frac{\hat{\beta}_1^*+\hat{\beta}_2^*+\hat{\beta}_3^*}3+\frac{\hat{\gamma}_{11}^*+\hat{\gamma}_{12}^*+\hat{\gamma}_{13}^*}3  \\
&=5.500+2.333+\frac{-2.000-0.500+2.500}3+\frac{-1.333+1.667-0.333}3 \\
&=7.833
\end{aligned}$$

当按预期使用置零限制的最小二乘解时，可以获得与 $\hat{\bar\mu}_{1\cdot}$ 相同的值，因为 ${\bar\mu}_{1\cdot}$ 在这个例子中是可估的。

当没有空单元格时，所有的总体边际均值都是可估的。如果有空单元格，那么任何涉及一个或多个缺失单元格的总体边际均值都是不可估的。例如，如果 $(2,2)$ 单元格在 2×2 处理结构中缺失，则没有关于 $\mu_{22}$ 的信息，因此 $\mu_{22}$ 是不可估的。第 2 列的总体边际均值为

$$\bar{\mu}_{\cdot 2}=\frac{\mu_{12}+\mu_{22}}2$$

由于 $\bar{\mu}_{\cdot 2}$ 依赖于 $\mu_{22}$，因此 $\bar{\mu}_{\cdot 2}$ 是不可估的（$\bar{\mu}_{2\cdot }$ 也不可估）。

**任何可估的总体边际均值都可以表示为重新参数化模型参数向量元素的线性组合**；也就是说，对于 $\boldsymbol a$ 的适当选择，总体边际均值可以表示为 $\boldsymbol a^\prime\boldsymbol \beta^*$. 总体边际均值估计的方差为

$$\mathrm{Var}(\boldsymbol a^{\prime}\hat{\boldsymbol{\beta}^*})=\sigma^2 \boldsymbol a^{\prime}(\boldsymbol X^{*\prime}\boldsymbol X^*)^{-1}\boldsymbol a$$

且 $\boldsymbol a^\prime\boldsymbol \beta^*$ 的标准误估计为

$$\widehat{s.e.}(\boldsymbol a^{\prime}\hat{\boldsymbol \beta}^*)=\hat{\sigma}\sqrt{\boldsymbol a^{\prime}(\boldsymbol X^{*\prime} \boldsymbol X^*)^{-1}\boldsymbol a}$$

对于置零限制的单向模型，总体边际均值估计为

$$\widehat{\mu+\tau_i}=\hat{\mu}^*+\hat{\tau}_i^*,\quad i=1,2,\ldots,t-1$$

且

$$\widehat{{\mu}+{\tau}_t}=\hat{{\mu}}^*$$

这些总体边际均值估计的方差为

$$\mathrm{Var}(\widehat{\mu+\tau_i})=\mathrm{Var}(\hat{\mu}^*)+2\mathrm{Cov}(\hat{\mu}^*,\hat{\tau}_i^*)+\mathrm{Var}(\hat{\tau}_i^*)\quad i=1,2,\ldots,t-1$$

且

$$\mathrm{Var}(\widehat{\mu+\tau_t})=\mathrm{Var}(\hat{\mu}^*)$$

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table6-3)双向处理结构的平均绩点数据</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%206.3.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table6-3) 中的数据是一所学校学生样本的平均绩点 (GPA)，学生按学年和性别进行分类。双向模型

$$y_{ijk}=\mu+\tau_i+\beta_j+\gamma_{ij}+\varepsilon_{ijk}\quad i=1,2,3,4,\quad j=1,2,\quad k=1,2,\ldots,n_{ij}$$

用于描述数据，其中 $y_{ijk}$ 是 GPA 值，$\tau_i$ 是第 $i$ 年的效应，$\beta_j$ 是第 $j$ 个性别的效应，$\gamma_{ij}$ 是交互作用效应。参数的置零最小二乘解如表 \@ref(tab:table6-4) 所示。对于那些置零的参数，标准误估计为零。第 $1$ 学年和女性的单元格均值估计为

$$\hat{\mu}_{1f}=\hat{\mu}+\hat{\tau}_1+\hat{\beta}_f+\hat{\gamma}_{1f}=3.333-0.133-0.033+0.205=3.372$$

其余单元格均值可类似计算，如表 \@ref(tab:table6-6) 所示。表 \@ref(tab:table6-6) 和 \@ref(tab:table6-7) 提供了边际年份和边际性别的均值估计。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table6-4)双向处理结构中 GPA 数据的最小二乘解</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%206.4.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table6-5)单元格均值及标准误估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%206.5.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table6-6)跨性别水平的年份原始均值、最小二乘均值和加权平均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%206.6.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table6-7)跨年份水平的女性和男性平均原始均值、最小二乘均值和加权平均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%206.7.png">
</td>
  </tr>
</tbody>
</table>

**原始均值** (raw measn) 是各个单元格均值的加权平均值，使用单元格中的观测数量作为权重。例如，第 1 年的原始均值或**未调整的均值** (unadjusted mean) 计算为

$$\hat{\bar{\mu}}_{1\cdot}=\frac{7^*3.371+5^*3.200}{7+5}=3.300$$

年份和性别效应的**最小二乘均值** (least square means) 是相应单元格均值的未加权平均值，第 1 年的最小二乘均值为

$$\hat{\bar{\mu}}_{1\cdot}=\frac{3.371+3.200}2=3.286$$

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table6-8)平均成绩研究中的学生在班级中的分布</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%206.8.png">
</td>
  </tr>
</tbody>
</table>

在论述边际年份或边际性别效应时，人们可能想要使用原始均值（根据观察到的单元格大小进行加权）或最小二乘均值（每个单元格均值具有相同的权重，就好像单元格大小相同一样）。对于设计的实验，最小二乘均值可能是人们感兴趣的均值，因为您可能设计了一个每个单元格中观测数量相等的实验，因此，即使某些数据可能缺失，提供假设单元格大小相等的边际均值估计也是一个合理的解决方案。但是，如果数据来自观察性研究，则未加权的均值可能不是人们感兴趣的边际均值。如果数据来自总体的简单随机样本，那么单元格大小可能反映了总体的比例成员结构 (proportional membership structure). 如果情况确实如此，那么原始均值或按样本量加权的均值就是感兴趣的边际均值。但如果单元格中的样本量不代表总体结构，那么最小二乘均值和原始均值都不是感兴趣的指标。当总体中的某些部分因设计或偶然原因被过度采样或采样不足时，就会出现这种现象。当已知总体结构，即每个单元格的总体比例已知时，通过使用已知的总体比例作为权重，可以获得感兴趣的边际均值。例如，假设研究每年中男性和女性的比例如表 \@ref(tab:table6-8) 所示。使用权重估计的第 1 年的边际均值为：

$$\hat{\tilde{\mu}}_{1\cdot}=\frac{11^*3.371+22^*3.200}{33}=3.257$$

使用权重对女性所有年份的边际均值估计为

$$\hat{\tilde{\mu}}_{\cdot f}=\frac{11^*3.371+12^*3.180+10^*3.475+8^*3.300}{41}=3.327$$

使用权重对男性所有年份的边际均值估计为

$$\hat{\tilde{\mu}}_{\cdot m}=\frac{22^*3.200+14^*3.300+13^*3.267+10^*3.333}{59}=3.261$$

计算几种调整的均值并从中选择你喜欢的值并不是问题，但有必要尽可能准确地指定总体结构中的比例，然后使用这些权重计算的边际平均估计。这可以通过允许使用 `estimate` 语句的软件包来完成，其中提供了所需边际均值的估计和标准误的估计。表 \@ref(tab:table6-9) 包含使用 SAS-GLM 程序的代码，以表 \@ref(tab:table6-8) 中的比例作为权重提供边际均值的估计。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table6-9)带 `estimate` 语句的SAS^®^-GLM代码，使用表 \@ref(tab:table6-8) 中的权重提供性别和年份水平的总体边际均值估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%206.9.png">
</td>
  </tr>
</tbody>
</table>

## 结束语 {#sec6-6}

本章仅供对理论背景感兴趣的人阅读，介绍了最小二乘估计程序，并讨论了可估性的重要概念。还讨论了总体边际均值的定义和估计。本章为那些想开发供自己使用的统计软件的人提供了通用公式。

## 练习 {#sec6-7}

