
# 模型的矩阵形式 {#chap6}

当使用不均衡 (unbalanced) 的固定效应模型、随机效应模型或混合效应模型时，求和符号变得非常费力，有时甚至几乎无法使用。这个问题可以通过使用模型的矩阵形式表达来解决。本章讨论了模型矩阵形式的构建，并描述了如何使用矩阵获得最小二乘估计、检验假设、计算最小二乘或总体边际均值以及构建置信区间。在 \@ref(sec6-3) 节中讨论了可估性 (estimability) 的概念。

## 基本符号 {#sec6-1}

模型的矩阵形式可以表示为

\begin{equation}
\underset{n\times1}{\boldsymbol y}=\underset{n × p}{\vphantom{\boldsymbol y}\boldsymbol X} \,\,\underset{p × 1}{\vphantom{\boldsymbol y}\boldsymbol \beta}+\underset{n × 1}{\vphantom{\boldsymbol y}\boldsymbol \varepsilon}
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

**均值模型** (means model) 可表示为

$$y_{ij}={\mu}_1x_{1ij}+{\mu}_2x_{2ij}+\cdots+{\mu}_tx_{tij}+{\varepsilon}_{ij}$$

其中 $i=1,2,\cdots,t$，$j=1,2,\cdots,n_i$. 或者以矩阵表示法

$$\begin{bmatrix}{y}_{11}\\{y}_{12}\\\vdots\\{y}_{1n_1}\\{y}_{22}\\\vdots\\{y}_{2{n}_2}\\\vdots\\{y}_{t1}\\\vdots\\{y}_{t{n}_t}\end{bmatrix}=\begin{bmatrix}1&0&\cdots&0\\1&0&\cdots&0\\\vdots&\vdots&\vdots&\vdots\\1&0&\cdots&0\\0&1&\cdots&0\\0&1&\cdots&0\\\vdots&\vdots&\vdots&\vdots\\0&1&\cdots&0\\\vdots&\vdots&\vdots&\vdots\\0&0&\cdots&1\\\vdots&\vdots&\ddots&\vdots\\0&0&\cdots&1\end{bmatrix}\begin{bmatrix}\mu_1\\\mu_2\\\vdots\\\mu_t\end{bmatrix}+\boldsymbol \varepsilon$$

均值模型通常表示为 $y_{ij}=\mu_i+\varepsilon_{ij}$，其中 $i=1,2,\cdots,t$，$j=1,2,\cdots,n_i$.

**效应模型** (effect model) 可表示为 $y_{ij}=\mu+\tau_1x_{1ij}+\tau_2x_{2ij}+\cdots+\tau_tx_{tij}+\varepsilon_{ij}$，其中 $i=1,2,\cdots,t$，$j=1,2,\cdots,n_i$. 或者以矩阵表示法


$$\begin{bmatrix}{y}_{11}\\{y}_{12}\\\vdots\\{y}_{1{n}_1}\\{y}_{21}\\{y}_{22}\\\vdots\\{y}_{2{n}_2}\\\vdots\\{y}_{t1}\\\vdots\\{y}_{t{n}_t}\end{bmatrix}=\begin{bmatrix}1&1&0&\cdots&0\\1&1&0&\cdots&0\\\vdots&\vdots&\vdots&\vdots&\vdots\\1&1&0&\cdots&0\\1&0&1&\cdots&0\\1&0&1&\cdots&0\\\vdots&\vdots&\vdots&\vdots&\vdots\\1&0&1&\cdots&0\\\vdots&\vdots&\vdots&\vdots&\vdots\\1&0&0&\cdots&1\\\vdots&\vdots&\vdots&\ddots&\vdots\\1&0&0&\cdots&1\end{bmatrix}\begin{bmatrix}\mu\\\tau_1\\\tau_2\\\vdots\\\tau_t\end{bmatrix}+\boldsymbol \varepsilon $$

效应模型通常表示为 $y_{ij}=\mu+\tau_i+\varepsilon_{ij}$，其中 $i=1,2,\cdots,t$，$j=1,2,\cdots,n_i$. 

均值模型和效应模型之间的区别在于，效应模型的设计矩阵包含模型截距 $\mu$ 的全 1 列，而均值模型不包含全 1 列。

### 双向处理结构模型 {#sec6-1-3}



