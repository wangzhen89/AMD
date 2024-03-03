
# 混合模型的分析 {#chap22}

混合模型用于刻画来自实验或研究的数据，这些数据需要多个方差-协方差参数，并涉及一些固定效应参数。第 \@ref(chap2) 章的不等方差模型是混合模型，因为它们涉及多个方差分量。第 \@ref(chap18) - \@ref(chap21) 章描述的模型称为随机效应模型，但每个模型都有一个未知的均值参数，因此这些模型本质上是混合模型。第 \@ref(chap18) 章中使用的混合模型的定义围绕着具有固定效应和多个方差分量的一些处理结构。因此，**混合模型的一般定义是具有一些固定效应参数和协方差结构中多于一个参数的模型**。符合此定义的模型包括随机完全区组模型 (randomized complete blocks models)、不完全区组模型 (incomplete blocks models)、裂区型模型 (split-plot-type models)、裂条区型模型 (strip-plot-type models)、重复测量型模型 (repeated measures type models)、随机系数模型 (random coefficients models)、多水平模型 (multilevel models) 和分层模型 (hierarchical models).

混合模型包含三部分：1）模型的固定效应部分；2）模型的随机效应部分；3）模型的残差部分。因此，混合模型的分析包括两种类型的分析，一种是模型随机效应和残差部分的分析，另一种是模型固定效应部分的分析。本章讨论混合模型的构建以及分析模型随机效应部分和固定效应部分所需的步骤。本章的结果为理论结果和概念理解之间提供了桥梁。

## 混合模型简介 {#sec22-1}

描述同时具有随机效应因子和固定效应因子的实验的模型称为混合模型。由于存在两种类型的因子，因此所得到的模型具有两个部分，随机效应部分和固定效应部分。为了构造这样的模型，使用以下规则来确定特定交互作用是随机效应还是固定效应。

:::: {.rmdnote data-latex="{规则}"}
**规则**：如果主效应是随机效应，则涉及该主效应的任何交互作用也是随机效应。只有那些所有对应的主效应都是固定效应的交互作用才是固定效应。
::::

例如，三因素处理结构，其中 A 和 B 的水平为固定效应，C 的水平为随机效应，模型为

$$y_{ijkm}=\mu+\alpha_i+\beta_j+\gamma_{ij}+c_k+d_{ik}+f_{jk}+g_{ijk}+\varepsilon_{ijkm}$$

其中 $\mu$ 表示平均响应，$\alpha_i$ 表示固定因子 A 第 i 个水平的效应，$\beta_j$ 表示固定因子 B 第 j 个水平的效应，$\gamma_{ij}$ 表示 A 水平和 B 水平之间交互作用，$c_k$ 表示随机因子 C 第 k 个水平的效应，$d_{ik}$ 表示 A 水平和 C 水平之间的交互作用，$f_{jk}$ 表示 B 水平和 C 水平之间的交互作用，$g_{ijk}$ 表示 A, B, C 水平之间的三向交互作用，$\varepsilon_{ijkm}$ 表示残差效应。

模型的固定效应部分为 $\mu+\alpha_i+\beta_j+\gamma_{ij}$，随机效应部分为 $c_k+d_{ik}+f_{jk}+g_{ijk}$，残差效应部分为 $\varepsilon_{ijkm}$. 关于随机效应项和模型残差部分的简单假设如下 $c_k\thicksim i.i.d.N(0,\sigma_c^2),d_{ik}\thicksim i.i.d.N(0,\sigma_d^2),f_{jk}\sim i.i.d.N(0,\sigma_f^2),g_{ijk}\sim i.i.d.N(0,\sigma_g^2),\varepsilon_{ijkm}\sim i.i.d.N(0,\sigma_\varepsilon^2)$ 并且 $c_k,d_{ik},f_{jk},g_{ijk},\varepsilon_{ijkm}$ 为独立随机变量。A 和 B 之间的交互作用（由 $\gamma_{ij}$ 表示）是固定效应，而所有其他交互作用都是随机效应，因为它们都涉及由下标 k 表示的随机因子 C. 通常，混合模型还将涉及与设计结构对应的项，但这些项不包括在上述模型中。

一般线性混合模型的矩阵表示是

$$\boldsymbol y=\boldsymbol X\boldsymbol{\beta}+\boldsymbol Z_1\boldsymbol u_1+\boldsymbol Z_2\boldsymbol u_2+\cdots+\boldsymbol Z_k\boldsymbol u_k+\boldsymbol{\varepsilon}$$

其中 $\boldsymbol y$ 为 N × 1 数据向量，$\boldsymbol X\boldsymbol{\beta}$ 为模型的固定效应部分，$\boldsymbol Z_1\boldsymbol u_1+\boldsymbol Z_2\boldsymbol u_2+\cdots+\boldsymbol Z_k\boldsymbol u_k$ 为模型的随机效应部分，以及 $\boldsymbol{\varepsilon}$ 为模型的残差部分。理想条件为 $\boldsymbol {u}_{i}\thicksim N(\boldsymbol0,{\sigma}_{i}^{2}\boldsymbol{I}_{n_{i}}),i=1,2,\ldots,k;\boldsymbol{\varepsilon}\thicksim N(\boldsymbol{0},{\sigma}_{\varepsilon}^{2}\boldsymbol{I}_{N})$ 且 $\boldsymbol{u}_i\left(i=1,2,\ldots,k\right),\boldsymbol \varepsilon$ 为独立随机变量。

给定 $\boldsymbol u_1,\boldsymbol u_2,\ldots,\boldsymbol u_k$ 的 $\boldsymbol y$ 的条件分布由固定效应模型表示：$\boldsymbol y=\boldsymbol X\boldsymbol{\beta}+\boldsymbol Z_1\boldsymbol u_1+\boldsymbol Z_2\boldsymbol u_2+\cdots+\boldsymbol Z_k\boldsymbol u_k+\boldsymbol{\varepsilon}$ 其中 $\boldsymbol{\varepsilon}\thicksim N(\boldsymbol{0},{\sigma}_{\varepsilon}^{2}\boldsymbol{I}_{N})$ 且其中的 $\boldsymbol u_1,\boldsymbol u_2,\ldots,\boldsymbol u_k$ 为固定效应（因为取了条件）。$\boldsymbol{y}$ 的边际分布为 $\boldsymbol{y}\thicksim N(\boldsymbol{X\beta},\boldsymbol{\Sigma})$，或 $\boldsymbol{\Sigma}={\sigma}_1^2\boldsymbol{Z}_1\boldsymbol{Z}_1'+{\sigma}_2^2\boldsymbol{Z}_2\boldsymbol{Z}_2'+\cdots+{\sigma}_k^2\boldsymbol{Z}_k\boldsymbol{Z}_k'+{\sigma}_\varepsilon^2\boldsymbol{I}_N$，或等价地 $\boldsymbol{y}=\boldsymbol{X\beta}+\boldsymbol{e}$ 其中 $\boldsymbol{e}\thicksim N(\boldsymbol{0},\boldsymbol{\Sigma})$. 

混合模型的总体参数为 $\boldsymbol{\beta},{\sigma}_1^2,{\sigma}_2^2,\ldots,{\sigma}_k^2$ 和 ${\sigma}_\varepsilon^2$. 模型随机效应的分析部分包括方差分量 ${\sigma}_1^2,{\sigma}_2^2,\ldots,{\sigma}_k^2$ 和 ${\sigma}_\varepsilon^2$ 的估计、检验假设和置信区间的构建。模型固定效应部分的分析包括 $\boldsymbol \beta$ 的可估函数的 估计、检验假设和置信区间的构建。 $\boldsymbol y$ 的这种边际分布基于这样的假设：随机效应的协方矩阵是单位阵的标量倍数。混合效应模型更一般的形式是 $\boldsymbol y=\boldsymbol X\boldsymbol{\beta}+\boldsymbol Z_1\boldsymbol u_1+\boldsymbol Z_2\boldsymbol u_2+\cdots+\boldsymbol Z_k\boldsymbol u_k+\boldsymbol{\varepsilon}$，其中 $\boldsymbol u_i\sim N(\mathbf{0},\boldsymbol{G}_i),i=1,2,\ldots,k,\boldsymbol{\varepsilon}\thicksim N(\mathbf{0},\boldsymbol{R})$，且 $\boldsymbol u_i\,(i=1,2,\ldots,k)$ 和 $\boldsymbol{\varepsilon}$ 为独立随机变量。$\boldsymbol y$ 的边际分布可以用更一般的术语表示为 $\boldsymbol{y}\thicksim N(\boldsymbol{X\beta},\boldsymbol{\Sigma})$ 其中 
$\boldsymbol{\Sigma}=\boldsymbol{Z}_1\boldsymbol{G}_1\boldsymbol{Z}_1^{\prime}+\boldsymbol{Z}_2\boldsymbol{G}_2\boldsymbol{Z}_2^{\prime}+\cdots+\boldsymbol{Z}_k\boldsymbol{G}_k\boldsymbol{Z}_k^{\prime}+\boldsymbol{R}$，或等价地 $\boldsymbol{y}=\boldsymbol{X\beta}+\boldsymbol{e}$ 其中 $\boldsymbol{e}\thicksim N(\boldsymbol{0},\boldsymbol{\Sigma})$. 

理想情况下，矩阵 $\boldsymbol G_1,\boldsymbol G_1,\ldots,\boldsymbol G_k$ 和 $\boldsymbol R$ 是正定的，并且包含对随机效应和模型残差部分的协方差结构进行建模所需的参数。混合模型的分析将在接下来的两节中介绍，每节对应一个分析部分。首先考虑模型的随机效应和残差部分的分析。

## 混合模型随机效应部分的分析 {#sec22-2}

通过首先拟合模型的固定效应部分，然后计算该模型的残差，可以从一般线性混合模型构建随机效应模型。生成的模型不依赖于模型的固定效应部分。表示为 $\boldsymbol y$ 边际分布的一般线性混合模型是

$$\boldsymbol y=\boldsymbol X\boldsymbol{\beta}+\boldsymbol e\quad\mathrm{where~}\boldsymbol e\thicksim N(\boldsymbol{0},\boldsymbol{\Sigma})\mathrm{~and~}\boldsymbol{\Sigma}=\sigma_1^2\boldsymbol{Z}_1\boldsymbol{Z}_1^{\prime}+\sigma_2^2\boldsymbol{Z}_2\boldsymbol{Z}_2^{\prime}+\cdots+\sigma_k^2\boldsymbol{Z}_k\boldsymbol{Z}_k^{\prime}+\sigma_\varepsilon^2\boldsymbol{I}_N$$

$\boldsymbol \beta$ 的（普通）最小二乘估计为 $\boldsymbol{\hat{\beta}}=(\boldsymbol X^{\prime}\boldsymbol{X})^-\boldsymbol{X^{\prime}}\boldsymbol{y}$ 其中 $(\boldsymbol X^{\prime}\boldsymbol{X})^-$ 表示 $\boldsymbol X^{\prime}\boldsymbol{X}$ 的 Moore–Penrose 广义逆。残差向量为

$$\boldsymbol r=\boldsymbol y-\boldsymbol X\hat{\boldsymbol{\beta}}=(\boldsymbol I-\boldsymbol X\boldsymbol X^-)\boldsymbol y$$

残差向量的模型（称为残差模型）为

$$\boldsymbol r=(\boldsymbol I-\boldsymbol X\boldsymbol X^-)\boldsymbol Z_1\boldsymbol u_1+(\boldsymbol I-\boldsymbol X\boldsymbol X^-)\boldsymbol Z_2\boldsymbol u_2+\cdots+(\boldsymbol I-\boldsymbol X\boldsymbol X^-)\boldsymbol Z_k\boldsymbol u_k+(\boldsymbol I-\boldsymbol X\boldsymbol X^-)\boldsymbol{\varepsilon}$$

残差模型不依赖于固定效应参数或模型的固定效应部分 $\boldsymbol {Xb}$，因此残差模型是随机效应模型。第 \@ref(chap18) - \@ref(chap20) 章中讨论的随机效应模型分析方法可用于分析残差随机效应模型。矩法、最大似然法、REML 法和 MINQUE 法这四种技术是接下来四小节的主题。第 \@ref(chap23) 章通过两个例子演示了这些方法，一个是均衡数据集，一个是不均衡数据集。

### 矩法 {#sec22-2-1}

正如第 \@ref(chap19) 章中描述的随机效应模型分析，矩法技术需要计算平方和，确定其期望，然后从方程组中估计方差分量，该方程组通过将观测到的均方与其期望值相等而获得。为了估计方差分量，必须获得其期望值不依赖于模型的混合效应部分的平方和。第 \@ref(sec18-3) 节讨论的拟合常数的方法提供了这样的平方和：首先拟合模型的固定效应部分，然后拟合随机效应部分。

为了演示混合模型随机效应部分的分析，考虑完全随机设计结构中的双向处理结构，其中一个混合因子用 B 表示，另一个随机因子用 T 表示。得到的双向混合模型为

$$y_{ijk}=\mu+\beta_i+t_j+g_{ij}+\varepsilon_{ijk}\quad i=1,2,\ldots,b\quad j=1,2,\ldots,t\quad k=1,2,\ldots,n_{ij}$$

其中假设 $t_i \thicksim i.i.d.N(0,\sigma_t^2),g_{ij} \thicksim i.i.d.N(0,\sigma_g^2)$ 以及 $\varepsilon_{ijk}\thicksim i.i.d.N(0,\sigma_\varepsilon^2)$. 从拟合常数方法获得的平方和为 $R(\beta|\mu),R(t|\mu,\beta),R(g|\mu,\beta,t)$ 以及 $SSERROR$. 对于某些常数 $k_1,k_2,k_3$，最后三个均方的期望值分别具有以下形式

$$\begin{aligned}E[MSR(t|\mu,\beta)]&=\sigma_\varepsilon^2+k_1\sigma_g^2+k_2\sigma_\varepsilon^2\\E[MSR(g|\mu,\beta,t)]&=\sigma_\varepsilon^2+k_3\sigma_g^2\\E[MSERROR]&=\sigma_\varepsilon^2\end{aligned}$$

$k_1,k_2,k_3$ 的值将取决于样本量和数据结构。这些均方期望不涉及混合效应参数，因此均方可用于估计方差分量，也可以用于检验关于它们的假设。

方程组是通过将观测到的均方等于期望均方来构建的，其中 $\tilde{{\sigma}}_{{\varepsilon}}^2,\tilde{{\sigma}}_t^2$ 和 $\tilde{{\sigma}}_g^2$ 表示解；即，

$$\begin{aligned}
MSR(t|\mu,\beta)& =\tilde{{\sigma}}_\varepsilon^2+k_1\tilde{{\sigma}}_g^2+k_2\tilde{{\sigma}}_\varepsilon^2  \\
MSR(g|\mu,\beta,t)& =\tilde{\sigma}_\varepsilon^2+k_3\tilde{\sigma}_g^2  \\
MSERROR& =\tilde{{\sigma}}_\varepsilon^2 
\end{aligned}$$

解为

$$\begin{aligned}
\tilde{{\sigma}}_\varepsilon^2 & =MSERROR  \\
\tilde{{\sigma}}_g^2& =[MSR(g|\mu,\beta,t)-\tilde{\sigma}_\varepsilon^2]/k_3  \\
\tilde{{\sigma}}_t^2&=[MSR(t|\mu,\beta)-k_1\tilde{\sigma}_g^2-\tilde{\sigma}_\varepsilon^2]/k_2
\end{aligned}$$

方差分量的估计为

$$\begin{aligned}
&\hat{{\sigma}}_\varepsilon^2 =\tilde{\sigma}_\varepsilon^2  \\
&\hat{{\sigma}}_g^2 =\begin{cases}\tilde{{\sigma}}_g^2&\mathrm{~if~}\tilde{{\sigma}}_g^2\geq0\\0&\mathrm{~if~}\tilde{{\sigma}}_g^2<0&\end{cases}  \\
&\hat{{\sigma}}_t^2 =\begin{cases}\tilde{{\sigma}}_t^2&\mathrm{~if~}\tilde{{\sigma}}_t^2\geq0\\0&\mathrm{~if~}\tilde{{\sigma}}_t^2<0&\end{cases} 
\end{aligned}$$

方差分量的近似 $(1-\alpha)100\%$ 置信区间为

$$\frac{df_{\hat{\sigma}_r^2}\hat{\sigma}_r^2}{\chi^2_{\alpha/2,df_{\hat{\sigma}_r^2}}}\leq\sigma_r^2\leq\frac{df_{\hat{\sigma}_r^2}\hat{\sigma}_r^2}{\chi^2_{1-(\alpha/2),df_{\hat{\sigma}_r^2}}},\quad r=t,g,\varepsilon $$

其中

$$df_{\hat{\sigma}_\varepsilon^2}=df_{MSERROR},\quad df_{\hat{\sigma}_s^2}=\frac{(\hat{\sigma}_g^2)^2}{\frac{\left[(1/k_3)MSR(g|\mu,\beta,t)\right]^2}{df_{MSR(g|\mu,\beta,t)}}+\frac{[(1/k_3)MSERROR]^2}{df_{MSERROR}}}$$

以及

$$df_{\hat{\sigma}_t^2}=\frac{(\hat{\sigma}_t^2)^2}{\frac{\left[(1/k_3)MSR(t|\mu,\beta)\right]^2}{df_{MSR(t|\mu,\beta)}}+\frac{\left[(k_1/k_3k_2)MSR(g|\mu,\beta,t)\right]^2}{df_{MSR(g|\mu,\beta,t)}}+\frac{\left[(1/k_2)(1-(k_1/k_3))MSERROR\right]^2}{df_{MSERROR}}}$$

期望均方可用于构建关于模型中方差分量假设的检验。检验 $H_{0}\colon{\sigma}_{t}^{2}=0\operatorname{vs}H_{a}\colon{\sigma}_{t}^{2}>0$ 的统计量为

$$\begin{aligned}F_{tc}&=\frac{MSR(t\mid\mu,\beta)}{Q}\quad\mathrm{where~}Q=\frac{k_1}{k_3}MSR(g|\mu,\beta,t)+\left(1-\frac{k_1}{k_3}\right)MSERROR\end{aligned}$$

其具有近似 $F$ 的抽样分布，分子自由度为 $df_{MSR(t|\mu,\beta)}$，分母自由度为 v，其中

$$v=\frac{Q^2}{\frac{\left[(k_1/k_3)MSR(g|\mu,\beta,t)\right]^2}{df_{MSR(g|\mu,\beta,t)}}+\frac{\left[\left(1-(k_1/k_3)\right)MSERROR\right]^2}{df_{MSERROR}}}$$

由 Satterthwaite 近似确定。

用于检验 $H_{0}\colon{\sigma}_{g}^{2}=0\operatorname{vs}H_{a}\colon{\sigma}_{g}^{2}>0$ 的统计量为

$$F_{gc}=\frac{MSR(g|\mu,\beta,t)}{MSERROR}$$

其具有近似 $F$ 的抽样分布，分子自由度为 $df_{MSR(g|\mu,\beta,t)}$，分母自由度为 $df_{MSERROR}$.

### 最大似然方法 {#sec22-2-2} 

最大似然方法可应用于完全似然函数，记作 $L(\boldsymbol{\beta},\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2,\sigma_\varepsilon^2)=(2\pi)^{-n/2}|\boldsymbol{\Sigma}|^{-1/2}\exp[-\frac12(\boldsymbol{y}-\boldsymbol X\boldsymbol{\beta})^{\prime}\boldsymbol{\Sigma}^{-1}(\boldsymbol{y}-\boldsymbol     X\boldsymbol{\beta})]$. 关于所有参数最大化该似然函数或等价地最小化 $-2\log[L(\boldsymbol{\beta},\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2,\sigma_\varepsilon^2)]$，提供了用于同时估计固定参数和随机效应参数的方程组 (Hartley and Rao, 1967). 所需的方程组为

$$\left.\frac{\partial l}{\partial\boldsymbol\beta}\right|_{\boldsymbol\beta=\hat{\boldsymbol\beta},\sigma=\hat{\sigma}}=0,\quad\left.\frac{\partial l}{\partial\sigma_1^2}\right|_{\boldsymbol\beta=\hat{\boldsymbol\beta},\sigma=\hat{\sigma}}=0,\quad\left.\frac{\partial l}{\partial\sigma_2^2}\right|_{\boldsymbol\beta=\hat{\boldsymbol\beta},\sigma=\hat{\sigma}}=0,\ldots,\quad\left.\frac{\partial l}{\partial\sigma_k^2}\right|_{\boldsymbol\beta=\hat{\boldsymbol\beta},\sigma=\hat{\sigma}}=0,\quad\left.\frac{\partial l}{\partial\sigma_\varepsilon^2}\right|_{\boldsymbol\beta=\hat{\boldsymbol\beta},\sigma=\hat{\sigma}}=0$$

其中 $l=-2\log[L(\boldsymbol{\beta},{\sigma}_1^2,{\sigma}_2^2,\ldots,{\sigma}_k^2,{\sigma}_\varepsilon^2)],{\sigma}^{\prime}=[{\sigma}_1^2,{\sigma}_2^2,\ldots,{\sigma}_k^2,{\sigma}_\varepsilon^2]$ 以及 $\hat{{\sigma}}^{\prime}=[\hat{{\sigma}}_1^2,\hat{{\sigma}}_2^2,\ldots,\hat{{\sigma}}_k^2,\hat{{\sigma}}_\varepsilon^2]$.

具有相等观测次数的双向混合模型用于演示最大似然估计的计算。该模型可以表示为

$$y_{ijk}=\mu_i+a_j+g_{ij}+\varepsilon_{ijk},\quad i=1,2,\ldots,t,\quad j=1,2,\ldots,a,\quad k=1,2,\ldots,n$$

其中 $\mu_i$ 表示混合效应因子 T 第 i 个水平的均值，$\alpha_j$ 表示随机效应因子 A 的第 j 个水平的效应，$\gamma_{ij}$ 表示随机相互作用效应，$\varepsilon_{ijk}$ 表示残差效应。在理想条件下，$a_j\thicksim i.i.d.N(0,\sigma_a^2),g_{ij}\thicksim i.i.d.N(0,\sigma_g^2),\varepsilon_{ijk}\thicksim i.i.d.N(0,\sigma_e^2)$，且 $a_j,g_{ij},\varepsilon_{ijk}$ 为独立随机变量。数据向量的方差为

$$\begin{aligned}
\boldsymbol{\Sigma}=&\,{\sigma}_a^2[\boldsymbol{J}_n\otimes\boldsymbol{I}_a\otimes\boldsymbol{J}_t]+{\sigma}_g^2[\boldsymbol{J}_n\otimes\boldsymbol{I}_a\otimes\boldsymbol{I}_t]+{\sigma}_\varepsilon^2[\boldsymbol{I}_n\otimes\boldsymbol{I}_a\otimes\boldsymbol{I}_t] \\
=&\,(\sigma_\varepsilon^2+n\sigma_g^2+nt\sigma_a^2)\left[\frac1{n}\boldsymbol J_n\otimes \boldsymbol I_a\otimes\frac1t\boldsymbol J_t\right]+(\sigma_\varepsilon^2+n\sigma_g^2)\left[\frac1{n}\boldsymbol J_n\otimes \boldsymbol I_a\otimes(\boldsymbol I_t-\frac1t\boldsymbol J_t)\right] \\
&+\sigma_\varepsilon^2\left[(\boldsymbol I_n-\frac1n\boldsymbol J_n)\otimes \boldsymbol I_a\otimes \boldsymbol I_t\right] \\
=&\,\lambda_1\boldsymbol V_1+\lambda_2\boldsymbol V_2+\lambda_3\boldsymbol V_3\text{ (say)}
\end{aligned}$$

其中 $\lambda_1=\sigma_\varepsilon^2+n\sigma_g^2+nt\sigma_a^2,\lambda_2=\sigma_\varepsilon^2+n\sigma_g^2,\lambda_3=\sigma_\varepsilon^2$ 以及 $\boldsymbol V_1,\boldsymbol V_2,\boldsymbol V_3$ 为幂等阵且两两正交。记号 $\boldsymbol A\otimes \boldsymbol B$ 表示矩阵 $\boldsymbol A,\boldsymbol B$ 的直积 (direct product) (Graybill 1976). 因此，$|\boldsymbol{\Sigma}|=\lambda_1^{a}\lambda_2^{a(t-1)}\lambda_3^{at(n-1)}$ 且 $\boldsymbol{\Sigma}^{-1}=(1/\lambda_1)\boldsymbol{V}_1+(1/\lambda_2)\boldsymbol{V}_2+(1/\lambda_3)\boldsymbol{V}_3$. $l=-2\log (L)$ 的值可表示为 $l=nta\log(2\pi)+a\log(\lambda_1)+a(t-1)\log(\lambda_2)+at(n-1)\log(\lambda_3)+Q$ 其中 $Q=[\boldsymbol{y}-(\boldsymbol{j}_n\otimes\boldsymbol{j}_a\otimes\boldsymbol{I}_t)\boldsymbol{\mu}]'\boldsymbol{\Sigma}^{-1}[\boldsymbol{y}-(\boldsymbol{j}_n\otimes\boldsymbol{j}_a\otimes\boldsymbol{I}_t)\boldsymbol{\mu}]$. 使用一些代数计算，并利用上述关系，$Q$ 的值可以表示为

$$Q=\frac1{\lambda_1}\boldsymbol y^{\prime}\boldsymbol A_1\boldsymbol y+\frac1{\lambda_2}\boldsymbol y^{\prime}\boldsymbol A_2\boldsymbol y+\frac1{\lambda_3}\boldsymbol y^{\prime}\boldsymbol A_3\boldsymbol y+\frac1{\lambda_1}nat(\bar{y}_{\cdot\cdot\cdot}-\bar{\mu}_{\cdot})^2+\frac{na}{\lambda_2}\sum_{i=1}^t(\bar{y}_{i\cdot\cdot}-\bar{y}_{\cdot\cdot\cdot}+\bar{\mu}_{\cdot}-\mu_{i})^2$$

其中

$$\begin{aligned}
&\boldsymbol{y^{\prime}}\boldsymbol A_1\boldsymbol{y} \left.=\boldsymbol y^{\prime}{\left[\frac1n\boldsymbol J_n\right.}\otimes\left(\boldsymbol I_a-\frac1a\boldsymbol J_a\right)\otimes\frac1t\boldsymbol J_t\right]\boldsymbol{y}=SSA  \\
&\boldsymbol{y'A}_2\boldsymbol{y} \left.=\boldsymbol y^{\prime}{\left[\frac1n\boldsymbol J_n\right.}\otimes\left(\boldsymbol I_a-\frac1a\boldsymbol J_a\right)\otimes\left(\boldsymbol I_t-\frac1t\boldsymbol J_t\right)\right]\boldsymbol{y}=SSTA  \\
&\boldsymbol{y'A_3}\boldsymbol{y} =\boldsymbol y^{\prime}\bigg[\left(\boldsymbol I_n-\frac1n\boldsymbol J_n\right)\otimes \boldsymbol I_a\otimes \boldsymbol I_t\bigg]\boldsymbol{y}=SSERROR 
\end{aligned}$$

$-2 \log(\text{likelihood})$ 的值可以表示为

$$\begin{aligned}
l =&\,nta\log(2\pi)+a\log(\lambda_1)+a(t-1)\log(\lambda_2)+at(n-1)\log(\lambda_3)  \\
&+\frac1{\lambda_3}SSERROR+\frac1{\lambda_1}SSA+\frac1{\lambda_2}SSTA+\frac1{\lambda_1}nat(\bar{y}_{\cdot\cdot\cdot}-\bar{\mu}_{\cdot})^2 \\
&+\frac{na}{\lambda_2}\sum_{i=1}^{t}(\bar{y}_{i\cdot\cdot}-\bar{y}_{\cdot\cdot\cdot}+\bar{\mu}_{\cdot}-\mu_i)^2
\end{aligned}$$

接下来，关于 $\lambda_1,\lambda_2,\lambda_3,\bar\mu_{\cdot}$ 和 $(\mu_i-\bar\mu_{\cdot})$ 对 $l$ 进行微分，并将导数设置为零，以提供方程组：

$$\begin{aligned}
&\left.\frac{\partial l}{\partial\lambda_3}\right|_{{\mu}=\hat{{\mu}},{\lambda}=\hat{{\lambda}}} =-\frac{at(n-1)}{\hat{\lambda}_3}+\frac{SSERROR}{\hat{\lambda}_3^2}=0\Longrightarrow\hat{\lambda}_3^2=\frac{SSERROR}{at(n-1)}=\hat{\sigma}_\varepsilon^2  \\
&\left.\frac{\partial l}{\partial\lambda_2}\right|_{{\mu}=\hat{{\mu}},{\lambda}=\hat{{\lambda}}} =-\frac{a(n-1)}{\hat{\lambda}_2}+\frac{SSTA}{\hat{\lambda}_2^2}=0\Longrightarrow\hat{\lambda}_2^2=\frac{SSTA}{a(n-1)}  \\
&\left.\frac{\partial l}{\partial\lambda_1}\right|_{{\mu}=\hat{{\mu}},{\lambda}=\hat{{\lambda}}} =-\frac a{\hat{\lambda}_1}+\frac{SSA}{\hat{\lambda}_1^2}=0\Longrightarrow\hat{\lambda}_1^2=\frac{SSA}a  \\
&\left.\frac{\partial l}{\partial\bar{\mu}.}\right|_{{\mu}=\hat{{\mu}},{\lambda}=\hat{{\lambda}}} =\frac{-2nat(\bar{y}_{\cdot\cdot\cdot}-\hat{\bar{\mu}}_{\cdot\cdot\cdot})}{\hat{\lambda}_1^2}=0\Longrightarrow\hat{\bar{\mu}}_{\cdot\cdot\cdot}=\bar{y}_{\cdot\cdot\cdot} 
\end{aligned}$$

以及

$$\left.\frac{\partial l}{\partial\mu_i-\bar{\mu}_{\cdot}}\right|_{{\mu}=\hat{{\mu}},{\lambda}=\hat{{\lambda}}}=\frac{-2na(\bar{y}_{i\cdot\cdot}-\bar{y}_{\cdot\cdot\cdot}+(\hat{\mu}_i-\hat{\bar{\mu}}_\cdot))}{\hat{\lambda}_2}=0\Longrightarrow(\hat{\mu}_i-\hat{\bar{\mu}}_{\cdot})=\bar{y}_{i\cdot\cdot}-\bar{y}_{\cdot\cdot\cdot}$$

或 $\hat{\bar{{\mu}}_i}=\bar{y}_{i\cdot\cdot}$.

利用 $\lambda_1=\sigma_\varepsilon^2+n\sigma_g^2+nt\sigma_a^2,\lambda_2=\sigma_\varepsilon^2+n\sigma_g^2$ 以及 $\lambda_3=\sigma_\varepsilon^2$ 的事实，方差分量的最大似然估计为：

$$\hat{\sigma}_\varepsilon^2=\hat{\lambda}_g,\quad\hat{\sigma}_g^2=\begin{cases}\dfrac{\hat{\lambda}_2-\hat{\lambda}_3}n&\mathrm{if~}\hat{\lambda}_2\geq\hat{\lambda}_3\\0&\mathrm{if~}\hat{\lambda}_2<\hat{\lambda}_3\end{cases},\quad\mathrm{and}\quad\hat{\sigma}_a^2=\begin{cases}\dfrac{\hat{\lambda}_1-\hat{\lambda}_2}{nt}&\mathrm{if~}\hat{\lambda}_1\geq\hat{\lambda}_2\\0&\mathrm{if~}\hat{\lambda}_1<\hat{\lambda}_2\end{cases}$$

请注意，$\hat{\sigma}_\varepsilon^2$ 是 ${\sigma}_\varepsilon^2$ 的无偏估计，而 $\hat{\sigma}_a^2,\hat{\sigma}_g^2$ 分别为 ${\sigma}_a^2,{\sigma}_g^2$ 的有偏估计。

### 残差最大似然法 {#sec22-2-3}

Corbeil and Searle (1976) 将似然函数分为两部分，一部分涉及固定效应，另一部分不包含固定效应。然后，他们从没有固定效应的模型部分获得了方差分量的最大似然估计，他们将其称为**受限最大似然估计** (restricted maximum likelihood estimators) 或 REML 估计。为了说明这一点，请再次考虑一般线性混合模型 $\boldsymbol y=\boldsymbol X\boldsymbol{\beta}+\boldsymbol Z_1\boldsymbol u_1+\boldsymbol Z_2\boldsymbol u_2+\cdots+\boldsymbol Z_k\boldsymbol u_k+\boldsymbol{\varepsilon}$ 以及第 \@ref(sec22-1) 节中给出的通常假设。假设 $\boldsymbol X$ 的秩为 p. 为了构造受限似然函数，令 $H$ 为秩为 N-p 的 N × (N - p) 矩阵，使得 $\boldsymbol {HH}^\prime = (\boldsymbol I - \boldsymbol X\boldsymbol X^- )$。定义变换

$$\boldsymbol z=\begin{bmatrix}\boldsymbol z_1\\\boldsymbol z_2\end{bmatrix}=\begin{bmatrix}\boldsymbol X^{\prime}\\\boldsymbol H^{\prime}\end{bmatrix}\boldsymbol y$$

其中

$$\begin{bmatrix}\boldsymbol X^{\prime}\\\boldsymbol H^{\prime}\end{bmatrix}$$

为 N × N 非奇异矩阵。因此从 $\boldsymbol y$ 到 $\boldsymbol z$ 的变换是一对一的变换。

$$\boldsymbol z=\begin{bmatrix}\boldsymbol{z}_1\\\boldsymbol{z}_2\end{bmatrix}$$

的分布为

$$\left.\begin{bmatrix}\boldsymbol z_1\\\boldsymbol z_2\end{bmatrix}\thicksim N\left[\begin{bmatrix}\boldsymbol X^{\prime}\boldsymbol X\boldsymbol{\beta}\\\boldsymbol 0\end{bmatrix}\right.,\begin{bmatrix}\boldsymbol X^{\prime}\boldsymbol{\Sigma}\boldsymbol X&\boldsymbol X^{\prime}\boldsymbol{\Sigma}\boldsymbol H\\\boldsymbol H^{\prime}\boldsymbol{\Sigma}\boldsymbol X&\boldsymbol H^{\prime}\boldsymbol{\Sigma}\boldsymbol H\end{bmatrix}\right]$$

$$\boldsymbol z=\begin{bmatrix}\boldsymbol{z}_1\\\boldsymbol{z}_2\end{bmatrix}$$

可以分解为给定 $\boldsymbol z_2$ 下 $\boldsymbol z_1$ 的似然乘以 $\boldsymbol z_2$ 的边际似然。$\boldsymbol z_2$ 的边际似然不取决于模型的固定效应参数，而是取决于方差分量。$\boldsymbol z_2$ 的边际似然为

$$L_{\boldsymbol H}(\boldsymbol z_2)=L_{\boldsymbol H}(\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2,\sigma_\varepsilon^2)=(2\pi)^{-(n-p)/2}|\boldsymbol{H^{\prime}\Sigma H}|^{-1/2}\exp\left[-\frac12\boldsymbol{yH}(\boldsymbol H^{\prime}\boldsymbol \Sigma \boldsymbol H)^{-1}\boldsymbol H^{\prime}\boldsymbol y\right]$$

方差分量的 REML 估计是使 $L_{\boldsymbol H}(\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2,\sigma_\varepsilon^2)$ 最大化的 $\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2,\sigma_\varepsilon^2$ 的值。该受限似然函数适用于随机向量 $\boldsymbol z_2$. 但固定效应模型的残差（参见第 \@ref(sec22-2) 节开头）可以表示为 $\boldsymbol r = \boldsymbol H \boldsymbol z_2$，即 $\boldsymbol z_2$ 的变换。因此，受限似然函数利用了残差模型中不依赖于固定效应的信息，或者说    受限似然是模型残差的函数。因此，为此方法选择的名称是残差最大似然。最大化 $l_R=\log[L_{\boldsymbol H}(\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2,\sigma_\varepsilon^2)]$ 的过程产生了一组需要求解的方程。不能保证方程组的解会提供一个位于参数空间中的解集。再次强调，获得方差分量的非负估计是很重要的，因为这些估计用于分析模型的固定效应部分。

在重新审视第 \@ref(sec22-2-2) 节中描述的均衡双向混合模型时，似然函数可表示为

$$\begin{aligned}
L(\mu,\lambda_1,\lambda_2,\lambda_3) =&\,\left[\left((2\pi)^{1/2}\lambda_1^{1/2}\right)^{-1}\exp\biggl(-\frac{nat}{2\lambda_1}(\bar{y}_{\cdot\cdot\cdot}-\bar{\mu}_{\cdot})^2\biggr)\biggr]\times\biggl[(2\pi)^{(t-1)/2}\lambda_2^{(t-1)/2}\biggr]^{-1}\right.   \\
&\times\exp\left[-\frac{na}{2\lambda_2}\sum_{i=1}^t(\bar{y}_{i\cdot\cdot}-\bar{y}_{\cdot\cdot\cdot}-\mu_i+\bar{\mu}_{\cdot})^2\right]\times\left[(2\pi)^{(nat-t)/2}\lambda_1^{(a-1)/2}\lambda_2^{(a-1)(t-1)/2}\lambda_3^{at(n-1)/2}\right]^{-1} \\
&\times\exp\left[-\frac12{\left(\frac{SSERROR}{\lambda_3}+\frac{SSA}{\lambda_1}+\frac{SSTA}{\lambda_2}\right)}\right] \\
=&\,L(\bar{\mu}_\cdot)\times L(\mu_i-\bar{\mu}_\cdot)\times L(\lambda_1,\lambda_2,\lambda_3)\mathrm{~(say)}
\end{aligned}$$

$\lambda_1,\lambda_2,\lambda_3$ 的 REML 估计通过最小化 $-2\log[L(\lambda_1,\lambda_2,\lambda_3)]$ 获得，这是通过对 $l_R$ 关于 $\lambda_1,\lambda_2,\lambda_3$ 求微分并将导数设置为零来实现的。微分和解是

$$\begin{aligned}
&\left.\frac{\partial l_R}{\partial\lambda_1}\right|_{\lambda=\hat{\lambda}} =\frac{(a-1)}{\hat{\lambda}_1}-\frac{SSA}{\hat{\lambda}_1^2}=0\Rightarrow\hat{\lambda}_1=\frac{SSA}{a-1}  \\
&\left.\frac{\partial l_R}{\partial\lambda_2}\right|_{\lambda=\hat{\lambda}} =\frac{(a-1)(t-1)}{\hat{\lambda}_2}-\frac{SSTA}{\hat{\lambda}_2^2}=0\Rightarrow\hat{\lambda}_2=\frac{SSTA}{(a-1)(t-1)}  \\
&\left.\frac{\partial l_R}{\partial\lambda_3}\right|_{\lambda=\hat{\lambda}} =\frac{at(n-1)}{\hat{\lambda}_3}-\frac{SSERROR}{\hat{\lambda}_3^2}=0\Rightarrow\hat{\lambda}_3=\frac{SSERROR}{at(n-1)} 
\end{aligned}$$

方差分量的估计计算如下：

$$\hat{\sigma}_\varepsilon^2=\hat{\lambda}_3,\quad\hat{\sigma}_g^2=\begin{cases}\frac{\hat{\lambda}_2-\hat{\lambda}_3}n&\mathrm{~if~}\hat{\lambda}_2\geq\hat{\lambda}_3\\0&\mathrm{~if~}\hat{\lambda}_2<\hat{\lambda}_3\end{cases},\quad\mathrm{~and~}\quad\hat{\sigma}_a^2=\begin{cases}\frac{\hat{\lambda}_1-\hat{\lambda}_2}n&\mathrm{~if~}\hat{\lambda}_1\geq\hat{\lambda}_2\\0&\mathrm{~if~}\hat{\lambda}_1<\hat{\lambda}_2&\end{cases}$$

REML 估计是方差分量的无偏估计，它们与使用矩法获得的估计相同。

### MINQUE 法 {#sec22-2-4}

随机模型的 MINQUE 法在第 \@ref(sec19-3) 节中描述。在该应用中，模型的均值为 $\boldsymbol j_n\mu$，而该混合模型的均值为 $\boldsymbol{Xb}$. 该估计可以推广到一般线性混合模型 (Swallow and Searle, 1978)，其中 $\boldsymbol \sigma=(\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2,\sigma_\varepsilon^2)^{\prime}$ 的 MINQUE 估计为 $\boldsymbol{\hat{\sigma}}=\boldsymbol{S}^{-1}\boldsymbol{q}$ 其中矩阵 $\boldsymbol S$ 的元素为 $s_{ii^{\prime}}=\mathrm{tr}[\boldsymbol{Z}_i\boldsymbol{Z}_i^{\prime}\boldsymbol{R}\boldsymbol{Z}_{i^{\prime}}\boldsymbol{Z}_i^{\prime}],i,i^{\prime}=1,2,\ldots,k+1$，其中 $i=k+1$ 对应于 $\boldsymbol \varepsilon$ 以及 $\boldsymbol Z_{k+1}=\boldsymbol I_n$，

$$\boldsymbol R=\boldsymbol{\Sigma}^{-1}[\boldsymbol I_n-\boldsymbol{X}(\boldsymbol X^{\prime}\boldsymbol{\Sigma}^{-1}\boldsymbol{X})^-\boldsymbol{X}^{\prime}]\boldsymbol{\Sigma}^{-1}$$
以及 $\boldsymbol q$ 的元素为 $\boldsymbol{q}_i=\boldsymbol{y}^{\prime}\boldsymbol{R}\boldsymbol{Z}_i\boldsymbol{Z}_i^{\prime}\boldsymbol{R}\boldsymbol{y},{i}=1,{2},\ldots,{k}+1$.

解取决于 $\boldsymbol \Sigma$ 的元素，对于 MINQUE 解，通常选择方差为 1，协方差为 0. 如果模型是均衡的，解通常不依赖于为 $\boldsymbol \Sigma$     选择的值，并且过程在一次迭代中收敛。不能保证系统的解提供参数空间中的值。请确保使用提供非负估计值的解，因为需要在分析模型的固定部分时使用这些估计。一次迭代 (One-iteration) MINQUE 估计在一些软件包中实现（称为零迭代 (zero iteration)），第 \@ref(chap23) 章中讨论的两个例子用此估计进行了计算。

## 混合模型固定效应部分的分析 {#sec22-3}

混合模型固定效应部分的分析包括如下所述的模型固定效应部分分析的所有方面。

### 估计 {#sec22-3-1}

有多种方法可以估计混合模型中 $\boldsymbol \beta$ 的可估函数。一般线性混合模型可表示为

$$\boldsymbol y=\boldsymbol X\boldsymbol{\beta}+\boldsymbol \varepsilon\quad\mathrm{~where~Var}(\boldsymbol \varepsilon)=\boldsymbol{\Sigma}$$

以及

$$\boldsymbol{\Sigma}={\sigma}_1^2\boldsymbol{Z}_1\boldsymbol{Z}_2^{\prime}+{\sigma}_2^2\boldsymbol{Z}_2\boldsymbol{Z}_2^{\prime}+\cdots+{\sigma}_k^2\boldsymbol{Z}_k\boldsymbol{Z}_k^{\prime}+{\sigma}_\varepsilon^2\boldsymbol{I}_N$$

线性组合 $\boldsymbol a^\prime \boldsymbol \beta$ 对于该混合模型是可估的，当且仅当存在一个向量 $\boldsymbol c$ 使得 $E(\boldsymbol c^\prime \boldsymbol y) = \boldsymbol a^\prime \boldsymbol \beta$ 时。该定义与一般线性模型的定义相同（参见第 \@ref(chap6) 章）。$\boldsymbol a^\prime \boldsymbol \beta$ 的估计是 $\boldsymbol a^\prime \hat {\boldsymbol \beta}$，为 $\boldsymbol \beta$ 的可估函数，其中 $\hat {\boldsymbol \beta}$ 是 $\boldsymbol \beta$ 的任意估计。

[^othersolutionequations]: 原文：The ordinary least squares estimator of $\boldsymbol a^\prime \boldsymbol\beta$ is $\boldsymbol a^\prime \hat{\boldsymbol\beta}_{LS}$ where $\boldsymbol{\hat{\beta}}_{LS}=(\boldsymbol{X}^{\prime}\boldsymbol{X})^-\boldsymbol{X}^{\prime}\boldsymbol{y}$ or some other solution for $\hat{\boldsymbol \beta}$ to the normal equations $\boldsymbol{X'X\hat{\boldsymbol{\beta}}}=\boldsymbol{X'y}$.

$\boldsymbol a^\prime \boldsymbol\beta$ 的普通最小二乘估计 (ordinary least squares estimator) 为 $\boldsymbol a^\prime \hat{\boldsymbol\beta}_{LS}$ ，其中 $\boldsymbol{\hat{\beta}}_{LS}=(\boldsymbol{X}^{\prime}\boldsymbol{X})^-\boldsymbol{X}^{\prime}\boldsymbol{y}$ 或者正规方程 $\boldsymbol{X'X\hat{\boldsymbol{\beta}}}=\boldsymbol{X'y}$ 的其他解 $\hat{\boldsymbol \beta}$[^othersolutionequations]。$\boldsymbol \beta$ 的最小二乘估计不依赖于 $\boldsymbol y$ 的协方差阵；也就是说，它不依赖于 $\boldsymbol \Sigma$.

如果 $\boldsymbol \Sigma$ 的元素已知（即，如果 $\sigma_1^2,\sigma_2^2,\ldots,\sigma_k^2$ 和 $\sigma_\varepsilon^2$ 已知），则 $\boldsymbol a^\prime \boldsymbol\beta$ 的**最佳线性无偏估计** (best linear unbiased estimator, **BLUE**) 为 $\boldsymbol a^\prime \hat{\boldsymbol\beta}_{BLUE}$，其中 $\hat{\boldsymbol\beta}_{BLUE} = \left(\boldsymbol{X'\Sigma}^{-1}\boldsymbol{X}\right)^-\boldsymbol{X'\Sigma}^{-1}\boldsymbol{y}$ 或者 $\boldsymbol{X'\Sigma^{-1}}\boldsymbol{X\hat{\beta}}_{BLUE}\boldsymbol{=X'\Sigma^{-1}y}$ 中的任何其他解 $\hat{\boldsymbol \beta}_{BLUE}$.

对于大多数均衡设计和一些简单的不均衡设计，$\boldsymbol{\hat{\beta}}_{BLUE}=\boldsymbol{\hat{\beta}}_{LS}$. 因此，对于这些设计，$\boldsymbol a^\prime \boldsymbol\beta$ 的  BLUE 是 $\boldsymbol a^\prime \hat{\boldsymbol\beta}_{LS}$，其中 $\hat{\boldsymbol{\beta}}_{LS}=(\boldsymbol{X}^{\prime}\boldsymbol{X})^{-}\boldsymbol{X}^{\prime}\boldsymbol{y}$，它不依赖于方差分量。

当设计不均衡并且方差分量未知时，生活就不那么容易了 (life is 
not so easy). 由于 BLUE 不存在（因为它取决于未知方差），因此必须获得**加权最小二乘估计** (**weighted least squares estimator**)，其中 $\hat{\boldsymbol\Sigma}$ 用作加权矩阵。协方差阵估计为

$$\hat{\boldsymbol{\Sigma}}=\hat{{\sigma}}_1^2\boldsymbol Z_1\boldsymbol Z_1^{\prime}+\hat{{\sigma}}_2^2\boldsymbol Z_2\boldsymbol Z_2^{\prime}+\cdots+\hat{{\sigma}}_k^2\boldsymbol Z_k\boldsymbol Z_k^{\prime}+\hat{{\sigma}}_\varepsilon^2\boldsymbol I_N$$

其中，$\hat{{\sigma}}_1^2,\hat{{\sigma}}_2^2,\ldots,\hat{{\sigma}}_k^2$ 和 $\hat{{\sigma}}_\varepsilon^2$ 是使用第 \@ref(sec22-2) 节中讨论的方法之一获得的方差分量估计。$\boldsymbol a^\prime \boldsymbol\beta$  的加权最小二乘估计或 $\boldsymbol a^\prime \boldsymbol\beta$ 的**估计的最佳线性无偏估计** (estimated BLUE, **EBLUE**) 是 $\boldsymbol a^\prime \hat{\boldsymbol\beta}_W$，其中 $\boldsymbol{\hat{\beta}}_W=\left(\boldsymbol{X}^{\prime}{\hat{\boldsymbol{\Sigma}}}^{-1}\boldsymbol{X}\right)^{-}\boldsymbol{X'}{\hat{\boldsymbol{\Sigma}}}^{-1}\boldsymbol{y}$
或者 ${\boldsymbol X'\hat{\boldsymbol\Sigma}^{-1}}\boldsymbol{X\hat{\beta}}_{W}={\boldsymbol X'\hat{\boldsymbol\Sigma}^{-1}\boldsymbol y}$ 的其他解  $\boldsymbol{\hat{\beta}}_W$. 对于大多数设计，随着样本量的增加，$\boldsymbol a^\prime \hat{\boldsymbol\beta}_W$ 收敛于 $\boldsymbol a^\prime {\boldsymbol\beta}$. 为了确保收敛，必须注意，随着样本量的增加，参数的数量不能趋于无穷大。$\boldsymbol a^\prime \hat{\boldsymbol\beta}_W$ 的大样本方差等于 $\mathrm{Var}(\boldsymbol a^{\prime}\boldsymbol{\hat{\beta}}_W)=\boldsymbol  a^{\prime}\left(\boldsymbol{X'\hat{\Sigma}}^{-1}\boldsymbol{X}\right)^-\boldsymbol{a}$. 这种对 $\boldsymbol a^\prime \hat{\boldsymbol\beta}_W$ 方差的近似没有考虑方差分量估计的变异性。Kackar and Harville (1984) 表明，$\boldsymbol  a^{\prime}\left(\boldsymbol{X'\hat{\Sigma}}^{-1}\boldsymbol{X}\right)^-\boldsymbol{a}$ 太小，需要进行调整，因为方差分量的真实值被各自的估计值取代。Kackar and Harville (1984) 和 Kenward and Roger (1997) 使用关于未知方差分量的泰勒级数展开来调整固定效应的标准误估计（该近似超出了本文的范围）。使用 SAS-Mixed 中的 `DDFM = KR` 选项可为固定效应参数的估计提供调整后的标准误估计。该选项还通过将 Satterthwaite 方法应用于 KR 调整的标准误估计来估计与标准误关联的自由度。

[^wherepartition]: 原文：where $\boldsymbol W$ is the partition of the generalized inverse of the information matrix corresponding to $\hat{\boldsymbol\beta}_{ML}$.

最大似然法也可以用来获得 $\boldsymbol a^\prime \boldsymbol\beta$ 的估计。$\boldsymbol a^\prime \boldsymbol\beta$ 的最大似然估计为 $\boldsymbol a^\prime \hat{\boldsymbol\beta}_{ML}$，其中 $\hat{\boldsymbol\beta}_{ML}$ 是未受限似然方程 (unrestricted likelihood equations) 的解。$\boldsymbol a^\prime \hat{\boldsymbol\beta}_{ML}$ 的方差为 $\boldsymbol a^{\prime} \boldsymbol W \boldsymbol a$，其中 $\boldsymbol W$ 是 $\hat{\boldsymbol\beta}_{ML}$ 对应的信息矩阵的广义逆的划分[^wherepartition]。第 \@ref(chap23) 章中的示例演示了均衡和不均衡设计的上述估计。

### 置信区间的构建 {#sec22-3-2}

利用 $\boldsymbol a^\prime \hat{\boldsymbol\beta}_W$ 的渐近抽样分布 $\boldsymbol a^{\prime}\boldsymbol{\hat{\beta}}_W\thicksim N\left[\boldsymbol a^{\prime}\boldsymbol{\beta},\boldsymbol{a}^{\prime}\left(\boldsymbol X^{\prime}\boldsymbol{\hat{\Sigma}}^{-1}\boldsymbol X\right)^-\boldsymbol{a}\right]$，可以获得关于 $\boldsymbol a^\prime \boldsymbol\beta$ 可估函数的 $(1 - \alpha)100\%$ 置信区间，其中 $\boldsymbol a^\prime \boldsymbol\beta$ 是 $\boldsymbol\beta$ 的可估计函数。$\boldsymbol a^\prime \hat{\boldsymbol\beta}_W$ 的标准误估计 $\widehat{s.e.}\left(\boldsymbol a^\prime \hat{\boldsymbol\beta}_W\right)=\sqrt{\boldsymbol{a}^{\prime}\left(\boldsymbol X^{\prime}\boldsymbol{\hat{\Sigma}}^{-1}\boldsymbol X\right)^-\boldsymbol{a}}$. 此时应使用标准误估计的 Kackar-Harville 型调整。通过使用 Satterthwaite 近似的推广 (Geisbrecht and Burns, 1985) 获得近似自由度：

$$\hat{v}=\frac{2\left[\boldsymbol a^{\prime}\left(\boldsymbol X^{\prime}\hat{\boldsymbol \Sigma}^{-1}\boldsymbol X\right)^-\boldsymbol a\right]^2}{\widehat{\mathrm{Var}}\left[\boldsymbol a^{\prime}\left(\boldsymbol X^{\prime}\hat{\boldsymbol \Sigma}^{-1}\boldsymbol X\right)^-\boldsymbol a\right]}$$

因此 $\boldsymbol a^\prime \boldsymbol\beta$ 的 $(1 - \alpha)100\%$ 近似置信区间为

$${\boldsymbol a'\hat{\boldsymbol{\beta}}}_W\pm\left(t_{\alpha/2,\hat{v}}\right)\sqrt[]{\boldsymbol{a'}\left(\boldsymbol{X'\hat{\boldsymbol{\Sigma}}}^{-1}\boldsymbol{X}\right)^-\boldsymbol{a}}$$

可以关于一组可估线性组合 $\boldsymbol a_1^\prime \hat{\boldsymbol\beta},\boldsymbol a_2^\prime \hat{\boldsymbol\beta},\ldots,\boldsymbol a_m^\prime \hat{\boldsymbol\beta}$ 构造一组同时置信区间，但需要近似每个线性组合的自由度。然后可以使用第 \@ref(chap3) 章的技术对多重比较进行调整。

### 假设检验 {#sec22-3-3}

通常，假设检验为固定效应的线性组合提供适当的推断过程。为检验 $H_0{:\boldsymbol{H\beta}}=\boldsymbol{b}\mathrm{~vs~}H_a{:\boldsymbol{H\beta}}\neq\boldsymbol{b}$ 计算统计量

$$Q=(\boldsymbol H\hat{\boldsymbol \beta}_W-\boldsymbol b)^{\prime}\left[\boldsymbol H\left(\boldsymbol X^{\prime}\hat{\boldsymbol \Sigma}^{-1}\boldsymbol X\right)^{-}\boldsymbol H^{\prime}\right]^{-1}(\boldsymbol H\hat{\boldsymbol \beta}_W-\boldsymbol b)$$

在原假设条件下，$Q$ 的渐进抽样分布为 $\chi^2_q$，其中 $q=\text{rank}(\boldsymbol H)$. 即 $q$ 是 $\boldsymbol H$ 中 $\boldsymbol \beta$ 的线性独立线性组合的个数。小样本检验统计量为 $F_c = Q/v$，其中 v 是与 $\boldsymbol H\left(\boldsymbol X^{\prime}\hat{\boldsymbol \Sigma}^{-1}\boldsymbol X\right)^{-}\boldsymbol H^{\prime}$ 相关的近似自由度。近似自由度 (SAS Institute, Inc., 1999) 是通过在 $\boldsymbol H\left(\boldsymbol X^{\prime}\hat{\boldsymbol \Sigma}^{-1}\boldsymbol X\right)^{-}\boldsymbol H^{\prime}$ 上进行谱分解来计算的：$\boldsymbol H\left(\boldsymbol X^{\prime}\hat{\boldsymbol \Sigma}^{-1}\boldsymbol X\right)^{-}\boldsymbol H^{\prime}=\boldsymbol P^{\prime}\boldsymbol \Delta \boldsymbol P$，其中 $\boldsymbol \Delta$ 是 q × q 对角    阵，对角元为 $\boldsymbol H\left(\boldsymbol X^{\prime}\hat{\boldsymbol \Sigma}^{-1}\boldsymbol X\right)^{-}\boldsymbol H^{\prime}$ 的特征根，$\boldsymbol P$ 是包含相应的特征向量的 q × q 矩阵。令 $\boldsymbol h_s$ 为 $\boldsymbol{PH}$ 的第 s 行，则 $v_s=2d_s^2/(\boldsymbol \xi_s^{\prime}\boldsymbol \Omega\boldsymbol \xi_s)$ 其中

$$\boldsymbol \xi_s=\frac{\partial\left[\boldsymbol{h'}_s\left(\boldsymbol{X'\Sigma}^{-1}\boldsymbol{X}\right)^-\boldsymbol{h}_s\right]}{\partial\boldsymbol{\sigma}}$$

为关于每个方差分量的导数向量，以及 $\boldsymbol\Omega$ 为方差分量估计 $\hat{\boldsymbol \sigma}$ 的渐近协方差阵。那么

$$\zeta=\sum_{s=1}^q\frac{V_s}{V_s-2}\mathrm{In}(v_s>2)$$

其中指示函数 $\mathrm{In}(v_s > 2)$ 删除 $v_s \le 2$ 的项。那么与 $F_c$  相关的分母自由度的近似数为

$$v=\begin{cases}\frac{2\zeta}{\zeta-q}&\mathrm{if~}\zeta>q\\0&\mathrm{if~}\zeta>q&\end{cases}$$

第 \@ref(chap23) 章讨论了示例。

## 最佳线性无偏预测 {#sec22-4}

在某些情况下，预测随机变量的值或预测研究中使用的随机效应的值是有意义的。假设有线性模型 $\boldsymbol y=\boldsymbol X\boldsymbol{\beta}+\boldsymbol \varepsilon$，其中 $\mathrm{Var}(\varepsilon) = \boldsymbol \Sigma$. 假设有一个随机变量 $\boldsymbol\omega$，其值未知，但假设 $\boldsymbol\omega\thicksim N(\boldsymbol k^\prime\boldsymbol\beta,\boldsymbol \sigma_\omega^2)$ 和 $\operatorname{Cov}(\boldsymbol y,\boldsymbol\omega)=\boldsymbol c$。目的是预测 $\boldsymbol\omega$ 的值。$\boldsymbol\omega$ 的预测值称为**最佳线性无偏预测** (best linear unbiased predictor, **BLUP**)，当

1. $\tilde{\boldsymbol{\omega}}=\boldsymbol{a'y}+b$
2. $E(\tilde{\boldsymbol\omega})=\boldsymbol k^{\prime}{\boldsymbol\beta}$ 以及
3. $E(\tilde{\boldsymbol\omega}-\boldsymbol\omega)^2$ 被最小化

所得的 $\boldsymbol\omega$ 的 BLUP 为 $\tilde{\boldsymbol{\omega}}=\boldsymbol{c'\Sigma}^{-1}(\boldsymbol{y}-\boldsymbol{X}\boldsymbol{\hat{\beta}}_{BLUE})+\boldsymbol{k'}\boldsymbol{\hat{\beta}}_{BLUE}$. 当 $\boldsymbol \Sigma$ 的元素未知且需要估计时，$\boldsymbol\omega$ 的**最佳线性无偏预测估计** (estimated BLUP, **EBLUP**) 为

$$\hat{\tilde{\boldsymbol \omega}}=c^{\prime}\boldsymbol{\hat{\Sigma}}^{-1}(\boldsymbol{y}-\boldsymbol{X}\boldsymbol{\hat{\beta}}_W)+\boldsymbol{k}^{\prime}\boldsymbol{\hat{\beta}}_W$$

如果 $\boldsymbol{k}^{\prime}\boldsymbol{\beta}=0$ 那么 $\hat{\boldsymbol{\tilde{\omega}}}=\boldsymbol{c}^{\prime}\boldsymbol{\hat{\Sigma}}^{-1}(\boldsymbol{y}-\boldsymbol{X}\boldsymbol{\hat{\beta}}_W)$

对于一般线性混合模型，$\boldsymbol{y}=\boldsymbol{X\beta}+\boldsymbol{Z}\boldsymbol{u}+\boldsymbol{\varepsilon}$ 其中 $\boldsymbol u \thicksim N(\boldsymbol 0,\boldsymbol G)$ 和 $\boldsymbol \varepsilon \thicksim N(\boldsymbol 0, \boldsymbol R)$，$\boldsymbol y$ 和 $\boldsymbol u$ 之间的协方差为 $\operatorname{Cov}(\boldsymbol y,\boldsymbol     u)=\boldsymbol G\boldsymbol{Z}^{\prime}$ 且 $\boldsymbol u$ 的 BLUP 为 $\tilde{\boldsymbol{u}}=\boldsymbol G\boldsymbol{Z}^{\prime}\boldsymbol{\Sigma}^{-1}(\boldsymbol{y}-\boldsymbol{X}\boldsymbol{\hat{\beta}})$.

## 混合模型方程组 {#sec22-5}

Henderson (1984) 开发了混合模型方程组，其解可以同时得出 $\boldsymbol \beta$ 可估计函数 $\boldsymbol a^\prime \boldsymbol\beta$ 的 BLUE 以及随机效应 $\boldsymbol u=(\boldsymbol u_1^{\prime},\boldsymbol u_2^{\prime},\ldots,\boldsymbol u_k^{\prime})^{\prime}$ 的 BLUPs. 将混合模型表示为 $\boldsymbol{y}=\boldsymbol{X\beta}+\boldsymbol{Z}\boldsymbol{u}+\boldsymbol{\varepsilon}$，其中

$$\begin{aligned}
 \boldsymbol{Z}\boldsymbol u&=\boldsymbol{Z}_1\boldsymbol u_1+\boldsymbol{Z}_2\boldsymbol u_2+\cdots+\boldsymbol{Z}_k\boldsymbol u_k  \\
\boldsymbol u_i&\thicksim N(\boldsymbol 0,\sigma_i^2\boldsymbol I_{n_i}),\quad i=1,2,\ldots,k \\
\boldsymbol{\varepsilon}&\thicksim N(\boldsymbol{0},{\sigma}_\varepsilon^2\boldsymbol{I}_N)
\end{aligned}$$

其中 $\boldsymbol{u}_1,\boldsymbol{u}_2,\ldots,\boldsymbol{u}_k,\boldsymbol{\varepsilon}$ 为独立随机变量。

注意到

$$\boldsymbol{u}=\begin{bmatrix}\boldsymbol{u}_1\\\boldsymbol{u}_2\\\vdots\\\boldsymbol{u}_k\end{bmatrix}\thicksim N(\boldsymbol{0},\boldsymbol{G})$$

其中

$$\boldsymbol{G}=\begin{bmatrix}{\sigma}_1^2\boldsymbol{I}_{n_1}&0&\cdots&0\\0&{\sigma}_2^2\boldsymbol{I}_{n_2}&\cdots&0\\\vdots&\vdots&\ddots&\vdots\\0&0&\cdots&{\sigma}_k^2\boldsymbol{I}_{n_k}\end{bmatrix}$$

使用这些假设，$\boldsymbol y$ 的边际分布可以表示为

$\boldsymbol y\thicksim N(\boldsymbol X\boldsymbol{\beta},\boldsymbol{\Sigma})$

其中

$$\mathrm{Var}(\boldsymbol y)=\boldsymbol{\Sigma}=\sum_{i=1}^k{\sigma}_i^2\boldsymbol{Z}_i\boldsymbol{Z}_i^{\prime}+\boldsymbol{\sigma}_\varepsilon^2\boldsymbol{I}=\boldsymbol{Z}\boldsymbol{G}\boldsymbol{Z}^{\prime}={\sigma}_\varepsilon^2\boldsymbol{I}_N$$

该模型及其假设意味着 $\operatorname{Cov}(\boldsymbol{y},\boldsymbol{u})=\boldsymbol{G}\boldsymbol{Z}^{\prime}$.

使用上述信息，$\boldsymbol y$ 和 $\boldsymbol u$ 的联合分布为

$$\begin{bmatrix}\boldsymbol{u}\\\boldsymbol{y}\end{bmatrix}\thicksim N\left(\begin{bmatrix}\boldsymbol{0}\\\boldsymbol{X}\boldsymbol{\beta}\end{bmatrix},\begin{bmatrix}\boldsymbol{G}&\boldsymbol{G}\boldsymbol{Z}^{\prime}\\\boldsymbol{Z}\boldsymbol{G}^{\prime}&\boldsymbol{\Sigma}\end{bmatrix}\right)$$

给定 $\boldsymbol u$ 的 $\boldsymbol y$ 的条件分布是 $\boldsymbol{y}|\boldsymbol{u}\sim N(\boldsymbol{X\beta}+\boldsymbol{Zu},\sigma_\varepsilon^2\boldsymbol{I})$ 以及 $\boldsymbol u$ 的边际分布为 $\boldsymbol{u}\thicksim N(\boldsymbol{0},\boldsymbol{G})$. 因此，$\boldsymbol y$ 和 $\boldsymbol u$ 的联合分布可以表示为给定 $\boldsymbol u$ 的 $\boldsymbol y$ 的条件分布与 $\boldsymbol u$ 的边际分布的乘积，如

$$\begin{aligned}
h(\boldsymbol y,\boldsymbol u)&=f(\boldsymbol y|\boldsymbol u)g(\boldsymbol u) \\
&=[2\pi\sigma_\varepsilon^2]^{-(N/2)}\exp\left[-\frac1{2\sigma_\varepsilon^2}\left(\boldsymbol y-\boldsymbol X\boldsymbol{\beta}-\boldsymbol Z\boldsymbol{u}\right)^{\prime}\left(\boldsymbol y-\boldsymbol X\boldsymbol{\beta}-\boldsymbol Z\boldsymbol{u}\right)\right][2\pi]^{-(q/2)}|\boldsymbol{G}|^{-1/2}\exp(-\frac12\boldsymbol{u}\boldsymbol G^{-1}\boldsymbol{u}) 
\end{aligned}$$

其中 $q=n_1+n_2+\ldots,n_k$.

[^henderson1984]: 原文：Henderson (1984) differentiated $-2\log[f(\boldsymbol y|\boldsymbol u)g(\boldsymbol u)]$ with respect to $\boldsymbol \beta$ and $\boldsymbol u$ to derive the mixed model equations whose solution provides the BLUE of estimable functions of $\boldsymbol \beta$ and BLUPs of $\boldsymbol u$ as

为导出混合模型方程（其解提供了 $\boldsymbol \beta$ 可估函数的 BLUE 和 $\boldsymbol u$ 的 BLUPs），Henderson (1984) 关于 $\boldsymbol \beta$ 和 $\boldsymbol u$ 对 $-2\log[f(\boldsymbol y|\boldsymbol u)g(\boldsymbol u)]$ 进行微分[^henderson1984]：

$$\begin{aligned}
-2\log[f(\boldsymbol y|\boldsymbol u)g(\boldsymbol u)] =&\,(N+q)\mathrm{log}(2\pi)+N\mathrm{log}(\sigma_\varepsilon^2)+\mathrm{log}|\boldsymbol G|  \\
&+\frac1{\sigma_\varepsilon^2}(\boldsymbol y-\boldsymbol X\boldsymbol{\beta}-\boldsymbol Z\boldsymbol{u})^{\prime}(\boldsymbol y-\boldsymbol X\boldsymbol{\beta}-\boldsymbol Z\boldsymbol{u})+\boldsymbol u^{\prime}\boldsymbol G^{-1}\boldsymbol{u} \\
=&\,h\mathrm{~(say)}
\end{aligned}\\
\frac{\partial h}{\partial\boldsymbol{\beta}}=\frac2{\sigma_\varepsilon^2}\boldsymbol X^{\prime}(\boldsymbol y-\boldsymbol X\boldsymbol{\beta}-\boldsymbol Z\boldsymbol{u}),\quad\frac{\partial h}{\partial\boldsymbol{u}}=\frac2{\sigma_\varepsilon^2}\boldsymbol Z^{\prime}(\boldsymbol y-\boldsymbol X\boldsymbol{\beta}-\boldsymbol Z\boldsymbol{u})+2\boldsymbol G^{-1}\boldsymbol{u}$$

将导数设置为零提供混合模型方程组

$$\begin{bmatrix}\boldsymbol X^{\prime}\boldsymbol X&\boldsymbol X^{\prime}\boldsymbol{Z}\\\boldsymbol{Z^{\prime}}\boldsymbol X&(\boldsymbol{Z^{\prime}}\boldsymbol{Z}+\sigma_\varepsilon^2\boldsymbol{G^{-1}})\end{bmatrix}{\begin{bmatrix}\hat{\boldsymbol \beta}_{BLUE}\\\tilde{\boldsymbol u}\end{bmatrix}}={\begin{bmatrix}\boldsymbol X^{\prime}\boldsymbol y\\\boldsymbol Z^{\prime}\boldsymbol{u}\end{bmatrix}}$$

如果 $\boldsymbol \varepsilon$（残差）的方差是 $\boldsymbol R$，则混合模型方程组变为

$$\begin{bmatrix}\boldsymbol X^{\prime}\boldsymbol R^{-1}\boldsymbol X&\boldsymbol X^{\prime}\boldsymbol R^{-1}\boldsymbol Z\\\boldsymbol Z^{\prime}\boldsymbol R^{-1}\boldsymbol X&(\boldsymbol Z^{\prime}\boldsymbol R^{-1}\boldsymbol Z+\boldsymbol G^{-1})\end{bmatrix}\begin{bmatrix}\hat{\boldsymbol\beta}_{BLUE}\\\tilde{\boldsymbol u}\end{bmatrix}=\begin{bmatrix}\boldsymbol X^{\prime}\boldsymbol R^{-1}\boldsymbol y\\\boldsymbol Z^{\prime}\boldsymbol R^{-1}\boldsymbol u\end{bmatrix}$$

混合模型方程组的解为

$$\hat{\boldsymbol{\beta}}=\left(\boldsymbol X^{\prime}{\boldsymbol{\Sigma}}^{-1}\boldsymbol X\right)^{-}\boldsymbol X^{\prime}{\boldsymbol{\Sigma}}^{-1}\boldsymbol y\quad\mathrm{and}\quad\tilde{\boldsymbol{u}}=\boldsymbol G\boldsymbol Z^{\prime}{\boldsymbol{\Sigma}}^{-1}(\boldsymbol y-\boldsymbol X\hat{\boldsymbol{\beta}})$$

通常，$\boldsymbol{\Sigma}$ 的元素是未知的，需要进行估计，因此估计的 BLUE（EBLUE）和估计的 BLUP（EBLUP）为

$$\boldsymbol{\hat{\beta}_W}=\left(\boldsymbol X^{\prime}\hat{\boldsymbol{\Sigma}}^{-1}\boldsymbol X\right)^-\boldsymbol{X}^{\prime}\boldsymbol{\hat{\Sigma}}^{-1}\boldsymbol{y}\quad\mathrm{and}\quad\hat{\tilde{\boldsymbol u}}=\boldsymbol{\hat{G}}\boldsymbol{Z}^{\prime}\boldsymbol{\hat{\Sigma}}^{-1}(\boldsymbol{y}-\boldsymbol{X}\boldsymbol{\hat{\beta}})$$

## 结束语 {#sec22-6}

本章介绍了一个混合模型的理论发展，用于描述实验中的处理和设计结构，其中涉及固定效应和随机效应。混合模型的分析包括方差分量以及固定效应参数函数的估计和推断。对于模型的随机效应部分，可以使用平方和方法而不是方差分量估计的渐近抽样分布来进行假设检验。为了对固定效应进行推断，对模型的固定效应部分进行了分析。可以使用可估函数的估计的渐近抽样分布或近似 $F$ 统计量（其中分母自由度通过 Satterthwaite 型近似确定）来执行置信区间和假设检验。还讨论了随机效应的最佳线性无偏预测以及混合模型方程组。

## 练习 {#sec22-7}
