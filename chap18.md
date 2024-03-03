# 随机效应模型和方差分量 {#chap18}

具有一个以上**随机分量** (random component) 的模型适用于多种情况，包括随机效应和混合效应模型，其中处理结构中的部分或全部因子是随机的，或者如裂区和重复测量设计那样有多种尺寸的实验单元。对于此类模型，感兴趣的参数包括与随机分量的分布相关的方差（通常称为**方差分量**, variance components）。重要的是能够识别模型的随机分量，并能够在模型分析中使用它们。当对给定模型进行方差分析时，需要均方的期望值（均方是方差分量的函数），以构造适当的检验统计量，并确定用于比较固定效应参数的标准误。还需要能够获得方差分量的估计，并检验假设和构造方差分的函数的置信区间。随机效应模型及其分析方法的讨论分为四章。本章定义了随机效应模型，并描述了一种计算平方和期望的一般程序。该程序可以很容易地被计算机软件用于评估平方和的期望值。估计问题在第 \@ref(chap19) 章中讨论，检验假设和构造置信区间的方法在第 \@ref(chap20) 章中介绍，第 \@ref(chap21) 章提供了一个详细的示例分析。

## 介绍 {#sec18-1}

在抽样方案和感兴趣的参数中，使用随机效应模型背后的哲学与使用固定效应模型（在前几章中讨论的）背后的哲学截然不同。在讨论这些差异之前，先给出随机效应和固定效应的定义。

::: {.definition #18-1}
♦

如果一个因子的水平由可能水平总体中的随机样本组成，则该因子是一个**随机效应因子** (random effect factor).
:::

::: {.definition #18-2}
♦

如果一个因子的水平是通过非随机过程选择的，或者如果它的水平由可能水平的整个总体组成，那么它就是一个**固定效应因子** (fixed effect factor).
:::

因此，为了确定一个因子是固定效应还是随机效应，我们需要知道实验者是如何选择该因子的水平的。如果实验中包括了该因子所有可能水平或一组选定的因子水平，则该因子被视为混合效应。如果使用某种形式的随机化来选择实验中包括的水平，那么该因子就是随机效应。

:::: {.rmdnote data-latex="{规则}"}
**规则**：因子的水平在被证明为随机之前都是固定的。
::::

为了确定一个因子的水平是随机的，必须描述可能水平的总体或概念总体，并且必须指定随机选择该因子水平的方法。要对水平总体进行推断，因此如果无法描述该总体，则推断可能没有意义。

例如，假设一个植物育种学家想要研究小麦品种的某个特性（比如说产量）。存在许多可能的小麦品种（品种的总体），但是如果他想要研究某一组特定的品种，那么他就会选择这些品种进行实验。在这种情况下，因子 “品种” 被称为固定效应，因为品种的水平是选定的 (chosen) 或固定的 (fixed). 然而，如果植物育种学家感兴趣的是特性在总体中的品种之间是如何分布的，那么他就对哪些品种包含在实验中不感兴趣。在这种情况下，植物育种学家可以从品种的总体中随机选择包含在实验中的品种。因此在第二个实验中，因子品种是一个随机效应。在构造模型来描述给定的实验情况时，必须指出一个因子是随机效应还是固定效应。前几章所考虑的模型是在假设处理结构中的所有因子都是固定效应，并且实验单元只有一个尺寸的假设下构造的。然而，在第四章引入区组时，已经暗示了随机效应的概念，其中假设因子 “区组” 是一个随机效应。此外，一个基本假设是，与模型的设计结构相关的那些因子是随机效应。一些因子既有一组固定效应水平，也有一组随机效应水平 (Njuho and Milliken, 2005)，但这里没有详细阐述这一主题。根据设计结构中实验单元的尺寸和对处理结构中因子的假设，可以构造三种基本类型的模型。这些类型的模型定义如下。

::: {.definition #18-3}
♦

如果处理结构中的所有因子都是固定效应，并且研究中只有一个尺寸的实验单元，没有区组，并且方差都相等，则该模型称为固定 (fixed) 或**固定效应** (fixed effects) 模型。
:::

::: {.definition #18-4}
♦

如果处理结构中的所有因子都是随机效应（设计结构中的全部因子已经被假设为随机效应），则该模型称为随机 (random) 或**随机效应** (random effects) 模型。
:::

::: {.definition #18-5}
♦

如果处理结构中的一些因子是固定效应，一些因子是随机效应，或者处理结构中所有因子都是固定效应并且设计结构中有一个以上尺寸的实验单元（或者模型中至少有两个方差分量），则该模型称为混合 (mixed) 或**混合效应** (mixed effects) 模型。
:::

本章讨论的模型都是随机效应模型。第 \@ref(chap22) 章和第 \@ref(chap23) 章介绍了对混合模型的讨论。下面的例子有助于激发随机效应模型的应用和分析。

### 示例 18.1：随机效应嵌套处理结构 {#sec18-1-1}

一个消费者团体研究了美国至少拥有 20,000 人口的城市中咖啡价格的变化。该团体希望调查的三个因子是州、州内的城市以及州内城市内的商店。处理结构或抽样设计是一个涉及州、城市和商店的三向两水平的嵌套系统，其中城市嵌套在州内，商店嵌套在城市内。抽样程序是从所有可能的州中随机选择 $r$ 个州（$r < 50$）。接下来，从第 $i$ 个州中的 $C_i$ 个城市中随机选择 $t_i$ 个城市（$t_i < C_i$），这些城市的人口至少为 20,000. 最后，从第 $i$ 个州中的第 $j$ 个城市的 $S_{ij}$ 家商店中随机选择 $n_{ij}$ 家商店（$n_{ij} < S_{ij}$），并确定每家随机选择的商店中特定等级的咖啡价格。可以用来描述咖啡价格变化的模型是

$$y_{ijk}=\mu+s_i+c_{j(i)}+a_{k(ij)}\quad i=1,2,\ldots,r,j=1,2,\ldots,t_i,\quad k=1,2,\ldots,n_{ij}$$

其中 $\mu$ 表示美国咖啡的平均价格，$s_i$ 表示第 $i$ 个随机选择的州的效应，$c_{j(i)}$ 表示第 $i$ 个州中第 $j$ 个随机选择的城市的效应，$a_{k(ij)}$ 表示从第 $i$ 个州中的第 $j$ 个城市随机选择的第 $k$ 家商店的效应。假设

1. $s_i\sim i.i.d.\,N(0,\sigma_{\mathrm{State}}^2)$
2. $c_{j(i)}\sim i.i.d.\,N(0,\sigma_{\mathrm{City}}^2)$
3. $a_{k(ij)}\sim i.i.d.\,N(0,\sigma_{\mathrm{Store}}^2)$

还假设所有随机效应都是彼此独立分布的。该随机效应模型中的参数为 $\mu,\sigma_{\mathrm{State}}^2,\sigma_{\mathrm{City}}^2,\sigma_{\mathrm{Store}}^2$. 项 $s_i,c_{j(i)},a_{k(ij)}$ 是随机变量，不是模型中的参数。大多数应用只对参数估计感兴趣，对预测随机变量的值不感兴趣。但是，当对随机变量的预测感兴趣时，可以获得随机效应的**最佳线性无偏预测估计** (estimated best linear unbiased predictors, EBLUP) (Littell et al., 2006; Milliken and Johnson, 2002).

咖啡价格的方差可以通过模型参数表示为

$$\mathrm{Var}(y_{ijk})=\sigma_{\mathrm{Price}}^2=\sigma_{\mathrm{State}}^2+\sigma_{\mathrm{City}}^2+\sigma_{\mathrm{Store}}^2$$

同一城市（和州）内商店的两种咖啡价格之间的协方差为

$$\mathrm{Cov}(y_{ij1},~y_{ij2})=\sigma_{y_{ij1}y_{ij2}}=\sigma_{\mathrm{State}}^2+\sigma_{\mathrm{City}}^2$$

这两个价格之间的相关性为

$$\rho_{y_{ij1}y_{ij2}}=\frac{\sigma_{\mathrm{State}}^2+\sigma_{\mathrm{City}}^2}{\sigma_{\mathrm{State}}^2+\sigma_{\mathrm{City}}^2+\sigma_{\mathrm{Store}}^2}$$

州内不同城市商店的两种咖啡价格之间的协方差为 $\mathrm{Cov}(y_{i11},y_{i22})=\sigma_{y_{i11}y_{i22}}=\sigma_{\mathrm{State}}^2$. 同一州内两个不同城市的两个价格之间的相关性为

$$\rho_{y_{i11}y_{i22}}=\frac{\sigma_{\mathrm{State}}^2}{\sigma_{\mathrm{State}}^2+\sigma_{\mathrm{City}}^2+\sigma_{\mathrm{Store}}^2}$$

因此，当模型包含随机效应时，会在结果数据上施加与咖啡价格相似的**相关结构** (correlation structure). 能够构造模型并识别各种观测分组之间的变异来源以及相关性是有意义的。以下示例用于演示构造一个或多个随机效应的模型，并评估数据集中观测结果之间的隐含协方差或相关结构。

## 矩阵表示法中的一般随机效应模型 {#sec18-2}

为了描述用于评估平方和期望的方法，有必要使用一些通用符号来描述一般随机效应模型 (general random effects model). 本节介绍一般随机效应模型的矩阵表示，后面的部分将使用该模型来演示计算期望均方的一般方法。为了帮助可视化一般随机效应模型及其在矩阵方面的表达，我们研究了完全随机设计结构中单向处理结构的随机效应模型。

### 示例 18.2：单向随机效应模型 {#sec18-2-1}

描述完全随机设计结构中的单向随机效应处理结构的模型是

\begin{equation}
y_{ij}={\mu}+{u}_i+{\varepsilon}_{ij},\quad i=1,2,\ldots,t,j=1,2,\ldots,n_i
(#eq:18-1)
\end{equation}

其中 $\mu$ 是响应的总体平均值，$\mu_i$ 表示第 $i$ 个随机选择的处理的效应，并假设其分布为 $i.i.d.\,N(0, \sigma^2_u)$，$\varepsilon_{ij}$ 表示第 $i$ 个处理的第 $j$ 个观测值的随机误差，并假设其分布为 $i.i.d.\,N(0, \sigma^2_\varepsilon)$. 还假设 $\mu_i$ 和 $\varepsilon_{ij}$ 是独立的随机变量。这些假设允许评估观测值的方差和协方差。观测值的方差是

$$\begin{aligned}
\operatorname{Var}(y_{ij})& ={\sigma}_y^2=\mathrm{Var}({\mu}+{\mu}_i+{\varepsilon}_{ij}){=}\mathrm{Var}({u}_i)+\mathrm{Var}({\varepsilon}_{ij})  \\
&=\sigma_u^2+\sigma_\varepsilon^2
\end{aligned}$$

$y_{ij}$ 的方差有两个分量，包括处理总体的方差或 $u$ 水平的方差和实验单元的方差，因此称为方差分量 (variance components or components of variance). 从 $u_i$ 的相同可能性获得的两个观测值的协方差为

$$\begin{aligned}\mathrm{Cov}(y_{ij},y_{ij^{\prime}})&=\mathrm{Cov}(\mu+u_i+\varepsilon_{ij},\mu+u_i+\varepsilon_{ij^{\prime}})\\&=\mathrm{Cov}(u_i,u_i)=\mathrm{Var}(u_i)=\sigma_u^2\end{aligned}$$

从不同 $i$ 值获得的两个观测值之间的协方差为零。因此，从相同 $i$ 获得的观测值是相关的，这种相关性称为**组内相关性** (intraclass correlation)，定义为

$$\rho=\frac{\mathrm{Cov}(y_{ij},y_{ij},)}{\sqrt{\mathrm{Var}(y_{ij})\mathrm{Var}(y_{ij},)}}=\frac{\sigma_u^2}{\sigma_u^2+\sigma_\varepsilon^2}$$

不同 $i$ 值的观察结果是不相关的。方程 \@ref(eq:18-1) 中的模型可以用矩阵表示法描述为

\begin{equation}
\boldsymbol y=\boldsymbol j \mu+\boldsymbol Z_1 \boldsymbol u+\boldsymbol \varepsilon 
(#eq:18-2)
\end{equation}

其中 $\boldsymbol j$ 是 N × 1 全一向量（$N=\sum_{i=1}^tn_i$），$\boldsymbol Z_1$ 是 N × t 设计矩阵，$\boldsymbol u$ 是假定服从多元正态分布 $N_t(\boldsymbol 0,\sigma^2_u\boldsymbol I_t)$ 的 t × 1 随机向量变量，且 $\boldsymbol \varepsilon$ 是假定服从多元正态分布 $N_N(\boldsymbol 0,\sigma^2_\varepsilon\boldsymbol I_t)$ 的 N × 1 随机向量变量。观测值向量 $\boldsymbol y$ 的协方差阵为

$$\begin{aligned}
\boldsymbol\Sigma& =\mathrm{Var}(\boldsymbol{y})=\mathrm{Var}(\boldsymbol{j} \mu+\boldsymbol{Z}_1\boldsymbol{u}+\boldsymbol{\varepsilon})  \\
&=\mathbf{Z}_1\operatorname{Var}(\boldsymbol{u})\mathbf{Z}_1^{\prime}+\operatorname{Var}(\boldsymbol{\varepsilon}) \\
&=\sigma_u^2\mathbf{Z}_1\mathbf{Z}_1^{\prime}+\boldsymbol{\sigma}_\boldsymbol{\varepsilon}^2\mathbf{I}_N
\end{aligned}$$

$yij$ 的方差是 $\boldsymbol\Sigma$ 的对角元，$yij$ 对子之间的协方差是 $\boldsymbol\Sigma$ 的非对角元。

等价地，模型可以写为

$$
\begin{bmatrix}y_{11}\\y_{12}\\\vdots\\y_{1n_1}\\y_{21}\\y_{22}\\\vdots\\y_{2n_2}\\\vdots\\y_{t1}\\y_{t2}\\\vdots\\y_{tn_t}\end{bmatrix}=
\begin{bmatrix}1\\1\\\vdots\\1\\1\\1\\\vdots\\1\\\vdots\\1\\1\\\vdots\\1\end{bmatrix}\mu+
\begin{bmatrix}
1&0&\cdots&0\\
1&0&\cdots&0\\
\vdots&\vdots&\cdots&\vdots\\
1&0&\cdots&0\\
0&1&\cdots&0\\
0&1&\cdots&0\\
\vdots&\vdots&\cdots&\vdots\\
0&1&\cdots&0\\
\vdots&\vdots&\cdots&\vdots\\
0&0&\cdots&1\\
0&0&\cdots&1\\
\vdots&\vdots&\cdots&\vdots\\
0&0&\cdots&1\\\end{bmatrix}\begin{bmatrix} \mu_1\\\mu_2\\\vdots\\\mu_t\end{bmatrix}+
\begin{bmatrix}\varepsilon_{11}\\\varepsilon_{12}\\\vdots\\\varepsilon_{1n_1}\\\varepsilon_{21}\\\varepsilon_{22}\\\vdots\\\varepsilon_{2n_2}\\\vdots\\\varepsilon_{t1}\\\varepsilon_{t2}\\\vdots\\\varepsilon_{tn_t}\end{bmatrix}
$$

[^generalrandomeffectmodel]: 原文：The general random effects model will have r random components representing the main effects and interactions for the random effect factors of the treatment structure and for those factors used to describe the design structure as well as possible interactions between components of the design and treatment structures used to describe the necessary error terms.

一般随机效应模型将具有 $r$ 个随机分量，代表了处理结构中随机效应因子的主效应和交互作用，同时也代表了用于描述设计结构的因子以及设计与处理结构分量之间可能存在的交互作用，此交互作用用于描述必要误差项[^generalrandomeffectmodel]。一般随机效应模型包括由 $\mu$ 表示的总体均值参数由向量 $\boldsymbol \varepsilon$ 表示的残差或最小尺寸实验单元的误差。一般随机效应模型可以用矩阵表示为

\begin{equation}
\boldsymbol y=\boldsymbol j \mu+\boldsymbol Z_1\boldsymbol u_1+\boldsymbol Z_2\boldsymbol u_2+\cdots+\boldsymbol Z_r\boldsymbol u_r+\boldsymbol\varepsilon
(#eq:18-3)
\end{equation}

其中 $\mu_s(s=1,2,\cdots,r)$ 表示随机效应以及 $\boldsymbol\varepsilon$ 表示残差，其中假设所有这些随机变量独立分布，且假设它们的边际分布为

$$\boldsymbol u_1\thicksim N(0,\sigma_1^2\boldsymbol I_{t_1}),\quad \boldsymbol u_2\thicksim N(0,\sigma_2^2\boldsymbol I_{t_2}),\ldots,\quad \boldsymbol u_r\thicksim N(0,\sigma_r^2\boldsymbol I_{t_r}),\quad\mathrm{and}\quad\boldsymbol\varepsilon\thicksim N(0,\sigma_\varepsilon^2\boldsymbol I_N)$$

$N$ 为数据向量 $\boldsymbol y$ 中的观测总数，$\boldsymbol Z_i$ 为对应于第 $i$ 个随机效应向量的 $N × t_i$ 设计矩阵。这种一般随机效应模型可以表示为

$$\begin{aligned}&\boldsymbol{y}=\boldsymbol{j}_{\Tiny N} \mu+\boldsymbol{Z}\boldsymbol{u}+\boldsymbol{\varepsilon}\quad\mathrm{where~}\boldsymbol{Z}=\boldsymbol{[Z}_1,\boldsymbol{Z}_2,\ldots,\boldsymbol{Z}_r],\\\\&\boldsymbol{u}=\begin{bmatrix}\boldsymbol{u}_1\\\boldsymbol{u}_2\\\vdots\\\boldsymbol{u}_r\end{bmatrix}\mathrm{~and~}\quad\mathrm{Var}(\boldsymbol{u})=\begin{bmatrix}{\sigma}_1^2\boldsymbol{I}_{t_1}&0&0&0\\0&{\sigma}_2^2\boldsymbol{I}_{t_2}&0&0\\0&0&\ddots&0\\0&0&0&{\sigma}_r^2\boldsymbol{I}_{t_r}\end{bmatrix}\end{aligned}$$

因此

$$\operatorname{Var}(\boldsymbol y)=\boldsymbol Z^{\prime}\operatorname{Var}(\boldsymbol u)\boldsymbol Z+\sigma_\varepsilon^2\boldsymbol I_N$$

数据向量 $\boldsymbol y$ 的协方差阵为

$$\begin{aligned}
\boldsymbol \Sigma& =\operatorname{Var}(\boldsymbol y)=\operatorname{Var}(\boldsymbol j \mu+\boldsymbol Z_1\boldsymbol u_1+\boldsymbol Z_2\boldsymbol u_2+\cdots+\boldsymbol Z_r\boldsymbol u_r+\boldsymbol\varepsilon)  \\
&=\mathbf{Z}_1\mathrm{Var}(\boldsymbol{u}_1)\mathbf{Z}_1^{\prime}+\mathbf{Z}_2\mathrm{Var}(\boldsymbol{u}_2)\mathbf{Z}_2^{\prime}+\cdots+\mathbf{Z}_r\mathrm{Var}(\boldsymbol{u}_r)\mathbf{Z}_r^{\prime}+\mathrm{Var}(\boldsymbol{\varepsilon}) \\
&=\sigma_1^2\mathbf{Z}_1\mathbf{Z}_1^{\prime}+\sigma_2^2\mathbf{Z}_2\mathbf{Z}_2^{\prime}+\cdots+\sigma_r^2\mathbf{Z}_r\mathbf{Z}_r^{\prime}+\sigma_\varepsilon^2\mathbf{I}_N
\end{aligned}$$

这种一般形式的随机效应模型可以在许多情况下使用，以帮助识别数据收集系统中的变异源。分析这些变异来源的一种方法是计算一组与每个随机效应对应的平方和，然后确定每个平方和估计的参数函数。通过获取效应平方和的期望值，可以评估由每个效应平方和所估计的方差分量函数。接下来的章节将使用这种矩阵形式的一般随机效应模型来描述涉及观测值的平方和期望值的评估方法。

## 计算期望均方 {#sec18-3}

在对随机效应模型进行方差分析时，平方和的期望值涉及到方差分量。对于给定的模型，至少有两种方法可用于评估期望均方（记住，均方是一个平方和除以其自由度）。第一种方法是利用模型假设代数地求解期望值，第二种方法是通过计算机算法来评估期望值。代数方法是通过将其应用于从具有单向随机效应处理结构的模型分析中获得的平方和来呈现的。计算机算法方法则从一般性角度进行了讨论，并通过更复杂的例子进行了演示。

### 代数方法 {#sec18-3-1}

方程 \@ref(eq:18-1) 中的单向随机效应模型中有两个方差分量，因此在分析中使用两个平方和来描述响应的变异性。通常计算的两个平方和是随机效应水平内的平方和，记作 $Q_1$ ，以及随机效应水平之间的平方和，记作 $Q_2$。 对于单向随机效应模型，这些平方和由下式给出

$$\begin{aligned}Q_1&=\sum_{i=1}^t\sum_{j=1}^{n_i}(y_{ij}-\bar{y}_{i\cdot})^2=\sum_{i=1}^t\sum_{j=1}^{n_i}y_{ij}^2-\sum_{i=1}^tn_i\bar{y}_{i\cdot}^2\\\\Q_2&=\sum_{i=1}^tn_i(\bar{y}_{i\cdot}-\bar{y}_{\cdot\cdot})^2=\sum_{i=1}^tn_i\bar{y}_{i\cdot}^2-n_{\cdot\cdot}\bar{y}_{\cdot\cdot}^2\end{aligned}$$

根据模型 \@ref(eq:18-1) 的随机变量，$Q_1,Q_2$ 中的量可以表示为

\begin{align}
y_{ij}&=\mu+\mu_i+\varepsilon_{ij}\\
\bar y_{i\cdot}&=\mu+\mu_i+\bar\varepsilon_{i\cdot}\\
\bar y_{\cdot\cdot}&=\mu+\tilde \mu_{\cdot}+\bar\varepsilon_{\cdot\cdot}\quad\text{where}\,\,\tilde\mu_{\cdot}=\frac1N\sum_{i=1}^tn_i\mu_i
(#eq:18-4)
\end{align}

将方程 \@ref(eq:18-4) 中的项代入 $Q_1$，表达式变为

$$Q_1=\sum_{i=1}^t\sum_{j=1}^{n_i}[(\mu+u_i+\varepsilon_{ij})-(\mu+u_i+\bar{\varepsilon}_{i\cdot})]^2=\sum_{i=1}^t\sum_{j=1}^{n_i}(\varepsilon_{ij}-\bar{\varepsilon}_{i\cdot})^2$$

使用 $\boldsymbol \varepsilon$ 的分布性质可以评估$Q_1$ 的期望值

$$\begin{aligned}
E(Q_1)& =\sum_{i=1}^t\sum_{j=1}^{n_i}E({\varepsilon}_{ij}-\bar{{\varepsilon}}_i.)^2=\sum_{i=1}^t\sum_{j=1}^{n_i}[E({\varepsilon}_{ij}^2)+E(\bar{{\varepsilon}}_{i\cdot}^2)-2E({\varepsilon}_{ij}\bar{{\varepsilon}}_{i\cdot})]\quad\text{(by squaring)}  \\
&=\sum_{i=1}^t\sum_{j=1}^{n_i}\left(\sigma_\varepsilon^2+\frac{\sigma_\varepsilon^2}{n_i}-2\frac{\sigma_\varepsilon^2}{n_i}\right)\quad\mathrm{using~}E(\varepsilon_{ij}^2)=\sigma_\varepsilon^2\quad\mathrm{and~}\quad E(\bar{\varepsilon}_{i\cdot}^2)=\frac{\sigma_\varepsilon^2}{N_i} \\
&=\sum_{i=1}^t\sum_{j=1}^{n_i}\frac{(n_i-1)\sigma_\varepsilon^2}{n_i}=\sigma_\varepsilon^2\sum_{i=1}^t\left(n_i-1\right)=\left(N-t\right)\sigma_\varepsilon^2
\end{aligned}$$

使用此表达式作为平方和的期望值，期望均方为

$$E(\text{mean square of Q}_1)=E{\left(\frac{Q_1}{N-t}\right)}=\sigma_\varepsilon^2$$

将方程 \@ref(eq:18-4) 中的表达式代入定义 $Q_2$ 的方程得到

$$\begin{gathered}
Q_2 \begin{aligned}=\sum_{i=1}^tn_i(\mu+\mu_i+\varepsilon_{ij}-\mu-\tilde{u}_{\cdot}+\bar{\varepsilon}_{i\cdot})^2\end{aligned} \\
=\sum_{i=1}^tn_i[(u_i-\tilde{u}_{\cdot})+(\varepsilon_{ij}-\bar{\varepsilon}_{i\cdot})]^2 
\end{gathered}$$

$Q_2$ 的期望为

$$\begin{aligned}
E(Q_{2})& \begin{aligned}=\sum_{i=1}^tn_i[E(u_i-\tilde{u}_{\cdot})^2+E(\bar{{\varepsilon}}_{i\cdot}-\bar{{\varepsilon}}_{\cdot\cdot})^2]\end{aligned}  \\
&\begin{aligned}=\sum_{i=1}^tn_i[E(u_i)^2+E(\tilde{u}_{\cdot})^2-2E(u_i\tilde{u}_{\cdot})+E(\bar{\varepsilon}_{i\cdot})^2+E(\bar{\varepsilon}_{\cdot\cdot})^2-2E(\bar{\varepsilon}_{i\cdot}\bar{\varepsilon}_{\cdot\cdot})]\end{aligned}
\end{aligned}$$

为了简化这一点，$Q_2$ 的期望分为两部分进行评估。第一部分评估了涉及 $\boldsymbol \varepsilon$ 的期望。与均值 $\bar \varepsilon_{i\cdot}$ 和 $\bar \varepsilon_{\cdot\cdot}$ 相关的分布为

$$\bar{{\varepsilon}}_{i\cdot}\sim N{\left(0,\frac{{\sigma}_\varepsilon^2}{n_i}\right)}\quad\mathrm{and}\quad\bar\varepsilon_{\cdot}\sim N{\left(0,\frac{{\sigma}_\varepsilon^2}N\right)}$$

那么 $Q_2$ 中涉及 $\boldsymbol \varepsilon$ 部分的期望是

$$\begin{aligned} &\sum_{i=1}^{t}n_{i}\left[E(\bar{{\varepsilon}}_{i\cdot})^{2}+E(\bar{{\varepsilon}}_{\cdot\cdot})^{2}-2E(\bar{{\varepsilon}}_{i\cdot}\bar{{\varepsilon}}_{\cdot\cdot})\right]  \\
&=\sum_{i=1}^tn_i\left\{\sigma_\varepsilon^2{\left(\frac1{n_i}\right)}+\frac{\sigma_\varepsilon^2}N-2E{\left[\bar{{\varepsilon}}_{i\cdot}\frac{\sum_{i^{\prime}=1}^tn_{i^{\prime}}\bar{{\varepsilon}}_{i^{\prime}\cdot}}N\right]}\right\} \\
&=\sum_{i=1}^tn_i\left\{\sigma_\varepsilon^2{\left(\frac1{n_i}\right)}+\frac{\sigma_\varepsilon^2}N-2\frac{n_i}NE(\bar{\varepsilon}_i)^2\right\} \\
&=\sum_{i=1}^tn_i\left\{\sigma_\varepsilon^2{\left(\frac1{n_i}\right)}+\frac{\sigma_\varepsilon^2}N-2\frac{n_i}N{\left(\frac{\sigma_\varepsilon^2}{n_i}\right)}\right\} \\
&=\sum_{i=1}^tn_i\left(\frac1{n_i}-\frac1N\right){\sigma}_\varepsilon^2 \\
&=\sum_{i=1}^t\left(1-\frac{n_i}N\right)\sigma_\varepsilon^2=\left(t-1\right)\sigma_\varepsilon^2
\end{aligned}$$

为了评估 $Q_2$ 期望中涉及 $u_i$ 的部分，设 $\tilde{u}_{\cdot}=\sum_{i=1}^t(n_iu_i/N)$. 由于 $u_i$ 独立且均值为零

$$\mathrm{Var}(\tilde{u}_\cdot)=E(\tilde{u}_\cdot)^2=\sum_{i=1}^t\left(\frac{n_i^2}{N^2}\right){\sigma}_u^2=\left(\frac1{N^2}\right){\sum}_{i=1}^t(n_i^2){\sigma}_u^2$$

由于对于 $i\ne i^\prime$ 有 $E(u_iu_{i^\prime})$,$u_i$ 和 $\tilde{u}_\cdot$ 之间的协方差为

$$\mathrm{Cov}(\tilde{u}_\cdot,u_i)=E(\tilde{u}_\cdot u_i)=\left(\frac{n_i}N\right)E(u_i)^2=\left(\frac{n_i}N\right)\sigma_u^2$$

把这些部分放在一起，$Q_2$ 的期望中涉及 $u_i$ 的部分是

$$\begin{aligned}
\sum_{i=1}^{t} &n_{i}\left[E(u_{i})^{2}+E(\tilde{u}_{\cdot})^{2}-2E(u_{i}\tilde{u}_{\cdot})\right]  \\
&=\sum_{i=1}^tn_i\sigma_u^2+\sum_{i=1}^t\frac{n_i\sigma_u^2}{N^2}\left(\sum_{i=1}^tn_i^2\right)-2\sum_{i=1}^t\left(\frac{n_i^2}N\right)\sigma_u^2 \\
&\left.=\left(\sum_{i=1}^tn_i{\sigma}_u^2\right.-\frac{\sum_{i=1}^tn_i^2}N{\sigma}_u^2\right)=\left(N-\frac{\sum_{i=1}^tn_i^2}N\right){\sigma}_u^2
\end{aligned}$$

最后，将这两个部分合在一起，$Q_2$ 的期望为

$$\begin{aligned}E(Q_2)=(t-1)\sigma_\varepsilon^2+\left(N-\frac{\sum_{i=1}^tn_i^2}N\right)\sigma_u^2\end{aligned}$$

$Q_2$ 的期望均方为

$$E{\left(\frac{Q_2}{t-1}\right)}=\sigma_\varepsilon^2+{\left(\frac1{t-1}\right)}{\left(N-\frac{\sum_{i=1}^tn_i^2}N\right)}\sigma_u^2$$

表 \@ref(tab:table18-1) 中给出了单向随机效应模型的方差分析表，包括来源、自由度、平方和以及期望均方。

<table>
<caption>(\#tab:table18-1)单向随机效应模型的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2018.1.png">
</td>
  </tr>
</tbody>
</table>

### Hartley 综合法的计算 {#sec18-3-2}

Hartley (1967) 提出了一种技术，用于评估特定设计和平方和的均方期望值，他称之为综合法 (Hartley's method of synthesis). 为了描述并随后应用这一技术，我们首先用矩阵表示从一般随机效应模型 \@ref(eq:18-3) 计算出的平方和的期望值。任何通过公式计算或通过计算机软件提取的平方和，总是可以表示为数据的**二次型** (quadratic form)

\begin{equation}
Q=\boldsymbol{y'Ay}
(#eq:18-5)
\end{equation}

其中 $\boldsymbol y$ 是观测向量，$\boldsymbol A$ 是适当选择的对称常数矩阵，称为**二次型的矩阵** (the matrix of the quadratic form) (Graybill, 1976). 例如，n 个观测值的向量的样本方差为

$$s^2=\sum_{i=1}^n\frac{(\boldsymbol y_i-\bar{\boldsymbol y}_{\cdot})^2}{n-1}=\boldsymbol y^{\prime}\biggl[\frac1{n-1}{\biggl(\boldsymbol I_n-\frac1n\boldsymbol J_n\biggr)}\biggr]\boldsymbol{y}$$
其中 $\boldsymbol I_n$ 为 n × n 单位阵，$\boldsymbol J_n$ 为 n × n 全一阵，二次型 \boldsymbol{y'Ay} 的矩阵为

$$\boldsymbol A=\frac1{n-1}{\biggl(\boldsymbol I_n-\frac1n\boldsymbol J_n\biggr)}$$

针对不同的模型，总存在某些特定的选择使得矩阵 $\boldsymbol A$ 产生所需的平方和，幸运的是，正如我们很快将看到的那样，无需了解矩阵 $\boldsymbol A$ 中的元素具体是什么，甚至不必知道如何确定矩阵 $\boldsymbol A$ 的元素。你只需要知道矩阵 $\boldsymbol A$ 是存在的。

对于一般随机效应模型 \@ref(eq:18-3) 及其相应的协方差阵 $\boldsymbol \Sigma$，二次型的期望 (Graybill, 1976) 为

\begin{equation}
E(\boldsymbol{y^{\prime}Ay})=\mathrm{Tr}[\boldsymbol{A\Sigma}]+\frac12\mu^2\boldsymbol{j}_n^{\prime}\boldsymbol A\boldsymbol{j}_n
(#eq:18-6)
\end{equation}

其中 $\operatorname{Tr}[\boldsymbol{B}]=\sum_{i=1}^nb_{ii}$ 且 $b_{ii},i=1,2,\cdots,n$ 表示方阵 $\boldsymbol B$ 的对角元。方差分析中的平方和被构造为使得 $\mu^2\boldsymbol{j}_n^{\prime}\boldsymbol A\boldsymbol{j}_n$. 因此，平方和的期望不取决于 $\mu^2\boldsymbol{j}_n^{\prime}\boldsymbol A\boldsymbol{j}_n$，并且由下式给出

\begin{equation}
E(\boldsymbol{y^{\prime}Ay})=\mathrm{Tr}[\boldsymbol{A\Sigma}]
(#eq:18-7)
\end{equation}

方程 \@ref(eq:18-3) 中的一般随机效应模型的协方差阵为

$$\boldsymbol{\Sigma}={\sigma}_1^2\boldsymbol{Z}_1\boldsymbol{Z}_1^{\prime}+{\sigma}_2^2\boldsymbol{Z}_2\boldsymbol{Z}_2^{\prime}+\cdots+{\sigma}_r^2\boldsymbol{Z}_r\boldsymbol{Z}_r^{\prime}+{\sigma}_\varepsilon^2\boldsymbol{I}_N$$

因此，二次型 $\boldsymbol{y^{\prime}Ay}$ 的期望是

$$\begin{aligned}E(\boldsymbol{y^{\prime}Ay})&=\mathrm{Tr}[\boldsymbol A\boldsymbol{\Sigma}]=\mathrm{Tr}[\boldsymbol A(\sigma_1^2\boldsymbol Z_1\boldsymbol Z_1^{\prime}+\sigma_2^2\boldsymbol Z_2\boldsymbol Z_2^{\prime}+\cdots+\sigma_r^2\boldsymbol Z_r\boldsymbol Z_r^{\prime}+\sigma_\varepsilon^2\boldsymbol I_N)]\\&=\sigma_1^2\mathrm{~Tr}[\boldsymbol A\boldsymbol Z_1\boldsymbol Z_1^{\prime}]+\sigma_2^2\mathrm{~Tr}[\boldsymbol A\boldsymbol Z_2\boldsymbol Z_2^{\prime}]+\cdots+\sigma_r^2\mathrm{~Tr}[\boldsymbol A\boldsymbol Z_r\boldsymbol Z_r^{\prime}]+\sigma_\varepsilon^2\mathrm{~Tr}[\boldsymbol A]\end{aligned}$$

因此，$\sigma^2_{\varepsilon}$ 的系数为 $\mathrm{Tr}[\boldsymbol A]$，这等于与平方和 $\boldsymbol{y^{\prime}Ay}$ 相关的自由度。对于 $s=1,2,\cdots,r$，$\sigma^2_{\varepsilon}$ 的系数为 $\mathrm{Tr}[\boldsymbol A\boldsymbol Z_s\boldsymbol Z_s^{\prime}]$.

迹运算符的一个性质是，$\mathrm{Tr}[\boldsymbol A\boldsymbol Z_s\boldsymbol Z_s^{\prime}]=\mathrm{Tr}[\boldsymbol Z_s^{\prime}\boldsymbol A \boldsymbol Z_s]$ 其中 $\mathrm{Tr}[\boldsymbol Z_s^{\prime}\boldsymbol A \boldsymbol Z_s]$ 是 $\boldsymbol Z_s^{\prime}\boldsymbol A \boldsymbol Z_s$ 对角元之和，即 $\mathrm{Tr}[\boldsymbol Z_s^{\prime}\boldsymbol A \boldsymbol Z_s]=\sum_{j=1}^{t_s}\boldsymbol z^\prime_{sj}\boldsymbol A \boldsymbol z_{sj}$ 因为 $\boldsymbol Z_s$ 中有 $t_s$ 列。但是，$\boldsymbol z^\prime_{sj}\boldsymbol A \boldsymbol z_{sj}$ 与 $\boldsymbol{y^{\prime}Ay}$ 是相同的平方和，只不过这里是以列向量 $\boldsymbol z^\prime _{sj}$ 代替 $\boldsymbol y$ 作为数据使用。因此，如果你有一个计算 $\boldsymbol{y^{\prime}Ay}$ 的方程或计算机程序，它也可以用来计算 $\boldsymbol z^\prime_{sj}\boldsymbol A \boldsymbol z_{sj}$. 因此，在 $\boldsymbol{y^{\prime}Ay}$ 的期望中，$\sigma^2_s$ 的系数是 $\boldsymbol{z'}_{s1}\boldsymbol{Az}_{s1}+\boldsymbol{z'}_{s2}\boldsymbol{Az}_{s2}+\cdots+\boldsymbol{z'}_{st_s}\boldsymbol{Az}_{st_s}$.

如果 $\boldsymbol A$ 的元素是已知的，那么可以明确地评估上述平方和。然而，由于通常会使用计算机代码来计算 $\boldsymbol{y^{\prime}Ay}$，此时 $\boldsymbol A$ 未知，那么每个 $\boldsymbol z^\prime_{sj}\boldsymbol A \boldsymbol z_{sj}$ 可以通过让计算机计算以列向量 $\boldsymbol z_{sj}$ 作为数据向量（替代 $\boldsymbol y$）所得的平方和来进行计算。因此，必须根据模型中的所有随机效应计算设计矩阵每列的平方和 $[\boldsymbol Z_1,\boldsymbol Z_2,\cdots,\boldsymbol Z_r]$，就好像它们是数据一样，然后平方和 $\boldsymbol{y^{\prime}Ay}$ 的期望被评估为

$$E(\boldsymbol{y^{\prime}Ay})=v\sigma_\varepsilon^2+\left(\sum_{q=1}^{t_1}\boldsymbol z_{1q}^{\prime}\boldsymbol A\boldsymbol z_{1q}\right)\sigma_1^2+\left(\sum_{q=1}^{t_2}\boldsymbol z_{2q}^{\prime}\boldsymbol A\boldsymbol z_{2q}\right)\sigma_2^2+\cdots+\left(\sum_{q=1}^{t_r}\boldsymbol z_{rq}^{\prime}\boldsymbol A\boldsymbol z_{rq}\right)\sigma_r^2$$

其中 $v$ 是与 $\boldsymbol{y^{\prime}Ay}$ 相关的自由度。

为了更好地说明综合法的概念，我们重新计算单向随机效应模型的平方和 $Q_1$ 和 $Q_2$ 的期望值。首先，我们会采用一个特定模型，其中 $t=4$，并且对每个 $i$ 都有 $n_i=4$，以此来展示如何计算期望值；接着为一般单向随机效应模型计算期望值。

在具有四个重复的完全随机设计结构中，用矩阵形式描述四个随机选取的小麦品种产量的模型为：

$$\begin{bmatrix}y_{11}\\y_{12}\\y_{13}\\y_{14}\\y_{21}\\y_{22}\\y_{23}\\y_{24}\\y_{31}\\y_{32}\\y_{33}\\y_{34}\\y_{41}\\y_{42}\\y_{43}\\y_{44}\end{bmatrix}=\begin{bmatrix}1\\1\\1\\1\\1\\1\\1\\1\\1\\1\\1\\1\\1\\1\\1\\1\end{bmatrix}\mu+\begin{bmatrix}1&0&0&0\\1&0&0&0\\1&0&0&0\\1&0&0&0\\0&1&0&0\\0&1&0&0\\0&1&0&0\\0&1&0&0\\0&0&1&0\\0&0&1&0\\0&0&1&0\\0&0&1&0\\0&0&0&1\\0&0&0&1\\0&0&0&1\\0&0&0&1\end{bmatrix}+\begin{bmatrix}\boldsymbol u_1\\\boldsymbol u_2\\\boldsymbol u_3\\\boldsymbol u_4\end{bmatrix}+\begin{bmatrix}
\varepsilon_{11}\\\varepsilon_{12}\\\varepsilon_{13}\\\varepsilon_{14}\\
\varepsilon_{21}\\\varepsilon_{22}\\\varepsilon_{23}\\\varepsilon_{24}\\
\varepsilon_{31}\\\varepsilon_{32}\\\varepsilon_{33}\\\varepsilon_{34}\\
\varepsilon_{41}\\\varepsilon_{42}\\\varepsilon_{43}\\\varepsilon_{44}\end{bmatrix}$$

或 $\boldsymbol{y}=\boldsymbol{j}_{16} \mu+[\boldsymbol z_1\,\,\boldsymbol z_2\,\,\boldsymbol z_3\,\,\boldsymbol z_4]\boldsymbol{u}+\boldsymbol{\varepsilon}$.

[^withinsumofsquares]: 是前文 "the sum of squares within levels of the random effect" 的简写。

组内平方和 (within sum of squares[^withinsumofsquares]) $Q_1$ 的期望值为

$$E(Q_1)=E(\boldsymbol{y'}A_w\boldsymbol{y})=\boldsymbol{\sigma}_u^2\left[\boldsymbol{\sum}_{j=1}^4\boldsymbol{z'}_j\boldsymbol{A}_w\boldsymbol{z}_j\right]+12{\sigma}_\varepsilon^2$$

其中 12 是与 $Q_1$ 相关的自由度，$\boldsymbol{A}_w$ 表示组内平方和 $Q_1$ 二次型的矩阵。为了获得 $\sigma^2_u$ 的系数，使用 $\boldsymbol z_1$ 作为数据计算组内平方和 $Q_1$，使用 $\boldsymbol z_2$ 作为数据计算 $Q_1$，使用 $\boldsymbol z_3$ 作为数据计算 $Q_1$，使用 $\boldsymbol z_4$ 作为数据计算 $Q_1$. 对于列 $\boldsymbol z_1$ 的组内平方和

$$\begin{aligned}Q_1(\boldsymbol z_1)=\sum_{i=1}^4\sum_{j=1}^4z_{1ij}^2-4\sum_{i=1}^4\bar{z}_{1i\cdot}^2&=4-4(1)=0\end{aligned}$$

同样，$Q_1(\boldsymbol z_2),Q_1(\boldsymbol z_3),Q_1(\boldsymbol z_4)$ 也为 0. 这意味着 $\sigma^2_u$ 在 $E(Q_1)$ 中的系数为零，因此 $E(Q_1)=12{\sigma}_\varepsilon^2$.

[^betweensumofsquares]: 是前文 "the sum of squares between the levels of the random effect" 的简写。

组间平方和 (between sum of squares[^betweensumofsquares]) $Q_2$ 的期望值为

$$E(Q_2)=E(\boldsymbol{y'}\boldsymbol A_B\boldsymbol{y})={\sigma}_u^2\biggl[\boldsymbol{\sum}_{j=1}^4\boldsymbol{z'}_{j}\boldsymbol A_B\boldsymbol{z}_j\biggr]+3{\sigma}_\varepsilon^2$$

其中

$$Q_2=\sum_{i=1}^tn_i(\bar{y}_{i\cdot}-\bar{y}_{\cdot\cdot})^2=\sum_{i=1}^tn_i\bar{y}_{i\cdot}^2-n_{\cdot\cdot}\bar{y}_{\cdot\cdot}$$

其中 3 是与 $Q_2$ 相关的自由度，$\boldsymbol{A}_B$ 表示组间平方和 $Q_2$ 二次型的矩阵。要计算 $\sigma^2_u$ 的系数，请使用每列 $\boldsymbol z_1,\boldsymbol z_2,\boldsymbol z_3,\boldsymbol z_4$ 计算组间平方和 $Q_2$ 的值，并将这四个值相加。使用 $\boldsymbol z_1$ 作为数据计算的平方和之间的值为

$$Q_2(\boldsymbol z_1){=}4\sum_{i=1}^4\bar{z}_{1i\cdot}^2-16\bar{z}_{\cdot\cdot\cdot}^2=4(1^2+0^2+0^2+0^2){-}16(0.25)^2{=}3$$

$Q_2(\boldsymbol z_2),Q_2(\boldsymbol z_3),Q_2(\boldsymbol z_4)$ 的值也等于 3；因此，$\sigma^2_u$ 的系数为 $3+3+3+3=12$. 使用这些值，$E(Q_2)$ 的表达式为

$$E(Q_2)=12\sigma_u^2+3\sigma_\varepsilon^2$$

接下来，使用 Hartley 的综合法来计算方程 \@ref(eq:18-1) 中的一般单向随机效应模型的组内和组间平方和的期望值。模型的矩阵形式为

$$\begin{bmatrix}y_{11}\\\vdots\\y_{1n_1}\\y_{21}\\\vdots\\y_{2n_2}\\\vdots\\y_{t1}\\\vdots\\y_{tn_l}\end{bmatrix}=\begin{bmatrix}1\\\vdots\\1\\1\\\vdots\\1\\\vdots\\1\\\vdots\\1\end{bmatrix} \mu+\begin{bmatrix}1&0&\cdots&0\\\vdots&\vdots&\cdots&\vdots\\1&0&\cdots&0\\0&1&\cdots&0\\\vdots&\vdots&\cdots&\vdots\\0&1&\cdots&0\\\vdots&\vdots&\cdots&\vdots\\0&0&\cdots&1\\\vdots&\vdots&\cdots&\vdots\\0&0&\cdots&1\end{bmatrix}\begin{bmatrix}\boldsymbol{u}_1\\\boldsymbol{u}_2\\\vdots\\\boldsymbol{u}_t\end{bmatrix}+\boldsymbol{\varepsilon}$$

或

$$\boldsymbol y=\boldsymbol j_{n} \mu+[\boldsymbol z_1,\boldsymbol z_2,\ldots,\boldsymbol z_t]\boldsymbol u+\boldsymbol \varepsilon $$

组内或误差平方和为

$$Q_1=\sum_{i=1}^t\sum_{j=1}^{n_i}y_{ij}^2-\sum_{i=1}^tn_i\bar{y}_{i\cdot}^2=\boldsymbol{y'A}_w\boldsymbol{y}$$

其期望为

$$\begin{aligned}E(Q_1)=\sigma_u^2\sum_{i=1}^t\boldsymbol z_i'\boldsymbol A_w\boldsymbol z_i+(N-t)\sigma_\varepsilon^2\end{aligned}$$

其中 $N=\sum_{i=1}^tn_i$ 以及 $\boldsymbol A_w$ 是组内平方和二次型的矩阵。当使用第一列 $\boldsymbol z_1$ 作为数据时，$Q_1$ 的值为

$$Q_1(\boldsymbol z_1)=\sum_{i=1}^t\sum_{j=1}^{n_i} z_{1ij}^2-\sum_{i=1}^tn_i\bar{ z}_{1i\cdot}^2=n_1-\left[n_1(1)+n_2(0)+\cdots+n_t(0)\right]=0$$

同样，$Q_1(\boldsymbol z_2),Q_1(\boldsymbol z_3),Q_1(\boldsymbol z_4)$ 也为 0. 因此，在 $\sigma^2_u$ 在 $E(Q_1)$ 中的系数为零，这意味着 $E(Q_1)=(N-t)\sigma^2_{\varepsilon}$，并且其中存在与 $Q_1$ 相关的 N-t 个自由度。

组间平方和为

$$Q_2=\sum_{i=1}^tn_i\bar{y}_{i\cdot}^2-n_{\cdot\cdot}\bar{y}_{\cdot\cdot}^2=\boldsymbol y^{\prime}\boldsymbol A_B\boldsymbol y$$

其期望为

$$E(Q_2)=\sigma_u^2\sum_{i=1}^t\boldsymbol z_i \boldsymbol A_B\boldsymbol z_i+(t-1)\sigma_\varepsilon^2$$

其中 $\boldsymbol A_B$ 是组间平方和二次型的矩阵。当使用第一列 $\boldsymbol z_1$ 作为数据时，$Q_2$ 的值为

$$Q_2(\boldsymbol z_1){=}\sum_{i=1}^tn_i\bar{z}_{1i\cdot}^2-N\bar{z}_{1\cdot\cdot}^2$$

对于列 $\boldsymbol z_1$ 有 $\bar z_{11\cdot}=1,\bar z_{12\cdot}=\bar z_{13\cdot}=\cdots \bar z_{1t\cdot}=0$ 以及 $\bar z_{1\cdot\cdot}=n_1/N$. 因此


$$Q_2(z_1)=n_1(1)^2+n_2(0)^2+\cdots+n_t(0)^2-N(n_1/N)^2=\frac{n_1-n_1^2}N$$

类似地，使用 $\boldsymbol Z$ 的其他列作为数据的 $Q_2$ 的值为

$$Q_2(\boldsymbol z_2)=\frac{n_2-n_2^2}N,\quad Q_2(\boldsymbol z_3)=\frac{n_3-n_3^2}N,\ldots,Q_2(\boldsymbol z_t)=\frac{n_t-n_t^2}N$$

合并这些结果，组间平方和的期望值为

$$E(Q_2)=(t-1)\sigma_\varepsilon^2+\left[\sum_{i=1}^tQ_2(\boldsymbol z_i)\right]\sigma_u^2=(t-1)\sigma_\varepsilon^2+\left[\sum_{i=1}^tn_i^2\right]\sigma_u^2$$

通过 Hartley 综合法获得的 $Q_1$ 和 $Q_2$ 的期望值与使用代数技术获得的期望值相等。

接下来，使用 Hartley 综合法来评估用于描述在完全随机设计结构中具有交互作用的双向随机效应处理结构的模型的平方和期望值。在这种情况下，行处理和列处理都是随机效应。该实验的数据可建模为

$$y_{ijk}=\mu+a_i+b_j+c_{ij}+\varepsilon_{ijk},\quad i=1,2,\ldots,s,j=1,2,\ldots,t,\quad k=1,2,\ldots,n_{ij}$$

其中 $a_i,i=1,2,\cdots,s$ 表示与行对应的随机效应，其具有 $i.i.d.\,N(0,\sigma^2_a)$ 分布；$b_i,i=1,2,\cdots,t$ 表示与列对应的随机效应，其具有 $i.i.d.\,N(0,\sigma^2_b)$ 分布；$c_{ij}$ 表示行-列组合的随机效应，其具有 $i.i.d.\,N(0,\sigma^2_c)$ 分布，$\varepsilon_{ijk}$ 实验单元误差，其具有 $i.i.d.\,N(0,\sigma^2_\varepsilon)$ 分布。图 \@ref(fig:figure18-1) 中的示意图代表来自完全随机设计结构中的不均衡双向处理结构的数据，其中每个单元格有两到三个观测结果。

<div class="figure" style="text-align: center">
<img src="figure/figure 18.1.png" alt="双向随机效应处理结构示例" width="781" />
<p class="caption">(\#fig:figure18-1)双向随机效应处理结构示例</p>
</div>

在图 \@ref(fig:figure18-1) 所示数据结构的完全随机设计结构中，这种双向随机效应处理结构的矩阵形式为

$$\begin{bmatrix}y_{111}\\y_{112}\\y_{113}\\y_{121}\\y_{122}\\y_{131}\\y_{132}\\y_{211}\\y_{212}\\y_{221}\\y_{222}\\y_{223}\\y_{231}\\y_{232}\end{bmatrix}=\begin{bmatrix}1\\1\\1\\1\\1\\1\\1\\1\\1\\1\\1\\1\\1\\1\end{bmatrix}\mu +\begin{bmatrix}1&0\\1&0\\1&0\\1&0\\1&0\\1&0\\1&0\\0&1\\0&1\\0&1\\0&1\\0&1\\0&1\\0&1\end{bmatrix}\begin{bmatrix}a_1\\a_2\end{bmatrix}+\begin{bmatrix}1&0&0\\1&0&0\\1&0&0\\0&1&0\\0&1&0\\0&0&1\\0&0&1\\1&0&0\\1&0&0\\0&1&0\\0&1&0\\0&1&0\\0&0&1\\0&0&1\end{bmatrix}\begin{bmatrix}b_1\\b_2\\b_3\end{bmatrix}+\begin{bmatrix}1&0&0&0&0&0\\1&0&0&0&0&0\\1&0&0&0&0&0\\0&1&0&0&0&0\\0&1&0&0&0&0\\0&0&1&0&0&0\\0&0&1&0&0&0\\0&0&0&1&0&0\\0&0&0&1&0&0\\0&0&0&0&1&0\\0&0&0&0&1&0\\0&0&0&0&1&0\\0&0&0&0&0&1\\0&0&0&0&0&1\end{bmatrix}\begin{bmatrix}c_{11}\\c_{12}\\c_{13}\\c_{21}\\c_{22}\\c_{23}\end{bmatrix}+\boldsymbol \varepsilon$$

或 $\boldsymbol{y}=\boldsymbol{j}_{14} \mu+\boldsymbol{Z}_1\boldsymbol{a}+\boldsymbol{Z}_2\boldsymbol{b}+\boldsymbol{Z}_3\boldsymbol{c}+\boldsymbol{\varepsilon}$ 其中模型参数为 $\mu,\sigma_a^2,\sigma_b^2,\sigma_c^2$ 和 $\sigma_\varepsilon^2$

这种不均衡双向随机效应模型的平方和可以用几种不同的方法计算（见第 \@ref(chap9) 章和第 \@ref(chap10) 章）。在该模型的分析中，有四个方差分量，因此需要四个不同的平方和。为了演示综合方法，选择了四个平方和，它们对应于 SSROWS, SSCOLUMNS, SSINTERACTION 和 SSERROR 的均衡情况，但针对不等的样本量进行了修改（对应于 Henderson 的 I 型平方和，在本节稍后讨论）。四个平方和是

$$\begin{aligned}
&Q_1 =\sum_{i=1}^s\frac{y_{i\cdot\cdot}^2}{n_{i\cdot}}-\frac{y_{i\cdot\cdot}^2}{n_{\cdot\cdot}}\mathrm{~(SSROWS)}  \\
&Q_2 =\sum_{j=1}^t\frac{y_{\cdot j\cdot}^2}{n_{\cdot j}}-\frac{y_{\cdot\cdot\cdot}^2}{n_{\cdot\cdot}}\mathrm{~(SSCOLUMNS)}  \\
&Q_3 =\sum_{i=1}^\text{s}{ \sum _ { j = 1 }^{t}\frac{y_{ij\cdot}^2}{n_{ij}}-Q_1-Q_2+\frac{y^2_{\cdot\cdot\cdot}}{n_{\cdot\cdot}}\quad\text{(SSINTERACTION)}}  \\
&Q_4 =\sum_{i=1}^s\sum_{j=1}^t\sum_{k=1}^{n_{ij}}{(y_{ijk}-\bar{y}_{ij\cdot})^2}\quad\text{(SSERROR)} 
\end{aligned}$$

对于图 \@ref(fig:figure18-1) 中的数据结构，

$$Q_2=\frac{y_{\cdot1\cdot}^2}5+\frac{y_{\cdot2\cdot}^2}5+\frac{y_{\cdot3\cdot}^2}4-\frac{y_{\cdot\cdot\cdot}^2}{14}$$

以及，对于 $k_1,k_2,k_3$ 的一些值，$E(Q_2)$ 具有形式 $E(Q_2)＝k_1\sigma_a^2+k_2\sigma_b^2+k_3\sigma_c^2+2\sigma_\varepsilon^2$，其中 2 是与 $Q_2$ 相关的自由度。下一步是使用 Hartley 的综合法来确定 $k_1,k_2,k_3$. 要确定 $k_1$ 的值，请为 $\boldsymbol Z_1$ 的两列中的每一列计算 $Q_2$. 使用 $\boldsymbol Z_1$ 的第一列作为数据的 $Q_2$ 的值为

$$\begin{aligned}Q_2( z_{11})=\frac{3^2}5+\frac{2^2}5+\frac{2^2}4-\frac{7^2}{14}=0.1\end{aligned}$$

并且使用 $\boldsymbol Z_1$ 的第二列作为数据的 $Q_2$ 的值为

$$Q_2(\boldsymbol z_{12})=\frac{2^2}5+\frac{3^2}5+\frac{2^2}4-\frac{5^2}{14}=0.1$$

因此，$k_1＝Q_2( z_{11})+Q_2( z_{11})＝0.1+0.1＝0.2$. 为确定 $k_2$ 的值，计算 $\boldsymbol Z_2$ 每列的 $Q_2$ 为

$$
\begin{aligned}
Q_2(z_{21})=\frac{5^2}5+\frac{0^2}5+\frac{0^2}4-\frac{5^2}{14}=3.214 \\
Q_2(z_{22})=\frac{0^2}5+\frac{5^2}5+\frac{0^2}4-\frac{5^2}{14}=3.214 \\
Q_{2}(z_{23})=\frac{0^{2}}{5}+\frac{0^{2}}{5}+\frac{4^{2}}{4}-\frac{4^{2}}{14}=2.857 
\end{aligned}
$$

$k_2$ 的值为 $k_2=Q_2( z_{21})+Q_2( z_{22})+Q_2( z_{23})=3.214+3.214+2.857=9.285$.

$k_3$ 的值是通过计算 $\boldsymbol Z_3$ 的每列的 $Q_2$ 来获得的。$Q_2$ 的值为

$$\begin{aligned}
Q_2({z}_{311}) =\frac{3^2}5+\frac{0^2}5+\frac{0^2}4-\frac{3^2}{14}=1.157 \\
Q_2(z_{312}) =\frac{0^2}5+\frac{2^2}5+\frac{0^2}4-\frac{2^2}{14}=0.514 \\
Q_2(z_{313}) =\frac{0^2}5+\frac{0^2}5+\frac{2^2}4-\frac{2^2}{14}=0.714 \\
Q_{2}({z}_{321}) =\frac{2^2}5+\frac{0^2}5+\frac{0^2}4-\frac{2^2}{14}=0.514 \\
Q_2(z_{322}) =\frac{0^2}5+\frac{3^2}5+\frac{0^2}4-\frac{3^2}{14}=1.157 \\
Q_2({z}_{323})=\frac{0^2}5+\frac{0^2}5+\frac{2^2}4-\frac{2^2}{14}=0.714
\end{aligned}$$

因此，$k_3=1.157+0.514+0.714+0.514+1.157+0.714=4.770$. 使用上述 $k_1,k_2,k_3$ 的值，SSROWS 的期望值为

$$E(Q_2)=0.200\sigma_a^2+9.285\sigma_b^2+4.770\sigma_c^2+2\sigma_\varepsilon^2$$

类似地，应用 Hartley 的综合法，$Q_1$ 和 $Q_3$ 的期望值被确定为

$$E(Q_1)=7.000{\sigma}_a^2+0.143{\sigma}_b^2+2.429{\sigma}_c^2+(1){\sigma}_\varepsilon^2$$

以及

$$E(Q_3)=4.371{\sigma}_c^2+2{\sigma}_\varepsilon^2$$

一般来说，SSERROR 的期望值等于与 SSERROR 相关的自由度乘以 ${\sigma}_\varepsilon^2$，在这种情况下，它提供

$$E(Q_4)=E(\text{SSERROR})=8\sigma_\varepsilon^2$$

对于不均衡处理结构，计算平方和的方法有很多种，包括 SAS^®^-GLM 的 I-IV 型平方和，SAS-MIXED 的 I-III 型平方和，以及 Henderson 提出的四种方法。下面展示了如何使用 SAS 计算 SAS-GLM 的 I-III 型期望值，并演示了如何使用 SAS 来提供平方和期望的计算。表 \@ref(tab:table18-2) 包含 SAS 代码，用于生成列 $\boldsymbol Z_1$（由 $\boldsymbol a_1$ 和 $\boldsymbol a_2$ 表示）、$\boldsymbol Z_1$（由 $\boldsymbol b_1,\boldsymbol b_2$ 和 $\boldsymbol b_3$ 表示）和 $\boldsymbol Z_3$（由 $\boldsymbol c_{11},\boldsymbol c_{12},\boldsymbol c_{13},\boldsymbol c_{21},\boldsymbol c_{22}$ 和 $\boldsymbol c_{23}$ 表示），并使用 SAS-GLM 计算每列的三种类型的平方和。每列的平方和如表 \@ref(tab:table18-3) 所示。由行导致的 I 型平方和期望值的计算为

$$\begin{aligned}
E[\mathrm{SSROW(I)}]=&(1)\sigma_\varepsilon^2+(0.6429+0.2857+0.2857+0.2857+0.6429+0.2857)\sigma_c^2\\
&+(0.0714+0.0714+0.000)\sigma_b^2+(3.500+3.500)\sigma_a^2 \\
=&\sigma_\varepsilon^2+2.4284\sigma_c^2+0.1428\sigma_b^2+7\sigma_a^2 
\end{aligned}$$

<table>
<caption>(\#tab:table18-2)图 \@ref(fig:figure18-1) 中数据结构和表 \@ref(tab:table19-2) 中的数据的 $[\boldsymbol Z_1,\boldsymbol Z_2,\boldsymbol Z_3]$ 和 $\boldsymbol y$ 的每列的 I-III 型平方和的 SAS 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2018.2.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table18-3)图 \@ref(fig:figure18-1) 中数据结构和表 \@ref(tab:table19-2) 中的数据的 $[\boldsymbol Z_1,\boldsymbol Z_2,\boldsymbol Z_3]$ 和 $\boldsymbol y$ 的每列的 I-III 型平方和</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2018.3.png">
</td>
  </tr>
</tbody>
</table>

各效应的平方和总和显示在表 \@ref(tab:table18-4) 中。因此，表 \@ref(tab:table18-5) 包含使用表 \@ref(tab:table18-4) 中显示的摘要计算的剩余平方和的期望。可以计算每个预期平方和的系数，其中 $\sigma^2_\varepsilon$ 的系数是与相应平方和相关的自由度。

<table>
<caption>(\#tab:table18-4)使用表 \@ref(tab:table18-3) 中的数据，针对每个效应，$[\boldsymbol Z_1,\boldsymbol Z_2,\boldsymbol Z_3]$ 的每列的 I-III 型平方和之和</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2018.4.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table18-5)使用 Hartley 综合法计算的 SAS I-III 型平方和的期望值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2018.5.png">
</td>
  </tr>
</tbody>
</table>

在接下来讨论的统计文献中，有几种计算平方和的替代方法（附有名称）。讨论包含关于何时可以将这些技术应用于混合效应模型以及随机效应模型的指示。Henderson (1953) 介绍了四种计算平方和的方法，称为 Henderson's methods I, II, III 和 IV （另见 Searle, 1987; Henderson, 1984). 以下讨论考虑了 Henderson's methods I 和 II. 方差分析法或 Henderson's methods I 仅适用于随机效应模型，是一种由计算平方和组成的技术，类似于为均衡数据集计算的平方和，不同之处在于它们被更改以考虑每个处理组合或不均衡数据集的不等数量的观测值。上面使用 Henderson's methods I 平方和时未给出理由，下面给出理由。使用双向分类随机效应模型来演示该方法。

可用于描述均衡双向数据集的模型是

$$\begin{aligned}y_{ijk}&=\mu+a_i+b_j+c_{ij}+\varepsilon_{ijk}&&i=1,2,\ldots,s,&j=1,2,\ldots,t,&k=1,2,\ldots,n\end{aligned}$$

其中

$$a_i\thicksim i.i.d.\mathrm{~N}(0,\sigma_a^2),\quad b_j\thicksim i.i.d.\mathrm{~N}(0,\sigma_b^2),\quad c_{ij}\thicksim i.i.d.\mathrm{~N}(0,\sigma_c^2),\quad\mathrm{and~}\quad\varepsilon_{ijk}\thicksim i.i.d.\mathrm{~N}(0,\sigma_\varepsilon^2)$$

方差分析中使用的平方和为

$$\begin{aligned}
\text{SSA}& =nt\sum_{i=1}^{s}(\bar{y}_{i\cdot\cdot}-\bar{y}_{\cdot\cdot\cdot})^2=\sum_{i=1}^{s}\frac{y_{i\cdot\cdot}^2}{nt}-\frac{y_{\cdot\cdot\cdot}^2}{nst}  \\
\text{SSB}& =ns\sum_{j=1}^t(\bar{y}_{\cdot j\cdot}-\bar{y}_{\cdot\cdot\cdot})^2=\sum_{j=1}^t\frac{y_{\cdot j\cdot}^2}{ns}-\frac{y_{\cdot\cdot\cdot}^2}{nst},  \\
\text{SSAB}& =n\sum_{i=1}^s\sum_{j=1}^t(\bar{y}_{ij\cdot}-\bar{y}_{i\cdot\cdot}-\bar{y}_{\cdot j\cdot}+\bar{y}_{\cdot\cdot\cdot})^2  \\
&=\sum_{i=1}^s\sum_{j=1}^t\frac{y_{ij\cdot}^2}n-\sum_{i=1}^s\frac{y_{i\cdot\cdot}^2}{nt}-\sum_{j=1}^t\frac{y_{\cdot j\cdot}^2}{ns}+\frac{y_{\cdot\cdot\cdot}^2}{nst}
\end{aligned}$$

要将这些平方和从等样本量转换为不等样本量，请将乘积 $nst$ 替换为 $n_{\cdot\cdot}$，将 $ns$ 替换为 $n_{\cdot j}$，将 $nt$ 替换为 $n_{i\cdot }$. 因此，不等样本量的平方和变为


$$
\begin{aligned}
\text{SSA} &=\sum_{i=1}^s\frac{y_{i\cdot\cdot}^2}{n_{i\cdot}}-\frac{y^2_{\cdot\cdot\cdot}}{n_{\cdot\cdot\cdot}}, \\
\text{SSB} &=\sum_{j=1}^t\frac{y_{\cdot j\cdot}^2}{n_{\cdot j}}-\frac{y_{\cdot\cdot\cdot}^2}{n_{\cdot\cdot\cdot}}, \\
\text{S} SAB&=\sum_{i=1}^s\sum_{j=1}^t\frac{y_{ij\cdot}^2}{n_{ij}}-\sum_{i=1}^s\frac{y_{i\cdot\cdot}^2}{n_{i\cdot}}-\sum_{j=1}^t\frac{y_{\cdot j \cdot}^2}{n_{\cdot j}}+\frac{y_{\cdot\cdot\cdot}^2}{n_{\cdot\cdot}}
\end{aligned}
$$

这些平方和是在讨论图 \@ref(fig:figure18-1) 中数据结构之初所使用的。

拟合常数法，或 Henderson's methods III，涉及将数据拟合到各种线性模型中，然后计算相应的平方和。该方法使用所谓的因拟合完全模型而导致的平方和减少以及因拟合各种子模型而导致的平方和减少（参见第 \@ref(chap10) 章）。为了设定符号，考虑模型

$$\boldsymbol y=\boldsymbol X_1\boldsymbol b_1+\boldsymbol X_2\boldsymbol b_2+\boldsymbol X_3\boldsymbol b_3+\boldsymbol\varepsilon $$

由于拟合完全模型而导致的总平方和的减少为

$$R(\boldsymbol b_1,\boldsymbol b_2,\boldsymbol b_3)=\boldsymbol y^{\prime}\boldsymbol y-\mathrm{SSERROR}(b_1,b_2,b_3)$$

其中 $\mathrm{SSERROR}(\boldsymbol b_1,\boldsymbol b_2,\boldsymbol b_3)$ 是拟合完全模型后的残差平方和。由于拟合 $\boldsymbol b_1$ 和 $\boldsymbol b_2$ 导致的减少为 $R(\boldsymbol b_1,\boldsymbol b_2)=\boldsymbol y^{\prime}\boldsymbol y-{\text{SSERROR}(\boldsymbol b_1,\boldsymbol b_2)}$ 其中 ${\text{SSERROR}(\boldsymbol b_1,\boldsymbol b_2)}$ 是模型 $\boldsymbol y=\boldsymbol X_1 \boldsymbol b_1+\boldsymbol X_2 \boldsymbol b_2 + \boldsymbol \varepsilon$ 的残差平方和。

拟合 $\boldsymbol b_1$ 和 $\boldsymbol b_2$ 后由于 $\boldsymbol b_3$ 造成的减少由 $R(\boldsymbol b_3 |\boldsymbol b_1 , \boldsymbol b_2)$ 表示，并由下式给出： $R(\boldsymbol b_3 |\boldsymbol b_l , \boldsymbol b_2) = R(\boldsymbol b_1 , \boldsymbol b_2, \boldsymbol b_3) - R(\boldsymbol b_1 , \boldsymbol b_2) = \text{SSERROR}(\boldsymbol b_1, \boldsymbol b_2) - \text{SSERROR}(\boldsymbol b_1, \boldsymbol b_2, \boldsymbol b_3)$。同样，拟合 $\boldsymbol b_2$ 后由于 $\boldsymbol b_1$ 导致的减少为 $R(\boldsymbol b_1 |\boldsymbol b_2) = R(\boldsymbol b_1, \boldsymbol b_2) - R(\boldsymbol b_2) = \text{SSERROR}(\boldsymbol b_2)-\text{SSERROR}(\boldsymbol b_1, \boldsymbol b_2)$。最后，拟合 $\boldsymbol b_3$ 后由于 $\boldsymbol b_1$ 和 $\boldsymbol b_2$ 导致的减少为 $R(\boldsymbol b_1 , \boldsymbol b_2 |\boldsymbol b_3) = R(\boldsymbol b_1 , \boldsymbol b_2, \boldsymbol b_3) - R(\boldsymbol b_3) = \text{SSERROR}(\boldsymbol b_3) - \text{SSERROR}(\boldsymbol b_1, \boldsymbol b_2, \boldsymbol b_3)$.

[^b1b2b3interaction]: 原文：unless $\boldsymbol b_3$ denotes an interaction between $\boldsymbol b_1$ and $\boldsymbol b_2$ or an interaction with $\boldsymbol b_1$ or $\boldsymbol b_2$.

使用此技术的优点之一是，$E[R(\boldsymbol b_1 , \boldsymbol b_2 |\boldsymbol b_3)]$ 不依赖于 $\boldsymbol b_3$（如果它是固定效应）或 $\sigma^2_\varepsilon$（如果它是随机效应），除非 $\boldsymbol b_3$ 表示 $\boldsymbol b_1$ 和 $\boldsymbol b_2$ 之间的交互作用或者与 $\boldsymbol b_1$ 或 $\boldsymbol b_2$ 的交互作用[^b1b2b3interaction]。

对于模型

$$\boldsymbol y=\boldsymbol j  \mu+\boldsymbol X_1\boldsymbol b+\boldsymbol X_2\boldsymbol t+\boldsymbol X_3\boldsymbol g+\boldsymbol \varepsilon $$

其中 $\boldsymbol{b}\sim N(0,\sigma_b^2),\boldsymbol{t}\sim N(0,\sigma_t^2),\boldsymbol{g}\sim N(0,\sigma_g^2)$ 以及 $\boldsymbol{\varepsilon}\sim N(0,\sigma_\varepsilon^2)$，一组可能的平方和（来自 SAS-GLM 和 SAS-MIXED 的 I 型平方和）为

$$\begin{aligned}
&R(\boldsymbol b| \mu)=R( \mu,\boldsymbol b)-R( \mu) \\
&R(\boldsymbol t| \mu,\boldsymbol b)=R( \mu,\boldsymbol b,\boldsymbol t)-R( \mu,\boldsymbol b) \\
&R(\boldsymbol g| \mu,\boldsymbol b,\boldsymbol t)=R( \mu,\boldsymbol b,\boldsymbol t,\boldsymbol g)-R( \mu,\boldsymbol b,\boldsymbol t) \\
&\mathrm{SSERROR}( \mu,\boldsymbol b,\boldsymbol t,\boldsymbol g)]=\boldsymbol y^{\prime}\boldsymbol y-R( \mu,\boldsymbol b,\boldsymbol t,\boldsymbol g)
\end{aligned}$$

这些平方和的期望具有以下形式

$$\begin{aligned}
&E[\text{SSERROR}( \mu,\boldsymbol b,\boldsymbol t,\boldsymbol g)]=(n-p)\sigma_\varepsilon^2 \\
&E[R(\boldsymbol g| \mu,\boldsymbol b,\boldsymbol t)]=k_1\sigma_\varepsilon^2+k_2\sigma_g^2 \\
&E[R(\boldsymbol t| \mu,\boldsymbol b)]=k_3\sigma_\varepsilon^2+k_4\sigma_g^2+k_5\sigma_t^2, \\
&E[R(\boldsymbol b| \mu)]=k_6\sigma_\varepsilon^2+k_7\sigma_g^2+k_8\sigma_t^2+k_9\sigma_b^2
\end{aligned}$$

也可以使用由 $R(\boldsymbol t| \mu),R(\boldsymbol b| \mu, \boldsymbol t),R(\boldsymbol g| \mu, \boldsymbol b, \boldsymbol t)$ 和 $\text{SSERROR}( \mu,\boldsymbol b,\boldsymbol t, \boldsymbol g)$ 给出的 I 型平方和，拟合常数的方法可以用于随机效应和混合效应模型，只要固定效应先于随机效应拟合即可。

综合法可用于评估涉及随机效应和/或多个误差项的任何模型的平方和的期望。一组平方和及其期望（或均方）可用于估计方差分量，以开发假设检验，并构造关于各个方差分量和/或方差分量函数的置信区间。第 \@ref(chap19) 章讨论了估计方法，第 \@ref(chap20) 章介绍了推断技术。许多统计软件程序自动使用 Hartley 综合法来计算期望均方 (expected mean squares).

## 结束语 {#sec18-4}

本章定义了随机效应和固定效应的概念以及随机效应模型的概念。随机效应模型以矩阵形式表示，以描述计算期望均方的方法。采用完全随机设计结构中不均衡单向处理结构来演示计算期望均方的代数方法和 Hartley's 综合法。描述了计算平方和的不同方法，并使用完全随机设计结构中不均衡双向处理结构来演示各自期望均方的计算。提供 SAS 代码来演示 Hartley's 综合法。

## 练习 {#sec18-5}
