
# 案例研究：混合模型 {#chap23}

> "Natural selection is a mechanism for generating an exceedingly high degree of improbability." - R. A. Fisher

第 \@ref(chap22) 章讨论了均衡和不均衡混合模型的数据分析方法。本章针对每种情况都通过实例进行了详细分析。不均衡情况的数据是通过从均衡情况的数据中随机删除一些观测值而获得的。该研究涉及一家公司，该公司希望更换其一家工厂中用于制造某种部件的机器。有三种不同品牌的机器可供选择，因此管理层设计了一项实验，以评估当由公司自己的员工操作时，这三台机器的生产率。从受过操作此类机器培训的员工中随机选择了六名员工。每个被选中的员工都需要在三个不同的班次中操作每台机器。记录的数据是整体生产率得分，考虑了生产部件的数量和质量。数据在表 \@ref(tab:table23-1) 中给出。

<table>
<caption>(\#tab:table23-1)机器-人示例的生产率得分</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.1.png">
</td>
  </tr>
</tbody>
</table>


## 双向混合模型 {#sec23-1}

本实验的处理结构是双向的，机器 (machine) 是固定效应，人 (person) 是随机效应。设计结构为完全随机设计。用于描述在 n 个不同轮班期间每个操作 t 台机器的 b 个人的数据的双向混合模型为

$$y_{ijk}=\mu+\tau_i+p_j+(\tau p)_{ij}+\varepsilon_{ijk},\quad i=1,2,\ldots,t,\quad j=1,2,\ldots,b,\quad k=1,2,\ldots,n$$

其中，$\mu$ 表示由工人总体操作的一组品牌机器的平均生产率得分，$\tau_i$ 表示第 i 台机器对生产率得分的效应，$p_j$ 表示第 j 个人对生产率得分的随机效应，$(\tau p)_{ij}$ 表示特定于第 j 个人操作第 i 台机器时对生产率得分的随机交互效应，$\varepsilon_{ijk}$ 表示与第 j 个人第 k 次操作第 i 台机器相关的随机误差项，以及 $p_j,(\tau p)_{ij},\varepsilon_{ijk}$ 为独立随机变量。

关于该模型中随机变量的额外假设为

$$\begin{aligned}p_j&\sim i.i.d.N(0,\sigma_{person}^2)\\(\tau p)_{ij}&\sim i.i.d.N(0,\sigma_{m\times p}^2)\\\varepsilon_{ijk}&\sim i.i.d.N(0,\sigma_{\varepsilon}^2)\end{aligned}$$

假设 ${\varepsilon}_{ijk}\thicksim i.i.d.N(0,{\sigma}_{\varepsilon}^2)$ 意味着，当对 machine × person 的组合进行测量时，不同班次之间的时间间隔足够长，使得误差项不相关。

第一步是分析模型的随机效应部分。SAS^®^-Mixed 代码用于获得三个方差分量的矩法、REML、最大似然和 MINQUE0 估计。III 型平方和、均方及其相应的期望均方如表 \@ref(tab:table23-2) 所示。表 \@ref(tab:table23-3) 给出了方差分量的最大似然、REML、矩法和 MINQUE0 估计量。表 \@ref(tab:table23-4) 中的 SAS-Mixed 代码用于拟合双向混合模型，并获得方差分量的 REML 估计。通过指定 `Method=ML`, `Method=MIVQUE0` 和 `Method=type3` 获得其他三种方法的估计。

<table>
<caption>(\#tab:table23-2)使用 III 型平方和的均衡数据方差表分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.2.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table23-3)使用 REML, ML, MIVQUE0 和矩法（III 型）来估计方差分量</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.3.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table23-4)用于获取均衡数据集方差分量 REML 估计的 SAS Mixed 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.4.png">
</td>
  </tr>
</tbody>
</table>

第一个要检验的假设是 $H_0\colon{\sigma}_{m\times p}^2=0\text{ vs }H_a\colon{\sigma}_{m\times p}^2>0$. 通过使用表 \@ref(tab:table23-2) 中的期望均方构建的用于检验该假设的统计量为

$$F_{m\times p}=\frac{MSPerson\times Machine}{MSResidual}=46.13$$

与该 $F$ 检验相关的显著性水平小于 0.0001，说明有强有力的证据表明 ${\sigma}_{m\times p}^2$ 非零，是数据变异性的重要因素。如果 machine × person 交互方差分量等于零，则表明三台机器的所有员工（而不仅仅是样本中的员工）之间的生产率差异相似。这也表明，三台机器之间的生产率差异对所有员工（不仅仅是样本中的员工）来说是相似的。换句话说，推断是针对所有可能的员工的总体，而不仅仅是从公司所有员工中随机选择的员工。

对显著的 machine × person 交互方差分量的解释是，机器之间的生产率差异取决于使用机器的人。对大的 machine × person 交互方差分量的一种解释是，有些人更适合某些机器，而不是其他机器。在数据中，可以看到 pearson 6 在 machine 1 和 2 上的表现不如样本中的其他人，但 pearson 6 在 machine 3 上的表现与样本中的其它人差不多。要检验的第二个假设是 $H_0\colon{\sigma}_{p}^2=0\text{ vs }H_a\colon{\sigma}_{p}^2>0$. 同样，检验这一假设的统计量是通过使用表 \@ref(tab:table23-2) 中的期望均方构建的，并由下式给出

$$F_p=\frac{MSPerson}{MSPerson\times Machine}$$

与 $H_0\colon{\sigma}_{p}^2=0\text{ vs }H_a\colon{\sigma}_{p}^2>0$ 相关的显著性水平为 0.0089，表明工厂员工的生产率得分存在相当大的差异。培训项目可以减少员工之间的差异。

表 \@ref(tab:table23-3) 显示了四种方法的方差分量估计。对于 REML、MIVQUE0 和 III 型，方差分量的估计是相同的，但对于 person 和 machine × person 方差分量，ML 的估计要小一些。

Satterthwaite 近似可用于确定在构建关于方差分量置信区间时使用的自由度。第一步是将方差分量的估计表示为表 \@ref(tab:table23-2) 中的均方的线性组合。${{\sigma}}_p^2$ 的矩法估计为 $\hat{{\sigma}}_p^2=\frac19{ MSPerson }-\frac19{ MSPerson }\times Machine=22.8584$. 通过 Satterthwaite 近似获得的近似卡方分布的自由度为

$$v=\frac{(\hat{\sigma}_p^2)^2}{\frac{[\frac19MSPerson]^2}5+\frac{[\frac19MSPerson\times Machine]^2}{10}}=3.38035$$

95% 置信区间的分位数为

$$\chi_{.025,3.38}^2=10.0467\quad\mathrm{~and~}\quad\chi_{.975,3.38}^2=0.30725$$

因此，${{\sigma}}_p^2$ 的 95% 置信区间为

$$\frac{3.38035(22.8584)}{10.0467}\leq\sigma_p^2\leq\frac{3.38035(22.8584)}{0.30725}$$

或

$$7.69102\leq\sigma_p^2\leq251.486$$

${{\sigma}}_p^2$ 的 95% 置信区间为 $2.773\leq\sigma_p^2\leq15.858$.

表 \@ref(tab:table23-5) 给出了 III 型分析的协方差参数估计。方差分量的估计和标准误的估计在前两列中。Z-value 列给出了每个估计与其相应的标准误估计的比值，Pr Z 列是与 Z-value 相对应的双侧显著性水平。Lower 和 Upper 列是使用 Wald 区间计算的置信区间；即，$\hat{{\sigma}}^2\pm Z_{0.025}[\widehat{s.e.}(\hat{{\sigma}}^2)]$. 当与待估方差分量相关的自由度数量较大时，Wald 置信区间才适用。关于 ${\sigma}_{{\varepsilon}}^{{2}}$ （残差）的置信区间基于 36 个自由度的卡方分布。表 \@ref(tab:table23-6) 包含来自 REML 分析的协方差参数估计。

<table>
<caption>(\#tab:table23-5)根据 III 型分析的协方差参数估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.5.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table23-6)来自 REML 分析的协方差参数估计，包括用于计算置信区间的计算自由度</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.6.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table23-5) 和 \@ref(tab:table23-6) 中方差分量估计及其标准误估计是相同的，但 Pr Z 列不同。表 \@ref(tab:table23-6) 中的显著性水平是单侧的（而且应该是单侧的）。计算出的自由度 $df=2(\text{Z-value})^2$ 与使用矩法估计计算出的自由度相同。残差方差分量的置信区间基于 36 个自由度的卡方分布。因此，对于均衡设计，REML、MIVQUE0 和矩法估计得到的方差分量的信息在与单个估计相关的近似自由度方面是相同的，但是 Wald 法被用于构建 person 和 machine × person 的置信区间，这两个方差分量是从矩法估计得来的，而对于 REML 和 MIVQUE0，则使用了更合适的 Satterthwaite 近似法。

估计和检验混合模型中固定效应假设所需的方法取决于数据是均衡的还是不均衡的。第 \@ref(sec22-3) 节讨论了不均衡设计的方法。

一般混合模型可以表示为

$$\boldsymbol y=\boldsymbol X\boldsymbol{\beta}+\boldsymbol Z_1\boldsymbol u_1+\boldsymbol Z_2\boldsymbol u_2+\cdots+\boldsymbol Z_k\boldsymbol u_k+\boldsymbol{\varepsilon}$$

其中 $\boldsymbol{X}=[\boldsymbol{j},\boldsymbol{X}_1],\boldsymbol{\beta}=({\mu},\tau_1,\tau_2,\ldots,\tau_t)'$

如果协方差阵

$$\boldsymbol{\Sigma}={\sigma}_1^2\boldsymbol{Z}_1\boldsymbol{Z}_1'+{\sigma}_2^2\boldsymbol{Z}_2\boldsymbol{Z}_2'+\cdots+{\sigma}_k^2\boldsymbol{Z}_k\boldsymbol{Z}_k'+{\sigma}_\varepsilon^2\boldsymbol{I}$$

的元素是已知的；即，方差分量是已知的，可估函数 $\boldsymbol a'\boldsymbol \beta$ 的 BLUE (best linear unbiased estimate) 为

$$\boldsymbol a'\hat{\boldsymbol\beta}_{BLUE} = \boldsymbol a'\left(\boldsymbol{X'\Sigma}^{-1}\boldsymbol{X}\right)^-\boldsymbol{X'\Sigma}^{-1}\boldsymbol{y}$$

对于大多数均衡混合模型，$\boldsymbol a'\boldsymbol \beta$ 的估计简化为

$$\boldsymbol a'\hat{\boldsymbol\beta}_{BLUE} = \boldsymbol a'\left(\boldsymbol{X'}\boldsymbol{X}\right)^-\boldsymbol{X'}\boldsymbol{y}$$

下一个例子展示了均衡双向混合模型的这种简化。模型可以重新参数化为

$$y_{ijk}=\mu_i+p_j+(\tau p)_{ij}+\varepsilon_{ijk}$$

其中 $\mu_i=\mu+\tau_i$ 且 $\boldsymbol X_1$ 为模型固定效应部分 $\boldsymbol{\mu}=(\mu_1,\mu_2,\ldots,\mu_t)^{\prime}$ 相应的 nbt × t 设计矩阵。

注意到

$$\boldsymbol{X}_1=\begin{bmatrix}\boldsymbol{j}_{nb}&0&0&\cdots&0\\0&\boldsymbol{j}_{nb}&0&\cdots&0\\0&0&\boldsymbol{j}_{nb}&\cdots&0\\\vdots&\vdots&\vdots&\ddots&\vdots\\0&0&0&\cdots&\boldsymbol{j}_{nb}\end{bmatrix}=\boldsymbol{j}_n\otimes\boldsymbol{j}_b\otimes\boldsymbol{I}_t$$

均衡模型的协方差阵为

$$\boldsymbol{\Sigma}=\sigma_p^2(\boldsymbol{J}_n\otimes\boldsymbol{I}_b\otimes\boldsymbol{J}_t)+\sigma_{m\times p}^2(\boldsymbol{J}_n\otimes\boldsymbol{I}_b\otimes\boldsymbol{I}_t)+\sigma_\varepsilon^2(\boldsymbol{I}_n\otimes\boldsymbol{I}_b\otimes\boldsymbol{I}_t)$$

可以表示为

$$\boldsymbol{\Sigma}=\lambda_1\left(\frac1n\boldsymbol{J}_n\otimes\boldsymbol{I}_b\otimes\frac1t\boldsymbol{J}_t\right)+\lambda_2\left[\frac1n\boldsymbol{J}_n\otimes\boldsymbol{I}_b\otimes\left(\boldsymbol{I}_t-\frac1t\boldsymbol{J}_t\right)\right]+\lambda_3\left[\left(\boldsymbol{I}_n-\frac1n\boldsymbol{J}_n\right)\otimes\boldsymbol{I}_b\otimes\boldsymbol{I}_t\right]$$

其中 $\lambda_1=nt\sigma_p^2+n\sigma_{m\times p}^2+\sigma_\varepsilon^2,\lambda_2=n\sigma_{m\times p}^2+\sigma_\varepsilon^2$ 以及 $\lambda_3=\sigma_\varepsilon^2$. 可以证明

$$\boldsymbol{\Sigma}^{-1}=\frac1{\lambda_1}\bigg(\frac1n\boldsymbol{J}_n\otimes\boldsymbol{I}_b\otimes\frac1t\boldsymbol{J}_t\bigg)+\frac1{\lambda_2}\bigg[\frac1n\boldsymbol{J}_n\otimes\boldsymbol{I}_b\otimes\bigg(\boldsymbol{I}_t-\frac1t\boldsymbol{J}_t\bigg)\bigg]+\frac1{\lambda_3}\bigg[\bigg(\boldsymbol{I}_n-\frac1n\boldsymbol{J}_n\bigg)\boldsymbol{\otimes}\boldsymbol{I}_b\otimes\boldsymbol{I}_t\bigg]$$

然后 $\boldsymbol{X'\Sigma}^{-1}\boldsymbol{X}_1$ 和 $\boldsymbol{X'\Sigma}^{-1}\boldsymbol{y}$
简化为

$$\begin{aligned}\boldsymbol{X}_1^{\prime}\boldsymbol{\Sigma}^{-1}\boldsymbol{X}_1&=\frac{n{b}}{{\lambda}_1}\bigg(\frac1t\boldsymbol{J}_t\bigg)+\frac{n{b}}{{\lambda}_2}{\bigg(\boldsymbol{I}_t-\frac1t\boldsymbol{J}_t\bigg)}\\\boldsymbol{X}_1^{\prime}\boldsymbol{\Sigma}^{-1}&=\frac1{{\lambda}_1}\boldsymbol{j}_n^{\prime}\otimes\boldsymbol{j}_b^{\prime}\otimes\frac1t\boldsymbol{J}_t+\frac1{{\lambda}_2}\boldsymbol{j}_n^{\prime}\otimes\boldsymbol{j}_b^{\prime}\otimes\bigg(\boldsymbol{I}_t-\frac1t\boldsymbol{J}_t\bigg)\end{aligned}$$

接着可以证明

$$(\boldsymbol X_1^{\prime}\boldsymbol \Sigma^{-1}\boldsymbol X_1)^{-1}=\frac{\lambda_1}{nb}{\left(\frac1t\boldsymbol J_t\right)}+\frac{\lambda_2}{nb}{\left(\boldsymbol I_t-\frac1t\boldsymbol J_t\right)}$$

[^messymu]: 原书在此处上下文中关于 $\boldsymbol\mu$ 和 $\mu$ 的加粗以及角标较为混乱，此处译者做了修改。

因此 $\boldsymbol \mu$ 的估计为[^messymu]

$$\hat{\boldsymbol\mu}=(\boldsymbol X_1^{\prime}\boldsymbol\Sigma^{-1}\boldsymbol X_1)^{-1}\boldsymbol X_1^{\prime}\boldsymbol\Sigma^{-1}\boldsymbol y=\left(\frac1n\boldsymbol j^{\prime}_n\otimes\frac1b\boldsymbol{j}_b^{\prime}\otimes\boldsymbol{I}_t\right)\boldsymbol{y}=(\boldsymbol X_1^{\prime}\boldsymbol{X}_1)^{-1}\boldsymbol{X}_1^{\prime}\boldsymbol{y}$$

因此，${\mu}$ 的估计值为 $\hat{\mu}_i=\bar{y}_{i\cdot\cdot},i=1,2,\ldots,t$.

$\hat{\boldsymbol\mu}$ 的方差为

$$\mathrm{Var}(\hat{\boldsymbol\mu})=(\boldsymbol X_1^{\prime}\boldsymbol{\Sigma}^{-1}\boldsymbol{X}_1)^{-1}=\frac{{\sigma}_\varepsilon^2+n{\sigma}_{m\times p}^2+nt{\sigma}_p^2}{nb}\biggl(\frac1t\boldsymbol{J}_t\biggr)+\frac{{\sigma}_\varepsilon^2+n{\sigma}_{m\times p}^2}{nb}\biggl(\boldsymbol{I}_t-\frac1t\boldsymbol{J}_t\biggr)$$

因此 $\hat{\mu}_i$ 的方差为

$$\mathrm{Var}(\hat{\mu}_i)=\frac{\sigma_\varepsilon^2+n\sigma_{m\times p}^2+n\sigma_p^2}{nb}$$

对于机器-人的例子，每个机器均值的方差为

$$\mathrm{Var}(\hat{\mu}_i)=\frac{\sigma_\varepsilon^2+3\sigma_{m\times p}^2+3\sigma_p^2}{18}$$

对比 $\boldsymbol{a}^{\prime}\boldsymbol{\mu}$（其中 $\boldsymbol{a'j}_t=\boldsymbol{0}$）的估计为 $\boldsymbol{a}^{\prime}\hat{\boldsymbol{\mu}}$，其方差为

$$\operatorname{Var}(\boldsymbol a^{\prime}\hat{\boldsymbol \mu})=\frac{\sigma_\varepsilon^2+n\sigma_{m\times p}^2}{nb}\boldsymbol{a'a}$$

差值 $\mu_i-\mu_{i^{\prime}}\left(i\neq i^{\prime}\right)$ 的估计为 $\hat{\mu_i}-\hat{\mu}_{i^{\prime}}$，其方差为

$$\mathrm{Var}(\hat{\mu}_i-\hat{\mu}_{i^{\prime}})=2{\left[\frac{\sigma_\varepsilon^2+n\sigma_{m\times p}^2}{nb}\right]}$$

两个机器均值之差的标准误估计为

$$\widehat{s.e.}_{\mu_i-\mu_{i^\prime}}=\sqrt{\frac{2(MSPerson\times Machine)}{18}}=2.177$$

表 \@ref(tab:table23-7) 包含机器均值估计、标准误估计、$t$ 值和显著性水平，用于检验单个均值是否等于零。机器均值之间的成对差值如表 \@ref(tab:table23-8) 所示，其中包括差值估计、标准误估计、$t$ 值和显著水平。用于比较两个均值的 $LSD_{0.05}$ 值计算为

$$LSD_{0.05}=(t_{0.025,10})(\widehat{s.e}._{\hat{\mu}_i-\hat{\mu}_{i^\prime}})=2.228(2.177)=4.85$$

<table>
<caption>(\#tab:table23-7)均衡数据集的机器均值和标准误估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.7.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table23-8)均衡数据集的机器均值之间的成对差值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.8.png">
</td>
  </tr>
</tbody>
</table>

所有的差值都大于 4.85，因此 LSD 显示机器均值彼此之间都有显著差异。

为了检验等均值假设，$H_0\colon\mu_1=\mu_2=\mu_3\mathrm{~vs~ }H_a\colon(\text{not }H_0)$，使用表 \@ref(tab:table23-2) 中的期望均方来构建检验统计量

$$F_{mc}=\frac{MSMachine}{MSPerson\times Machine}=20.58$$

与等均值假设对应的显著性水平为 0.0004。这个例子表明，对于均衡的情况，分析是非常简单的。然而，对不均衡情况的分析并不是那么容易。

## 不均衡双向混合模型 {#sed23-2}

此例中表 \@ref(tab:table23-1) 的数据与均衡示例中的数据相同，只是随机删除了一些观测值。这样做是为了演示在分析不均衡数据集时出现的问题，并比较均衡和不均衡情况下的估计程序。描述此数据的模型与第 \@ref(sec23-1) 节中均衡数据集的模型相同，或者

$$y_{ijk}=\mu+\tau_i+p_j+(\tau p)_{ij}+\varepsilon_{ijk},\quad i=1,2,\ldots,t,\quad j=1,2,\ldots,b,\quad k=1,2,\ldots,n_{ij}$$

不均衡数据的方差分析表如表 \@ref(tab:table23-9) 所示。平方和是通过拟合常数法获得的，或 I 型平方和；即，

$$\begin{aligned}SSMachine&=R(\tau|\mu)\\SSPerson&=R(p|\mu,\tau)\\SSPerson\times Machine&=R((\tau p)|\mu,\tau,p)\end{aligned}$$

表 \@ref(tab:table23-9) 包含了相应的期望均方。

<table>
<caption>(\#tab:table23-9)不均衡数据集基于 I 型平方和的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.9.png">
</td>
  </tr>
</tbody>
</table>

首先，考虑估计 $\sigma_p^2$. 在不均衡的情况下，MSPerson 和 MSPerson × Machine 的期望均方中的 $\sigma_{m\times p}^2$ 系数不相同。因此为了估计 $\sigma_p^2$，我们需要找到一个关于所有三个期望均方的函数，该函数等于 $\sigma_p^2$. 可以证明

$$\begin{aligned}
\tilde{\mathcal{\sigma}}_p^2 &=\frac1{7.219}{\left[(\tilde{\sigma}_\varepsilon^2+2.5866\tilde{\sigma}_{m\times p}^2+7.2190\tilde{\sigma}_p^2)-\frac{2.5866}{2.3162}(\tilde{\sigma}_\varepsilon^2+2.3162\tilde{\sigma}_{m\times p}^2)+\left(\frac{2.5866}{2.3162}-1\right)\tilde{\sigma}_\varepsilon^2\right]} \\
&=\frac1{7.219}{\left[(\tilde{\sigma}_\varepsilon^2+2.5866\tilde{\sigma}_{m\times p}^2+7.2190\tilde{\sigma}_p^2)-1.1167(\tilde{\sigma}_\varepsilon^2+2.3162\tilde{\sigma}_{m\times p}^2)+(0.1167)\tilde{\sigma}_\varepsilon^2\right]} \\
&=\frac1{7.219}[{E}(MSPerson)-1.1167{E}(MSPerson\times\text{М}SMachine)+(0.1167){E}(MSResidual)] 
\end{aligned}$$

I 型平方和的矩法方程组为

$$\begin{aligned}
\text{MSPerson}& =\tilde{\sigma}_\varepsilon^2+2.5866\tilde{\sigma}_{m\times p}^2+7.2190\tilde{\sigma}_p^2  \\
\text{MSPerson × Machine}& =\tilde{\sigma}_\varepsilon^2+2.3162\tilde{\sigma}_{m\times p}^2  \\
MSResidual&=\tilde{\sigma}_\varepsilon^2
\end{aligned}$$

同时求解上述方程组给出了每个方差分量的 I 型矩法估计为

$$\begin{aligned}
\tilde{{\sigma}}_p^2& =\frac1{7.219}[MSPerson-1.1167 MSPerson\times Machine+0.1167MSResidual]  \\
&=\frac1{7.219}{\left[201.7527-1.1167(40.4315)+0.1167(0.8726)\right]}=21.7073 \\
\tilde{{\sigma}}_{m\times p}^2& =\frac1{2.3162}[MSPerson\times Machine-MSResidual]  \\
&=\frac1{2.3162}{\left[40.4315-0.8726\right]}{=17.0792} \\
 \tilde{{\sigma}}_\varepsilon^2&=MSResidual=0.8726 
\end{aligned}$$

表 \@ref(tab:table23-10) – \@ref(tab:table23-13) 给出了从 REML, ML, MIVQUE0, I 型, II 型和 III 型方法获得的方差分量估计。使用 MIVQUE0 法的残差方差分量的解是 0.（如果在 SAS-Mixed 中使用 `nobound` 选项，则解为负。如果有一个我们应该始终能够估计的方差分量，那就是残差）表 \@ref(tab:table23-13) 中的总结表明，REML 方法似乎处于其他方法估计的中间。

<table>
<caption>(\#tab:table23-10)不均衡数据集残差方差分量的估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.10.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table23-11)不均衡数据集 Person 方差分量的估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.11.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table23-12)不均衡数据集 Person × machine 方差分量的估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.12.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table23-13)不均衡数据集方差分量估计的总结</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.13.png">
</td>
  </tr>
</tbody>
</table>

用于获得这些结果的 SAS-Mixed 代码如表 \@ref(tab:table23-14) 所示，其中 `Method=REML` 可替换为 ML, MIVQUE0, type1, type2 或 type3，以产生其他结果。这些表格包含方差分量估计、标准误估计、重新计算的自由度、置信区间和重新计算的 Satterthwaite 置信区间。回想一下，REML, ML 和 MIVQUE0 的方差分量置信区间是基于使用 $df=2(\text{Z-value})^2$ 的卡方分布的。lower 和 newlower 以及 upper 和 newupper 对于这些方法是相同的。当 Method 为 I 型、II 型或 III 型时，使用 Wald 区间为方差分量（残差除外）提供区间。newlower 和 newupper 是重新计算的置信区间，将是所选的区间。

<table>
<caption>(\#tab:table23-14)用于生成不均衡数据集 REML 估计和分析的 SAS Mixed 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.14.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table23-9) 中的 $F$ 统计量可用于检验关于方差分量的假设。用于检验 $H_0\colon{\sigma}_{p}^2=0\text{ vs }H_a\colon{\sigma}_{p}^2>0$ 的统计量是 $F=4.48$，分子自由度为 5，分母自由度为 9.9549. 显著性水平为 0.0003. 用于检验 $H_0\colon{\sigma}_{m\times p}^2=0\text{ vs }H_a\colon{\sigma}_{m\times p}^2>0$ 的统计量是 $F=46.34$，分子自由度为 10，分母自由度为 26. 显著性水平小于 0.0001. 似然比检验也可用于检验关于方差分量的假设。该过程包括获取完全模型方差分量的 ML 估计，并计算 -2 log(full likelihood function)，如表 \@ref(tab:table23-15) 所示。为了检验 $H_0\colon{\sigma}_{m\times p}^2=0\text{ vs }H_a\colon{\sigma}_{m\times p}^2>0$，我们需要拟合一个不含 person × machine 的模型，并计算 -2 log(reduced likelihood function). 似然比检验统计量计算为 LR test = -2 log(reduced likelihood function) - [-2 log(full likelihood function)]，等于 55.5665，并且渐近分布于自由度为 1 的卡方分布。此检验的显著性水平小于 0.0001. 为了检验 $H_0\colon{\sigma}_{m\times p}^2=0\text{ vs }H_a\colon{\sigma}_{m\times p}^2>0$，我们需要拟合一个不含 person 的模型，并计算 -2 log(reduced likelihood function). 似然比检验统计量计算为 LR test = -2 log(reduced likelihood function) - [-2 log(full likelihood function)]，等于 6.4226，并且渐近分布于自由度为 1 的卡方分布。此检验的显著性水平为 0.0113.

<table>
<caption>(\#tab:table23-15)关于 ${\sigma}_{p}^2$ 和 ${\sigma}_{g}^2$ 假设的似然比检验的三种模型的方差分量的 ML 估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.15.png">
</td>
  </tr>
</tbody>
</table>

固定效应参数的估计计算为

$$\hat{\boldsymbol{\beta}}_W=\left(\boldsymbol X^{\prime}\hat{\boldsymbol{\Sigma}}^{-1}\boldsymbol X\right)^{-}\boldsymbol X^{\prime}\hat{\boldsymbol{\Sigma}}^{-1}\boldsymbol y$$

[^evaluatedcovariance]: 原文：where $\hat{\boldsymbol{\Sigma}}$ is the estimated covariance matrix evaluated at the estimates of the variance components from the specified method of estimating the variance components.

其中 $\hat{\boldsymbol{\Sigma}}$ 是在通过特定方法估计方差分量时，在方差分量估计处评估得到的协方差阵估计[^evaluatedcovariance]，以及 $\hat{\boldsymbol{\beta}}$ 的方差-协方差阵为 $\mathrm{Var}({\hat{\boldsymbol\beta}})=\left({\boldsymbol{X'\Sigma}}^{-1}\boldsymbol{X}\right)^-$. 对于每种估计方差分量的方法，都获得了机器等均值的检验。表 \@ref(tab:table23-16) 总结了固定效应的 III 型检验。MIVQUE0 结果没有用处。REML 的分母自由度为 10.1，接近于均衡数据集的 10。ML 方法分母自由度太大，而 I–III 型方法太少。

<table>
<caption>(\#tab:table23-16)不均衡数据集固定效应 III 型检验的总结</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.16.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table23-17) - \@ref(tab:table23-19) 包含了机器均值的估计。REML 与 I 型、II 型或 III 型结果之间没有太大差异。表 \@ref(tab:table23-20) – \@ref(tab:table23-22) 包含了在每种估计方差分量方法下，每对机器均值之差的估计。由于 REML 的分母自由度大约等于 10（如果数据集是均衡的，则这是预期结果），因此这些方法的结果之间存在差异。REML 方法似乎提供了最全面的结果。

<table>
<caption>(\#tab:table23-17)Machine 1 均值的估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.17.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table23-18)Machine 2 均值的估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.18.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table23-19)Machine 3 均值的估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.19.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table23-20)Machine 1 和 2 均值之差的估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.20.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table23-21)Machine 1 和 3 均值之差的估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.21.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table23-22)Machine 2 和 3 均值之差的估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.22.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table23-23) 显示了基于 REML 估计的等均值假设 $H_0{:}\mu_1=\mu_2=\mu_3\mathrm{~vs~}H_a{:}(\text{not }H_a)$ 的检验，其中基于 2 个分子自由度和 10.1 个分母自由度计算的 $F$ 值为 19.96. 均衡数据集分析的分母自由度为 10（对应于 machine × person 交互的自由度），因此 10.1 是一个近似值。

<table>
<caption>(\#tab:table23-23)不均衡数据集等机器均值的 III 型检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2023.23.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table23-20) – \@ref(tab:table23-22) 中的结果是每对机器均值之间的成对差值。同样，基于方差分量估计的 REML 方法的结果似乎是最好的，因为成对比较的自由度是 10.2, 10.1 和 10，这些值接近 10，就好像均衡分析的情况一样。

## 不均衡双向数据集的 JMP 分析 {#sec23-3}

使用 JMP^®^ 对不均衡数据集的分析包括构建如图 \@ref(fig:figure23-1) 所示的数据表。变量 machine, person 和 rep 被定义为名义的，并且 rating 是连续的。fit model 界面如图 \@ref(fig:figure23-2) 所示，其中 rating 已被选为 “Y” 变量，machine 是一种固定效应，person 以及 person 与 machine 的交互是随机效应。选择 REML 方法作为估计方差分量的方法。方差分量的 REML 估计以及标准误估计和 Wald 置信区间如图 \@ref(fig:figure23-3) 所示。等均值假设 $H_0{:}\mu_1=\mu_2=\mu_3\mathrm{~vs~}H_a{:}(\text{not }H_a)$ 的检验如图 \@ref(fig:figure23-4) 所示，其中基于 2 个分子自由度和 10.11 个分母自由度计算的 $F$ 值为 19.9639.

<div class="figure" style="text-align: center">
<img src="figure/figure 23.1.png" alt="JMP 数据表屏幕" width="775" />
<p class="caption">(\#fig:figure23-1)JMP 数据表屏幕</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 23.2.png" alt="JMP fit model 界面，包含响应变量和模型效应" width="980" />
<p class="caption">(\#fig:figure23-2)JMP fit model 界面，包含响应变量和模型效应</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 23.3.png" alt="JMP 方差分量的 REML 估计" width="978" />
<p class="caption">(\#fig:figure23-3)JMP 方差分量的 REML 估计</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 23.4.png" alt="JMP 机器效应检验" width="554" />
<p class="caption">(\#fig:figure23-4)JMP 机器效应检验</p>
</div>

最小二乘均值及其标准误估计如图 \@ref(fig:figure23-5) 所示。使用 LSD 多重比较方法对均值进行成对比较，结果包括差值、标准误估计和 95% 置信区间，如图 \@ref(fig:figure23-6) 所示。使用 JMP 的 estimate 过程的结果与使用 REML 方法的 SAS Mixed 的结果相同，同之处在于置信区间是使用 Wald 法构建的，而不是使用 Satterthwaite 近似法。

<div class="figure" style="text-align: center">
<img src="figure/figure 23.5.png" alt="机器最小二乘均值和标准误估计" width="393" />
<p class="caption">(\#fig:figure23-5)机器最小二乘均值和标准误估计</p>
</div>

<div class="figure" style="text-align: center">
<img src="figure/figure 23.6.png" alt="机器均值的 LSD 多重比较" width="574" />
<p class="caption">(\#fig:figure23-6)机器均值的 LSD 多重比较</p>
</div>

## 结束语 {#sec23-4}

[^mainindicator]: 原文：The main indicator is that the error degrees of freedom for comparing the machine means using the REML estimates of the variance components were similar to those for the balanced data set.

**对于均衡数据集，只要方差分量的解为正，REML, MIVQUE0 和矩法这三种估计方差分量的方法就能提供相同的结果**。这些方法对于固定效应的分析也提供相同的结果。**当数据集不均衡衡时，REML 方法似乎为固定效应的分析提供了最佳结果**。主要指标是，使用方差分量的 REML 估计来比较机器均值的误差自由度与均衡数据集的误差自由度相似[^mainindicator]。使用 REML 选项的 JMP 提供的结果与使用 REML 方法的 SAS-Mixed 获得的方差分量和固定效应的估计相同。不幸的是，关于方差分量构建的置信区间使用了 Wald 法，这仅当与每个方差分量的估计相关联的自由度数量很大时才适用。**对于小样本量数据集，使用 Satterthwaite 自由度近似的卡方置信区间优于 Wald 法提供的置信区间**。

## 练习 {#sec23-5}
