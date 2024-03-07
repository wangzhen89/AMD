
# 分析重复测量实验的方法 {#chap26}

[^subject]: 本书将 "subject" 译作个体，原因见介绍页的[术语规范](#术语规范)。
[^timeinterval]: 原文：The smaller unit is the interval of time during which the subject is exposed to a treatment or an interval just between time measurements.

与使用裂区设计的实验一样，使用**重复测量设计** (repeated measures designs) 的实验具有涉及多个尺寸的实验单元的结构。例如，可以随时间测量个体 (subject[^subject])，其中时间是实验处理结构中的因素之一。通过在几个不同的时间测量个体，个体本质上被 “分裂” 成多个部分（时间间隔），并测量每个部分的响应。较大的实验单元是个体或时间间隔的集合。较小的单元是个体接受处理的时间间隔，或只是时间测量之间的间隔[^timeinterval]。

重复测量设计不同于裂区设计，因为一个或多个因素的水平不能随机分配到实验中一个或多个实验单元的尺寸。在这种情况下，**时间的水平不能随机分配到时间间隔，因此像裂区实验一样分析重复测量实验可能是无效的**。由于这种非随机分配方式，相应实验单元误差的协方差阵的结构可能不符合常规裂区实验所对应协方差阵结构。将重复测量实验视为裂区实验进行分析通常称为**时间的裂区分析** (split-plot in time analysis).

第 \@ref(sec26-1) 节描述了重复测量模型，给出了必要的假设，这些假设使**<u>时间的裂区方差分析</u>** (split-plot in time analysis of variance) 有效。第 \@ref(sec26-2) 节给出了三个示例，演示了用于<u>时间的裂区方差分析</u>的计算，包括用于在均值之间进行各种比较的标准误的计算。第 \@ref(chap27) 章提出了当不满足<u>时间的裂区分析</u>的假设时分析重复测量实验的方法。此外，还提出了一些方法，可以检验<u>时间的裂区分析</u>的假设是否满足。

## 模型指定和理想条件 {#sec26-1}

重复测量设计可以应用于多种情况。下一节将分析的示例 [26.1](#sec26-2-1) 研究了三种药物对心率的影响，其中每种药物对八个人进行了给药。然后在给药后的 5, 10, 15 和 20 分钟测量每个人的心率。

在简单重复测量实验的一般模型中，$n_i$ 名个体被随机分配处理 i，每个个体在 p 个时间点进行测量。表 \@ref(tab:table26-1) 说明了简单重复测量实验的布局。

<table>
<caption>(\#tab:table26-1)一个简单重复测量实验的布局</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.1.png">
</td>
  </tr>
</tbody>
</table>

当使用时间的裂区这种表述时，较大尺寸的实验单元是个体，较小尺寸的实验单元是时间间隔。描述每个测量响应的模型类似于完全随机设计中的裂区模型，由

\begin{equation}
y_{ijk}=\mu+\alpha_i+\delta_{ik}+\tau_j+(\alpha\tau)_{ij}+\varepsilon_{ijk}
(#eq:26-1)
\end{equation}

或

$$y_{ijk}=\mu_{ij}+\delta_{ik}+\varepsilon_{ijk}$$

给出。其中 $\mu+\alpha_i+\delta_{ik}$ 表示模型的个体部分，$\tau_j+(\alpha\tau)_{ij}+\varepsilon_{ijk}$ 是模型的个体内（时间间隔）部分。在时间 j 处的处理 i 的均值为

$$\mu_{ij}=\mu+\alpha_i+\tau_j+(\alpha\tau)_{ij}$$

$\delta_{ik}$ 表示个体误差分量，$\varepsilon_{ijk}$ 表示个体间 (within-subject)（时间间隔）误差。<u>时间的裂区分析</u>的理想条件为

1. $\delta_{ik}$ 独立同分布于 $N(0,\sigma_{\delta}^{2})$
2. $\varepsilon_{ijk}$ 独立同分布于 $N(0,\sigma_\varepsilon^2)$
3. 所有 $\delta_{ik}$ 和 $\varepsilon_{ijk}$ 相互独立

请注意，这些假设与在完全随机设计中对裂区实验进行分析时所做的假设相同（见第 \@ref(chap24) 章）。这种假设可能并不总是适用于重复测量设计。然而，在更一般的条件下，裂区分析也是正确的分析。更一般的条件将在第 \@ref(chap27) 章中给出。

在示例 [26.3](#sec26-2-3) 中，研究了家庭成员的态度是否与其居住环境（乡村或城市）有关。研究所需的家庭类型包括三位成员：儿子、父亲和母亲。随后，从城市环境中随机选择了 10 个这样的家庭，从乡村环境中随机选择了 7 个这样的家庭。确定了每个家庭成员对道德问题的态度。此外，研究还关注态度是否会随着时间发生变化；因此，在 0 个月、6 个月和 12 个月的时间点对每个人进行了测量。该设计中共有三种尺寸的实验单元：对应于环境处理效应的实验单元是家庭；对应于家庭成员效应的实验单元是个人；而对应于时间的实验单元是 6 个月的时间间隔。描述态度的模型具有三个误差分量，分别对应实验中的每种实验单元尺寸。模型为

$$\begin{aligned}
Y_{ijkm}& =\mu+\alpha_{i}+f_{im}  &&\text{family experimental unit}\\
&+\beta_j+(\alpha\beta)_{ij}+p_{ijm} &&\text{person experimental unit}\\
&+\tau_k+(\alpha\tau)_{ik}+(\beta\tau)_{jk}+(\alpha\beta\tau)_{ijk}+\varepsilon_{ijkm}&&\text{time interval experimental unit}
\end{aligned}$$

或

\begin{equation}
y_{ijkm}=\mu_{ijk}+f_{im}+p_{ijm}+\varepsilon_{ijkm}
(#eq:26-2)
\end{equation}

在上述模型中，$i=1,2;j=1,2,3;k=1,2,3;m=1,2,\ldots,n_i$，其中 $n_1=10,n_2=7$.

[^familymember]: 原文：family member designation cannot be randomly assigned to family members. 译者注：这句话的意思是家庭成员的身份不可交换。如父母和孩子虽然都是三个 “人” ，但 “人” 的属性不同，例如是先有的爸爸再有的儿子，而不能反过来。

误差项为，$f_{im}$ 表示环境内家庭误差、$p_{ijm}$ 表示家庭内个人误差，$\varepsilon_{ijkm}$ 表示个人内时间误差。请注意，本实验有两组重复测量，三个家庭成员形成一组重复测量值，因为家庭成员的指定不能随机分配给家庭成员[^familymember]，时间形成第二组重复测量。如果两组重复测量都满足传统的裂区型假设，那么上述模型中误差项的理想条件是

1. $f_{im}$ 独立同分布于 $N(0,\sigma_f^2)$.
2. $p_{ijm}$ 独立同分布于 $N(0,\sigma_p^2)$.
3. $\varepsilon_{ijkm}$ 独立同分布于 $N(0,\sigma_\varepsilon^2)$.
4. 所有 $f_{im},p_{ijm}$ 和 $\varepsilon_{ijkm}$ 相互独立

第 \@ref(chap27) 章将探讨家庭成员示例的一种分析方法，该方法不需要满足理想条件。

## 时间的裂区分析 {#sec26-2}

通常的<u>时间的裂区方差分析</u>是指第 \@ref(chap24) 章和第 \@ref(chap25) 章中讨论的裂区和条区设计的分析。这种分析可由大多数计算机程序包提供。如果满足理想条件，则这些分析提供了有效的 $F$ 检验。

本节中使用了三个例子来演示重复测量设计的分析，并展示如何确定感兴趣效应的估计及其标准误估计（进行各种多重比较所必需的），并提供研究各种均值之间对比的方法。这里使用的例子也在第 \@ref(chap24) 章中用于演示当模型误差项不满足理想条件时所需的技术。

### 示例 26.1：药物对心率的影响 {#sec26-2-1}

进行了一项涉及 t 种药物的实验，以研究每种药物对人类心率的影响。给药后，每隔五分钟测量一次每个人的心率，总共 p 次。在研究开始时，n 名女性受试者被随机分配每种药物。使用类似于等式 \@ref(eq:26-1) 的模型来描述数据。模型为

$$y_{ijk}=\mu+\alpha_i+\delta_{ik}+\tau_j+(\alpha\tau)_{ij}+\varepsilon_{ijk},~i=1,2,\ldots,t;~j=1,2,\ldots,p;~k=1,2,\ldots,n$$

该模型有两个误差项：$\delta_{ik}$ 表示个体误差分量，$\varepsilon_{ijk}$ 表示时间误差分量。<u>时间的裂区分析</u>的理想条件为

1. $\delta_{ik}$ 独立同分布于 $N(0,\sigma_\delta^2)$.
2. $\varepsilon_{ijk}$ 独立同分布于 $N(0,\sigma_\varepsilon^2)$.
3. 所有 $\delta_{ik}$ 和 $\varepsilon_{ijk}$ 相互独立

表 \@ref(tab:table26-2) 给出了示例 [26.1](#sec26-2-1) 的<u>时间的裂区方差分析</u>表。

<table>
<caption>(\#tab:table26-2)示例 [26.1](#sec26-2-1) 的方差表分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.2.png">
</td>
  </tr>
</tbody>
</table>

在表 \@ref(tab:table26-2) 中，$Q_{Drug},Q_{Time}$ 和 $Q_{Drug\times Time}$ 分别是测量药物效应、时间效应和药物 × 时间交互效应的非中心参数。当且仅当相应的效应等于零时，这些 $Q$ 值为零。

为检验 $H_{01}\colon Q_{Drug}=0$，拒绝 $H_{01}$ 如果

$$F=\frac{MSDrug}{MSError(Subject)}>F_{\alpha,t-1,t(n-1)}$$

为检验 $H_{02}\colon Q_{Time}=0$，拒绝 $H_{02}$ 如果

$$F=\frac{MSTime}{MSError(Time)}>F_{\alpha,p-1,t(p-1)(n-1)}$$

为检验 $H_{03}\colon Q_{Drug\times Time}=0$，拒绝 $H_{03}$ 如果

$$F=\frac{MSDrug\times Time}{MSError(Time)}>F_{\alpha,(t-1)(p-1),t(p-1)(n-1)}$$

与裂区实验类似，两个方差分量的矩法解由下式给出

$$\begin{aligned}\tilde{{\sigma}}_\varepsilon^2&=MSError(Time)\\\tilde{{\sigma}}_\delta^2&=\frac{MSError(Subject)-\tilde{{\sigma}}_\varepsilon^2}p\end{aligned}$$

那么两个方差分量的估计为

$$\begin{aligned}\hat{\sigma}_\varepsilon^2&=\tilde{\sigma}_\varepsilon^2\\\hat{\sigma}_\delta^2&=\begin{cases}\tilde{\sigma}_\delta^2&\mathrm{if~}\tilde{\sigma}_\delta^2>0\\0&\mathrm{if~}\tilde{\sigma}_\delta^2\leq0&\end{cases}\end{aligned}$$

当在数据中发现显著效应时，通常希望对各种均值进行相互比较。如果没有药物与时间的交互作用，那么我们通常会想要比较药物主效应的均值和时间主效应的均值。如果交互效应显著，那么可能希望在同一时间点上对不同药物之间的均值进行比较，或者对每种药物在不同时间点上的均值进行比较。$\mu_{ij},\bar{\mu}_{i\cdot},\bar{\mu}_j$ 和 $\bar{\mu}_{\cdot\cdot}$ 的估计分别为 $\hat{{\mu}}_{ij}=\bar{{y}}_{ij\cdot},\hat{\bar{{\mu}}}_{i\cdot}=\bar{{y}}_{i\cdot\cdot},\hat{\bar{{\mu}}}_{\cdot j}=\bar{{y}}_{\cdot j\cdot}$ 和 $\hat{{\bar{\mu}}}_{\cdot\cdot}=\bar{{y}}_{\cdot\cdot\cdot}$. 由于重复测量设计涉及多种尺寸的实验单元和多个误差项，因此每次比较的方差可能涉及各误差项方差的不同函数，因此需要确定。表 \@ref(tab:table26-3) 给出了重复测量实验中均值参数各种函数的最佳估计及其各自的标准误估计。

<table>
<caption>(\#tab:table26-3)各种参数函数、最佳估计和的标准误估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.3.png">
</td>
  </tr>
</tbody>
</table>

最后，$\mu_{ij}-\mu_{i'j}-\mu_{ij'}+\mu_{i'j'}$ 的估计为 $\hat{\mu}_{ij}-\hat{\mu}_{i'j}-\hat{\mu}_{ij'}+\hat{\mu}_{i'j'}$，当 $i\ne i',j\ne j'$ 时，标准误估计为 $\sqrt{(4\hat{\sigma}_{\varepsilon}^{2}/n)}$.

可以使用最佳估计与其标准误估计的比值来对上述每个参数函数进行推断。每个比值要么具有 $t$ 分布，要么可以通过 $t$ 分布来近似。对于那些标准误估计仅涉及 $\hat \sigma^2_\varepsilon$ 的比值，该比率具有精确的 $t$ 分布，自由度为 t(n - 1)(p - 1). 对于那些标准误估计为 $\sqrt{(\hat{\sigma}_{\varepsilon}^{2}+p\hat{\sigma}_{\delta}^{2})}$ 倍数的比值，该比值也具有精确的 $t$ 分布，但自由度为 t(n - 1). 最后，对于那些标准误估计为 $\sqrt{(\hat{\sigma}_{\varepsilon}^{2}+\hat{\sigma}_{\delta}^{2})}$ 倍数的比值，该比率具有近似 $t$ 分布，其自由度必须通过 Satterthwaite 法进行近似。

一般来说，令 ${\varphi}=\sum_{i=1}^{t}{\sum}_{j=1}^{p}{c}_{i{j}}{\mu}_{i{j}}$ 为 $\mu_{ij}$ 的任何线性组合，那么 $\hat \mu_{ij}$ 的相同线性组合 ${\hat\varphi}=\sum_{i=1}^{t}{\sum}_{j=1}^{p}{c}_{i{j}}{\hat\mu}_{i{j}}$ 是其最佳估计。对于 $k,c_1,c_2$ 的某种取值，$\hat \phi$ 的方差估计等于 $k[c_1{\hat{\sigma}}_1^2+c_2{\hat{\sigma}}_2^2]$，其中 $v_i\hat{\sigma}_i^2/\sigma_i^2$ 独立分布于 $\chi^2(v_i)$，其中 $i=1,2$. 那么 ${(\widehat{\phi}-\phi)}\widehat{s.e.}(\widehat{\phi})$ 近似为 $t(\hat v)$，其中

$$\hat{v}=\frac{(c_1\hat{\sigma}_1^2+c_2\hat{\sigma}_2^2)^2}{(c_1^2\hat{\sigma}_1^2/v_1)+(c_2^2\hat{\sigma}_2^2/v_2)}$$

对于正在考虑的例子，$\hat{{\sigma}}_1^2$ 是自由度为 $v_1=t(n-1)$ 的 MSError(Subject)，$\hat{{\sigma}}_2^2$ 是自由度为 $v_2=t(n-1)(p-1)$ 的 MSError(Time). 此外，对于那些标准误估计为 $\hat{{\sigma}}_\varepsilon^2+\hat{{\sigma}}_\delta^2$ 倍数的比值，则 $c_1 = 1/p$ 且 $c_2 = (p - 1)/p$.

接下来考虑药物主效应均值的线性对比 $\sum_{i=1}^tc_i\bar{\mu}_{i\cdot}$。这种对比的标准误估计为

$$\sqrt{\frac{(\hat{\sigma}_\varepsilon^2+p\hat{\sigma}_\delta^2)\sum_{i=1}^tc_i^2}{np}}$$

因此

$$\frac{\sum_{i=1}^tc_i\hat{\bar{{\mu}}}_{i\cdot}-\sum_{i=1}^tc_i\bar{{\mu}}_{i\cdot}}{\sqrt{\frac{(\hat{{\sigma}}_\varepsilon^2+p\hat{{\sigma}}_\delta^2)\sum_{i=1}^tc_i^2}{np}}}\sim t[t(n-1)]$$

考虑时间主效应均值的线性对比 $\sum_{i=1}^td_j\bar{\mu}_{\cdot j}$。这种对比的标准误估计是

$$\sqrt{\frac{\hat{{\sigma}}_\varepsilon^2\sum_{j=1}^pd_j^2}{nt}}$$

因此

$$\frac{\sum_{j=1}^pd_j\hat{{\bar{\mu}}}_{\cdot j}-\sum_{j=1}^pd_j{{\bar{\mu}}}_{\cdot j}}{\sqrt{\frac{{\hat{\sigma}}_\varepsilon^2\sum_{j=1}^pd_j^2}{nt}}}\sim t[t(n-1)(p-1)]$$

这种对比可用于检查随时间变化的线性、二次和类似趋势。

如果存在药物 × 时间交互作用，那么我们需要比较每个时间段的药物和/或比较每种药物的时间段。也就是说，可以考虑 j 的每种可能性的 $\mu_{ij}$ 中的对比，和/或 i 的每种可能性的 $\mu_{ij}$ 中的对比。对于前者，可以证明，对于每个 j

$$\frac{\sum_{i=1}^tc_i\hat{{\mu}}_{ij}-\sum_{i=1}^tc_i\mu_{ij}}{\sqrt{\frac{(\hat{\sigma}_\varepsilon^2+\hat{\sigma}_\delta^2)\sum_{i=1}^tc_i^2}n}}$$

近似分布于自由度为 $\hat v$ 的 $t$ 分布。对于后者，可以证明，对于每个 i

$$\frac{\sum_{i=1}^td_j{\hat{\mu}}_{ij}-\sum_{i=1}^td_j{\mu}_{ij}}{\sqrt{\frac{{\hat{\sigma}}_\varepsilon^2\sum_{i=1}^td_j^2}n}}$$

正好是自由度为 t(n - 1)(p - 1) 的 t 分布。这两个结果可以用来检验假设，并为行内或列内的 $\mu_{ij}$ 的对比构造置信区间。

表 \@ref(tab:table26-4) 中的数据用于演示上述分析。在这个实验中，有三种药物，每种药物八个人，以及四个时间段。方差分析表见表 \@ref(tab:table26-5). 存在显著的时间与药物交互作用；因此，我们需要对每种药物的不同时间进行比较，并在每个时间点对不同药物进行比较。

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table26-4)示例 [26.1](#sec26-2-1) 的心率数据</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.4.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> $T_i$ 表示第 i 个时间段</td></tr></tfoot>
</table>

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table26-5)数据表 \@ref(tab:table26-4) 的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.5.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> Q 表示相应的非中心参数</td></tr></tfoot>
</table>

为比较每种药物的不同时间，两均值之差的标准误估计（见表 \@ref(tab:table26-3)）为

$$\widehat{s.e.}(\hat{\mu}_{ij}-\hat{\mu}_{íj})=\sqrt{2\hat{\sigma}_\varepsilon^2/n}=\sqrt{2(7.45)/8}=1.365$$

每种药物中用于比较时间均值的 5% LSD 值为

$$\mathrm{LSD}_{0.05}=t_{0.025,63}[\widehat{s.e.}(\hat{\mu}_{ij}-\hat{\mu}_{i^{\prime}j})]=(2.00)(1.365)=2.730$$

表 \@ref(tab:table26-6) 中给出了相同药物中时间均值的比较。

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table26-6)表 \@ref(tab:table26-6) 中数据在相同药物中时间均值比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.6.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> 在 5% 的显著水平上，在一列内具有相同字母的均值没有显著差异。$\text{LSD}_{0.05}=2.730$.</td></tr></tfoot>
</table>

由于时间水平是定量的且等距的，因此可以使用正交多项式来检查对每种药物的响应的线性和二次趋势。第一种药物的线性时间趋势的测量值为 $\theta_{LT_1}=-3\mu_{11}-1\mu_{12}+1\mu_{13}+3\mu_{14}$，其估计为

$$\hat{\theta}_{LT_1}=-3(70.50)-1(80.50)+1(81.00)+3(73.13)=8.39$$

相应的标准误估计为

$$\widehat{s.e.}(\hat{{\theta}}_{LT_1})=\sqrt{\frac{{\hat{\sigma}}_\varepsilon^2\sum_{i=1}^td_j^2}n}=\sqrt{\frac{7.45(9+1+1+9)}8}=4.316$$

相应的 $t$ 统计量为 $t_c=8.39/4.316=1.94$. 第一种药物的时间二次趋势的测量为

$$\theta_{QT_1}=1{\mu}_{11}-1{\mu}_{12}-1{\mu}_{13}+1{\mu}_{14}$$

其估计为

$$\hat{\theta}_{QT_1}=1\hat{\mu}_{11}-1\hat{\mu}_{12}-1\hat{\mu}_{13}+1\hat{\mu}_{14}=70.50-80.50-81.00+73.13=-17.87$$

其标准误估计为

$$\widehat{s.e.}(\hat{\theta}_{LT_1})=\sqrt{\frac{\hat{{\sigma}}_\varepsilon^2\sum_{i=1}^td_j^2}n}=\sqrt{\frac{7.45(1+1+1+1)}8}=1.930$$

相应的 $t$ 统计量为 $t_c=-17.87/1.930=9.259$.

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table26-7)表 \@ref(tab:table26-4) 中数据每种药物的线性和二次趋势估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.7.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> 括号中的值是相应的 $t$ 统计量</td></tr></tfoot>
</table>

<div class="figure" style="text-align: center">
<img src="figure/figure 26.1.png" alt="心率数据中每种药物各时间点的均值" width="682" />
<p class="caption">(\#fig:figure26-1)心率数据中每种药物各时间点的均值</p>
</div>

表 \@ref(tab:table26-7) 总结了所有药物在时间上的线性和二次趋势。药物 BWW9 呈负线性趋势，药物 AX23 呈强二次趋势。图 \@ref(fig:figure26-1) 中的图表显示了这些关系。为了在每个时间点相互比较药物，标准误估计为 $\widehat{s.e.}(\hat{{\mu}}_{ij}-\hat{{\mu}}_{i'j})=\sqrt{[2(\hat{{\sigma}}_\varepsilon^2+\hat{{\sigma}}_\delta^2)/n]}$. 为了评估这个量，必须估计每个方差分量。可以得到

$$\hat{{\sigma}}_\varepsilon^2=MSError(Time)=7.45$$

和

$$\hat{\sigma}_\delta^2=\frac{MSError(Subject)-\tilde{\sigma}_\varepsilon^2}p=\frac{111.33-7.45}4=25.97$$

那么

$$\widehat{s.e.}(\hat{\mu}_{ij}-\hat{\mu}_{i^{\prime}j})=\sqrt{2(\hat{\sigma}_\varepsilon^2+\hat{\sigma}_\delta^2)/n}=\sqrt{\frac{2(7.45+25.97)}8}=2.891$$

该标准误估计的 Satterthwaite 自由度估计为

$$\hat{v}=\frac{[c_1\hat{\sigma}_1^2+c_2\hat{\sigma}_2^2]^2}{\frac{c_2^2\hat{\sigma}_1^4}{v_1}+\frac{c_2^2\hat{\sigma}_2^4}{v_2}}=\frac{[\frac14(111.33)+\frac34(7.45)]^2}{\frac{\frac1{16}(111.33)^2}{21}+\frac{\frac9{16}(7.45)^2}{63}}=\frac{(33.42)^2}{36.888+0.496}=29.9$$

在每个时间点将对药物均值彼此比较的 5% 显著性水平下的近似 LSD 为

$$\mathrm{LSD}_{0.05}=t_{0.025,29.9}\times\widehat{s.e.}(\hat{\mu}_{ij}-\hat{\mu}_{i^{\prime}j})=(2.042)(2.891)=5.903$$

表 \@ref(tab:table26-8) 给出了各时间点药物之间的比较。

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table26-8)表 \@ref(tab:table26-4) 中数据各时间点药物之间的比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.8.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> 在 5% 的显著性水平上，在一行中具有相同字母的均值没有显著差异。$\text{LSD}_{0.05}=5.903$</td></tr></tfoot>
</table>

假设药物 AX23 和 BWW9 是实验性的，我们想在每个时间点将它们的平均值与对照进行比较。在时间 1，与该比较相关的对比为 $\theta=\mu_{11}+\mu_{21}-2\mu_{31}$. 其估计为 ${\hat{\theta}}=70.5+81.75-2\times72.75=6.75$，标准误估计为

$$\widehat{s.e.}\left(\hat{\theta}\right)=\sqrt{\frac{(\hat{\sigma}_{\varepsilon}^{2}+\hat{\sigma}_{\delta}^{2})\sum_{i=1}^{t}c_{i}^{2}}n}=\sqrt{\frac{(7.45+25.97)(1+1+4)}8}=5.006$$

相应的 $t$ 统计量为 $t_c=6.75/4.088=1.651$. 我们无法拒绝 $H_0{:{\theta}}=0$ 因为 $t_{0.025,29.9}=2.042$. 

表 \@ref(tab:table26-9) 给出了在每个时间段两种实验药物的平均值与对照药物的比较。结果表明，在 5% 显著水平下，在最后三个时间间隔，药物的平均值有显著性差异，但在时间 1 时没有差异。

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table26-9)在每个时间点药物 AX23 和 BWW9 的平均值与对照组的比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.9.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> * 表示在 5% 水平下显著，$t_{0.025,29.9}=2.042$</td></tr></tfoot>
</table>

### 示例 26.2：一个复杂的舒适度实验 {#sec26-2-2}

一位工程师有三种环境 (environments) 来检验两种类型的衣物 (clothing). 由于男性和女性对环境的反应也有所不同，因此性别 (sex) 也被视为一个因素。四人（两名男性和两名女性）被放入一个环境室 (chamber)（该室分配了三种环境之一）。一名男性和一名女性着衣物类型 1，另一名男性和女性则着衣物类型 2. 每个人的舒适度评分在 1 小时、2 小时和 3 小时结束时都被记录下来。该实验的数据如表 \@ref(tab:table26-10) 所示。实验单元有三种大小。最大的实验单元是一个室，或等价地，四个人的组。室实验单元的实验设计是完全随机设计结构中的单向处理结构（环境是处理），每个环境水平有三次重复。中等大小的实验单元是一个人 (person). 一个人的实验设计是随机完全区组设计结构中的双向处理结构（性别 × 衣物），共有九个区组（每个区组包含四个实验单元（人）。最小的实验单元是一个 1 小时的时间间隔，我们将其称为小时 (hour). 小时的实验设计是随机完全区组设计结构中的单向处理结构（时间），共有 36 个区组【每个区组包含三个实验单元（1 小时时间间隔）】。

<table>
<caption>(\#tab:table26-10)示例 [26.2](#sec26-2-2) 舒适性实验的数据</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.10.png">
</td>
  </tr>
</tbody>
</table>

该实验的模型为（其中模型被分为与三种尺寸的实验单元相对应的部分，并假设小时测量满足时间的裂区假设）

$$\begin{aligned}
Y_{ijkmn} =&\,{\mu}+E_{i}+{r}_{in} && \text{Chamber part}  \\
&+S_j+C_k+(SC)_{jk}+(ES)_{ij}+(EC)_{ik}+(ESC)_{ijk}+p_{ijkn}&& \text{Person part}  \\
&+T_m+(ET)_{im}+(ST)_{jm}+(CT)_{km}+(SCT)_{jkm}+(EST)_{ijm}+(ECT)_{ikm}&&\text{Hour part} \\
&+(ESCT)_{ijkm}+E_{ijkmn}
\end{aligned}$$

其中，E 表示环境，S 表示性别，C 表示衣服类型，T 表示时间，$r_{in}$ 表示随机室效应并假定其分布为 $i.i.d. N(0, \sigma_r^2)$，$p_{ijkn}$ 表示假定在随机人效应并假定其分布为 $i.i.d. N(0, \sigma_p^2)$，$\varepsilon_{ijkmn}$ 表示给定小时的测量误差并假定其分布为 $i.i.d. N(0, \sigma_\varepsilon^2)$。此外，所有 $r_{in},p_{ijkn},\varepsilon_{ijkmn}$ 都是独立分布的误差项。

表 \@ref(tab:table26-11) 给出了在上述误差项假设下的方差分析表。$F$ 统计量是在期望均方的指导下计算的。

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table26-11)示例 [26.2](#sec26-2-2) 方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.11.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> 表中的所有数值都四舍五入为两位有效数字，但计算是以双精度进行的。$Q_\ell$ 表示对应于给定效果的非中心参数</td></tr></tfoot>
</table>

表 \@ref(tab:table26-11) 中的误差平方和计算如下

$$
\begin{aligned}
SSERROR(PERSON) =&\,Replication \times Sex \times Clothing(Environment) SS\\
&+ Replication \times Sex(Environment) SS\\
&+ Replication \times Clothing(Environment) SS
\end{aligned}
$$

以及

$$SSERROR(CHAMBER)=Chamber(Environment)SS$$

分析的下一步是进行必要的比较。如果选择 $\alpha = 0.01$ 作为 I 类错误概率，则存在两个显著的交互作用：Environment × Time 和 Sex × Clothing. 对于 Sex × Clothing 交互作用，我们需要比较四种 Sex × Clothing 的均值。由于两种处理均应用于相同尺寸的实验单元（人），因此只需计算一个标准误。

令

$$\begin{aligned}\mu_{ijkm}=&\,\mu+E_i+S_j+C_k+(SC)_{jk}+(ES)_{ij}+(EC)_{ik}+(ESC)_{ijk}+T_m+(ET)_{im}\\&+(ST)_{jm}+(CT)_{km}+(SCT)_{jkm}+(EST)_{ijm}+(ECT)_{ikm}+(ESCT)_{ijkm}\end{aligned}$$

差值 $\bar{{\mu}}_{\cdot jk{\cdot}}-\bar{{\mu}}_{\cdot j'{k'\cdot}}$ 可被估计为 $\bar{Y}_{\cdot jk\cdot\cdot}-\bar{Y}_{\cdot j^{\prime}k^{\prime}\cdot\cdot}$，$\bar{Y}_{\cdot jk\cdot\cdot}-\bar{Y}_{\cdot j^{\prime}k^{\prime}\cdot\cdot}$ 的标准误估计为

$$\widehat{s.e.}(\bar{Y}_{\cdot jk\cdot\cdot}-\bar{Y}_{\cdot j^{\prime}k^{\prime}\cdot\cdot})=\sqrt{\frac{2MSERROR(PERSON)}{3\cdot3\cdot3}}=\sqrt{\frac{2\cdot(0.58)}{27}}=0.207$$

自由度为 18. 1% 显著性水平下的 Fisher's LSD 为

$$(t_{0.005,18})\mathrm{~}\widehat{s.e.}(\bar{Y}_{\cdot jk\cdot\cdot}-\bar{Y}_{\cdot j^{\prime}k^{\prime}\cdot\cdot})=2.878\cdot0.207=0.596$$

表 \@ref(tab:table26-12) 列出了 Sex × Clothing 的均值，所有这些都彼此显著不同。

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table26-12)示例 [26.2](#sec26-2-2) 带 LSD 的 Sex × Clothing 均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.12.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> $\text{LSD}_{0.01} = 0.595$. 所有均值基于 27 个观测</td></tr></tfoot>
</table>

Environment × Time 交互作用涉及两个不同尺寸的实验单元；因此，有两种类型的比较。首先，可以比较每个环境的两个时间平均值，并且可以检查每个环境的时间平均值之间的线性和二次趋势。其次，人们可以在每个时间点对环境进行比较。

相同环境下两个时间均值比较的标准误为

$$s.e.(\bar{Y}_{i\cdot\cdot m\cdot}-\bar{Y}_{i\cdot\cdot m^{\prime}\cdot})=s.e.\left[(\bar{r}_{i\cdot}-\bar{r}_{i\cdot})+(\bar{p}_{i\cdot\cdot\cdot}-\bar{p}_{i\cdot\cdot\cdot})+(\bar{\varepsilon}_{i\cdot\cdot m\cdot}-\bar{\varepsilon}_{i\cdot\cdot{m^{\prime}}\cdot})\right]=s.e.\left[(\bar{\varepsilon}_{i\cdot\cdot m\cdot}-\bar{\varepsilon}_{i\cdot\cdot{m^{\prime}}\cdot})\right]=\sqrt{\frac{2\sigma_\varepsilon^2}{2\cdot2\cdot3}}$$

以及

$$\widehat{s.e.}\left(\bar{Y}_{i\cdot\cdot m\cdot}-\bar{Y}_{i\cdot\cdot m^{\prime}\cdot}\right)=\sqrt{\frac{2\hat{\sigma}_\varepsilon^2}{2\cdot2\cdot3}}=\sqrt{\frac{2(0.06)}{12}}=0.01$$

同一环境下两时间点比较的 1% LSD 值为 $t_{0.005,48}(0.01)=2.682(0.01) = 0.268$.

表 \@ref(tab:table26-13) 包含 Environment × Time 均值以及每个环境的时间之间的比较。

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table26-13)示例 [26.2](#sec26-2-2) 在每个环境下比较时间的 Environment × Time 均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.13.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> 在给定的环境（列）中，相同字母的时间均值没有显著差异；$\text{LSD}_{0.01}=0.268$</td></tr></tfoot>
</table>

由于时间的三个水平的间隔相等，因此可以使用正交多项式来测量每个环境随时间的线性和二次趋势。环境 1 的线性趋势通过 $\hat{\theta}_{LT_1}=-\bar{Y}_{1\cdot\cdot1\cdot}+0\bar{Y}_{1\cdot\cdot2\cdot}+\bar{Y}_{1\cdot\cdot3\cdot}=-15.11+12.01=-3.10$ 来测量，其标准误估计为

$$s.e.\left(\hat{\theta}_{LT_1}\right)=\sqrt{\frac{\hat{\sigma}_\varepsilon^2[(-1)^2+(1)^2]}{12}}=\sqrt{\frac{(0.06)(2)}{12}}=0.1$$

相应的 $t$ 统计量为 $t_c = -3.10/0.1 = -31.0$，自由度为 48. 环境 1 中时间点的二次趋势通过 $\hat{{\theta}}_{QT_{1}}=1\bar{Y}_{1{\cdot\cdot}1{\cdot}}-2\bar{Y}_{1{\cdot\cdot}2{\cdot}}+1\bar{Y}_{1{\cdot\cdot}3{\cdot}}=15.11 - 2(9.11) + 12.01 = 8.90$ 来测量。其标准误估计由下式给出

$$s.e.\left(\hat{\theta}_{{Q}T_1}\right)=\sqrt{\frac{\hat{\sigma}_{\varepsilon}^2[(1)^2+(-2)^2+(1)^2]}{12}}=\sqrt{\frac{(0.06)(6)}{12}}=0.173$$

相应的 $t$ 统计量为 $t_c = 8.90/0.173 = -51.45$，自由度为 48.

表 \@ref(tab:table26-14) 给出了每种环境的线性和二次趋势。环境 1 和 2 存在显著的线性和二次趋势，但环境 3 则没有。图 \@ref(fig:figure26-2) 展示了这些关系。

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table26-14)示例 [26.2](#sec26-2-2) 种环境的线性和二次趋势的度量</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.14.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> $t$ 值在括号中</td></tr></tfoot>
</table>

<div class="figure" style="text-align: center">
<img src="figure/figure 26.2.png" alt="舒适度数据每种环境下各时间点的响应" width="684" />
<p class="caption">(\#fig:figure26-2)舒适度数据每种环境下各时间点的响应</p>
</div>

第二种 Environment × Time 比较是对同一时间或不同时间的不同环境进行比较。这种比较的标准误是

$$\begin{aligned}
s.e.\left(\bar{Y}_{i\cdot\cdot m\cdot}-\bar{Y}_{i^{\prime}\cdot\cdot m\cdot}\right)& =s.e.\left[(\bar{r}_{i\cdot}-\bar{r}_{i'\cdot})+(\bar{p}_{i\cdot\cdot m}-\bar{p}_{i'\cdot\cdot m})+(\bar{\varepsilon}_{i\cdot\cdot m\cdot}-\bar{\varepsilon}_{i'\cdot\cdot m\cdot})\right]  \\
&=\sqrt{\frac{2\sigma_r^2}3+\frac{2\sigma_p^2}{2\cdot2\cdot3}+\frac{2\sigma_\varepsilon^2}{2\cdot2\cdot3}}=\sqrt{\frac2{12}(\sigma_\varepsilon^2+\sigma_p^2+4\sigma_r^2)}
\end{aligned}$$

量 $\sigma_\varepsilon^2+\sigma_p^2+4\sigma_r^2$ 可以通过 $(1/3) MSERROR(CHAMBER) + (2/3) MSERROR(HOUR)$ 来估计。因此 Satterthwaite 自由度估计为

$$\hat{v}=\frac{[c_1\hat{\sigma}_1^2+c_2\hat{\sigma}_2^2]^2}{\frac{c_1^2\hat{\sigma}_1^4}{v_1}+\frac{c_2^2\hat{\sigma}_2^4}{v_2}}=\frac{[\frac13(29.21)+\frac23(0.06)]^2}{\frac{\frac19(29.21)^2}{6}+\frac{\frac49(0.06)^2}{48}}=\frac{(9.777)^2}{15.800+0.006}=6.05$$

因为

$$\sqrt{\frac{2\left[\frac13MSERROR(Chamber)+\frac23MSERROR(Hour)\right]}{12}}=\sqrt{\frac{2\left[\frac13(29.21)+\frac23(0.06)\right]}{12}}=1.276$$

用于比较每个时间点内的环境的 1% LSD 为 $\mathrm{LSD}_{0.01}=t_{0.005,6}(1.276)=(3.707)(1.276) = 4.73$. 表 \@ref(tab:table26-15) 包含在同一时间不同环境之间进行的Fisher LSD比较的 Environment × Time 平均值。此示例说明了为什么在一组交互均值中进行比较时使用正确的误差项非常重要：两个 LSD 值可能差异极大。在此示例中，这两个值为 0.268 和 4.73. 因此，用于比较相同时间水平下环境的 LSD 0.01 值比用于比较相同环境下时间的 LSD0.01 值的 17 倍以上。

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table26-15)用于在每个时间点比较环境均值的 Environment × Time 均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.15.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> 在给定的时间（行）内，相同字母的环境均值没有显著差异；$\text{LSD}_{0.01}=4.73$</td></tr></tfoot>
</table>

为了演示如何使用第 \@ref(sec24-3) 中的信息来计算使用多于两种尺寸实验单元的情况下的标准误和相应的 $t$ 统计量，我们构建了一个 LSD 来比较同一时间点内的环境和性别。假设感兴趣的比较为 $\bar{{\mu}}_{11{\cdot}1}-\bar{{\mu}}_{22{\cdot}1}$. 对此的最佳估计为 $\bar{Y}_{11{\cdot}1{\cdot}}-\bar{Y}_{22{\cdot}1{\cdot}}=8.04$. 接下来注意到该估计的方差是

$$\begin{aligned}
\mathrm{Var}(\bar{Y}_{11{\cdot}1{\cdot}}-\bar{Y}_{22{\cdot}1{\cdot}})& =\mathrm{Var}[(\bar{r}_{1\cdot}-\bar{r}_{2\cdot})+(\bar{p}_{11\cdot\cdot}-\bar{p}_{22\cdot\cdot})+(\bar{\varepsilon}_{11\cdot\cdot}-\bar{\varepsilon}_{22\cdot\cdot})]  \\
&=\frac{2\sigma_r^2}3+\frac{2\sigma_p^2}{2\cdot3}+\frac{2\sigma_\varepsilon^2}{2\cdot3}=\frac26(\sigma_\varepsilon^2+\sigma_p^2+2\sigma_r^2)
\end{aligned}$$

括号内的量估计为

$$\begin{aligned}(1/6)&MSERROR(CHAMBER)+(1/6)MSERROR(PERSON)\\&+(2/3)MSERROR(HOUR)=(1/6)(29.21)+(1/6)(0.58)+(2/3)(0.06)=5.005\end{aligned}$$

$\bar{Y}_{11{\cdot}1{\cdot}}-\bar{Y}_{22{\cdot}1{\cdot}}$ 的标准误估计为 $\widehat{s.e.}(\bar{Y}_{11\cdot1\cdot}-\bar{Y}_{22\cdot1\cdot})=\sqrt{[2(5.005)/6]}=1.2916$. 与该标准误估计相关的 Satterthwaite 自由度估计为

$$\frac{\left[(\frac16)(29.21)+(\frac16)(0.58)+(\frac23)(0.06)\right]^2}{\frac{(\frac16)^2(29.21)^2}6+\frac{(\frac16)^2(0.58)^2}{18}+\frac{(\frac23)^2(0.06)^2}{48}}=\frac{(5.005)^2}{3.9501+0.0005+0.00003}=\frac{25.05}{3.9506}=6.34$$

用于计算上述两个标准误的方法可以应用于其他情况，其中给定的比较可以划分为分量的总和，其中每个分量是仅针对一种尺寸的实验单元的比较。

### 示例 26.3：家庭态度 {#sec26-2-3}

来自农村和城市环境的家庭的态度每六个月测量一次，为期三个时间段。数据来自 7 个农村家庭和 10 个城市家庭，每个家庭由一个儿子、父亲和母亲组成，如表 \@ref(tab:table26-16) 所示。用于描述数据的模型先前在等式 \@ref(eq:26-2) 中给出。本节中的分析假设满足第 \@ref(sec26-1) 节中给出的理想条件。第 \@ref(chap27) 章将考虑不要求满足理想条件的更一般的分析。

<table>
<caption>(\#tab:table26-16)示例 [26.3](#sec26-2-3) 的家庭态度研究数据</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.16.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table26-17) 给出了与方程 \@ref(eq:26-2) 中的模型相对应的方差分析表。如果我们在 $\alpha=0.05$ 下操作，则有四个显著的效应：Area, Person, Time 和 Area × Time.

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table26-17)示例 [26.3](#sec26-2-3) 的家庭态度数据的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.17.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> $Q_\ell$ 表示对应于各个平方和的非中心参数</td></tr></tfoot>
</table>

因此，我们希望进行三种比较：人与人之间、每个地区的时间之间以及每个时间的地区之间。

家庭成员之间比较的方差为

$$\begin{aligned}
\mathrm{Var}(\bar{y}_{\cdot j\cdot\cdot}-\bar{y}_{\cdot j'\cdot\cdot}) &=\mathrm{Var}[(\bar{f}_{\cdot\cdot}-\bar{f}_{\cdot\cdot})+(\bar{p}_{\cdot j\cdot}-\bar{p}_{\cdot j^{\prime}\cdot})+(\bar{\varepsilon}_{\cdot j\cdot\cdot}-\bar{\varepsilon}_{\cdot j^{\prime}\cdot\cdot})] \\
&=\frac2{17}\sigma_p^2+\frac2{3\cdot17}\sigma_\varepsilon^2=\frac2{51}(\sigma_\varepsilon^2+3\sigma_p^2) 
\end{aligned}$$

两个家庭成员均值差的标准误估计为

$$\widehat{s.e.}(\bar{y}_{\cdot j\cdot\cdot}-\bar{y}_{\cdot j'\cdot\cdot})=\sqrt{\frac2{51}MSError(\text{Person})}=\sqrt{\frac2{51}(25.512)}=0.9906$$

5% LSD 值为 $LSD_{0.05}=2.042(0.99)=2.022$. 表 \@ref(tab:table26-18) 总结了成对比较。指定地区内两个时间均值之差的方差为

$$\mathrm{Var}(\bar{y}_{i\cdot k\cdot}-\bar{y}_{i\cdot k^{\prime}\cdot})=\mathrm{Var}\left[(\bar{f}_{i\cdot}-\bar{f}_{i\cdot})+(\bar{p}_{i\cdot}-\bar{p}_{i\cdot})+(\bar{\varepsilon}_{i\cdot k}.-\bar{\varepsilon}_{i\cdot k^{\prime}\cdot})\right]=\frac2{3n_{i}}\sigma_{\varepsilon}^{2}$$

其中 $n_i$ 表示地区 i 的家庭数量。

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table26-18)示例 [26.3](#sec26-2-3) 人均值的比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.18.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> 在一列中具有相同字母的均值没有显著差异。$\text{LSD}_{0.05}=2.022$</td></tr></tfoot>
</table>

城市地区两时间比较的标准误估计为

$$\widehat{s.e.}(\bar{y}_{1\cdot k\cdot}-\bar{y}_{1\cdot k^{\prime}\cdot})=\sqrt{\frac2{3n_1}\hat{\sigma}_\varepsilon^2}=\sqrt{\frac2{30}(0.370)}=0.157$$

以及在农村地区两时间比较的标准误估计为

$$\widehat{\mathrm{s.e.}}(\bar{y}_{2\cdot k\cdot}-\bar{y}_{2\cdot k^{\prime}\cdot})=\sqrt{\frac2{3n_2}\hat{\sigma}_\varepsilon^2}=\sqrt{\frac2{21}(0.370)}=0.188$$

相应的 5% LSD 为 $LSD_{0.05}(\text{Urban})=1.987(0.157)=0.312$ 和 $LSD_{0.05}(\text{Rural})=1.987(0.188) = 0.374$。表 \@ref(tab:table26-19) 总结了多重比较。

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table26-19)示例 [26.3](#sec26-2-3) 在每个地区内时间均值的比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.19.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> 在一列中具有相同字母的均值没有显著差异。农村的 $\text{LSD}_{0.05}=2.022$，城市的 $\text{LSD}_{0.05}=0.312$</td></tr></tfoot>
</table>

接下来，在给定的时间点上，城市和农村均值之差的方差是

$$\begin{aligned}
\mathrm{Var}(\bar{y}_{1\cdot k\cdot}-\bar{y}_{2\cdot k\cdot})& =\text{Var}\left[(\bar{f}_{1\cdot}-\bar{f}_{2\cdot})+(\bar{p}_{1\cdot}-\bar{p}_{2\cdot})+(\bar{\varepsilon}_{1\cdot k\cdot}-\bar{\varepsilon}_{2\cdot k\cdot})\right]  \\
&=\sigma_f^2\left(\frac1{10}+\frac17\right)+\sigma_p^2\left(\frac1{3\cdot10}+\frac1{3\cdot7}\right)+\sigma_\varepsilon^2\left(\frac1{3\cdot10}+\frac1{3\cdot7}\right)=\left(\frac1{30}+\frac1{21}\right)(\sigma_\varepsilon^2+\sigma_p^2+3\sigma_f^2)
\end{aligned}$$

上述表达式中方差分量的函数可以通过下式来估计

$$(1/3)MSError(Family)+(2/3)MSError(Time)=(1/3)(54.155)+(2/3)(0.370)=18.2983$$

因此，标准误估计为

$$\widehat{s.e.}(\bar{y}_{1\cdot k\cdot}-\bar{y}_{2\cdot k\cdot})=\sqrt{\left(\frac1{30}+\frac1{21}\right)(18.2983)}=1.217$$

与该估计相对应的 Satterthwaite 自由度估计是

$$\frac{\left[(\frac13)(54.155)+(\frac23)(0.370)\right]^2}{\frac{(\frac13)^2(54.155)^2}{15}+\frac{(\frac23)^2(0.370)^2}{90}}=\frac{(18.2983)^2}{21.7242+0.0007}=\frac{334.8278}{21.7249}=15.41$$

因此，在给定时间点比较城市和农村的近似 5% LSD 值为 $LSD_{0.05}=2.131(1.217)=2.59$。表 \@ref(tab:table26-20) 总结了多重比较。可以使用示例 [26.1](#sec26-2-1) 中描述的方法来研究每个地区中的时间的线性和二次趋势。

<table style="NAborder-bottom: 0;">
<caption>(\#tab:table26-20)示例 [26.3](#sec26-2-3) 每个时间点下地区均值的比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.20.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> 在一行中具有相同字母的均值没有显著差异。$\text{LSD}_{0.05}=2.59$</td></tr></tfoot>
</table>

## 使用 SAS-Mixed 程序的数据分析 {#sec26-3}

本节说明对上一节中描述的三个示例中的每一个进行 SAS^®^-Mixed 分析。考虑表 \@ref(tab:table26-4) 中的数据。为了重现第 \@ref(sec26-2) 节所示的结果，可以使用表 \@ref(tab:table26-21) 中所示的 SAS 命令。

<table>
<caption>(\#tab:table26-21)用于分析表 \@ref(tab:table26-4) 中数据的 SAS-Mixed 代码</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.21.png">
</td>
  </tr>
</tbody>
</table>

### 示例 26.1 {#sec26-3-1}

应注意的一部分输出是提供两个方差分量 $\sigma^2_\varepsilon,\sigma^2_\delta$ 的估计的部分。表 \@ref(tab:table26-22) 给出了估计值。可以看出，$\hat\sigma^2_\varepsilon=7.45,\hat\sigma^2_\delta=25.97$. 

<table>
<caption>(\#tab:table26-22)表 \@ref(tab:table26-4) 中数据的协方差参数估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.22.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table26-23) 中给出了与表 \@ref(tab:table26-5) 中 ANOVA 相对应的检验。请注意，MIXED 程序没有提供误差平方和，但它确实给出了与表 \@ref(tab:table26-5) 中给出的相同的 $F$ 统计量。

<table>
<caption>(\#tab:table26-23)表 \@ref(tab:table26-4) 中数据药物、时间和药物 × 时间的统计检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.23.png">
</td>
  </tr>
</tbody>
</table>

药物主效应均值、时间主效应均值和药物时间双向均值如表 \@ref(tab:table26-24) 所示，这些均值的子集之间的成对比较如表 \@ref(tab:table26-25) - \@ref(tab:table26-27) 所示。表 \@ref(tab:table26-28) 中给出了可在 MIXED 程序中用于计算每种药物的线性和二次对比的 estimate 选项。这些选项的结果如表 \@ref(tab:table26-29) 所示。

<table>
<caption>(\#tab:table26-24)表 \@ref(tab:table26-4) 中数据药物、时间和药物 × 时间的均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.24.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table26-25)表 \@ref(tab:table26-4) 中数据药物主效应均值和时间主效应均值之间的成对比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.25.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table26-26)表 \@ref(tab:table26-4) 中数据在每种药物内时间均值的比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.26.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table26-27)表 \@ref(tab:table26-4) 中数据在每个时间下药物均值的比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.27.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table26-28)用于计算表 \@ref(tab:table26-4) 中数据中每种药物内时间的线性和二次对比的 Estimate 选项</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.28.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table26-29)表 \@ref(tab:table26-4) 中数据中每种药物内时间的线性和二次对比</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.29.png">
</td>
  </tr>
</tbody>
</table>

### 示例 26.2 {#sec26-3-2}

本节提供复杂舒适度实验的 Mixed 分析。使用表 \@ref(tab:table26-30) 中所示的 SAS 命令对示例 [26.2](#sec26-2-2) 中的数据进行分析。表 \@ref(tab:table26-31) 给出了与室、人和时间相对应的三个方差分量的估计。固定效应的 III 型检验如表 \@ref(tab:table26-32) 所示。表 \@ref(tab:table26-33) 给出了每种环境下时间的线性和二次对比。该表还包含时间 1 下的环境 1 与时间 2 下的环境 2 的比较，以说明如果有人对此类比较感兴趣，应该如何进行比较。表 \@ref(tab:table26-34) 给出了 Sex × Clothing 均值以及 Environment × Time 均值。表 \@ref(tab:table26-35) - \@ref(tab:table26-37) 提供了这些均值的不同子集之间的成对比较。

<table>
<caption>(\#tab:table26-30)用于分析复杂舒适度实验的 SAS 命令</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.30.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table26-31)舒适度研究的室、人和时间方差分量估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.31.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table26-32)复杂舒适度实验固定效应的 III 型检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.32.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table26-33)舒适度研究中每种环境下时间的线性和二次对比</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.33.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table26-34)舒适度研究中 Sex × Clothing 均值和 Environment × Time 均值</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.34.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table26-35)舒适度研究中 Sex × Clothing 均值的成对比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.35.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table26-36)舒适度研究中每种环境内时间的成对比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.36.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table26-37)舒适度研究中每种环境内时间的成对比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.37.png">
</td>
  </tr>
</tbody>
</table>

### 示例 26.3 {#sec26-3-3}

本节提供家庭成员实验的 Mixed 分析。使用表 \@ref(tab:table26-38) 中所示的 SAS 命令对示例 [26.3](#sec26-2-3) 中的数据进行分析。表 \@ref(tab:table26-39) 给出了与地区、家庭成员和时间相对应的三个方差分量的估计。固定效应的 III 型检验如表 \@ref(tab:table26-40) 所示。表 \@ref(tab:table26-41) 给出了 Area × Time 均值。表 \@ref(tab:table26-42) 提供了这些均值的不同子集之间的成对比较。

<table>
<caption>(\#tab:table26-38)用于分析家庭态度实验的 SAS 命令</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.38.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table26-39)家庭态度实验的协方差参数估计</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.39.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table26-40)家庭态度实验固定效应的 III 型检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.40.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table26-41)家庭态度实验地区 × 时间双向均值检验</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.41.png">
</td>
  </tr>
</tbody>
</table>

<table>
<caption>(\#tab:table26-42)家庭态度实验地区 × 时间双向均值的成对比较</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%2026.42.png">
</td>
  </tr>
</tbody>
</table>

## 结束语 {#sec26-4}

本章针对涉及重复测量的三个示例（在这些示例中，时间的裂区假设成立），描述了重复测量设计的分析。讨论中包括了每个示例的模型和假设。给出了计算多重比较标准误的公式以及研究均值间各种对比的方法。在第 \@ref(sec26-3) 节中，使用 SAS-Mixed 程序进行了分析，可以与第 \@ref(sec26-2) 节中给出的结果进行比较。

## 练习 {#sec26-5}
