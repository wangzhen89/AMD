
# 实验设计基础 {#chap4}

经过适当设计和分析的实验提供了有关所使用资源的调查条件的最大信息量。本章介绍了实验者设计和分析实验时使用的概念和方法。本章讨论的基本概念是**处理结构** (treatment structure) 和**设计结构** (design structure)，以及**重复** (replication)、**区组** (block) 和**实验单元** (experimental unit) 的思想。给出了设计和处理结构相结合的实例，以说明完全区组设计和不完全区组设计的概念。这些设计涉及一种尺寸的实验单元。第五章描述了实验单元尺寸的概念，并给出了涉及多于一种实验单元尺寸的各种设计。本章提出的设计结构包括**完全随机化** (completely randomized, CRD)、**随机完全区组** (randomized complete block, RCBD)、**不完全区组** (incomplete block, IBD)和**拉丁方** (Latin square, LSD)。所考虑的处理结构包括单向、双向、双向对照、部分析因和 $n$ 向结构。在本一章中，介绍了具有必要的变异源和自由度的模型和方差分析表。讨论提供了确定用于计算误差平方和的变异源的方法，以及用于计算结果自由度的算法。一般来说，误差平方和是通过对观测值或经过相同处理的观测值的线性组合 (observations that are treated alike) 的比较来获得的。其他平方和的计算将在后面的章节中讨论。本章的基本方法是通过示例来演示概念。裂区 (split-plot)、重复测量 (repeated measures)、裂条区 (strip-plot) 和交叉设计 (crossover designs ) 使用不同**尺寸** (size) 的实验单元的概念，并在第 \@ref(chap5) 章中进行了描述。涉及嵌套 (nesting) 的设计也在第 \@ref(chap5) 章中进行了讨论。

## 介绍基本概念 {#sec4-1}

实验设计涉及计划实验，以便从可用资源中获得最多的信息。实验的设计应该从要达到的目的开始。应指定一个或多个主要响应变量，并充分描述实验目的和响应变量之间的联系。应解释用作实验单元的实体，并将其与实验目的联系起来。通常，在设计实验时，实验者可以控制某些被称为处理、总体或处理组合的因素。实验者通常控制要使用的实验单元的选择以及这些实验单元是否可以分组（称为区组）。典型的实验涉及要进行比较或要研究其效应的处理（或处理组合）。

在进行实验之前，必须回答几个问题：

1) **实验的目的是什么？响应变量是什么？**这是设计的一个非常重要的方面，因为必须了解生成数据的过程以及数据与目的的关系。
2) **需研究多少种处理？**处理的数量可能已经被指定，但有时应该讨论处理的选择如何与研究的目的相关。例如，一名营养学专业的学生想要设计一项研究，以确定在饼干产品中可以使用的大豆粉的最大量来代替小麦粉，这样所得的饼干中就没有大豆味了。该学生选择的饼干配方涉及用大豆粉替代 0、10、20、30、40、50、60、70、80、90 和 100% 的小麦粉。第一个问题是：“如果使用100%大豆粉，你能尝到大豆的味道吗？”因为如果你不能品尝到100%大豆粉的大豆味，那么就没有必要进行实验，如果你能品尝到大豆味，那么就没有必要在一组处理中包括100%大豆粉。我们认为尚不清楚用 20%、30%、40% 和 50% 大豆粉制成的产品是否可以品尝到大豆味。即，确定大豆粉含量为10％时无法品尝到大豆味，而大豆粉含量超过50％时则可品尝到大豆味。确定了不能品尝到含有10%大豆粉的大豆味，而可以品尝到含有超过50%大豆粉的大豆味，在将处理的选择与研究目的联系起来后，只需要最初11个处理中的5个（0、20、30、40和50%大豆粉）。0%作为小麦粉对照。这一过程大大减少了研究所需的样本数量。
3) **每种处理需要观测多少次？**这个问题与实现特定目的所需的样本量有关。
4) **实验单元是什么？**许多实验只涉及一个尺寸的实验单元。但其中涉及的重要思想是独立的重复。一些实验涉及多个尺寸的实验单元，这一概念在第 \@ref(chap5) 章中有所描述。研究人员通常认为一个实验只涉及一个尺寸的实验单元，而无法识别涉及多种尺寸的实验单元的情况。这个问题需要认真解决，如第 \@ref(chap5) 章所述。
5) **实验者如何将处理应用于可用的实验单元，然后对响应进行观测？**这个问题涉及使用随机化将处理分配给实验单元，以及在提供数据过程的其他部分中使用随机化。重要的是，随机化用于将处理分配给实验单元。在医学研究中，以随机顺序在实验室中运行样本或使用随机顺序评估受试者也很重要。
6) **是否可以对最终设计进行分析或进行所需的比较？**这可能是研究人员最重要的问题，但本书的一个主要目的是使读者能够使用更复杂的设计，并且仍然能够进行分析，以估计重要参数并检验所需的假设。

这些问题的答案不一定是直截了当的，也不能笼统地回答这些问题。希望这里讨论的想法和概念能帮助实验者收集足够的信息，为他们的研究提供答案。

为了继续，考虑一个涉及 $t$ 个处理的实验，其中每个处理应用于 $r$ 个不同的实验单元。令 $y_{ij}$ 是分配给第 $i$ 个处理的第 $j$ 个实验单元观测到的响应，用于描述 $y_{ij}$ 的数学模型为

\begin{equation}
y_{ij}=\mu_i+\varepsilon_{ij}\quad\mathrm{for~}i=1,2,\ldots,t,\mathrm{~and~}j=1,2,\ldots,r
(#eq:4-1)
\end{equation}

其中 $\mu_i$ 是第 $i$ 个处理的响应均值的未知真实值，$\varepsilon_{ij}$ 是一个随机变量，表示由自然变异和其他可能的随机和非随机误差源引起的噪声。研究人员应尽力控制非随机误差源，包括模型误差、测量误差、观测误差和处理误差的错误指定。为了进行这个实验，研究人员必须选择 $rt$ 个实验单元，然后将 $r$ 个实验单元随机分配给每个处理。**该过程的随机化部分对于防止处理分配中出现偏差非常重要**。正是由于实验单元被随机分配给处理，随机化或**置换分析** (permutation analysis) 可用于开发适当分析的理论 (Kempthorne, 1952). 关于随机化的使用，至少可以说，它可以防止在实验中引入系统偏差。如果实验者不使用随机化，则无法判断观测到的差异是由于实验单元对处理的反应差异造成的，还是由于用于将实验单元分配给处理的系统方法造成的。

**实验的统计目的是比较在实验单元上观测到的处理响应**。例如，如果研究人员想比较高血压化合物对人类血压的影响，那么在研究中使用小白鼠作为实验单元将无法对人类做出推断。只能对实验单元的总体进行推断，研究中使用的实验单元是代表性样本。对想要进行推断的实验单元的总体进行表征是非常重要的。研究中使用的实验单元样本必须从总体中随机选择，以对该总体进行推断。通常不可能从要纳入研究的实验单元总体中随机选择实验单元。但至少实验单元的样本必须能够代表总体或一个**概念总体** (conceptual population)。例如，在对生产线进行抽样时，人们希望为研究选择的物品能够代表尚未生产的物品总体。如果生产过程没有改变，那么可以合理地假设所选物品将代表未来的物品。描述实验中使用的实验单元所代表的实验单元的总体是非常重要的。如果实验单元不代表要对其进行**推断** (inference) 的实验单元的总体，则不能进行推断，相反，人们可以**推测** (conjecture) 处理对未抽样总体的效应。例如，从一项涉及白鼠的研究中获得的关于高血压化合物的信息可以用来推测（而不是推断）高血压化合物对人类的效应。

当为研究选择实验单元时，当研究中使用的实验单元组是**同质的** (homogeneous) 或非常相似时，可以获得更好的处理之间的比较。在许多实验中，不可能选择 $rt$ 个相同的实验单元。不相同的实验单元通过 $\varepsilon_{ij}$ 导致数据中的噪声。当实验单元不是同质的或不相似时，有三种方法可以用来解释实验单元之间的差异，一种方法是测量描述实验单元之间差异的特征，如体重和年龄，并使用 Milliken and Johnson (2001) 中描述的**协方差分析** (analysis of covariance)。第二种方法是将实验单元分成几组几乎相同的实验单元。几乎相同的实验单元称为同质的。在这种情况下，处理可以在组内的相似实验单元上进行比较，其中组与组之间的变异可以在分析中被解释。相似实验单元的组被称为**区组** (block)。解释实验单元之间变异性的第三种方法是同时使用区组和协变量 (covariate) 的分析 。假设有 $r$ 个区组，每个区组有 $t$ 个实验单元，每个区组中每个处理发生一次。表示第 $j$ 个区组中第 $i$ 个处理的观测到的响应的模型为

\begin{equation}
y_{ij}=\mu_i+b_j+\varepsilon_{ij}^*\quad\mathrm{~for~}i=1,2,\ldots,t,\mathrm{~and~}j=1,2,\ldots,r
(#eq:4-2)
\end{equation}

对于模型 \@ref(eq:4-2)，模型 \@ref(eq:4-1) 中的 $\varepsilon_{ij}$s 替换为 $\varepsilon_{ij} = b_j + \varepsilon_{ij}^*$；也就是说，实验单元组或区组之间的变异已被识别并从 $\varepsilon_{ij}^*$ 中分离出来，它代表区组内实验单元的变异。通过将区组效应与实验单元分离，**区组内变异** (within-block variation) 可用于比较处理效应，这涉及计算处理对比 (contrasts) 的标准误估计。

两种处理（或任何处理对比）可以在没有区组效应的情况下进行比较，方法是将两种处理的响应的区组内差异 (within-block differences) 视为

$$\begin{aligned}
y_{ij}-y_{i^{\prime}j}& =(\mu_i+b_j+\varepsilon_{ij}^*)-(\mu_{i^{\prime}}+b_j+\varepsilon_{i^{\prime}j}^*)  \\
&=\mu_i-\mu_{i^{\prime}}+\varepsilon_{ij}^*-\varepsilon_{i^{\prime}j}^*
\end{aligned}$$

这不取决于区组效应 $b_j$。这种差异的结果是，一个区组内两种处理响应差异的变异取决于实验单元之间的区组内变异，而不是**区组间变异** (between-block variation)。

实验设计的目的是选择实验材料并进行分组，以便尽可能减少实验中各组内实验单元之间的噪声或实验误差。因此，要比较处理的实验单元应该尽可能相似，这样两种处理之间的较小差异就可以被检测为显著差异。

如果有 $t$ 个处理和 $t$ 个实验单元，则可以进行实验，并且可以根据观测结果估计每个处理的平均值。但无法获得误差方差的估计值。只有当部分或全部处理被重复时，才能获得误差方差的估计值（与 $\varepsilon_{ij}$ 或 $\varepsilon_{ij}^*$ 相关）。处理的**重复** (replications) 是对处理的独立观测。因此，一种处理的两次重复必须涉及两个实验单元。实验单元是应用了处理的实体。但处理的重复必须涉及过程，以便在每个实验单元上独立应用并观测处理。因此，假设一名研究人员想比较两种饮食对兔子生长速度的影响，并且他有 $10$ 只兔子可用于研究。这个过程是将五只兔子随机分配到每种处理中。但是，如果他将分配给其中一种饮食的五只兔子放在一个笼子里，而将分配给第二种饮食的五只兔子放在另一个笼子里，而同一个笼子里的兔子则从一个普通的碗里喂食这种饮食，结果会怎样呢？每只兔子被随机分配了饮食，但由于所有五只兔子都被放入一个笼子并用一个共同的碗喂食，因此不能对兔子进行独立观测。在这种情况下，实验单元是笼子而不是单个兔子。笼子之间是独立处理的，但笼子里的兔子却不是独立处理的。因此，笼子里的兔子不提供独立的饮食重复。了解进行研究的完整过程非常重要，以便能够了解何时发生独立重复以及何时不发生。必须遵守此定义以确定实验期间何时使用重复。研究人员经常使用**复制** (duplicate) 或分割样本 
(split samples) 来生成两个观测结果，并将它们称为**重复** (replications)，而实际上，它们实际上是**子样本** (sub-samples) 或**重复测量** (repeated measures)。**复制项必然不会提供与独立重复项相同的信息**。

假设一名研究人员想要研究某所大学男女学生的身高差异。假设大学有 $22,000$ 名学生。一个过程可以是从大学学生中随机选择 $100$ 名女学生和 $100$ 名男学生。下一步是找到随机样本中的所有学生并测量他们的身高，从而产生 $100$ 个女性身高测量值和 $100$ 个男性身高测量值。要找到所有 $200$ 名学生，工作量似乎很大。相反，假设研究人员选择一名女性和一名男性，并测量每人的身高 $100$ 次。第二个过程产生 $100$ 个女性身高测量值（本例中为一名女性）和 $100$ 个男性身高测量值（本例中也是一名男性）。每个数据集都有 $200$ 个观测值，但第一个过程中一个数据集内 $100$ 个测量值之间的变异性 (variability) 提供了大学女生和男生之间的变异性的测量。在第二种情况下，一组中的 $100$ 个测量值的变异性提供了对同一个人进行的 $100$ 个测量值之间的差异的度量。第二种情况提供了有关测量过程变异性的信息（对同一个人的重复测量），但没有提供关于特定大学的女性或男性身高之间的变异性的信息。对一个人身高的独立测量不能提供对总体身高的真实变异的测量。对 $100$ 名女性身高的测量提供了 $100$ 次重复，而对同一个人的 $100$ 次测量提供了对该人的重复测量。对同一人的这 $100$ 项测量被称为重复测量或子样本，但不是重复测量。

一位面包师进行了一项实验，比较了三种防腐剂在某种类型的蛋糕产品中抑制霉菌生长的能力。面包师把每种防腐剂混合在一起烤了一个蛋糕。每立方厘米蛋糕的霉菌孢子数量是在储存 $9$ 天后测量的。面包师想要 $10$ 次重复的分析，所以他把每个蛋糕分成 $10$ 片，并获得了每片上的孢子数量。然而，这 $10$ 次测量并不是由 $10$ 种单独使用的防腐剂产生的。他的子样品测得的变异是蛋糕内变异 (within-cake variation) 的指标，而不是防腐剂实验单元间变异 (experimental-unit-to-experimental-unit variation) 或蛋糕间变异 (cake-to-cake variation) 的指标。为了得到每种防腐剂的 $10$ 次重复，面包师需要用每种防腐剂烘焙 $10$ 个蛋糕。这些蛋糕需要相互混合并独立烘焙。将一大批蛋糕面团搅拌在一起，加入防腐剂，然后将混合物倒入 $10$ 个蛋糕平底锅中，可能会更容易。这一过程提供了 $10$ 个蛋糕，但这些蛋糕并不是相互独立地混合在一起的，因此它们不是独立的重复。面包师需要将 $10$ 批面团与防腐剂混合，然后在不同的烤箱 （或在不同的时间）从每批面团中烘焙一个蛋糕，以获得 $10$ 种防腐剂的重复。

[^row]: 也可译作 “排”，即图 \@ref(fig:figure4-1) 中的竖线。

另一个非重复的例子涉及一些研究人员所称的**条状试验** (strip trial)。在农学中，条状试验包括在一行 (rows)（或行组，a group of rows）中种植给定品种的所有种子，每行（或行组）种植不同的品种。然后将行（或行组）划分为八个部分，这些部分称为重复。图 \@ref(fig:figure4-1) 中的图表表示六个**裂条区** (strip plot)，每个裂条区有四行[^row]，并被分成八个部分表示为重复。用条状试验代替八次独立重复的优点是，研究人员不需要在条状试验中不断地将种子从一个播种箱位置换到另一个播种箱位置，而且种植计划非常简单。现在，如果她想进行精心设计的实验，她就必须根据特定的随机方案更换播种箱八次。进行这样的随机化将在**随机完全区组设计结构** (randomized complete block design structure) 中提供八个区组。如果实验者像分析随机完全区组设计一样分析条状试验数据，那么她的分析将是错误的。事实上，条状试验不能用来推断品种差异，因为每个品种只有一次独立的观测。四行或裂条区是应用该品种的实验单元。在条状试验中，研究人员可以轻松地将裂条区划分为 $20$ 个部分 (parts) 或 $100$ 个部分；毕竟，有了更多的“复制”，人们可以根据两种均值较小差异检测出显著的结果。但是，在这些**裂条** (strips) 中的每一条上进行的这些测量都不是真正的重复；相反，它们是裂条区的子样本或重复测量。因此，获得越来越多的部分无助于检验均值之间的差异。确定某个部分或观测是否是真实重复的一个检验如下：如果研究人员可以通过分裂同样轻松地获得更多“重复”，那么她就没有获得真正的重复，而是获得了子样本或重复测量。区分子样本和重复样本非常重要，因为从子样本之间估计的误差方差通常比从重复样本或实验单元之间估计的误差方差小得多。使用从子样本计算出的误差方差构建的 $F$ 统计量将比应有的值大得多，导致实验者确定更多的差异，因为其统计显著性超过了应有的差异。

<div class="figure" style="text-align: center">
<img src="figure/figure 4.1.png" alt="裂条区设计示意图，其中四个行和六个品种排列在八个伪重复中" width="464" />
<p class="caption">(\#fig:figure4-1)裂条区设计示意图，其中四个行和六个品种排列在八个伪重复中</p>
</div>

## 设计实验的结构 {#sec4-2}

设计实验 (designed experiment) 由两个基本结构组成，能够识别和区分每个结构至关重要。在为特定设计构建适当的模型和分析之前，实验设计中使用的所有因素必须分类为属于处理结构或设计结构。以下定义和讨论用于帮助描述两种结构之间的差异。

::: {.definition #4-1}
♦

设计实验的**处理结构** (treatment structure) 由实验者选择用于研究和/或比较的一组处理、因素、处理组合或总体组成。
:::

处理结构是根据那些要比较的因素或处理构建的，通过它们对给定响应变量的效应来测量。必须选择处理结构中的因素来实现既定的实验目的。处理结构可以是一组处理，称为单向处理结构，也可以是一组处理组合，例如双向因素安排 (arrangement) 或高阶因素安排，加上任何对照或标准处理。

::: {.definition #4-2}
♦

设计实验的**设计结构** (design structure) 包括将实验单元分组为同质组或区组。
:::

设计实验的设计结构涉及用于形成实验单元组的因素，以便观测处理的条件尽可能统一。如果所有实验单元都非常同质，则只需要一组或一组观测，并且实验单元可以完全随机地分配给处理。这样的设计结构称为**完全随机设计结构** (completely randomized design structure).

如果需要多于一组的实验单元，使得组内的单元比组间的实验单元更加同质，则设计结构涉及某种类型的区组设计。有多种因素可用于形成实验单元的区组，但是，如下所述，用于构建区组的一个或多个因素不得与处理结构中的因素发生交互 (interact). 一旦选择了处理结构和设计结构，就通过准确描述将处理结构的处理随机分配（随机化）到设计结构中的实验单元的方法来指定设计的实验。

因此，总的设计实验包括 1）处理结构的选择，2）设计结构的选择和 3）将处理结构中的处理或处理组合随机分配给设计结构中的实验单元的方法。图 \@ref(fig:figure4-2) 显示了设计实验的两个部分。

<div class="figure" style="text-align: center">
<img src="figure/figure 4.2.png" alt="使用随机化将处理结构与设计结构相结合以形成实验总体设计的图示" width="512" />
<p class="caption">(\#fig:figure4-2)使用随机化将处理结构与设计结构相结合以形成实验总体设计的图示</p>
</div>

总的设计实验规定了用于获得适当分析的适当模型。在构建描述设计的模型时，对设计和处理结构进行了两个基本假设。首先，假设设计结构的组成部分是随机效应，即设计中使用的区组是来自实验单元的可能区组的总体的随机样本。这一假设意味着研究人员希望对大量的实验单元进行推断。其次，假设设计结构的组成部分（或因素）与处理结构的组成部分（或因素）之间没有交互作用。换言之，假设处理之间存在的关系从一个区组到另一个区组是一致的（随机变异除外），或者换句话说，区组因素不会影响处理之间的关系。许多教科书将区组因素描述为**讨厌因素** (nuisance factors)，而没有解决讨厌因素与处理结构中的因素交互作用的可能性 (Cobb, 1997). 如果可以发生这种交互作用，则这些讨厌因素必须包括在处理结构中，或者讨厌因素的值可以被视为具有不等斜率可能性的协变量 (Milliken and Johnson, 2001). 必须仔细评估可能用于构建区组的因素的选择，以防止与处理结构中的因素可能发生交互。

设计结构的选择是利用实验单元的所有可用知识，并且是独立于处理结构选择的（不要让处理结构影响不良设计结构的选择）。同样，实验者应该在不了解设计结构的情况下选择处理结构（不要让设计结构妨碍选择必要的处理集合）。在指定了适当的设计结构并选择了所需的处理结构之后，可能需要在一个或两个结构中做出一些妥协，以使它们彼此兼容，并使实验者能够进行有效的实验。

### 设计结构类型 {#sec4-2-1}

设计结构由区组类型或将实验单元分组为同质组来确定，并由用于形成区组的因素来指定。以下是对一些常见设计结构的描述。有两种基本的设计结构，它们被描述为**完全区组** (complete block) 和**不完全区组** (incomplete block) 设计结构。\@ref(sec4-3) 节对几种特定设计结构进行了更详细的讨论。

1) **完全随机设计结构** (Completely randomized design structure). 在完全随机设计结构中，假设所有实验单元都是同质的，并且实验单元完全随机地分配给处理。通常，这些处理被分配相同数量的实验单元，尽管这不是必需的。当实验单元不是同质的，并且实验者不能确定允许将它们分组到更同质的组中的因素时，也可以使用这种设计结构。当区组因素的值被用作协变量而不是使用它们来形成区组时，可以使用协方差分析 (Milliken and Johnson, 2001).
2) **随机完全区组设计** (Randomized complete block design). 如果有 $t$ 个处理，那么随机完全区组设计结构包括具有实验单元的区组，每个区组中有 $t$ 个或多个实验单元。当区组尺寸 (block size) 等于或大于处理数时，有可能在每个区组内出现一组完整的处理，因此被称为随机完全区组。如果每个区组的区组尺寸恰好等于 $t$，则每个处理恰好被随机分配到每个区组内的一个实验单元。如果每个区组内有 $t$ 个以上的实验单元，那么每个处理可以分配到一个实验单元，一些处理（也可能是全部）可以分配到多个实验单元。如果 $t=5$ 并且区组的尺寸为 $8$，则每个处理将被分配到每个区组内的实验单元，并且其中三个处理将被分配到一个额外的实验单元以进行这些处理的总共两次重复。一种策略是制定一种区组的处理模式，其中每个处理的观测数量尽可能**均衡** (balanced)。图 \@ref(fig:figure4-3) 包含一种将五种处理分配给三个尺寸为 $8$ 的区组的方法。在每个区组内，将各自的处理随机分配给该区组内的实验单元。所示的安排是一组可能的随机化。如果每个区组由 $c\times t$ 个实验单元组成，其中 $c$ 是整数，则每个处理可以分配给每个区组内的 $c$ 个实验单元。这也是一个随机完全区组设计结构。随机完全区组设计结构是指区组内实验单元的数量大于或等于处理数量的任何区组方案，因此可以将一整套处理分配给每个区组中的实验单元。

<div class="figure" style="text-align: center">
<img src="figure/figure 4.3.png" alt="随机完整区组设计结构，在尺寸为 $8$ 的区组中包含五种处理" width="552" />
<p class="caption">(\#fig:figure4-3)随机完整区组设计结构，在尺寸为 $8$ 的区组中包含五种处理</p>
</div>
3) **拉丁方设计** (Latin square design). 拉丁方设计结构由两个方向的区组组成。对于涉及 $t$ 个处理的实验，$t^2$ 个实验单元被安排为 $t\times t$ 的**方块** (square)，其中行被称为**行区组** (row blocks)，列被称为**列区组** (column blocks)。因此，实验单元的这种 $t\times t$ 安排在两个方向上被分区 (blocked)：行区组和列区组。为了构建拉丁方设计结构，将处理随机分配到方块中的实验单元，使得每种处理在每一行区组中发生一次且仅发生一次，在每一行区组中发生一次且仅发生一次。参见 Cochran and Cox (1957) 以了解在行区组和列区组中的各种处理安排。两个或多个方向的区组在许多学科中很常见，特别是当实验单元出现在矩形中时（其中一个维度可能是时间），按行和列进行分区很有用。Graeco–Latin 方块可用于形成三个方向的区组 (Cochran and Cox, 1957).
4) **不完全区组设计** (Incomplete block design)。当一个或多个区组中的处理数超过实验单元数时，即出现不完全区组设计。当这种情况发生时，每个区组中不可能出现完整的处理集，因此被称为“不完全区组”。有几种特殊的不完全区组设计结构，如**均衡不完全区组** (balanced incomplete blocks) 和**部分均衡不完全区组** (partially balanced incomplete blocks). 均衡不完全区组设计结构是指将处理分配给区组，使得每对处理在同一区组中出现的次数相等。当一组处理在区组内一起出现相同的次数，而其他处理在一些区组内一起出现不同的次数时，就会出现部分均衡不完全区组设计结构。第 \@ref(chap5) 章中讨论的裂区设计结构是部分不完全区组设计结构的一个示例。例 [4.5](#sec4-3-5) 描述了一些不完全区组设计结构。
5) **各种组合和推广** (Various combinations and generalizations). 有多种方法可以对实验单元进行分组。有时分组不满足上述定义，但仍然提供了有效的设计结构。一个例子是区组尺寸因区组而异，其中一些区组是不完全的，而另一些区组是完全的。在任何情况下，这些其他的区组方案都可以为实验者提供非常可行的设计结构，以便进行有效的实验。

### 处理结构类型 {#sec4-2-2}

处理结构由实验者希望研究的各种处理或处理组合或因素和因素组合组成。处理结构组成部分的选择应与实验目的相关，实验目的应在研究方案或描述中具体说明。接下来描述了一些常见的处理结构类型，每一种都在 \@ref(sec4-3) 节中进行了更详细的讨论。

1) **单向处理结构** (One-way treatment structure). 单向处理结构由一组 $t$ 个处理或总体组成，其中处理之间没有假定的结构。处理之间可能存在关系，例如使用四个温度，120、130、150和160°C。如果处理是通过组合两个或多个因素构建的，则这些因素不会用于表示模型中的处理。单向处理结构可用于表示任何一组处理，当一些可能的处理组合缺失时，单向处理结构通常用于表示**析因** (factorial) 处理结构。第 \@ref(chap13) 章使用了这种方法。
2) **双向处理结构** (Two-way treatment structure). 双向处理结构由通过组合两个不同因素的水平 (levels) 或可能性 (possibilities) 而构建的处理集组成。所得到的处理集，称为处理组合，是通过将其中一个因素的每个可能性与另一个因素的每个可能性组合而生成的。如果第一个因素有 $s$ 个可能性，第二个因素有 $r$ 个可能性，该组合产生 $sr$ 种处理组合。图 \@ref(fig:figure4-4) 给出了一个双向处理结构的例子，其中因素 $A$ 有三种可能性，因素 $B$ 有四种可能性，交叉产生 $12$ 种处理组合。

<div class="figure" style="text-align: center">
<img src="figure/figure 4.4.png" alt="双向处理结构，其中因素 $A$ 有三个水平，因素 $B$ 有四个水平，生成 $12$ 种处理组合" width="528" />
<p class="caption">(\#fig:figure4-4)双向处理结构，其中因素 $A$ 有三个水平，因素 $B$ 有四个水平，生成 $12$ 种处理组合</p>
</div>

3) **析因安排处理结构** (Factorial arrangement treatment structure). 析因安排处理结构由通过组合两个或多个因素的水平而构建的处理组合的集合组成。双向处理结构是双向析因安排。三向到 $n$ 向处理结构也是析因安排。$n$ 向处理结构是通过组合 $n$ 个因素的可能性而生成的，其中，因素分别具有 $s_1,s_2,\dots,s_n$ 个可能性，从而生成 $s_1\times s_2\times\dots\times s_n$ 种处理组合。析因安排处理结构的例子散布在整本书中。
4) **部分析因安排处理结构** (Fractional factorial arrangement treatment structure). 部分析因安排处理结构仅由析因安排处理结构中可能的处理组合的一部分或一小部分组成。有许多系统的技术来选择适当的分数 (fraction)，其中大多数取决于实验者对处理结构中各种类型因素之间交互作用的假设。一个实验可能涉及 $8$ 个具有两水平的不同因素，因此总共有 $2^8$ 个处理组合。研究人员可能想寻找重要的主效应和两个因素的交互作用。在这种情况下，研究中可以使用四分之一分数 (one-fourth fraction) 或64个处理组合。这样的设计通常表示为 $2^{8-4}$ 部分析因安排，更多细节参见 Milliken and Johnson (1989). 拉丁方排列处理结构涉及具有 $n$ 行处理，$n$ 列处理和 $n$ 个单元格 (cells) 处理的三向析因排列。拉丁方排列由 $n^3$ 种可能的处理组合中的 $n^2$ 种组成，或者是$n^3$ 种可能的处理组合的 $1/n =(n^2/n^3)$ 部分。部分拉丁方安排的一个可能用途是，可以假设这三个因素之间没有双向或三向交互作用。
5) **最优设计处理结构** (Optimal design treatment structures)。对于许多实验情况，处理结构中的所有因素都是定量的，研究目的是收集数据以建立特定的线性或非线性模型。使用其中一个设计准则 (St. John and Draper, 1975) 选择的处理组合形成了最优设计。这组处理组合被称为最优设计处理结构。
6) **具有一个或多个对照的部分安排** (Factorial arrangement with one or more controls). 用于满足实验目的的期望处理结构可以包括将多种处理结构相结合。例如，实验的处理结构可以包括将 $c$ 个对照的单向处理结构与双向析因安排处理结构相结合。图 \@ref(fig:figure4-5) 包含一种这样的处理结构，其中因素安排由两个水平的因素 $A$ 和三个水平的因素 $B$ 以及三个对照组成。

<div class="figure" style="text-align: center">
<img src="figure/figure 4.5.png" alt="将单向和双向处理结构组合，形成九种处理组合的处理结构用于实验" width="574" />
<p class="caption">(\#fig:figure4-5)将单向和双向处理结构组合，形成九种处理组合的处理结构用于实验</p>
</div>
出于分析目的，所有上述处理结构始终可以被视为一种单向处理结构。特别是，当处理结构是两个或多个处理结构的复杂组合时（如图 \@ref(fig:figure4-5) 所示），在分析结果数据时通常最好将一组处理视为单向处理结构。

裂区和重复测量设计结构由不完全区组设计结构和涉及两个或多个因素的析因安排处理结构组成。实际上，裂区和重复测量设计的设计结构和处理结构的组合产生不同尺寸的实验单元，这是一个必须解决的问题，以便获得适当的分析。这样的设计将在第 \@ref(chap5) 章讨论。

## 不同设计实验的示例 {#sec4-3}

有大量关于用于多种类型实验的各种类型设计的已发表信息，例如，
Cochra and Cox (1957), Davies (1954), Federer  (1955), Hicks (1993), John (1971), Kirk (1968), Cornell (1990), Anderson and McLean (1974), Box et al. (1978), Cobb (1997), Kempthorne (1952), Laundsby and Weese (1993), Lentner and Bishop (1986), Meed (1988), Montgomery (1991) and Winer (1971). 本节包含几个示例，以演示 \@ref(sec4-2) 节中描述的设计结构和处理结构。希望这一讨论将帮助读者将这些概念应用到他们自己的实验中。在大多数情况下，通过指定设计结构的类型和处理结构的类型来命名所得到的设计实验。例如，设计的实验可以由随机完全区组设计结构中的双向处理结构组成。这种描述设计实验的方法不同于文献中通常使用的方法，但作者认为使用设计结构和处理结构是识别设计实验的最好方法。此外，这一描述还有助于构建适当的模型并进行适当的分析。针对每种实验情况，给出了相应的设计结构和处理结构，并给出了相应的模型，以及带有相应的变异源及相应自由度的的方差分析表。本节不给出平方和的计算公式，但在其他章节中提供了一些有计算的例子。本章的重点是确定误差平方和的计算方法，并确立相应的自由度。

### 示例 4.1: 饮食 {#sec4-3-1}

一位营养学家想要研究健康饮食对体重减轻的影响。本实验的处理结构是单向分类，涉及单个因素，称为饮食，具有五个水平或五个处理。许多不同的设计结构可以用来评估饮食。如果有 $20$ 个同质的人，那么可以使用完全随机设计结构，其中每种饮食被随机分配给四个人。完全随机设计结构中的单向处理结构的一个模型是

\begin{equation}
y_{ij}=\mu_i+\varepsilon_{ij}\quad i=1,2,\ldots,t,j=1,2,\ldots,n_i
(#eq:4-3)
\end{equation}

其中 $\mu_i$ 表示第 $i$ 种处理（饮食）的均值，$\varepsilon_{ij}$ 表示随机误差。模型 \@ref(eq:4-3) 的方差分析表见表 \@ref(tab:table4-1)，其中假设了 $\varepsilon_{ij}\sim \text{i.i.d.}\,N(0,\sigma^2)$. 表 \@ref(tab:table4-1) 包含公式 \@ref(eq:4-3) 中模型的不同变异源及其各自的自由度。实验误差的 $15$ 个自由度是从经过相同处理的人之间的变异中得到的。有四个人接受了饮食 1，提供了经过相同处理的四个人，这四个人体重减轻值之间的差异为误差提供了三个自由度。每种饮食有四个人经过相同处理。因此，五种饮食中的每种饮食的数据都有三个可供误差的自由度。如果方差相等，那么这五组三个自由度可以合并为一个涉及 $15$ 个自由度的误差项。第 \@ref(chap2) 章中的方法可用于评估等方差假设的合理性。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-1)完全随机设计结构中单向处理结构的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.1.png">
</td>
  </tr>
</tbody>
</table>

假设研究者在给定时间段内只能招募并采访五个人。让时间段成为一个区组因素，并在每个时间段内将五种饮食随机分配给一个人。该设计是随机完全区组设计结构中的单向处理结构，具有四个尺寸为五的区组。可用于描述此设计中的数据的模型是

\begin{equation}
y_{ij}=\mu_i+b_j+\varepsilon_{ij}\quad i=1,2,3,4,5,j=1,2,3,4
(#eq:4-4)
\end{equation}

其中 $\mu_i$ 表示第 $i$ 个处理（饮食）的均值，$b_j$ 表示第 $j$ 个区组的效应，$\varepsilon_{ij}$ 表示与分配给第 $j$ 个区组中第 $i$ 个处理的人相关的随机误差。模型 \@ref(eq:4-4) 的方差分析表如表 \@ref(tab:table4-2) 所示。误差自由度是通过使用每个区组内处理的四个正交对比来计算的，例如，

$$\begin{array}{ll}q_{1j}=y_{1j}-y_{2j}&j=1,2,3,4\\q_{2j}=y_{1j}+y_{2j}-2y_{3j}&j=1,2,3,4\\q_{3j}=y_{1j}+y_{2j}+y_{3j}-3y_{4j}&j=1,2,3,4\\q_{4j}=y_{1j}+y_{2j}+y_{3j}+y_{4j}-4y_{5j}&j=1,2,3,4\end{array}$$

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-2)随机完全区组设计结构中单向处理结构的方差分析表，每个区组中每种处理有一次重复</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.2.png">
</td>
  </tr>
</tbody>
</table>

相同 $i$ 值的 $q_{ij}$ 值都具有相同的均值，表明它们都经相同处理。因此，每组四个 $q_{ij}$ 值的方差为误差提供了三个自由度。因此，可以将四组三个自由度合并（如果方差相等）为具有 $12$ 个自由度的误差平方和。然而，$q_{ij}$ 的方差并不相等，必须在合并 (pooling) 之前重新缩放。对于 $j = 1, 2, 3, 4$，方差为 $\text{Var}(q_{1j}) = 2\sigma^2,\text{Var}(q_{2j})= 6\sigma^2,\text{Var}(q_{3j})= 12\sigma^2$ 以及 $\text{Var}(q_{4j})= 20\sigma^2$，因此方差为在合并之前需要将 $q_{ij}$ 除以 $\sigma^2$ 的系数。实际上，这些合并的平方和提供了区组与饮食交互作用的平方和。因此，获得相似处理的方差的第二种方法是使用区组与处理的交互作用 (block-by-treatmen interaction) 来获得随机完全区组设计结构的实验误差估计。

接下来，假设没有 $20$ 个同质的人，但有 $10$ 个同质的男性和 $10$ 个同质的女性。一种策略是使用人的性别作为区组因素，其中有两个尺寸为 $10$ 的区组。可以使用随机完全区组设计结构，其中每种饮食将随机分配给两名男性和两名女性，因此每个区组中每种饮食有两次重复。随机完全区组设计结构中单向处理结构的模型为

\begin{equation}
y_{ijk}=\mu_{i}+b_{j}+\boldsymbol{\varepsilon}_{ijk}\quad i=1,2,3,4,5,j=1,2,k=1,2
(#eq:4-5)
\end{equation}

其中 $\mu_i$ 表示第 $i$ 个处理（饮食）效应的均值，$b_j$ 表示第 $j$ 个区组的效应，$\varepsilon_{ijk}$ 表示与分配给第 $j$ 个区组中第 $i$ 个处理第 $k$ 个人相关的随机误差。模型 \@ref(eq:4-5) 的方差分析表如表 \@ref(tab:table4-3) 所示。有一个与设计结构或人的性别相关的自由度已从表 \@ref(tab:table4-1) 的误差项中删除。随机完全区组设计结构中的单向处理结构的误差项，其中每个区组中每个处理有一个观测结果，是根据区组与处理的交互作用来计算的。这种设计包括区组和每个区组中每种处理的多次观测。区组与处理的交互作用的自由度等于 $(2-1)×(5-1)=4$. 每个区组内每种饮食的两个观测结果的变异性提供了一个自由度。因此，在每个区组内对经相同处理的人进行比较时，有五个误差自由度，或者在两个区组内总共有 $10$ 个自由度。将区组与处理的交互作用与区组内变异合并起来共提供了 $14$ 个误差自由度。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-3)随机完全区组设计结构中单向处理结构的方差分析表，每个区组中每种处理有两次重复</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.3.png">
</td>
  </tr>
</tbody>
</table>

在大多数情况下，人的性别不是区组因素好的选择，因为处理（这个例子为饮食）可能会与人的性别发生交互 (interact)。如果选择用于构建区组的一个或多个因素可能与处理发生交互，则该因素需要作为分层因素或在具有可能不等斜率的协方差分析中包括在处理结构中，以进行适当的分析 (Milliken and Johnson, 2001). 在这个例子中，一个人的性别必须与五种饮食相结合，以形成双向因素安排或双向处理结构。这种双向处理结构由 $10$ 种处理组合组成，通过将人的两个性别水平与五个饮食水平相结合而产生。通过将人的性别从设计结构切换到处理结构，得到的设计结构是一个完全随机的设计，每个处理组合有两次重复。随机化方案是将每种饮食随机分配给两名男性和两名女性。这是与模型 \@ref(eq:4-5) 相同的随机化方案，但具有不同的处理和设计结构。完全随机化设计结构中的双向处理结构的一个模型是均值模型：

\begin{equation}
y_{ijk}=\mu_{ij}+\varepsilon_{ijk}\quad i=1,2,\ldots,5,j=1,2,k=1,2
(#eq:4-6)
\end{equation}

其中 $\mu_{ij}$ 表示第 $ij$ 种处理组合（人的性别与饮食）的平均值，$\varepsilon_{ijk}$ 表示随机误差。有时，均值 $\mu_{ij}$ 表示为效应模型：

$$\mu_{ij}=\mu+\tau_i+\beta_j+\gamma_{ij}$$

[^overallmean]: 这里的 $\mu$ 还不足以称为“总体平均值”，严谨地说是“截距”，如果对模型其他参数施加约束才可称为“总体平均值”。更多的讨论详见 Stroup (2012)，例如前言的 xxi 页以及 P.54 的 Note.

其中 $\mu$ 是总体平均值[^overallmean](overall mean)，$\tau_i$ 是第 $i$ 种饮食的效应，$\beta_j$ 是人的第 $j$ 种性别的效应，$\gamma_{ij}$ 是交互作用效应。表 \@ref(tab:table4-4) 给出了模型 \@ref(eq:4-6) 这两种 $\mu_{ij}$ 表达式的方差分析表。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-4)完全随机设计结构中双向处理结构方差分析表，包含均值模型和效应模型</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.4.png">
</td>
  </tr>
</tbody>
</table>

接下来，假设饮食结构由一个对照饮食和四个饮食组成，这四个饮食由两个蛋白质水平和两个碳水化合物水平的四个组合组成，如图 \@ref(fig:figure4-6) 所示。饮食处理结构是一种双向析因安排，当与人的性别交叉时，对照组会产生一种三向处理结构（蛋白质×碳水化合物×性别），其中有两个对照，一个对照用于男性，一个对照用于女性。设计结构是完全随机的，每个处理组合分配给两个人。可以用来描述这些数据的模型是

\begin{equation}
y_{ijk}=\mu_{ij}+\varepsilon_{ijk}\quad i=0,1,2,3,4,j=1,2,k=1,2
(#eq:4-7)
\end{equation}

其中 $\mu_{01}$ 和 $\mu_{02}$ 表示对照组的平均值， $\mu_{ij},i = 1, 2, 3, 4, j = 1, 2$ 表示饮食与性别处理组合的平均值。表 \@ref(tab:table4-5) 为模型 \@ref(eq:4-7) 的方差分析表。其中 "Control vs 2^2^" 表示对照饮食与四种蛋白质和碳水化合物处理组合的平均值之间的比较。完整的分析最有可能使用饮食与性别的双向处理结构，并具有适当的处理均值的对比，以提供涉及蛋白质、碳水化合物、人的性别和对照的比较。这种对比的构建将在第 \@ref(chap6) 章及以后进行讨论。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-5)完全随机设计结构中三向析因安排与两个对照相结合的处理结构的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.5.png">
</td>
  </tr>
</tbody>
</table>


<div class="figure" style="text-align: center">
<img src="figure/figure 4.6.png" alt="具有双向析因安排和一个对照的处理结构的饮食结构" width="557" />
<p class="caption">(\#fig:figure4-6)具有双向析因安排和一个对照的处理结构的饮食结构</p>
</div>

### 示例 4.2: 房屋油漆 {#sec4-3-2}

一家油漆公司想要比较四种白色房屋油漆承受环境条件的能力。实验中准备了四栋方形房屋，每栋房屋的一面朝北，因此房屋可以作为区组因素。房屋的每一面都可能暴露在不同类型的天气下，因此房屋的侧面（这里用北、南、东、西方向表示）也可以用作区组因素。由于处理数量（四种油漆）与两个区组因素的水平数相同，因此可以使用拉丁方设计结构来研究油漆。在这里，油漆可以分配到房屋的侧面，其中每种油漆可以在每个房屋上发生且仅发生一次，并且在每个方向上发生且仅发生一次。存在三种基本的拉丁方安排（参见 Cochran and Cox, 1957）。随机化过程是从三种可能的安排中随机选择一种，然后将房屋随机分配给行，将方向随机分配给列，并将油漆类型随机分配给方块中的字母。表 \@ref(tab:table4-6) 显示了一种将油漆分配给房屋和方向的安排。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-6)将一组处理从单向处理结构分配到拉丁方设计结构</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.6.png">
</td>
  </tr>
</tbody>
</table>

本实验的设计是一种单向处理结构中的拉丁方设计结构。可以用来描述该实验数据的模型是

\begin{equation}
y_{ijk}=\mu_i+h_j+d_k+\boldsymbol{\varepsilon}_{ijk}\quad i=1,2,3,4,j=1,2,3,4,k=1,2,3,4
(#eq:4-8)
\end{equation}

其中 $\mu_i$ 表示第 $i$ 种油漆的平均耐磨性得分，$h_j$ 表示第 $j$ 个房屋的效应，$d_k$ 表示第 $k$ 个方向的效应，$\varepsilon_{ijk}$ 表示实验单元误差。表 \@ref(tab:table4-7) 给出了模型 \@ref(eq:4-8) 的方差分析表。拉丁方设计结构的误差由对比组成，这些对比测量了油漆、房屋与方向的交互作用。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-7)拉丁方处理结构中单向处理结构的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.7.png">
</td>
  </tr>
</tbody>
</table>

接下来，假设油漆具有由 1）基础油漆，2）基础油漆 + 添加剂 I，3）基础油漆 + 添加剂 II，以及 4）基础油漆 + 添加剂 I +添加剂 II。油漆的这种结构提供了双向处理结构，其中一个因素是两个水平[无和添加一些]的添加剂 I，第二个因素是两个水平[无和添加一些]的添加剂 II。由此产生的四种处理组合见表 \@ref(tab:table4-8).

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-8)房屋油漆的双向处理结构示例</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.8.png">
</td>
  </tr>
</tbody>
</table>

拉丁方设计结构中的双向处理结构的一个模型是

\begin{equation}
y_{ijkm}=\mu+\gamma_i+\beta_j+(\gamma\beta)_{ij}+h_k+d_m+\varepsilon_{ijkm}\quad i=1,2,\quad j=1,2,\quad k=1,2,3,4,\quad m=1,2,3,4
(#eq:4-9)
\end{equation}

其中 $\gamma_i$ 表示添加剂 I 的效应，$b_j$ 表示添加剂 II 的效应，$(gb)_{ij}$ 表示两种添加剂之间的交互作用。表 \@ref(tab:table4-9) 给出了模型 \@ref(eq:4-9) 的方差分析表。分析模型 \@ref(eq:4-8) 和 \@ref(eq:4-9) 的唯一区别在于，在模型 \@ref(eq:4-9) 中，油漆具有一种结构，用于将油漆效应划分为添加剂 I、添加剂 II 以及添加剂 I 和 II 的交互作用所产生的效应。即使处理结构的分析发生了变化，与设计结构相对应的分析部分仍然不受影响。

最后，假设有八栋房屋可用，以便可以使用两种不同的拉丁方设计结构进行实验。表 \@ref(tab:table4-10) 显示了将油漆分配给房屋-方向组合的一种可能分配。如果油漆具有表 \@ref(tab:table4-8) 中的双向处理结构，则模型由下式给出

\begin{align}
y_{ijkmn}&=\mu+\gamma_i+\beta_j+(\gamma\beta)_{ij}+s_k+h_{n(k)}+d_m+\varepsilon_{ijkmn}\\
i&=1,2,\,j=1,2,\,k=1,2,\,m=1,2,3,4,\,n=1,2,3,4
(#eq:4-10)
\end{align}

其中 $s_k$ 表示方块 $k$ 的效应，$h_{n(k)}$ 表示方块 $k$ 中房屋 $n$ 的效应。模型 \@ref(eq:4-10) 的方差分析表如表 \@ref(tab:table4-11) 所示。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-9)拉丁方处理结构中双向处理结构的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.9.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;border-bottom: 0;">
<caption>(\#tab:table4-10)在重复拉丁方设计结构中显示单向处理结构的安排</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.10.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> Square 1 和 Square 2 为用于重复拉丁方设计的两个房屋分组</td></tr></tfoot>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;border-bottom: 0;">
<caption>(\#tab:table4-11)重复拉丁方设计结构中双向处理结构的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.11.png">
</td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<span style="font-style: italic;">注释: </span> <sup></sup> Square 为用于重复拉丁方设计的房屋分组</td></tr></tfoot>
</table>

### 示例 4.3: 钢板 {#sec4-3-3}

当需要在两个方向上分区 (block) 时，拉丁方设计结构非常有用，但是实验者使用的每个拉丁方安排都不是拉丁方设计结构。此示例用于演示使用拉丁方安排处理结构的后果。将两种油漆添加剂组合起来，对钢板涂漆。该实验的目的是研究油漆组合保护钢材免受高温影响的能力。每种油漆添加剂有五个水平和五个温度来检查保护能力。该实验适用于拉丁方，其中添加剂 I 的水平分配给行，添加剂 II 的水平分配给列，温度水平分配给方块内的单元格。这种排列产生 $25$ 种处理组合。实验单元为 $25$ 块厚 $0.2cm$、面积为 $1 m^2$ 的钢板。随机化过程是将 $25$ 种处理组合中的一种随机分配给 $25$ 张钢板中的每一张。在这种情况下，处理结构是 $5^3$ 析因安排的一部分（因为它由添加剂 I × 添加剂 II × 温度的 $125$ 种可能的处理组合中的 $25$ 种组成），称为拉丁方安排 (Latin square arrangement) 或拉丁方处理结构 (Latin square treatment structure). 该设计结构是完全随机设计，因为处理组合完全随机地分配给钢板。由于这是部分析因，因此每个主效应与其他两个因素的两因素交互作用和三因素交互作用**部分混叠** (partially aliased)（参见 Cochran and Cox, 1957, p. 245）。为了正确分析该实验设计，必须对模型中的参数做出一些假设。通常的假设是不存在双向交互作用，也不存在三向交互作用。然而，人们应该非常小心，如果没有一些先验信息（可以来自其他实验、现有文献等）表明交互作用实际上可以忽略不计，就不要做出这样的假设。表 \@ref(tab:table4-12) 给出了一种这样的拉丁方安排，完全随机设计结构中拉丁方处理结构的模型为

\begin{align}
y_{ijk}&=\mu+\mathrm{AI}_i+\mathrm{AII}_j+T_k+\varepsilon_{ijk},\quad(i,j,k)\in\mathrm{Index}
(#eq:4-11)
\end{align}

其中 $\mathrm{AI}_i$ 表示的添加剂 I 第 $i$ 个水平的效应，$\mathrm{AII}_j$ 表示添加剂 II 第 $j$ 个水平的效应，$T_k$ 表示温度第 $k$ 个水平的效应，$\mathrm{Index}$ 表示由实验中观测到的 $25$ 种处理组合组成的指数集。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-12)拉丁方安排处理结构，其中 $T_i$ 表示温度第 $i$ 个水平</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.12.png">
</td>
  </tr>
</tbody>
</table>

如果忽略温度水平的应用，则得到的数据表是完全随机设计结构中双向处理结构的数据表，如表 \@ref(tab:table4-13) 所示。该双向处理结构的方差分析表见表 \@ref(tab:table4-14). 该设计包括每个处理组合的一个观测结果，因此没有可用的自由度来估计误差，因为没有经过相同处理的钢板。AI 水平和 AII 水平之间的交互作用与 $16$ 个自由度有关。现在，当温度包括在结构中时，由于检验温度均值相等而产生的平方和是 AI 与 AII 交互作用平方和的一部分，如表 \@ref(tab:table4-15) 所示。这意味着 AI 与 AII 交互作用的一部分与温度效应相同，或者说，与温度相关的四个自由度与此交互作用相关的四种自由度**混叠** (aliased). 同样，与 AI 效应相关联的四个自由度与 AII 和温度交互作用的四个自由度混叠，与 AII 效应相关的四个自由度与 AI 和温度交互作用的四个自由度混叠。表 \@ref(tab:table4-16) 中给出了使用模型 \@ref(eq:4-11) 的拉丁方处理结构的方差分析表。使用术语“残差” (residual) 而不是“误差” (error)，因为相应的平方和涉及误差加上可能不为零的任何交互作用效应。如果零交互作用的假设不正确，那么残差均方将过大，由此产生的 $F$ 检验将过小。因此，如果实验中存在交互作用，就无法发现，并且任何其他可检测的处理效应都可能被掩盖。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-13)拉丁方安排中添加剂 I 水平和添加剂 II 水平的双向处理结构，忽略了温度水平</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.13.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-14)拉丁方安排双向处理结构部分的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.14.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-15)拉丁方安排双向处理结构部分的方差分析表，从两种添加剂的交互作用中划分出温度的变异</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.15.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-16)拉丁方安排处理结构方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.16.png">
</td>
  </tr>
</tbody>
</table>

### 示例 4.4: 氮和钾的水平 {#sec4-3-4}

模型和结果分析由三个基本组成部分组成：i) 处理结构，ii) 设计结构，以及 iii) 误差结构。此示例演示了如何使用三个基本组件来构建模型。一位植物育种家想要研究三种水平的氮 (N) 和四种水平的钾 (K) 的组合对其新玉米品种产量的影响。他的处理结构是双向析因安排，有 $12$ 种处理组合（$3$ 个水平的氮 × $4$ 个水平的钾）。他有三块土地用于进行实验，他将这些土地作为区组。每个区组被分为 $12$ 个部分，称为“区” (plots). 每种处理组合被随机分配到每个区组内的一个区。因此，设计结构是随机完全区组设计，因为每种处理组合在每个区组中出现一次。实验的总体设计称为随机完全区组设计结构中的双向处理结构（随机完全区组设计中的区组被一些作者称为重复 (replications)；然而，我们更喜欢将它们称为区组而不是重复，以便将它们与完全随机设计结构中的重复 (replications) 区分开来。示例 [4.5](#sec4-3-5) 中的讨论描述了区组和重复之间的重要区别）。此示例的模型是

\begin{equation}
y_{ijk}=\mu_{ij}+b_{k}+\boldsymbol{\varepsilon}_{ijk}\quad i=1,2,3,j=1,2,3,4,k=1,2,3
(#eq:4-12)
\end{equation}

其中 $\mu_{ij}$ 是氮的第 $i$ 个水平与钾的第 $j$ 个水平的平均值，$b_k$ 是第 $k$ 个区组的效应，$\varepsilon_{ijk}$ 表示与每个区组内的区相关的随机误差。一般来说，通过将三种结构中的每一种结构的模型相加来构建一般模型

\begin{equation}
Y= \mathrm{treatment structure} + \mathrm{design structure} + \mathrm{error structure(s)}
(#eq:4-13)
\end{equation}

同样，相应的方差分析表也有三个部分。表 \@ref(tab:table4-17) 给出了模型 \@ref(eq:4-13) 一般的方差分析表。表 \@ref(tab:table4-18) 给出了模型 \@ref(eq:4-12) 的方差分析表。一般来说，必须考虑到可能存在一个以上的误差项。例如，裂区和重复测量模型有一个以上的误差项（见第 \@ref(chap5) 章）。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-17)一般模型的方差表分析</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.17.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-18)随机完全区组设计结构中双向处理结构的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.18.png">
</td>
  </tr>
</tbody>
</table>

### 示例 4.5: 区组和重复 {#sec4-3-5}

在许多关于实验设计的教科书中，要么没有区分区组和重复，要么至少对这种区别存在混淆 (confusion)。包含这组示例是为了演示两个概念之间的差异，并说明处理结构与设计结构的组合可以产生实验的多个总体设计。假设研究人员想使用一个设计结构来研究单向处理结构中四种处理的效应，该设计结构有六个尺寸为二的区组（每个区组只有两个同质的实验单元）。在这种情况下，所需的设计结构是不完全区组设计。如果有足够的区组，使得每个处理对都可以在一个区组中同时进行相同次数，那么就可以使用均衡不完全区组设计结构 (Cochran and Cox, 1957). 例如，如表 \@ref(tab:table4-19) 所示，可以将四种处理方式分配给区组。在这种情况下，设计结构中有六个区组，每个处理重复三次。此示例用于指出**区组和重复的概念是不同的，并强调区组和重复应始终保持分离**。**区组和重复仅在随机完全区组设计结构的情况下是等价的**，其中每个处理在每个区组内观测一次且仅观测一次。不完全区组设计结构的随机化过程包括为区组编号分配区组，然后将分配给该区组的两个处理随机分配给该区组内的两个实验单元。对于表 \@ref(tab:table4-19) 中的示例，设计结构与六个区组相关联（见表 \@ref(tab:table4-19) ），而不与由于分配过程恰好发生的三次重复相关联。表 \@ref(tab:table4-19) 安排的模型为

\begin{equation}
y_{ij}=\mu_i+b_j+\boldsymbol{\varepsilon}_{ij},\quad\mathrm{for~}(i,j)\in\mathrm{Index}
(#eq:4-14)
\end{equation}

[^connectd]: 见 \@ref(sec6-3-2) 节。

[^connectedblocktreatmentdf]: 原文：The degrees of freedom associated with this connected block–treatment arrangement are computed from the degrees of freedom associated with the block by treatment interaction as if all combinations were observed minus the number of empty cells.

其中 $\mathrm{Index}=\{(A, 1), (B, 1), (A, 2), (C, 2), (A, 3), (D, 3), (B, 4), (C, 4), (B, 5), (D, 5), (C, 6), (D, 6)\}$，且对子 $(i,j)$ 只能取如索引集所示观测到的处理 × 区组组合的值。表 \@ref(tab:table4-20) 给出了模型 \@ref(eq:4-14) 的方差分析表。与这种**连通的[^connectd]区组-处理** (connected block-treatment) 安排相关的自由度是如此计算的：与区组处理交互作用相关的自由度（就好像观测到了所有组合）减去空单元格的数量[^connectedblocktreatmentdf]。表 \@ref(tab:table4-21) 显示了观测到的区组-处理组合（用 "X" 表示）。有六个区组和四个处理，因此如果所有组合都存在，则区组处理交互作用将基于 $(6-1)(4-1)=15$ 个自由度。有 $12$ 个单元格缺失，因此与误差项相关的自由度为 $15-12=3$.

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-19)将四种处理分配到两个实验单元的六个区组的第一种分配：提供均衡不完全区组设计结构</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.19.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-20)表 \@ref(tab:table4-19) 中的均衡不完全区组设计结构的方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.20.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-21)区组与处理的双向结构，其中 $X$ 表示观测到的组合，被填充和空缺的各有 $12$ 个单元格</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.21.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table4-22) 包含第二种分配模式，用于将四种处理分配给六个尺寸为二的区组。处理 $A$ 出现在所有区组中并重复六次。处理 $B,C,D$ 出现在两个区组中且每个区组提供两次重复。模型（4.14）可用于描述此数据，其中索引集为

$$\mathrm{Index}=\{(A, 1), (B, 1), (A, 2), (C, 2), (A, 3), (D, 3), (A, 4), (B, 4), (A, 5), (C, 5), (A, 6), (D, 6)\}$$

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-22)将四种处理分配给两个实验单元的六个区组，第二种分配：提供不均衡不完全区组设计结构</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.22.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table4-20) 中的方差分析表也适用于表 \@ref(tab:table4-22) 中的安排。如果希望将所有四种处理相互比较，即检验 $\mu_A = \mu_B = \mu_C = \mu_D$，则表 \@ref(tab:table4-19) 中的安排是最优的（使用 D-最优准则；St. John and Draper, 1975）。如果希望将处理 $A$ 与其他三种处理进行比较，即检验 $\mu_A = \mu_B,\mu_A = \mu_C,\mu_A = \mu_D$，则表 \@ref(tab:table4-22) 中的安排是最优的。
 
表 \@ref(tab:table4-23) 给出了将四种处理分配到六个区组的第三种方法。将处理分为两组，其中 $A$ 和 $B$ 一起出现在三个区组中，处理 $C$ 和 $D$ 一起出现在三个区组中。处理 $A$ 不与处理 $C$ 或 $D$ 的同时出现在一个区组中。处理和区组之间的结构**不连通** (not connected)，因为处理分为这两组。定义一个新的变量“组” (group) 来表示两组处理。两组之间的比较是区组 $1,3,5$ 的平均值与区组 $2,4,6$ 的平均值之间的比较，或者说是区组间比较 (between-block comparisons). 处理 $A$ 与处理 $B$ 之间的比较以及处理 $C$ 与处理 $D$ 之间的比较属于区组内比较 (within-block comparisons). 对于处理 $A$ 和处理 $B$，处理与区组交互作用的自由度为二。同样，对于处理 $C$ 和处理 $D$，处理与区组交互作用的自由度也为二。误差平方和是将这两个处理与区组交互作用的平方和（及其自由度）合并得到的。表 \@ref(tab:table4-23) 中的安排的方差分析表见表 \@ref(tab:table4-24)。区组的平方和分为组的平方和和组内嵌套的区组的平方和。各组的平方和作为误差以检验假设 $\mu_A + \mu_B = \mu_C + \mu_D$），该假设具有由处理产生的一个自由度。两组处理与区组混杂 (confounded) 在一起；也就是说，如果两组均值之间存在差异，则不知道是由于处理集之间的差异还是由于两组区组之间的差异造成的。混杂的概念与混叠 (aliasing) 的概念类似，但混叠涉及两个（或多个）项无法区分，其中两个项都来自处理结构，而混杂涉及两个项无法区分，其中一项来自处理结构，另一项来自设计结构。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-23)将四种处理分配给两个实验单元的六个区组的第三种分配：提供不连通不完全区组设计结构</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.23.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-24)表 \@ref(tab:table4-23) 中的不完全区组设计方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.24.png">
</td>
  </tr>
</tbody>
</table>

表 \@ref(tab:table4-19)、\@ref(tab:table4-22) 和 \@ref(tab:table4-23) 中三种设计之间的差异在于如何将处理分配给区组。请记住，所有三种设计的设计结构和处理结构都是相同的，**因此设计结构和处理结构并不描述实验的总体设计**。**还必须指定将处理结构中的处理随机分配到设计结构中的实验单元的方法**。

### 示例 4.6：行区组和列区组 {#sec4-3-6}

区组可以通过多种方式出现，拉丁方设计结构是一种同时存在行区组和列区组的设计。拉丁方设计结构有多种变化，其中行比处理或多或少（或者列比处理或多或少）（Cochran and Cox, 1957）。本示例由按行和列划分的实验单元组成，其中每行和列的交集包含多个实验单元。处理结构是双向的，因素 $A$ 有两个水平，因素 $B$ 有两个水平，从而产生四种处理组合。共有 $24$ 个实验单元，排列成三列两行，每行列组合包含 $4$ 个实验单元。将四种处理组合随机分配给每个行列组内的四个实验单元。表 \@ref(tab:table4-25) 显示了实验单元的处理组合（非随机）的分配。该设计结构本质上由六个尺寸为四的区组组成，但评估行区组之间和列区组之间的变异性可能很有趣，因此可以构建一个模型以包含这些因素：

\begin{align}
y_{ijkm}&=\mu+\alpha_i+\beta_j+(\alpha\beta)_{ij}+r_k+c_m+(rc)_{km}+\varepsilon_{ijkm}\\
i&=1,2,j=1,2,k=1,2,m=1,2,3
(#eq:4-15)
\end{align}

其中 \mu + \alpha_i + \beta_j + (\alpha\beta)_{ij} 表示因素 A 和 B 的主效应和交互作用，r_k 表示行区组效应，c_m 表示列区组效应，(rc)_{km} 表示行区组和列区组之间的交互作用，\varepsilon_{ijkm} 表示实验单元误差。模型 \@ref(eq:4-15) 对应的方差分析表见表 \@ref(tab:table4-26). 分析可能只涉及六个区组，并且可以合并由于行区组、列区组及其交互作用而产生的平方和，从而为设计结构提供五个自由度。通过设计结构与处理结构的交互作用来获得误差平方和；即通过合并 A × 行区组、A × 列区组、A × 行区组 × 列区组、B × 行区组、B × 列区组、B × 行区组 × 列区组、A × B × 行区组、A × B × 列区组，以及 A × B × 行区组 × 列区组的平方和。

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-25)在设计结构中使用行区组和列区组，在处理结构中使用双向安排的设计（非随机形式）</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.25.png">
</td>
  </tr>
</tbody>
</table>

<table class="table table-condensed" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:table4-26)双向处理结构的行-列设计结构方差分析表</caption>
 <thead>
  <tr>
   <th style="text-align:center;color: white !important;background-color: white !important;font-size: 0px;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;">  <img src="table/table%204.26.png">
</td>
  </tr>
</tbody>
</table>

通过结合各种设计结构和处理结构，还有许多其他方法可以构建设计实验。希望以上例子将使实验者能够构建所需的设计实验，构建适当的模型，并进行相应的分析。

## 结束语 {#sec4-4}

本章介绍了实验者在设计和分析实验时可以使用的概念和方法。还介绍了设计良好的实验的基本概念。所有设计的实验都包括两个基本特征：处理结构和设计结构。这些概念通常不会在许多其他统计分析书籍中使用。了解设计实验的这两个特征之间的差异将有助于数据分析师为他们的实验选择适当的分析。在必须不与处理结构中的因素发生交互的情况下讨论了区组因素的选择。最后，必须能够识别何时有真正的重复，何时只有子样本。

## 练习 {#sec4-5}

## 参考文献 {-}

Kempthorne, O. (1952). Design and Analysis of Experiments. John Wiley and Sons, New York, NY.

Milliken, G. A. and Johnson, D. E. (2002). Analysis of Messy Data, Vol 3: Analysis of Covariance. Chapman and Hall/CRC, Boca Raton, FL.

Cobb, G. W. (1997). Introduction to Design and Analysis of Experiments. Springer, New York, NY.

Cochran, W. G., and Cox, G. M. 1957. Experimental Design. 2nd ed. John Wiley and Sons, New York, NY.

Davies, O. L. 1954. Design and Analysis of Industrial Experiments. Oliver and Boyd, London, UK.

Federer, W. T. (1995). Experimental Design. Macmillan, New York, NY.

Hicks, Charles R. (1993). Fundamental Concepts in the Design of Experiments, 4th ed. W.B. Saunders Co., Philadelphia, PA.

John, P. W. M. (1971). Statistical Design and Analysis of Experiments. John Wiley and Sons, New York, NY.

Kirk, R. E. (1968). Experimental Design: Procedures for the Behavioral Sciences. Brooks/Cole, Belmont, CA.

Cornell, J. A. (1990). Experiments with Mixtures: Designs, Models and Analysis of Mixture Data. John Wiley and Sons, New York, NY.

Anderson, V. L. and McLean, R. A. (1974). Design of Experiments: A Realistic Approach. Marcel Dekker, New York, NY.

Box, G. E. P., Hunter, W. G., and Hunter, J. S. (1978). Statistics for Experimenters. John Wiley and Sons, New York, NY.

Laundsby, R. G. and Weese, D. L. (1993). Straight Talk on Designing Experiments: An Introductory Design of Experiments Reference Handbook. Launsby Consulting, Colorado Springs, CO.

Lentner, M. and Bishop, T. (1986). Experimental Design and Analysis. Valley Book Company, Blacksburg, VA.

Meed, R. (1988). The Design of Experiments: Statistical Principles for Practical Application. Cambridge University Press, Cambridge, UK.

Montgomery, D. C. (1991). Design and Analysis of Experiments, 3rd edn. John Wiley and Sons, New York, NY.

Winer, B. J. (1971). Statistical Principles in Experimental Design, 2nd ed. McGraw-Hill, New York, NY.

St John, R. C. and Draper, N. R. (1975). D-Optimality for regression designs: A review. Technometrics, 17, 15.
