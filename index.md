--- 
title: "混乱数据分析：设计的实验"
author: 
  - Wang Zhen
date: "2024-02-22"
site: bookdown::bookdown_site
documentclass: elegantbook
bibliography: 
 - book.bib
colorlinks: true
link-citations: true
mathspec: true
graphics: true
geometry: margin=2cm
lot: true
lof: true
github-repo: wangzhen89/AMD
classoption: 
 - lang=cn
 - titlestyle=hang
 - chinesefont=ctexfont
description: "Analysis of Messy Data Volume 1: Designed Experiments的翻译"
---




## 说明 {-}

[混乱数据分析：设计的实验](https://github.com/wangzhen89/AMD) 为 [Analysis of Messy Data Volume 1: Designed Experiments, 2nd edn](https://doi.org/10.1201/EBK1584883340) 的翻译。

本书/翻译基于 [Elegant Bookdown Template](https://bookdown.org/xiangyun/elegantbookdown/)，进度还处于早期。尚未包括每章的练习。

> 本书探讨了多重比较程序、随机效应模型、混合模型、裂区实验和重复测量设计的各种技术。作者使用几个统计软件包简化了这些技术，并强调了设计结构和处理结构之间的区别。他们用例子介绍每一个主题，接着进行理论讨论，最后进行案例研究。本书向您展示如何有效地分析现实世界中的非标准数据集。

本书需要读者对（数理）统计学基本概念有一定的了解，并且最好有线性模型和实验设计的基础知识。在本书阅读过程中常常需要交替往返地进行查阅（或向前或向后）以获得完整的学习过程，出于此特殊结构，本书不适合作为教科书，更适合作为学习完线性模型和实验设计后的一本进阶学习资料或案头参考书。

**由于公式较多，本书可能需要加载数秒**。

## 本书框架 {-}

作为开篇，第 \@ref(chap1) 章讨论了最简单的情况：具有同质误差的完全随机设计结构中的单向处理结构。在这类实验中，研究人员希望考察某（一）种因素对响应变量的影响，并假定响应变量在该因素的不同水平中的方差相等。

第 \@ref(chap2) 章接续第 \@ref(chap1) 章的基本结构，但将同质误差改为异质误差，即该因素可能会改变响应变量的方差。本章着重讨论了有关方差齐性检验的一系列技术。

第 \@ref(chap3) 章重点关注实验中常常出现的一个问题：如何控制多重比较中的各种错误率？本章回顾了多重比较的一系列技术并给出了建议。

第 \@ref(chap4) 章对实验设计基本概念进行了详尽的描述。并强调了重复和复制的区别，提醒研究者注意勿将子样本或重复测量误认为重复。此外本章还讨论了其他实验设计教科书较少着墨甚至忽略的内容：将设计实验的结构划分为处理结构和设计结构，并列举了这两种结构常用的类型。最后结合两种结构的划分讨论了一些示例。

第 \@ref(chap5) 章讨论了多水平设计，包括裂区设计、裂条设计、重复测量设计及这些设计的组合。

第 \@ref(chap29) 章讨论了交叉设计。
