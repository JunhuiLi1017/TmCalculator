# 末端起始项与拟合基线：issue #8 的原文证据

核查日期：2026-09-29。对应 [Harald 的问题](https://github.com/JunhuiLi1017/TmCalculator/issues/8#issuecomment-5866928627)：
TMM 增量是否以“把起始项算在原始外侧碱基上”的预测为基线，从而已经吸收了这部分差值？

**结论：支持以配对核心为基线；证据强于物理直觉，但需要区分原文定义和数值重建。**
Bommarito 2000 对 DE 明确给出实验核心差分定义；2001 专利的 TMM Tables 2–3
全部48项的ΔG37以及47项ΔH可从其双链和核心实验数据重建，支持同一基线。
未取得原拟合程序，因此不把数值重建等同于复现作者当年的完整软件或统计流程。
这次只补文献证据，没有修改计算程序或热力学参数，也未向 GitHub 发布评论。

## 1. 为什么“两项研究都不合理”不能单独作为证明

即使两个研究分别测量 DE 和 TMM，也可能沿用同一种参数定义或分配约定；两篇文章
还共享研究人员。物理直觉不能排除拟合参数吸收基线差值的可能性。应检查参数是
相对于哪个核心求出的，以及是否存在额外端点调整。

若原模型写为 `full = core + increments`，改用原始外侧残基来选择起始项就变成
`full = core + increments + (init_outer - init_core)`。
后一个式子要与前一个等价，增量必须另减这项差值。下面检查的就是这种补偿是否存在。

## 2. DE：原文直接定义为相对于配对核心的增量

Bommarito S, Peyret N, SantaLucia J Jr (2000), *Thermodynamic parameters for DNA
sequences with dangling ends*, NAR 28:1929–1934，
[DOI 10.1093/nar/28.9.1929](https://doi.org/10.1093/nar/28.9.1929)，
[PMC原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC103285/)。
定位：Results and Discussion → Determination of dangling-end contributions to
 duplex stability，Equation 4；Tables 1–2。

作者说明：核心有实测值时，从含悬挂端双链的热力学量中减去核心实测值；没有实测值
时才预测**核心**。Equation 4 使用两个悬挂端的实验设计，差值除以2。

论文示例为 `(AGTAGCTAC)₂`，去掉两个5′悬挂A后的配对核心是 `(GTAGCTAC)₂`：

```text
5′ A GTAGCTAC   3′
3′   CATCGATG A 5′
     └ 配对核心 ┘
```

Table 1 的 ΔH / ΔG37（均为 kcal/mol）：

| 对象 | ΔH | ΔG37 |
|---|---:|---:|
| 含悬挂端双链 | −59.0 | −8.16 |
| 实测核心 | −51.6 | −7.01 |
| 差值 / 2 | −3.70 | −0.575 |
| Table 2 的 AG/C 项（包内 `AG/.C`） | −3.7 | −0.58 |

核心两端均为GC；外侧悬挂的是A。这一设计直接区分两种取端点方法。
差分值已经直接对应表中增量，没有再减“外侧A代替GC核心端点”的起始项差额。
因此将该增量加回预测核心时，应保留**核心GC末端**的起始项。
原文Equation 6还把含悬挂端双链显式拆成起始项、sym、WC堆叠和DE项；这些不是
让外侧未配对残基取得一个新起始项的规则。

## 3. TMM：专利原始数据能够重建同样的核心基线

SantaLucia J Jr & Peyret N (2001), WO2001094611A2：
[专利页面](https://patents.google.com/patent/WO2001094611A2/en)，
[原始PDF](https://patentimages.storage.googleapis.com/e5/cf/4b/ecca1b55126224/WO2001094611A2.pdf)。
定位：Appendix Table 1，印刷页47–51；Table 2，p52；Table 3，pp53–54。
PDF第48页对应印刷页47，以印刷页码为准。

Table 1 p51专门列出四条 **Core sequences**；其脚注说明热力学数据来自
熔解曲线拟合和 `1/Tm vs ln Ct` 的误差加权平均。Tables 2–3脚注写明这些参数由
Table 1配合equations 4和5计算。但专利附录没有一并给出清晰独立的TMM方法正文；
主文的同号方程属于平衡浓度计算，不能把它们冒认为这里的实验差分公式。
因此下述差分关系标为**从原表重建的证据**，不是声称找到了原文的完整方法段落。

### 3.1 能区分端点约定的实例：A·A 错配，GC 核心

```text
5′ A GTAGCTAC A 3′
3′ A CATCGATG A 5′
     └ 配对核心 ┘
```

两个末端错配增量在表中均可写为 `CA/GA`（另一端按链方向等价变换）。

| 对象 | ΔH (kcal/mol) | ΔG37 (kcal/mol) |
|---|---:|---:|
| Table 1 p47：带A·A末端错配的双链 | −60.3 | −9.03 |
| Table 1 p51：`GTAGCTAC` 核心 | −51.6 | −7.01 |
| 两者之差 / 2 | −4.35 | −1.01 |
| Table 2 p52：`CA/GA` | −4.3 | −1.01 |

在表格报告精度内一致。特别是这里外侧为A，实际核心两端为GC；增量是相对于
**实测GC核心**求出的，没有额外的外侧AT起始项补偿。

换回预测核心并不会改变增量的定义。用本包DNA_NN4的参数说明：若仍以外侧A
而非GC核心选择起始项，两端总共会多加 `(ΔH, ΔS) = (4.4, 13.8)`；这不是上述
实验核心差分所要求的项。在37°C时相应ΔG差约+0.11993 kcal/mol。
这里DNA_NN4只是说明当前两种实现之间的差额，不声称它是2001实验使用的WC版本。

### 3.2 反向实例：C·C 错配，AT 核心

`CTGAGCTCAC` 与配对核心 `TGAGCTCA`：Table 1 给出 ΔH=−50.6/−50.5，
ΔG37=−8.15/−7.73。差值除以2为 −0.05 / −0.21；Table 2对应
`AC/TC` 为 −0.1 / −0.21。这里外侧C不能把AT核心原有的起始项消掉。

### 3.3 扩大检查到 Tables 2–3 的全部48项

逐项由 Table 1 的实测核心差分重算：

- ΔG37：**48/48**在报告精度内与Tables 2–3一致；
- ΔH：**47/48**在报告精度内一致；
- 唯一ΔH异常是已发现的 `GG/CG`：核心差分得到−6.9，表中打印−0.7。
  这给疑似排印错误增加了新线索，仍未擅自更改包中参数。

还直接重建了issue示例所用的混合错配项 `CC/GA`：Table 1中
`AGTAGCTACC`双链ΔH=−57.1、ΔG37=−8.71，与`GTAGCTAC`核心的差值除以2
为−2.75、−0.85，对应Table 3的−2.7、−0.85。

逐行数据见 [terminal_core_reconstruction.tsv](terminal_core_reconstruction.tsv)。
可复算脚本：`tools/parameter_audit/check_terminal_core_differences.py`（仅需Python标准库）。
此检查没有把各行独立报告的ΔS当作简单差分；原表部分H/S/G本身不完全自洽，
不能为了让所有列符合一个恒等式而悄悄修改原始数据。

## 4. 结论的适用范围

- 可以说：**DE有明确的原文定义，TMM Tables 2–3全部48项有实验数据重建证据，均支持
  “预测配对核心 + 末端增量”的约定。** 将起始项转移到外侧残基，需要一个额外的
  补偿转换；所核查的数据并未显示这种转换。
- 不宜说：仅凭末端残基没有配对方，就已经证明历史拟合程序没有吸收差值。
- 不宜说：已找到全部48项TMM推导所用的原程序、已完全复现原回归，或其他NN参数
  版本（尤其Breslauer）的所有组合都因此获得实验验证。
- SantaLucia & Hicks 2004 Table 1的末端AT定义是辅助证据；TMM段落没有重列48项
  及完整推导。最有力的新增证据是上述原始实验核心差分。
- Peyret 2000博士论文可能提供完整方法，但学校页面标明WSU访问限制，本轮未取得全文：
  [Wayne State原始目录](https://digitalcommons.wayne.edu/oa_dissertations/3088/)。

因此，Harald的疑问可以从“只有模型自洽性论证”推进为“已有原始数据支持核心约定”；
未闭合部分应限缩为原始统计流程和历史程序实现，以及个别原表数值矛盾，
而不是继续笼统说缺少拟合基线的文献证据。
