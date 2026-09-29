# tm_nn 参数及参考文献核查

核查日期：2026-09-28。范围：`R/zzz.R` 构建并由 `R/sysdata.rda` 提供给
`tm_nn()` 的全部 NN、IMM、TMM、DE 表，及 Zuber 末端伴随表；另检查盐与化学
校正的引文。两份运行时常量完全一致。本次只更正文献、说明与运行结果的引用文字，
没有改动热力学数值、表名或计算公式。核查过程中已有的工作区修改予以保留。

**不能把“DOI 与论文对应”当作“所有数字已获原文验证”。** 下表逐套列明证据等级；
其中仍有原文数值待确认项。NULL 占位表不是已实现的参数集，不计入覆盖数量；
`GC_VARTAB` 属于 `tm_gc()`，不属于本次 NN 数值范围。

## 1. 全量清单和证据等级

逐行数值、来源、核查状态和实现性扩展见
[tm_nn_parameter_sources.tsv](tm_nn_parameter_sources.tsv)。

- `primary`：已打开原文数值表或扫描页核对（不表示整个模型实现正确）。
- `primary_reprint`：已核研究论文重列的原参数，而非重新取得最早论文全文。
- `primary_database`：已核原始参数数据库的机器可读数据及换算。
- `author_data`：已核该参数作者发布的 VarGibbs `.par` 数据。
- `secondary`：引文身份已核；数值与 VarGibbs 收录版本一致，最早原文表仍待独立核对。
- `mixed`：逐行证据等级不同，见 TSV。
- `pending`：引文身份正确，但本轮没有取得原始补充表，不能宣称数字已复核。

注意：通常 23 个存储行不等于 23 个独立拟合参数。常规同质双链有 10 个独立
堆叠项、6 个链对称展开项，以及起始／对称项；未使用的起始项填零属于实现。
RNA/DNA 的 16 个堆叠项彼此独立，不能套用交换 RNA/DNA 链的对称关系。

| 运行时表名 | 行数 | 来源及位置 | 本轮证据 |
|---|---:|---|---|
| `DNA_NN_Breslauer_1986` | 23 | [10.1073/pnas.83.11.3746](https://doi.org/10.1073/pnas.83.11.3746); Table 2; p3749 initiation discussion | primary |
| `DNA_NN_Sugimoto_1996` | 23 | [10.1093/nar/24.22.4501](https://doi.org/10.1093/nar/24.22.4501); VarGibbs P-SG96.par | secondary |
| `DNA_NN_Allawi_1998` | 23 | [10.1021/bi962590c](https://doi.org/10.1021/bi962590c), [10.1073/pnas.95.4.1460](https://doi.org/10.1073/pnas.95.4.1460); 1997 Table 1; 1998 Table 2 | primary |
| `DNA_NN_SantaLucia_2004` | 23 | [10.1146/annurev.biophys.32.110601.141800](https://doi.org/10.1146/annurev.biophys.32.110601.141800); Table 1 | primary |
| `RNA_NN_Freier_1986` | 23 | [10.1073/pnas.83.24.9373](https://doi.org/10.1073/pnas.83.24.9373); Table 2 | primary |
| `RNA_NN_Xia_1998` | 23 | [10.1021/bi9809425](https://doi.org/10.1021/bi9809425), [10.1021/bi3002709](https://doi.org/10.1021/bi3002709); Xia parameters reprinted in Chen 2012 Table 3, outside parentheses | primary_reprint |
| `RNA_NN_Chen_2012` | 34 | [10.1021/bi3002709](https://doi.org/10.1021/bi3002709); Table 3 and footnote c | primary |
| `RNA_DNA_NN_Sugimoto_1995` | 23 | [10.1021/bi00035a029](https://doi.org/10.1021/bi00035a029); VarGibbs P-SG95.par | secondary |
| `DNA_IMM_Peyret_1999` | 87 | [10.1021/bi962590c](https://doi.org/10.1021/bi962590c), [10.1093/nar/26.11.2694](https://doi.org/10.1093/nar/26.11.2694), [10.1021/bi9803729](https://doi.org/10.1021/bi9803729), [10.1021/bi9724873](https://doi.org/10.1021/bi9724873), [10.1021/bi9825091](https://doi.org/10.1021/bi9825091), [10.1093/nar/gki918](https://doi.org/10.1093/nar/gki918); Six-source composite; see row-level inventory | mixed |
| `DNA_TMM_Bommarito_2000` | 48 | [WO2001094611A2](https://patents.google.com/patent/WO2001094611A2/en); Tables 2-3, printed pp52-54 | primary |
| `DNA_DE_Bommarito_2000` | 32 | [10.1093/nar/28.9.1929](https://doi.org/10.1093/nar/28.9.1929); Table 2 | primary |
| `RNA_DE_Turner_2010` | 48 | [10.1093/nar/gkp892](https://doi.org/10.1093/nar/gkp892); NNDB Turner 2004 dangle_dh.txt and dangle_dg.txt | primary_database |
| `DNA_NN_Weber_OW04_69` | 23 | [10.1093/bioinformatics/btu751](https://doi.org/10.1093/bioinformatics/btu751); Author release VarGibbs 5.0 data/AOP-OW04-69.par | author_data |
| `DNA_NN_Weber_OW04_119` | 23 | [10.1093/bioinformatics/btu751](https://doi.org/10.1093/bioinformatics/btu751); Author release VarGibbs 5.0 data/AOP-OW04-119.par | author_data |
| `DNA_NN_Weber_OW04_220` | 23 | [10.1093/bioinformatics/btu751](https://doi.org/10.1093/bioinformatics/btu751); Author release VarGibbs 5.0 data/AOP-OW04-220.par | author_data |
| `DNA_NN_Weber_OW04_621` | 23 | [10.1093/bioinformatics/btu751](https://doi.org/10.1093/bioinformatics/btu751); Author release VarGibbs 5.0 data/AOP-OW04-621.par | author_data |
| `DNA_NN_Weber_OW04_1020` | 23 | [10.1093/bioinformatics/btu751](https://doi.org/10.1093/bioinformatics/btu751); Author release VarGibbs 5.0 data/AOP-OW04-1020.par | author_data |
| `DNA_NN_Weber_2015` | 23 | [10.1093/bioinformatics/btu751](https://doi.org/10.1093/bioinformatics/btu751); Author release VarGibbs 5.0 data/AOP-CMB.par | author_data |
| `RNA_NN_Weber_VIF_71` | 23 | [10.1016/j.chemphys.2019.01.016](https://doi.org/10.1016/j.chemphys.2019.01.016); Author release VarGibbs 5.0 data/AOP-VIFRW-71.par | author_data |
| `RNA_NN_Weber_VIF_121` | 23 | [10.1016/j.chemphys.2019.01.016](https://doi.org/10.1016/j.chemphys.2019.01.016); Author release VarGibbs 5.0 data/AOP-VIFRW-121.par | author_data |
| `RNA_NN_Weber_VIF_221` | 23 | [10.1016/j.chemphys.2019.01.016](https://doi.org/10.1016/j.chemphys.2019.01.016); Author release VarGibbs 5.0 data/AOP-VIFRW-221.par | author_data |
| `RNA_NN_Weber_VIF_621` | 23 | [10.1016/j.chemphys.2019.01.016](https://doi.org/10.1016/j.chemphys.2019.01.016); Author release VarGibbs 5.0 data/AOP-VIFRW-621.par | author_data |
| `RNA_NN_Weber_VIF_1021` | 23 | [10.1016/j.chemphys.2019.01.016](https://doi.org/10.1016/j.chemphys.2019.01.016); Author release VarGibbs 5.0 data/AOP-VIFRW-1021.par | author_data |
| `RNA_NN_Weber_FIF_71` | 23 | [10.1016/j.chemphys.2019.01.016](https://doi.org/10.1016/j.chemphys.2019.01.016); Author release VarGibbs 5.0 data/AOP-FIFRW-71.par | author_data |
| `RNA_NN_Weber_FIF_121` | 23 | [10.1016/j.chemphys.2019.01.016](https://doi.org/10.1016/j.chemphys.2019.01.016); Author release VarGibbs 5.0 data/AOP-FIFRW-121.par | author_data |
| `RNA_NN_Weber_FIF_221` | 23 | [10.1016/j.chemphys.2019.01.016](https://doi.org/10.1016/j.chemphys.2019.01.016); Author release VarGibbs 5.0 data/AOP-FIFRW-221.par | author_data |
| `RNA_NN_Weber_FIF_621` | 23 | [10.1016/j.chemphys.2019.01.016](https://doi.org/10.1016/j.chemphys.2019.01.016); Author release VarGibbs 5.0 data/AOP-FIFRW-621.par | author_data |
| `RNA_NN_Weber_FIF_1021` | 23 | [10.1016/j.chemphys.2019.01.016](https://doi.org/10.1016/j.chemphys.2019.01.016); Author release VarGibbs 5.0 data/AOP-FIFRW-1021.par | author_data |
| `RNA_DNA_NN_Weber_2019_FT` | 23 | [10.1016/j.bpc.2019.106189](https://doi.org/10.1016/j.bpc.2019.106189); Author release VarGibbs 5.0 data/AOP-DRFT.par | author_data |
| `RNA_DNA_NN_Weber_2019_VH` | 23 | [10.1016/j.bpc.2019.106189](https://doi.org/10.1016/j.bpc.2019.106189); Author release VarGibbs 5.0 data/AOP-DRVH.par | author_data |
| `RNA_DNA_NN_Weber_2019_LS` | 23 | [10.1016/j.bpc.2019.106189](https://doi.org/10.1016/j.bpc.2019.106189); Author release VarGibbs 5.0 data/AOP-DRLS.par | author_data |
| `RNA_DNA_NN_Banerjee_2020` | 23 | [10.1093/nar/gkaa572](https://doi.org/10.1093/nar/gkaa572); Table 2; footnotes b/c | primary |
| `RNA_NN_Zuber_2022` | 34 | [10.1093/nar/gkac261](https://doi.org/10.1093/nar/gkac261); Tables 1A-1B | primary |
| `RNA_NN_Zuber_2022_END` | 24 | [10.1093/nar/gkac261](https://doi.org/10.1093/nar/gkac261); Tables 1A-1B, six end classes expanded to 24 keys | primary |
| `RNA_NN_Ghosh_2023_PEG200` | 23 | [10.1093/nar/gkad020](https://doi.org/10.1093/nar/gkad020); Table 1 | primary |
| `DNA_NN_Ghosh_2020_PEG200` | 23 | [10.1073/pnas.1920886117](https://doi.org/10.1073/pnas.1920886117); SI Tables S7-S8, as recorded in existing v1.1.1 provenance | pending |

合计 **36 套表、974 个存储行**；表名和行数按本次工作树导出。

## 2. Allawi 表：1997 原文确实含有全部 13 项

正确主引文是 Allawi & SantaLucia (1997)，Biochemistry 36:10581–10594，
[DOI 10.1021/bi962590c](https://doi.org/10.1021/bi962590c)，**Table 1，p10583**。
标题虽为内部 G·T 错配，文中也重新拟合并列出 Watson–Crick 参数。
因此不能凭标题排除它；Biopython 此处有原文依据。
[SantaLucia 1998，Table 2](https://doi.org/10.1073/pnas.95.4.1460) 也列出相同参数并引用1997论文。
1996 年 [SantaLucia, Allawi & Seneviratne](https://doi.org/10.1021/bi951907q)
是较早的不同参数版本，不能作为这组数值的唯一出处。

| key | ΔH (kcal/mol) | ΔS (cal/mol/K) | 与1997 Table 1 |
|---|---:|---:|---|
| `init_A/T` | 2.30 | 4.10 | 一致 |
| `init_G/C` | 0.10 | -2.80 | 一致 |
| `sym` | 0.00 | -1.40 | 一致 |
| `AA/TT` | -7.90 | -22.20 | 一致 |
| `AT/TA` | -7.20 | -20.40 | 一致 |
| `TA/AT` | -7.20 | -21.30 | 一致 |
| `CA/GT` | -8.50 | -22.70 | 一致 |
| `GT/CA` | -8.40 | -22.40 | 一致 |
| `CT/GA` | -7.80 | -21.00 | 一致 |
| `GA/CT` | -8.20 | -22.20 | 一致 |
| `CG/GC` | -10.60 | -27.20 | 一致 |
| `GC/CG` | -9.80 | -24.40 | 一致 |
| `GG/CC` | -8.00 | -19.90 | 一致 |

旧 DOI `10.1093/nar/26.11.2694` 是 **C·T 内部错配**来源，现从此 WC 表引用中移除，
保留在 IMM 引文中。旧表名 `DNA_NN_Allawi_1998` 为兼容保留，说明文字标明真实年份。
与2004表相比是 **9/10 个独立堆叠项相同**；`AA/TT` 与 `TT/AA` 是同一独立项的
对称写法，不能算作两个独立差别。起始项还涉及全局 init 与末端罚项的重新分配。

## 3. IMM：六篇论文，87 行，而不是90行

| 类别 | 行数 | 正确文献 | 数值证据 |
|---|---:|---|---|
| G·T | **11** | Allawi & SantaLucia 1997, 10.1021/bi962590c, Table 5 | VarGibbs交叉核对；含8个单错配项和3个串联G·T项 |
| G·A | 8 | Allawi & SantaLucia 1998, 10.1021/bi9724873, Table 4 | VarGibbs交叉核对 |
| A·C | 8 | Allawi & SantaLucia 1998, 10.1021/bi9803729, Table 4 | VarGibbs交叉核对；pH 7参数 |
| C·T | 8 | Allawi & SantaLucia 1998, 10.1093/nar/26.11.2694, Table 4 | VarGibbs交叉核对 |
| A·A、C·C、G·G、T·T | 16 | Peyret et al. 1999, 10.1021/bi9825091 | VarGibbs交叉核对 |
| I·A、I·C、I·G、I·T、I·I | 36 | Watkins NE Jr & SantaLucia J Jr 2005, 10.1093/nar/gki918, Table 2 | 原文HTML表逐行比对，36/36一致 |

合计 **11+8+8+8+16+36=87**。把G·T计为14时合计是90，不能用来说明87行全覆盖。
每一个 key 的归属已写入 TSV，不再把整个复合表只归给Peyret。
全部87行与VarGibbs收录数值一致，但这不能替代尚未逐行复查的早期原文表。
VarGibbs `P-AL98C.par` 的G·A文件头本身写错DOI，不能直接复制其元数据。

## 4. TMM、DE：来源不同，必须分开

**TMM 48行**：已目视核对[原专利 WO2001094611A2](https://patents.google.com/patent/WO2001094611A2/en)
Tables 2–3，印刷页52–54，数值全部与包一致。作者/发明人为 SantaLucia J Jr &
Peyret N，公开日期2001-12-13。Bommarito 2000不是这套TMM表的出处。

**原始来源疑点**：专利 `GG/CG` 列 ΔH=−0.7、ΔS=−19.2，却列 ΔG37=−0.96。
由前两项计算得到 **+5.25488 kcal/mol**，无法由舍入解释。
代码忠实保留专利值；不能擅自把−0.7改成−6.7或其他推测值。
应联系原作者或寻找勘误/后续可靠数据再决定改值。

**关于Harald问的拟合**：专利将末端项从其熔解实验热力学数据结合模型方程提取，
参数表给出误差；它不是Bommarito的32项悬垂端实验表。
SantaLucia & Hicks 2004在Terminal Mismatches部分讨论末端错配，但归于
S. Varma & J. SantaLucia “manuscript in preparation”，没有重列这48个数，
也不足以证明其每一步拟合与专利相同。2026-09-29的进一步核查发现：专利Table 1
确实提供实测核心；由双链减核心再除以2，可重建Tables 2–3全部48项ΔG37及47项ΔH。
这支持以配对核心为基线，详见[末端起始项证据](terminal_initiation_evidence.md)。
尚未取得历史拟合程序；数值重建不等于复现其完整软件与统计流程。内部G·T的1997工作所述SVD及
非唯一二聚体分解也不能直接挪给TMM。

**DNA DE 32行**：与[Bommarito, Peyret & SantaLucia 2000 Table 2](https://doi.org/10.1093/nar/28.9.1929)
的三幅原表扫描逐行一致，原引文正确。

**RNA DE 48行**：与[NNDB Turner 2004 dangling ends](https://rna.urmc.rochester.edu/NNDB/rna_2004/rna_2004_dangling_ends.html)
匹配。ΔH直接取表，ΔS=(ΔH−ΔG37)×1000/310.15，与存储精度一致（最大允许差0.055）。
2010是Turner & Mathews的NNDB说明论文年份，不是48项都在2010年测得；
数据库说明这些参数由Serra & Turner 1995汇编，原始实验来自多篇更早论文。
所以正确说明应同时给出数据库版本、数据库论文及其原始实验参考链。

## 5. 其他NN表的结论及模型问题

- **Breslauer1986**：原文Table 2列解链方向；包用成链方向，符号转换正确。
  p3749的起始自由能以25°C为基准，在ΔHinit=0约定下5/6 kcal转换为约−16.8/−20.1
  entropy units；不要拿外部文件的另一种init/sym约定直接替换。
- **SantaLucia2004、Freier1986**：相应原表数值一致。2004参考文献作者顺序应为
  **SantaLucia J Jr, Hicks D**，不是“Hicks LD, Santalucia J”。
- **Sugimoto1996 DNA、Sugimoto1995 RNA/DNA**：作者、年份和DOI对应；数值与
  VarGibbs所收录版本相符。本轮未取得最初原文数值表逐项复查，列为secondary。
- **Xia1998**：与Chen2012 Table3重列的旧WC数值及VarGibbs一致。
- **Chen2012**：存储数字可以在Table3找到，但包组合的是括号内重拟合WC项和
  GU项；Table3脚注c明确说GU项以Xia1998的WC项为基础拟合。
  因此“每个数字在论文里”不等于这组组合复现论文GU模型。
  此外旧模型的GGUC/CUGG特殊项没有实现。需单独决定兼容性与模型修复。
- **Weber/VarGibbs全部19套**：逐行核作者发布的VarGibbs5.0 `.par` 文件，
  在包的四位小数精度内一致。DNA6套对应Weber2015；RNA10套对应**Ferreira et al.2019**；
  RNA/DNA3套对应**Basilio Barbosa et al.2019**。`Weber`在后两类名称里不是第一作者。
  文件名逐表记录在上表；发布包来自
  [作者服务器](https://bioinf.fisica.ufmg.br/software/vargibbs-5.0/Debian_12/vargibbs_5.0.orig.tar.gz)。
- **Banerjee2020**：16个堆叠值及链方向映射与Table2一致，**起始项应用不一致**。
  脚注b/c定义：任一端GC选一次GC-init；两端均AT选一次AT-init。
  当前实现按每个端点计费。两端GC的ΔSinit当前−9.8而非−4.9；混合端−11.9而非−4.9；
  两端AT−14.0而非−7.0。论文init的ΔG与H/S还存在单独的不自洽，不能混淆两类问题。
- **Zuber2022及24行末端伴随表**：Tables1A/B一致；六种末端类别正确展开。
  **新版模型不需要额外GGUC/CUGG项**，原文明确如此。旧provenance把它写成遗漏
  必需项并声称必然低估稳定性，这段已纠正。
- **Ghosh2023 RNA**：Table1的堆叠、init、terminal AU及sym一致；40wt%PEG200、100mMNaCl。
- **Ghosh2020 DNA**：DOI和文章对应，包记录为SI S7/S8的参考项加拥挤增量。
  本轮文章可读，但PMC/出版社补充材料下载返回验证页，**未重新核到S7/S8原件**。
  既有provenance的数值与验证说明属于先前记录，不能当作本轮独立核实证据。

## 6. 盐、添加剂及其他参数

| 参数/方法 | 来源对应与处理 |
|---|---|
| `Schildkraut2010` | 实为Schildkraut & Lifson **1965**；保留兼容名称，更正说明 |
| `Wetmur1991` | Wetmur1991综述对应经典单价盐公式 |
| `SantaLucia1996` | SantaLucia, Allawi & Seneviratne1996，10.1021/bi951907q；补参考文献 |
| `SantaLucia1998-1/-2` | SantaLucia1998，10.1073/pnas.95.4.1460；Tm式与entropy式不是同一个数值变换 |
| `Owczarzy2004` | 10.1021/bi034621r，原论文研究钠离子；纠正“原式包括二价离子”的说明 |
| `Owczarzy2008` | 10.1021/bi702363u，镁/单价离子竞争模型；本轮不修改公式 |
| 其他盐式前的Mg钠当量调整 | von Ahsen, Wittwer & Schutz2001，10.1093/clinchem/47.11.1956，不能归到1965原式 |
| 默认DMSO 0.75 | von Ahsen et al.2001；其他可选系数的原始逐项来源本轮未全部确认，不再统称已核“published values” |
| molar formamide | Blake & Delcourt1996，10.1093/nar/24.11.2095；补引文 |
| percent formamide 0.65 | 软件采用的经验默认值，不能冒充Hutton1977的精确结果；该论文报告约0.60。0.72选项原始数值出处尚待独立确认 |
| 浓度、shift、self_comp、默认开关 | 用户输入/算法设置，不是作者年份热力学拟合表，无须伪造一一对应的实验出处 |

本节核的是来源与公式身份，不是所有盐/添加剂公式的全系数原文复算。

## 7. 已实施修正与后续事项

已同步修正 `tm_nn()` 输出引用、roxygen说明及生成帮助页，修正旧provenance
关于Banerjee与Zuber的断言。保留旧API表名，避免破坏调用代码；若未来新增规范别名，
建议Allawi1997、SantaLucia_Peyret2001、Ferreira2019、BasilioBarbosa2019。

后续必须分开处理：
1. 原件复核：Ghosh2020 SI S7/S8、Sugimoto两套，以及IMM中尚未逐行复查的早期原表。
2. 模型修复：Banerjee起始项、Chen WC/GU组合及特殊motif。
3. 原始来源勘误：专利GG/CG异常；没有证据前不臆造修正值。
4. 可选化学校正系数逐项追溯。

这些未闭环项已明确列出，故本次结论不是“所有模型都已无误”。

## 8. 本次验证

- 来源清单覆盖全部36套运行时表、974行，无重复/遗漏，逐行H/S与运行时常量一致。
- 参数构建结果与核查前快照及`R/sysdata.rda`均完全一致。
- 修改后的R文件及Rd帮助页解析通过；实际调用`tm_nn()`确认WC、TMM和六篇IMM引用正确。
- 现有`test_nn_rc_completion.R`测试无失败；运行环境产生locale警告，不涉及参数差异。
