# 复用前需要修正或澄清的行为

以下是源码观察及其适配含义，未对已有论文结果作结论，也未运行真实队列。路径相对于 `<参考项目根目录>`。行号对应 2026-10-21 读取的根目录版本；后续修改须重新定位。

## 标准化与评分不一致

`project007/r.03_model.R:107–117` 显式对 LASSO 输入使用 `scale()`；`r.04_evaluate.R:9–20` 与 `r.06_testmodel.R:20–31` 用表达值直接乘系数并取指数。脚本中未见评分时应用训练中心及尺度。`project001`、`project002`、`project004` 也观察到相同结构。

适配时追踪最终模型：若下游实际改用原尺度重新拟合的多变量 Cox，不能机械套用 LASSO 的缩放问题。例如 `project003` 有这种再拟合分支，应核对实际预测对象。对于仍使用显式缩放后系数的分支，保存训练中心/尺度并在评分使用，或在训练阶段采用明确一致的管线，不能悄悄改变旧模型。

## 生存 ROC 的事件字段

`project007/r.04_evaluate.R:123–127`、`project001/r.08_evaluate.R:123–130` 等将派生的 `riskgroup` 传给 `survivalROC(status=...)`。这不能回答模型对生存事件的预测能力。

应使用真实事件字段（如已核实编码的 `os_status`）并把连续分数传给 marker。官方接口将 status 定义为事件/死亡=1、否则=0，Stime 为事件或删失时间；参见 [survivalROC 官方手册](https://cran.r-project.org/web/packages/survivalROC/survivalROC.pdf)。

## 验证集再次选择分组阈值

`project007/r.06_testmodel.R:37–48` 在验证数据上用其时间和事件重新调用 `surv_cutpoint()`。`project001/r.10_testmodel.R:37–45` 与多个其他项目有同类行为。

冻结模型验证应保持预先确定的评分和分组规则。验证集重新优化的分组可以单列为探索或更新分析，不将其 KM 显著性当作冻结规则的验证表现。`project005` 的 train/test 来自同一 TCGA 队列，其合并分析脚本也不是外部独立验证。内部留出、外部队列及合并结果须分开描述。

## 训练切分与候选来源

不要假设所有模型基因都来自单细胞。`project012/r.03_model.R:230–281` 用 Alive/Dead 做 `svmRadial` RFE，并与 LASSO 结果取交集；该分类步骤不建模随访时间及删失。可以如实复现并解释它，但不能把其准确率作为生存性能。

`project013/r.04_venn.R:6–11` 整合 bulk DEGs、WGCNA 及单细胞 marker；其 WGCNA 脚本在当前可见部分未显式导出下游读取的 `08.GS_MM.gene.txt`，因此交接文件还需验证或补上明确导出。不是所有缺失交接都能通过按编号运行全部脚本解决。

## 依赖交互会话与版本差异

多个建模脚本把 `save.image()` 注释，而下游直接 `load()` 工作区。例如 `project007/r.03_model.R:222` 与 `r.04_evaluate.R:6`。服务器现有对应工作区文件确实存在，但仅按源码重跑仍不保证重建它。新实现要显式保存模型及所需对象。

历史根目录并不总比归档正确：`project001/r.00_rawdata.R:12` 原先有多余的 `8`，其日期归档同处没有，随包参考副本已修复；重注释聚类与 CellChat 内容也有实质变化。适配先语法检查，再比较对应阶段，不直接整套选择根目录。

代码库同时存在旧 Seurat slot、Assay5/layer 和不同 Monocle 使用方式；读取对象前检查实际版本、assay 与 layer。[Seurat v5 官方数据访问说明](https://satijalab.org/seurat/articles/seurat5_essential_commands) 区分 counts、归一化数据、assay 及细胞集合。不要只做字符串替换升级旧对象。

`glmnet` 默认行为也可能随版本变化。复现时记录版本、标准化和 ties 处理；支持的参数先查询实际安装版本，再参考 [glmnet 官方 Cox 文档](https://glmnet.stanford.edu/articles/Coxnet.html)，不能向旧版本盲目添加新参数。

## 表达尺度与对比方向

`project007/r.01_DEG.R:12` 对输入执行 `round(2^count-1)`；只有确定原输入为 log2(count+1) 才有对应含义。`:24` 的 contrast 是 Normal 对 Tumor，解释上下调时要跟随该方向，不能默认火山图的正方向表示肿瘤更高。

GEO 临床表的固定列号也属于特定数据版本假设。例如 `project007/r.05_testdata.R:14–15` 按位置取列；新队列应按字段含义和数据字典转换，不能照抄列号。

其他已观察问题：`project002/r.07_venn.R:12–13` 在赋值前写 hubgenes；project017 的 CopyKAT 根脚本只保存 Epithelium，但 CellChat 又期望多类细胞；project023 的部分 Scissor 对象在使用前未明确赋值；project024 的 phenotype 先 factor 再 as.numeric，编码可能从预期 0/1 变成 1/2。复用相应模块时先核查这些交接。
