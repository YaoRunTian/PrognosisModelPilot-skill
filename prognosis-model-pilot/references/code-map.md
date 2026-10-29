# 项目代码地图

项目编号已统一重编；示例中的历史项目文件名也已按新编号改写，原始路径与编号对应表只保存在维护者本地。

来源：用户指定服务器的 `<参考项目根目录>`，2026-10-21 只读盘点。现附五个代表项目的根目录参考脚本；完整代码树、表达矩阵及临床数据未打包。详见 [project-code/INDEX.md](project-code/INDEX.md)。

## 覆盖范围

24 个项目；3,296 个目录；19,316 个文件。候选脚本 772 个，其中 2 个为第三方 node_modules 内 Makefile，排除后扫描 770 个分析及辅助脚本，117,640 行，750 个精确不同的文件内容。包含 387 个根目录 R 脚本，另外有日期归档、结果目录重复脚本及少量辅助代码。候选类型为 R/r、Rmd/Rmarkdown、Python、Notebook、Shell 等；此统计不包含全部网页或演示文稿生成依赖代码。

全部 770 个文件读取用于结构、算法标记与去重扫描；这不等于逐行语义理解或验证。关键方法另做选段深入阅读：基础建模、单细胞和候选整合均有重点代码段，五个随包项目的具体文件见参考代码索引；绘图、药物附件和归档版本未全部逐行审阅。未做依赖安装验证或真实结果复现。

另用远端 R 4.3.3 的 `parse()` 静态解析了760个R/r文件；759个通过、1个失败。显式UTF-8复核后仍仅project001根目录的 `r.00_rawdata.R` 失败，其前12行检查确定为多余数字导致的语法错误。此检查不执行脚本，也不验证方法、依赖或运行时对象。

## 24 个项目

“代码数”排除第三方 Makefile，但包括日期归档及重复副本；不是独立流程数。队列标识从源码提取，不能单独证明患者独立性；复制残留的癌种名须用真实数据再核对。

| 文件夹 | 根脚本/代码总数 | 主要路线与复用边界 |
|---|---:|---|
| project001 | 20/40 | GSE148071单细胞→Scissor耐药/敏感关联→hdWGCNA∩NSCLC DEGs→Cox/LASSO，GSE41271 |
| project002 | 23/54 | BLCA、PRR11上皮细胞虚拟敲除+hdWGCNA+DEGs→Cox/LASSO，GSE13507 |
| project003 | 20/42 | TCGA-LIHC训练、LASSO后原尺度多变量Cox再拟合，GSE144269；单细胞在后续 |
| project004 | 18/38 | SIRT7虚拟敲除与表达分组DEGs，GSE62452训练、TCGA-PAAD验证；DEG分组含结局驱动步骤 |
| project005 | 20/47 | TGFB1成纤维细胞虚拟敲除，TCGA-UCEC 70:30内部留出；all脚本合并训练测试，非外部验证 |
| project006 | 24/48 | LUSC、AQP4分组差异与肿瘤差异交集，Cox/LASSO及多变量Cox分支，GSE157010验证 |
| project007 | 16/32 | LUAD、DEGs∩ERGs，Cox/LASSO，GSE72094验证；基础接口参考 |
| project008 | 16/32 | COAD/READ、DEGs∩MMRG，Cox/LASSO，GSE39582验证 |
| project009 | 19/36 | GEO双差异及先验集合整合，TCGA-LAML建模，GSE71014验证；bulk入口目录不同 |
| project010 | 20/40 | CESC、DEGs∩EGs，含DFS信息与多变量/step分支，GSE44001；单细胞用于解释 |
| project011 | 17/34 | BRCA、DEGs∩AMRGs，Cox/LASSO后多变量Cox分支，GSE42568 |
| project012 | 18/36 | ESCA、LASSO与svmRadial RFE交集，后多变量Cox；RFE用Alive/Dead分类 |
| project013 | 19/38 | BRCA、DEGs∩WGCNA∩单细胞marker，Cox/LASSO及多变量Cox；GSE20685 |
| project014 | 17/40 | GSE183795、DEGs∩金属离子基因集合，Cox/LASSO，单细胞解释及GraphBAN附件 |
| project015 | 7/14 | HRS糖尿病二分类：LASSO、RF、SVM、KNN、NB、glmnet；不是肿瘤预后 |
| project016 | 8/8 | 单细胞→AUCell→hdWGCNA，BayesPrism去卷积；没有完整预后建模链 |
| project017 | 24/47 | 肺癌单细胞→CopyKAT→CellChat配受体基因∩bulk DEGs→Cox/LASSO；含免疫治疗队列 |
| project018 | 3/3 | TNBC、TLR4单基因表达与生存，GSE58812；不是多基因机器学习模型 |
| project019 | 20/40 | CRC、PRNP分组及双重DEG交集→Cox/LASSO；单细胞机制解释、GSE39582 |
| project020 | 13/26 | DRESIS分子∩DEGs，CGGA693训练、TCGA验证；单因素阈值与通用模板不同 |
| project021 | 16/28 | 原发/肝转移分类：binomial LASSO、手写SVM-RFE、RF候选整合；不是生存终点 |
| project022 | 17/34 | STAD单因素及PH筛选→多变量Cox→predict(lp)，GSE84433；无活跃LASSO主线 |
| project023 | 9/9 | BRCA/TNBC、双Scissor→AUCell→hdWGCNA及先验交集→Cox/LASSO；根代码止于建模 |
| project024 | 3/4 | ESCC单细胞与放疗敏感性Scissor；未见生存建模/验证脚本 |

## 代表路线的文件链

### 基础 bulk / 生存：project007

`r.00_rawdata.R` → `r.01_DEG.R` → `r.02_venn.R` → `r.03_model.R` → `r.04_evaluate.R`；GEO验证由 `r.05_testdata.R` → `r.06_testmodel.R`。候选、模型、评分三个接口需同时修正和验证，不能只复制模型脚本。

### 单细胞与 bulk 整合：project001

`r.01_sc.R` → `r.02_reannotated.R` → `r.03_scissor.R` → `r.04_hdWGCNA.R`；结合 `r.05_DEG.R`，由 `r.06_venn.R` 取候选，再进入 `r.07_model.R` → `r.08_evaluate.R`，验证是 `r.09_GSE41271.R` → `r.10_testmodel.R`。

重要交接：单细胞最终对象 `RData/R.06_scRNA_done.Rds` → `pbmc_all_with_subcelltype.Rds` → `R.07_scRNA_Scissor.RData`（文件名虽为RData，实际用saveRDS）；hdWGCNA的 `06.hub_module_gene_df.xlsx` 与DEG基因表取交集。Scissor phenotype为耐药/敏感、binomial，不能把该调用直接换成生存。细胞亚群、模块ID、soft power必须重新确定。

### 虚拟敲除：project002、project004、project005

先确认目标细胞和目标基因、counts layer、HVG及强制保留敲除基因、运行实现和响应筛选。部分脚本使用 `scTenifoldNet.lite` 的 fast 实现；库中载入 `scTenifoldKnk` 不代表调用了同一实现。project002还整合hdWGCNA，project004/project005的原始p值阈值也不能默认当FDR。

### 其他入口

- 要沿用多变量Cox预测对象，先看project022；要研究LASSO后重新拟合多变量Cox，先看project003。
- 单细胞去卷积看project016；单基因描述/生存看project018。这两者都不自动包含完整机器学习预后流程。
- project015和project021只借用适合的分类实现；不能宣传其分类性能为预后性能。
- project023/project024只可复用确实存在的阶段，不能假定编号后续文件存在。
