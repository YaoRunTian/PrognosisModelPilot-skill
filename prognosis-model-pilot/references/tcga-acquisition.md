# TCGA 转录组及临床信息的默认获取方式

2026-10-23 用户指定默认使用project003项目中的 TCGAbiolinks/GDC 下载路线。参考原文件：`<参考项目根目录>/project003/r.00_rawdata.R:8–19`，另有 `r.00_rawdata1.R:10–22`。本技能的辅助脚本是对下载逻辑的参数化重新实现，不是原始项目代码的完整副本。

## 默认接口与运行

表达查询：`GDCquery(project=..., data.category="Transcriptome Profiling", data.type="Gene Expression Quantification", workflow.type="STAR - Counts", access="open")`。同一查询对象及下载目录用于 `GDCdownload()` 和 `GDCprepare()`。临床单独使用 `GDCquery_clinic(project=..., type="clinical")`。

复制 [download_tcga.R](../scripts/download_tcga.R) 到生成项目的数据获取阶段，或直接运行：

```bash
Rscript scripts/download_tcga.R --project TCGA-LIHC --out data/TCGA-LIHC
```

`--project` 和 `--out` 必填，癌种及目录不硬编码。默认 API 下载、每块20个文件；可通过 `--files-per-chunk` 调整。20只是可配置的传输分块大小。需要已安装的 TCGAbiolinks 和 SummarizedExperiment；脚本不会自动安装依赖。不要在仅要求生成或检查代码时启动完整队列下载。

## 输出契约

- `expression_SE.rds`：完整 SummarizedExperiment。
- `counts.rds`：`unstranded` 原始 counts，保留下载对象中的基因和样本标识。
- `tpm.rds`、`log2_tpm_plus1.rds`：TPM与 `log2(TPM+1)`，分别命名，不能与counts混用。
- `clinical_raw.rds`：临床API原始返回表；普通无列表列的表另输出 `clinical_raw.csv`，否则在清单中明确记录CSV未生成，完整内容仍在RDS中。
- `sample_metadata.rds`、`gene_metadata.rds`：样本及基因元数据。
- `query.rds`、`acquisition_manifest.rds`：查询对象、查询参数、GDC文件清单、下载时间、assay名称、输出列表及R/包版本。清单最后写出，缺清单不作为完整阶段交付。
- `GDCdata/`：下载缓存，不删除已下载的源文件。

脚本先核对实际 `assayNames()`，要求counts与TPM可用。基因ID转换、protein_coding过滤、symbol去重和样本筛选在后续阶段做；此处保留全部返回行列，不继承原项目的NT/TP排序去重、固定样本排除、固定低表达阈值。

原始临床信息不是现成的分析用OS表。患者与样本的连接、重复患者处理、事件映射、缺失随访与具体终点需单独核对；不直接沿用旧脚本的 `ifelse(days_to_death, days_to_last_follow_up)`。GDC随访字段可能随版本变化，实际构造终点时查当前数据字典。

## 已有数据与失败处理

优先使用用户已提供且来源、尺度、完整性可核实的本地数据。该默认值约定新获取数据的入口，不触发自动重下载或其他数据源替换。脚本不覆盖已有同名输出；已有完整工件应由项目预检后直接复用，部分失败则检查缺项或选择新的导出目录。下载失败时保留错误和缓存，不用空表或其他队列顶替。

## 依据与验证边界

- [TCGAbiolinks 官方手册](https://bioconductor.org/packages/release/bioc/manuals/TCGAbiolinks/man/TCGAbiolinks.pdf)
- [维护者的GDC数据整理实现](https://github.com/BioinformaticsFMRP/TCGAbiolinks/blob/master/R/prepare.R)
- [NCI GDC mRNA流程及counts/TPM定义](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)

本次核对了本地 TCGAbiolinks 2.34.1 的函数签名；脚本做离线替身与小型真实SE对象测试，不启动TCGA下载。实际网络可用性、完整队列下载及临床字段覆盖仍以运行时为准。
