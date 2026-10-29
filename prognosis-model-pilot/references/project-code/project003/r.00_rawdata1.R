rm(list = ls()); gc()
ORIGINAL_DIR <- "<REFERENCE_DATA_ROOT>/project/project003"
output <- file.path(ORIGINAL_DIR, "00_rawdata/TCGA_LIHC")
if (!dir.exists(output)) { dir.create(output, recursive = TRUE) }
setwd(ORIGINAL_DIR)

library(TCGAbiolinks)

# 下载count数据（用于差异分析）
query_count <- GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query_count, directory = file.path(output, "count"))
data_count <- GDCprepare(query_count, directory = file.path(output, "count"))
saveRDS(data_count, file.path(output, "TCGA_LIHC_count.Rds"))

# 下载临床信息
clinical <- GDCquery_clinic("TCGA-LIHC", type = "clinical")
write.csv(clinical, file.path(output, "TCGA_LIHC_clinical.csv"), row.names = FALSE)

library(SummarizedExperiment)
library(dplyr)
library(tibble)

# 读取已下载的 RDS 对象
tcga_raw <- readRDS(file.path(output, "TCGA_LIHC_count.Rds"))

# 1. 筛选样本：只保留 TP（肿瘤）和 NT（正常），排除 TR
sample_info <- as.data.frame(colData(tcga_raw))
sample_info <- sample_info[sample_info$shortLetterCode %in% c("TP", "NT"), ]
cat("筛选TP/NT后样本数：", nrow(sample_info), "\n")

# 2. 去除重复患者：每个 patient 只保留一个样本，优先保留 TP
# TP 字母序小于 NT，排序后 TP 排前面，duplicated 保留第一个即保留 TP
sample_info <- sample_info[order(sample_info$shortLetterCode), ]
sample_info <- sample_info[!duplicated(sample_info$patient), ]
cat("去除重复患者后样本数：", nrow(sample_info), "\n")
cat("  Tumor(TP):", sum(sample_info$shortLetterCode == "TP"),
    "  Normal(NT):", sum(sample_info$shortLetterCode == "NT"), "\n")

# 3. 筛选蛋白编码基因
gene_info <- as.data.frame(rowData(tcga_raw))
gene_info <- gene_info[
  gene_info$gene_type == "protein_coding" &
    !is.na(gene_info$gene_name) &
    gene_info$gene_name != "", ]
cat("蛋白编码基因数：", nrow(gene_info), "\n")

# 4. 提取 count 矩阵
count_raw <- assay(tcga_raw, "unstranded")
count_raw <- count_raw[rownames(count_raw) %in% rownames(gene_info),
                       rownames(sample_info)]
rownames(count_raw) <- gene_info[rownames(count_raw), "gene_name"]

# 处理一个 symbol 对应多行的情况：rowMeans 降序，保留第一个
count_df <- as.data.frame(count_raw)
count_df$symbol <- rownames(count_df)
rownames(count_df) <- NULL
count_df <- count_df %>%
  dplyr::select(symbol, everything()) %>%
  mutate(rowMean = rowMeans(dplyr::select(., starts_with("TCGA")))) %>%
  arrange(desc(rowMean)) %>%
  distinct(symbol, .keep_all = TRUE) %>%
  dplyr::select(-rowMean)
rownames(count_df) <- count_df$symbol
count_df$symbol <- NULL
cat("count symbol去重后基因数：", nrow(count_df), "\n")

# 5. 过滤低表达基因（至少 20% 样本 count > 0）
keep <- rowSums(count_df > 0) >= ncol(count_df) * 0.2
count_df <- count_df[keep, ]
cat("过滤低表达后基因数：", nrow(count_df), "\n")

# 6. 提取 TPM 矩阵，与 count 同步处理
tpm_raw <- assay(tcga_raw, "tpm_unstrand")
tpm_raw <- tpm_raw[rownames(tpm_raw) %in% rownames(gene_info),
                   rownames(sample_info)]
rownames(tpm_raw) <- gene_info[rownames(tpm_raw), "gene_name"]

tpm_df <- as.data.frame(tpm_raw)
tpm_df$symbol <- rownames(tpm_df)
rownames(tpm_df) <- NULL
tpm_df <- tpm_df %>%
  dplyr::select(symbol, everything()) %>%
  mutate(rowMean = rowMeans(dplyr::select(., starts_with("TCGA")))) %>%
  arrange(desc(rowMean)) %>%
  distinct(symbol, .keep_all = TRUE) %>%
  dplyr::select(-rowMean)
rownames(tpm_df) <- tpm_df$symbol
tpm_df$symbol <- NULL

# 与 count 矩阵保持相同基因集（低表达过滤后）
tpm_df <- tpm_df[rownames(count_df), ]
tpm_log2_df <- log2(tpm_df + 1)
cat("TPM矩阵基因数：", nrow(tpm_log2_df), "\n")

# 7. 构建分组和临床信息文件
group_df <- data.frame(
  sample                = rownames(sample_info),
  patient               = sample_info$patient,
  group                 = ifelse(sample_info$shortLetterCode == "TP", "Tumor", "Normal"),
  vital_status          = sample_info$vital_status,
  os_status             = ifelse(sample_info$vital_status == "Dead", 1, 0),
  os_time               = ifelse(!is.na(sample_info$days_to_death),
                                 sample_info$days_to_death,
                                 sample_info$days_to_last_follow_up),
  ajcc_pathologic_stage = sample_info$ajcc_pathologic_stage,
  gender                = sample_info$gender,
  age_at_index          = sample_info$age_at_index,
  stringsAsFactors      = FALSE
)

cat("\n分组统计：\n")
print(table(group_df$group))
cat("肿瘤样本中有生存时间的样本数：",
    sum(!is.na(group_df$OS.time) & group_df$group == "Tumor"), "\n")

# 8. 保存文件
write.csv(count_df,    file.path(output, "01.dat.TCGA_LIHC_count.csv"))
saveRDS(count_df,"./00_rawdata/count.Rds")
write.csv(tpm_log2_df, file.path(output, "01.dat.TCGA_LIHC_tpm_log2.csv"))
saveRDS(tpm_log2_df,"./00_rawdata/tpm.Rds")
write.csv(group_df,   file.path(output, "01.group.TCGA_LIHC.csv"), row.names = FALSE)
saveRDS(group_df,"./00_rawdata/clinical.Rds")
# 9. 验证检查
cat("\n=== 验证检查 ===\n")
cat("count矩阵维度：", nrow(count_df), "基因 ×", ncol(count_df), "样本\n")
cat("TPM矩阵维度：  ", nrow(tpm_log2_df), "基因 ×", ncol(tpm_log2_df), "样本\n")
cat("分组文件行数： ", nrow(group_df), "\n")
cat("列名一致性（应为TRUE）：",
    all(colnames(count_df) == colnames(tpm_log2_df)), "\n")
cat("列名与分组文件一致（应为TRUE）：",
    all(colnames(count_df) == group_df$sample), "\n")

message("\n00_rawdata 预处理完成，结果保存至：", output)
