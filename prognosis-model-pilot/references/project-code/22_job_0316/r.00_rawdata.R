rm(list = ls()); gc()
ORIGINAL_DIR <- "<REFERENCE_DATA_ROOT>/project/22_job_0316"
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
saveRDS(data_count, file.path(output, "TCGA_LIHC_raw.Rds"))
# 下载临床信息
clinical <- GDCquery_clinic("TCGA-LIHC", type = "clinical")
write.csv(clinical, file.path(output, "TCGA_LIHC_clinical.csv"), row.names = FALSE)
library(SummarizedExperiment)
library(dplyr)
library(tibble)
# 读取已下载的 RDS 对象
tcga_raw <- readRDS(file.path(output, "TCGA_LIHC_raw.Rds"))
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
#保留01A和11A样本
table(substring(sample_info$sample,14,16))
# 01A 01B 11A 
# 319   2  50 
sample_info <- sample_info[substring(sample_info$sample, 14, 16) != "01B", ]
table(substring(sample_info$sample,14,16))

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
  Stage                 = sample_info$ajcc_pathologic_stage,
  M_Stage               = sample_info$ajcc_pathologic_m,
  N_Stage               = sample_info$ajcc_pathologic_n,
  T_Stage               = sample_info$ajcc_pathologic_t,
  Gender                = sample_info$gender,
  Age                   = sample_info$age_at_index,
  stringsAsFactors      = FALSE
)
cat("\n分组统计：\n")
print(table(group_df$group))
cat("肿瘤样本中有生存时间的样本数：",
    sum(!is.na(group_df$os_time) & group_df$group == "Tumor"), "\n")
# 肿瘤样本中有生存时间的样本数： 318
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
##GSE144269----------------------------------------------------------
rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/22_job_0316/00_rawdata")
library(GEOquery)
library(tidyverse)
library(limma)
library(data.table)
#获取测试集数据
geo_number = "GSE144269"
geo_data <- getGEO(geo_number, destdir = "./", getGPL = TRUE,AnnotGPL = TRUE)
clinical <- pData(geo_data[[1]])
colnames(clinical)
table(clinical$`tumor/non-tumor:ch1`)
# non-tumor     tumor 
# 70        70 
samples <- rownames(subset(clinical,clinical$`tumor/non-tumor:ch1` == "tumor"))
#表达矩阵
# expr <- exprs(geo_data[[1]]) %>% as.data.frame()
expr <-read.table("GSE144269_voom_GeneExpression.txt",header = TRUE)
expr <- expr %>%
  mutate(gene_symbol = sub("^.*\\|", "", entrez_id)) %>%
  relocate(gene_symbol, .after = entrez_id)
# 1. 计算每一行的均值（排除非数值列）
row_means <- rowMeans(expr[, !(names(expr) %in% c("entrez_id", "gene_symbol"))], na.rm = TRUE)
# 2. 将均值作为新列，并按降序排列
expr$row_mean <- row_means
expr_sorted <- expr[order(expr$row_mean, decreasing = TRUE), ]
# 3. 按 gene_symbol 去重，保留每个 symbol 的第一行
expr_dedup <- expr_sorted[!duplicated(expr_sorted$gene_symbol), ]
# 4. 删除辅助列 row_mean
expr_dedup$row_mean <- NULL
expr_dedup$entrez_id <- NULL
expr <- expr_dedup
rownames(expr) <- expr$gene_symbol
expr$gene_symbol <- NULL
qx <- as.numeric(quantile(expr, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6] - qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) {
  # 确保 expr 是矩阵（如果原本是数据框）
  expr <- as.matrix(expr)
  expr[expr <= 0] <- 0          # 替换所有 ≤0 的值
  expr <- log2(expr + 1)        # log2变换
  print("log2 transform finished")
} else {
  print("log2 transform not needed")
}
exp <- as.data.frame(expr)
#注释探针
# GPL<-fData(geo_data[[1]])	## 获取平台信息
# gpl<-GPL[,c(1,3)]
# colnames(gpl) <- c("ID","Gene symbol") 
# gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]
# exp$ID <- rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息
# exp_symbol<-merge(exp,gpl,by="ID")
# exp_symbol<-na.omit(exp_symbol)
# table(duplicated(exp_symbol$`ID`))
# matrix0 <- order(rowMeans(exp_symbol[,c(-1,-ncol(exp_symbol))]),decreasing = T)
# exp_symbol=exp_symbol[matrix0,]
# keep=!duplicated(exp_symbol$`Gene symbol`)
# exp_symbol=exp_symbol[keep,]
# rownames(exp_symbol) <- exp_symbol$`Gene symbol`
# exp_unique<-exp_symbol[,c(-1,-ncol(exp_symbol))] %>% as.data.frame()
# exp_final <- exp_unique[,samples]
#生存信息
clinical1 <-read.xlsx("41467_2020_18186_MOESM22_ESM.xlsx")
clinical1 <- clinical1[!is.na(clinical1$RNASeq_T), ]
rownames(clinical1) <- clinical1$RNASeq_T
os <- clinical1[intersect(clinical1$RNASeq_T,colnames(exp)),c(19,20)]
colnames(os) <- c("os_status","os_time")
exp_final <- exp[,rownames(os)]
write.table(os,"GSE144269_os.txt",sep = "\t",quote = FALSE,row.names = TRUE)
saveRDS(os, file = "GSE144269_os.Rds")
#保存
write.table(exp_final,"GSE144269_expr.txt",quote = FALSE,sep = "\t")
write.csv(exp_final,"GSE144269_expr.csv",row.names = FALSE,quote = FALSE)
saveRDS(exp_final, file = "GSE144269_expr.Rds")
##scRNA--------------------------------
