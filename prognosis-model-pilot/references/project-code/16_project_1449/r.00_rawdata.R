rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/16_project_1449")
if (!dir.exists("00_rawdata")) {dir.create("00_rawdata")}
setwd("00_rawdata")
#single cell data--------------------------------
library(Seurat)
library(data.table)
library(Matrix)
# 获取所有文件
files <- list.files("GSE148071_RAW/", 
                    pattern = "\\.txt\\.gz$", 
                    full.names = TRUE)
# 函数：读取单个表达矩阵文件并添加样本信息
read_expression_matrix <- function(file) {
  cat("Processing:", basename(file), "\n")
  # 读取数据（假设第一列是基因名）
  dt <- fread(file)
  # 获取基因名（第一列）
  gene_names <- dt[[1]]
  # 获取表达矩阵（去除第一列）
  expr_matrix <- as.matrix(dt[, -1])
  rownames(expr_matrix) <- gene_names
  # 设置细胞名（列名）
  sample_id <- gsub("\\.txt\\.gz$", "", basename(file))
  # 如果文件没有合适的列名，用文件名+序号
  if (all(colnames(dt)[-1] == paste0("V", 2:ncol(dt)))) {
    colnames(expr_matrix) <- paste0(sample_id, "_Cell_", 1:ncol(expr_matrix))
  }
  cat("Dimensions:", nrow(expr_matrix), "genes ×", ncol(expr_matrix), "cells\n")
  # 返回矩阵和样本ID
  return(list(matrix = expr_matrix, sample_id = sample_id))
}
# 读取所有文件并处理
matrices <- list()
sample_ids <- character()
for (i in seq_along(files)) {
  result <- read_expression_matrix(files[i])
  matrices[[i]] <- result$matrix
  sample_ids <- c(sample_ids, rep(result$sample_id, ncol(result$matrix)))
}
# 检查所有矩阵的基因是否一致
all_genes <- lapply(matrices, rownames)
common_genes <- Reduce(intersect, all_genes)
# 只保留共有基因
matrices_filtered <- lapply(matrices, function(mat) {
  mat[common_genes, ]
})
# 合并所有矩阵
combined_matrix <- do.call(cbind, matrices_filtered)
# 创建Seurat对象
seurat_obj <- CreateSeuratObject(
  counts = combined_matrix,
  project = "GSE148071",
  min.cells = 3,
  min.features = 200
)
# 添加样本信息到meta.data
seurat_obj$sample_id <- sample_ids
# 验证添加的信息
cat("\n样本信息已添加到meta.data:\n")
print(table(seurat_obj$sample_id))
print(seurat_obj)
saveRDS(seurat_obj, "GSE148071_pbmc.Rds")
#GSE108214-------------------
library(GEOquery)
library(tidyverse)
library(limma)
library(data.table)
geo_number = "GSE108214"
geo_data <- getGEO(geo_number, destdir = "./", getGPL = TRUE,AnnotGPL = TRUE)
clinical <- pData(geo_data[[1]])
colnames(clinical)
table(clinical$source_name_ch1)
groupInf <- data.frame("Sample" = rownames(clinical),"group" = clinical$source_name_ch1)
index <- grep("control",groupInf$group)
groupInf <- groupInf[index,]
res_index <- grep("resistant",groupInf$group)
sen_index <- grep("sensitive",groupInf$group)
groupInf[res_index,"group"] <- 0
groupInf[sen_index,"group"] <- 1
saveRDS(groupInf, file = "GSE108214_group.Rds")
#expr
expr <- exprs(geo_data[[1]]) %>% as.data.frame()
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
exp <- exp[,groupInf$Sample]
#注释探针
GPL<-fData(geo_data[[1]])	## 获取平台信息
gpl<-GPL[,c(1,7)]
colnames(gpl) <- c("ID","GENE_SYMBOL") 
gpl$`GENE_SYMBOL`<-data.frame(sapply(gpl$`GENE_SYMBOL`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]
exp$ID <- rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息
exp_symbol<-merge(exp,gpl,by="ID")
exp_symbol<-na.omit(exp_symbol)
table(duplicated(exp_symbol$`ID`))
matrix0 <- order(rowMeans(exp_symbol[,c(-1,-ncol(exp_symbol))]),decreasing = T)
exp_symbol=exp_symbol[matrix0,]
keep=!duplicated(exp_symbol$`GENE_SYMBOL`)
exp_symbol=exp_symbol[keep,]
rownames(exp_symbol) <- exp_symbol$`GENE_SYMBOL`
exp_unique<-exp_symbol[,c(-1,-ncol(exp_symbol))] %>% as.data.frame()
exp_unique
write.table(exp_unique,"02.GSE108214_expr.txt",quote = FALSE,sep = "\t")
write.csv(exp_unique,"02.GSE108214_expr.csv",row.names = FALSE,quote = FALSE)
saveRDS(exp_unique, file = "GSE108214_expr.Rds")

#TCGA----------------------------------
LUAD_rawdata <- read_tsv("TCGA-LUAD.star_fpkm.tsv.gz")
LUSC_rawdata <- read_tsv("TCGA-LUSC.star_fpkm.tsv.gz")
TCGA_rawdata <- merge(LUAD_rawdata,LUSC_rawdata,by = 1)
#读取注释文件并转换为dataframe
library(dplyr)
library(rtracklayer)
gtf_v36 <- import('gencode.v36.annotation.gtf.gz') %>% as.data.frame()
#提取mRNA的基因
gtf_mrna_v36 <- subset(gtf_v36,gene_type == "protein_coding")
gtf_mrna_v36 <- dplyr::select(gtf_mrna_v36,c("gene_id", "gene_name")) %>% unique()
gtf_mrna_v36 <- gtf_mrna_v36[-11238,]
#根据注释文件，转换gene_ID
exprSet <-TCGA_rawdata %>%
  inner_join(gtf_mrna_v36, by = c("Ensembl_ID" = "gene_id")) %>%
  dplyr::select(gene_name, starts_with("TCGA") )
samples <- colnames(exprSet)[-1]
#生存信息
LUAD_survival <- read_tsv("TCGA-LUAD.survival.tsv.gz") %>% as.data.frame()
LUSC_survival <- read_tsv("TCGA-LUSC.survival.tsv.gz") %>% as.data.frame()
survival <- rbind(LUAD_survival,LUSC_survival)

table(substring(survival$sample,14,16))

group <- ifelse(substring(survival$sample,14,16) == "01A", "Tumor",
                ifelse(substring(survival$sample,14,16) == "11A", "Normal", NA))
survival$group <- group
rownames(survival) <- survival$sample
survival <- survival[samples,] %>% na.omit()
table(survival$group)

groupInf <- survival[,c("sample","group")]

#存储数据
expr <- as.data.frame(exprSet)
matrix0 <- order(rowMeans(expr[,-1]),decreasing = T)
expr=expr[matrix0,]
keep=!duplicated(expr$gene)
expr=expr[keep,]
rownames(expr) <- expr$gene
colnames(expr) <- gsub("-","_",colnames(expr))
rownames(expr) <- gsub("-","_",rownames(expr))
survival$sample <- gsub("-","_",survival$sample)
length(intersect(survival$sample,colnames(expr)))
#1098
expr <- expr[,intersect(survival$sample,colnames(expr))]
write.table(expr,"01.fpkm.txt",sep = "\t",quote = F,row.names = T)
write.table(survival,"02.survival.txt",sep = "\t",quote = F,row.names = T)
#count数据
LUSC_count <- read_tsv("TCGA-LUSC.star_counts.tsv.gz")
LUAD_count <- read_tsv("TCGA-LUAD.star_counts.tsv.gz")
count <- merge(LUSC_count,LUAD_count,by = 1)
#根据注释文件，转换gene_ID
exprSet <- count %>%
  inner_join(gtf_mrna_v36, by = c("Ensembl_ID" = "gene_id")) %>%
  dplyr::select(gene_name, starts_with("TCGA") )
#存储数据
count <- as.data.frame(exprSet)
matrix0 <- order(rowMeans(count[,-1]),decreasing = T)
count=count[matrix0,]
keep=!duplicated(count$gene)
count=count[keep,]
rownames(count) <- count$gene
colnames(count) <- gsub("-","_",colnames(count))
rownames(count) <- gsub("-","_",rownames(count))
count <- count[,intersect(survival$sample,colnames(count))]
write.table(count,"03.counts.txt",sep = "\t",quote = F,row.names = T)
saveRDS(expr, file = "fpkm.Rds")
saveRDS(survival,file = "os.Rds")
saveRDS(count,file = "count.Rds")
saveRDS(groupInf,file = "group.Rds")
