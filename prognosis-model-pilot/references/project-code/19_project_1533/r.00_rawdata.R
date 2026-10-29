rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/19_project_1533")
if (! dir.exists("00_rawdata")){dir.create("00_rawdata")}
setwd("00_rawdata")
#读取下载的表达矩阵
library(readr)
TCGA_rawdata <- read_tsv("TCGA-BLCA.star_fpkm.tsv.gz")
#读取注释文件并转换为dataframe
library(dplyr)
library(rtracklayer)
gtf_v36 <- import('gencode.v36.annotation.gtf.gz') %>% as.data.frame()
#提取mRNA的基因
gtf_mrna_v36 <- subset(gtf_v36,gene_type == "protein_coding")
gtf_mrna_v36 <- dplyr::select(gtf_mrna_v36,c("gene_id", "gene_name")) %>% unique()
#根据注释文件，转换gene_ID
exprSet <-TCGA_rawdata %>%
  inner_join(gtf_mrna_v36, by = c("Ensembl_ID" = "gene_id")) %>%
  dplyr::select(gene_name, starts_with("TCGA") )
samples <- colnames(exprSet)[-1]
#生存信息
survival <- read_tsv("TCGA-BLCA.survival.tsv.gz") %>% as.data.frame()
table(substring(survival$sample,14,16))
# 01A 01B 06A 11A 
# 404   6   1  22 
group <- ifelse(substring(survival$sample,14,16) == "01A", "Tumor",
                ifelse(substring(survival$sample,14,16) == "11A", "Normal", NA))

survival$group <- group
rownames(survival) <- survival$sample
survival <- survival[samples,] %>% na.omit()
table(substring(survival$sample,14,16))
# 01A 11A 
# 398  18 
table(survival$group)
# Normal  Tumor 
# 18    398 
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
# [1] 416
expr <- expr[,intersect(survival$sample,colnames(expr))]
write.table(expr,"01.fpkm.txt",sep = "\t",quote = F,row.names = T)
write.table(survival,"02.survival.txt",sep = "\t",quote = F,row.names = T)
#count数据
count <-read_tsv("TCGA-BLCA.star_counts.tsv.gz")
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
# save.image("../RData/00_rawdata.RData")


#单细胞数据读取-----------------------------
rm(list = ls()); gc()
# untar("GSE135337_RAW.tar", exdir = "GSE135337_RAW")
raw_dir <- "GSE135337_RAW"
all_files <- list.files(raw_dir, full.names = TRUE)
cat("共解压出", length(all_files), "个文件\n")
print(all_files)

library(Seurat)
library(readxl)

# 分离txt文件和excel文件
txt_files <- all_files[1:7]
xls_file <- all_files[8]

# 读取txt文件的函数
read_txt_file <- function(file_path) {
  data <- read.table(gzfile(file_path), header = TRUE, sep = "\t", row.names = 1)
  gene_symbols <- data$Symbol
  counts <- as.matrix(data[, -1])
  rownames(counts) <- make.unique(gene_symbols)
  return(counts)
}

# 读取excel文件的函数（需要先解压）
read_xls_file <- function(file_path) {
  temp_file <- tempfile()
  gunzip(filename = file_path, destname = temp_file, remove = FALSE, overwrite = TRUE)
  
  # 用read.table读取
  data <- read.table(temp_file, header = TRUE, sep = "\t", row.names = 1)
  unlink(temp_file)
  
  gene_symbols <- data$Symbol
  counts <- as.matrix(data[, -1])
  rownames(counts) <- make.unique(gene_symbols)
  
  return(counts)
}

# 创建Seurat对象列表
seurat_list <- list()

# 处理7个txt文件
for (i in 1:7) {
  cat("Processing:", basename(txt_files[i]), "\n")
  counts <- read_txt_file(txt_files[i])
  
  # 从文件名提取样本ID (BC1, BC2, etc.)
  sample_id <- gsub(".*_(BC[0-9]+)_.*", "\\1", basename(txt_files[i]))
  
  seurat_list[[sample_id]] <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )
}

# 处理excel文件
cat("Processing:", basename(xls_file), "\n")
counts_xls <- read_xls_file(xls_file)
dim(counts_xls)
seurat_list[["BCN"]] <- CreateSeuratObject(
  counts = counts_xls,
  project = "BCN",
  min.cells = 3,
  min.features = 200
)

# 查看各样本细胞数
sapply(seurat_list, ncol)
# BC1  BC2  BC3  BC4  BC5  BC6  BC7  BCN 
# 3218 5006 3100 5131 8643 6167 5523 6035
# 合并所有样本
merged_obj <- merge(
  seurat_list[[1]], 
  y = seurat_list[-1], 
  add.cell.ids = names(seurat_list),
  project = "GSE135337"
)

cat("\n合并后的对象:\n")
merged_obj

# 查看每个样本的细胞数
table(merged_obj$orig.ident)

saveRDS(merged_obj,"GSE135337_raw_pbmc.Rds")
