rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/24_job_0780")
if (! dir.exists("00_rawdata")){dir.create("00_rawdata")}
setwd("00_rawdata")
#读取下载的表达矩阵
library(readr)
TCGA_rawdata <- read_tsv("TCGA-UCEC.star_fpkm.tsv.gz")
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
survival <- read_tsv("TCGA-UCEC.survival.tsv.gz") %>% as.data.frame()
table(substring(survival$sample,14,16))
# 01A 01B 01C 02A 11A 
# 538   8   2   1  35 
group <- ifelse(substring(survival$sample,14,16) == "01A", "Tumor",
                ifelse(substring(survival$sample,14,16) == "11A", "Normal", NA))

survival$group <- group
rownames(survival) <- survival$sample
survival <- survival[samples,] %>% na.omit()
table(substring(survival$sample,14,16))
# 01A 11A 
# 535  23 
table(survival$group)
# Normal  Tumor 
# 23    535
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
# [1] 558
expr <- expr[,intersect(survival$sample,colnames(expr))]
write.table(expr,"01.fpkm.txt",sep = "\t",quote = F,row.names = T)
write.table(survival,"02.survival.txt",sep = "\t",quote = F,row.names = T)
#count数据
count <-read_tsv("TCGA-UCEC.star_counts.tsv.gz")
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