rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project004")
if (! dir.exists("00_rawdata")){dir.create("00_rawdata")}
setwd("00_rawdata")
#读取下载的表达矩阵
library(readr)
TCGA_rawdata <- read_tsv("TCGA-PAAD.star_fpkm.tsv.gz")
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
survival <- read_tsv("TCGA-PAAD.survival.tsv.gz") %>% as.data.frame()
table(substring(survival$sample,14,16))
# 01A 06A 11A 
# 184   1  10 
group <- ifelse(substring(survival$sample,14,16) == "01A", "Tumor",
                ifelse(substring(survival$sample,14,16) == "11A", "Normal", NA))

survival$group <- group
rownames(survival) <- survival$sample
survival <- survival[samples,] %>% na.omit()
table(substring(survival$sample,14,16))
# 01A 11A 
# 177   4 
table(survival$group)
# Normal  Tumor 
# 4    177 
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
count <-read_tsv("TCGA-PAAD.star_counts.tsv.gz")
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

##GSE62452----------------------------------------------------------
rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project004/00_rawdata")
library(GEOquery)
library(tidyverse)
library(limma)
library(data.table)
#获取测试集数据
geo_number = "GSE62452"
geo_data <- getGEO(geo_number, destdir = "./", getGPL = TRUE,AnnotGPL = TRUE)
clinical <- pData(geo_data[[1]])
colnames(clinical)
table(clinical$`tissue:ch1`)
# adjacent pancreatic non-tumor              Pancreatic tumor 
# 61                            69 
samples <- rownames(subset(clinical,clinical$`tissue:ch1` == "Pancreatic tumor"))
group <- data.frame(
  "Sample" = rownames(clinical),
  "Group" = clinical$`tissue:ch1`
)
table(group$Group)
grep("adjacent",group$Group)
group$Group[grep("adjacent",group$Group)] <- "Normal"
group$Group[grep("tumor",group$Group)] <- "Tumor"
table(group$Group)
# Normal  Tumor 
# 61     69 
saveRDS(group,file = "GSE62452_group.Rds")
#生存信息
os <- clinical[,c(40,39)]
colnames(os) <- c("os_status","os_time")
# os$os_status[os$os_status == "survival"] <- 0
# os$os_status[os$os_status == "death"] <- 1
tumor_sample <- subset(group,Group == 'Tumor')$Sample
os <- os[rownames(os) %in%tumor_sample,]
os$os_time <- as.numeric(os$os_time)/12
os$os_status <- as.numeric(os$os_status)
os <- na.omit(os)#65个样本
saveRDS(os,file = "GSE62452_os.Rds")
#表达矩阵
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
#注释探针
GPL<-fData(geo_data[[1]])	## 获取平台信息
gpl<-GPL[,c(1,3)]
colnames(gpl) <- c("ID","Gene symbol") 
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]
exp$ID <- rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息
exp_symbol<-merge(exp,gpl,by="ID")
exp_symbol<-na.omit(exp_symbol)
table(duplicated(exp_symbol$`ID`))
matrix0 <- order(rowMeans(exp_symbol[,c(-1,-ncol(exp_symbol))]),decreasing = T)
exp_symbol=exp_symbol[matrix0,]
keep=!duplicated(exp_symbol$`Gene symbol`)
exp_symbol=exp_symbol[keep,]
rownames(exp_symbol) <- exp_symbol$`Gene symbol`
exp_unique<-exp_symbol[,c(-1,-ncol(exp_symbol))] %>% as.data.frame()
exp_final <- exp_unique[,samples]
#保存
write.table(exp_final,"02.GSE62452_expr.txt",quote = FALSE,sep = "\t")
write.csv(exp_final,"02.GSE62452_expr.csv",row.names = FALSE,quote = FALSE)
saveRDS(exp_final, file = "GSE62452_expr.Rds")

#scRNA---------------------------------------
