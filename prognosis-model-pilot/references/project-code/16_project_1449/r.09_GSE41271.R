rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/16_project_1449")
if (!dir.exists("09_GSE41271")) {dir.create("09_GSE41271")}
setwd("09_GSE41271")
library(GEOquery)
library(tidyverse)
library(limma)
library(data.table)
#获取测试集数据
geo_number = "GSE41271"
geo_data <- getGEO(geo_number, destdir = "./", getGPL = TRUE,AnnotGPL = TRUE)
clinical <- pData(geo_data[[1]])
colnames(clinical)
table(clinical$characteristics_ch1)
index1 <- grep("SCLC & NSCLC",clinical$characteristics_ch1)
clinical <- clinical[-index1,]
index2 <- grep("LCC-NE",clinical$characteristics_ch1)
clinical <- clinical[-index2,]
os <- clinical[,c(45,50,54)]
os[,1] <- as.Date(os[,1])
os[,2] <- as.Date(os[,2])
os$os_status <- as.numeric(os[,2] - os[,1])
os <- os[,c(4,3)]
colnames(os) <- c("os_time","os_status")
dead_index <- grep("D",os$os_status)
live_index <- grep("A",os$os_status)
os[dead_index,"os_status"] <-1
os[live_index,"os_status"] <-0
os$os_time <- as.numeric(os$os_time)
os$os_status <- as.numeric(os$os_status)
os <- na.omit(os)#270
saveRDS(os,file = "os.Rds")
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
#已归一化
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
exp_unique
write.table(exp_unique,"02.expr.txt",quote = FALSE,sep = "\t")
write.csv(exp_unique,"02.expr.csv",row.names = FALSE,quote = FALSE)
saveRDS(exp_unique, file = "expr.Rds")

save.image("../RData/09_GSE41271.RData")

