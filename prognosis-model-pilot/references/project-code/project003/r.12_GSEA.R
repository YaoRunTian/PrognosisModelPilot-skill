rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project003")
if (!dir.exists("12_GSEA")) {dir.create("12_GSEA")}
setwd("12_GSEA")
library(org.Hs.eg.db)
library(clusterProfiler)
library(cowplot)
library(dplyr)
library(ggplot2)
library(enrichplot)
library(limma)
library(metagenomeSeq)
library(edgeR)
library(dplyr)
library(data.table)
library(tidyverse)
library(msigdbr)
library(DESeq2)
ggfonttheme <- theme(text = element_text(family = "Times"))
###高低风险组GSEA--------------------------------------------------
#00.参考基因集------------------------------------------------------------------
gmt1 <- read.gmt("h.all.v2025.1.Hs.symbols.gmt")
#01.差异分析
risk_table <-  read.table("../04_evaluate/04.risk_table.txt",sep = "\t",header = TRUE)#读取高低风险组信息
count <- readRDS("../00_rawdata/count.Rds")
#01.差异分析------------------------------------------------------------------
expr <- count[,rownames(risk_table)]#高低表达组样本
expr <- expr[rowMeans(expr) > 1, ]  # 保留行均值大于1的行（排除不重要的极值基因）
write.csv(expr,'01.cleaned_expr.csv',row.names = TRUE)#重新写入整理后的矩阵文件
#分组信息
groupInf <- data.frame(sample = rownames(risk_table),group = risk_table$risk_list_cut)
groupInf$group <- factor(groupInf$group)
Data <- groupInf
group <- groupInf$group
dds <- DESeqDataSetFromMatrix(countData = expr,
                              colData = Data,
                              design = ~ group)
dds2 <- DESeq(dds)
res <- results(dds2, contrast=c("group", "high", "low"))
res <- res[order(res$pvalue),]#按P值从小到大排序
summary(res)#检查
my_result <- as.data.frame(res)#转成容易查看的数据框
my_result <- na.omit(my_result)#删除倍数为0的值
#保存差异分析的结果
my_result$Gene_symbol<-rownames(my_result)
my_result <- my_result %>% dplyr::select('Gene_symbol',colnames(my_result)[1:dim(my_result)[2]-1],everything())
rownames(my_result) <- NULL
write.csv(my_result,file="02.my_result_deseq2.csv")#写入
#筛选差异表达基因
DEG<-my_result
DEG$regulate <- ifelse(DEG$pvalue > 0.05, "unchanged",
                       ifelse(DEG$log2FoldChange > 1, "up-regulated",
                              ifelse(DEG$log2FoldChange < -1, "down-regulated", "unchanged")))
table(DEG$regulate)
# down-regulated      unchanged   up-regulated 
# 324          16697            322 
df <- DEG[grep("up-regulated|down-regulated",DEG$regulate),]
df <- df[,c("Gene_symbol","log2FoldChange")]
FCgenelist <- df$log2FoldChange
names(FCgenelist) <- df$Gene_symbol
FCgenelist = sort(FCgenelist,decreasing = T) 
head(FCgenelist)
DEGsList <- df$Gene_symbol
write.table(DEGsList,"03.DEGsList.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)
#clusterprofiler-------------------------------
set.seed(42)
em_risk<- GSEA(FCgenelist, TERM2GENE = gmt1,pAdjustMethod = "none")
em_risk@result <- em_risk@result %>%
  filter(p.adjust < 0.05 & abs(NES) > 1 & qvalue < 0.25)
print(paste("富集到",dim(em_risk@result)[1],"条"))
write.csv(em_risk@result,file="04.riskgroup_gsea_res.csv")#2条
p <- gseaplot2(em_risk,1:2, pvalue_table = FALSE,ES_geom = "line",
               base_size = 15,title = "",color = ggsci::pal_npg()(2))
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical",
                       legend.text = element_text(size = 15))+ggfonttheme
p[[2]] <- p[[2]]+ggfonttheme
p[[3]] <- p[[3]]+ggfonttheme
p
ggsave("05.riskgroup_gsea.pdf",p,device = cairo_pdf,width = 12, height = 8)
ggsave("05.riskgroup_gsea.png",p,dpi = 300,width = 12, height = 8)


##预后基因GSEA--------------------------------------------------------
rm(list = ls()); gc()
ggfonttheme <- theme(text = element_text(family = "Times"))
gmt1 <- read.gmt("h.all.v2025.1.Hs.symbols.gmt")
fpkm <- readRDS("../00_rawdata/tpm.Rds")
fpkm <- t(fpkm)
genes <- read.table("../03_model/10.multiresult.txt") %>% rownames()
for(gene in genes){
  print(gene)
  #计算相关性
  gene_expression <- fpkm[, gene]
  other_genes_matrix <- fpkm[, setdiff(colnames(fpkm),gene)]
  # 计算所有相关性
  correlations <- apply(other_genes_matrix, 2, function(x) {
    cor(gene_expression, x, method = "spearman")
  })
  # 计算p值
  cor_p_values <- apply(other_genes_matrix, 2, function(x) {
    cor_test <- cor.test(gene_expression, x, method = "spearman")
    cor_test$p.value
  })
  # 创建结果数据框
  efficient_results <- data.frame(
    gene = names(correlations),
    correlation = correlations,
    p_value = cor_p_values,
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(abs(correlation)))
  head(efficient_results, 20)
  
  FCgenelist <- efficient_results$correlation
  names(FCgenelist) <- efficient_results$gene
  FCgenelist = sort(FCgenelist,decreasing = T)
  
  # df1 <- read.xlsx("../04_DEG2/04.DEG.xlsx")
  # FCgenelist <- df1$log2FoldChange
  # names(FCgenelist) <- df1$Gene_symbol
  # FCgenelist = sort(FCgenelist,decreasing = T) 
  # head(FCgenelist)
  set.seed(42)
  em_risk<- GSEA(FCgenelist, TERM2GENE = gmt1)
  em_risk@result <- em_risk@result %>%
    filter(p.adjust < 0.05 & abs(NES) > 1 & qvalue < 0.25)
  print(paste("富集到",dim(em_risk@result)[1],"条"))
  csvfilename <- paste0(gene,"_gsea_res.csv")
  write.csv(em_risk@result,file=csvfilename)
  
  p <- gseaplot2(em_risk,1:10, pvalue_table = FALSE,ES_geom = "line",
                 base_size = 15,title = gene,color = ggsci::pal_npg()(10))
  p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical",
                         legend.text = element_text(size = 15))+ggfonttheme
  p[[2]] <- p[[2]]+ggfonttheme
  p[[3]] <- p[[3]]+ggfonttheme
  p
  pdf_filename <- paste0(gene,"_gsea.pdf")
  png_filename <- paste0(gene,"_gsea.png")
  ggsave(pdf_filename,p,device = cairo_pdf,width = 12, height = 8)
  ggsave(png_filename,p,dpi = 300,width = 12, height = 8)
}
# [1] "FCN3"
# preparing geneSet collections...
# GSEA analysis...
# leading edge analysis...
# done...
# [1] "富集到 35 条"
# [1] "PGLYRP4"
# preparing geneSet collections...
# GSEA analysis...
# leading edge analysis...
# done...
# [1] "富集到 16 条"
# There were 50 or more warnings (use warnings() to see the first 50)

save.image("../RData/12_GSEA.RData")

