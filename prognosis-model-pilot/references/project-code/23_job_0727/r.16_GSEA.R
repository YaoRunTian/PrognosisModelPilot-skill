rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/23_job_0727")
if (!dir.exists("16_GSEA")) {dir.create("16_GSEA")}
setwd("16_GSEA")
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
##预后基因GSEA--------------------------------------------------------
gmt1 <- read.gmt("h.all.v2025.1.Hs.symbols.gmt")
fpkm <- readRDS("../00_rawdata/GSE62452_expr.Rds")
fpkm <- t(fpkm)
genes <- read.table("../06_model/06.lassogeneslist.txt")[,1]
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
# [1] "FGL2"
# preparing geneSet collections...
# GSEA analysis...
# leading edge analysis...
# done...
# [1] "富集到 28 条"
# [1] "HGF"
# preparing geneSet collections...
# GSEA analysis...
# leading edge analysis...
# done...
# [1] "富集到 31 条"
# [1] "ADM"
# preparing geneSet collections...
# GSEA analysis...
# leading edge analysis...
# done...
# [1] "富集到 40 条"
# [1] "SLC2A1"
# preparing geneSet collections...
# GSEA analysis...
# leading edge analysis...
# done...
# [1] "富集到 38 条"
# There were 50 or more warnings (use warnings() to see the first 50)

save.image("../RData/16_GSEA.RData")
