rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/24_job_0780")
if (!dir.exists("12_CellCycle")) {dir.create("12_CellCycle")}
setwd("12_CellCycle")
library(Seurat)
library(SeuratData)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(svglite)
library(openxlsx)
library(RColorBrewer)
ggfonttheme <- theme(text = element_text(family = "Times"))
scRNA <- readRDS("../RData/R.06_scRNA_done.Rds")
scRNA$celltype <- factor(scRNA$celltype)
pbmc <- subset(scRNA, celltype == "Fibroblasts")

# 获取S期marker基因 
s.genes <- cc.genes$s.genes 
# 获取G2M期marker基因 
g2m.genes <- cc.genes$g2m.genes

# 细胞周期评分（注意：该步骤无p值阈值，报告相关表述需删除）
fibro_cc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# 在细胞周期基因上运行PCA
fibro_cc <- RunPCA(fibro_cc, features = c(s.genes, g2m.genes))

# 可视化：显示细胞周期阶段分布
ncluster <- length(unique(fibro_cc[[]]$Phase))
mycol <- colorRampPalette(brewer.pal(8, "Set1"))(ncluster)

p <- DimPlot(fibro_cc, reduction = "umap", label = TRUE, repel = TRUE,
             label.size = 6,
             label.box = TRUE,
             cols = mycol,
             pt.size = 0.9,
             alpha = 0.3) +
  ggfonttheme +
  coord_fixed(ratio = 1)   # 保留图例，不再使用 NoLegend()
p
ggsave("01.CellCycle.pdf", plot = p, device = cairo_pdf, width = 8, height = 8)
ggsave("01.CellCycle.png", plot = p, dpi = 300, width = 8, height = 8)

save.image("../RData/12_CellCycle.RData")