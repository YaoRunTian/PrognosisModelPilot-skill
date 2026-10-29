rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/19_project_1533")
if (!dir.exists("21_CellCycle")) {dir.create("21_CellCycle")}
setwd("21_CellCycle")
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
pbmc <- subset(scRNA,celltype == "Epithelial cells")
# 获取S期marker基因 
s.genes <- cc.genes$s.genes 
# 获取G2M期marker基因 
g2m.genes <- cc.genes$g2m.genes
marrow <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by phase 
marrow <- RunPCA(marrow, features = c(s.genes, g2m.genes)) 
# 数据可视化，可以看到细胞按不同的细胞周期进行了分群 
ncluster <- length(unique(marrow[[]]$Phase))
mycol <- colorRampPalette(brewer.pal(8, "Set1"))(ncluster)
DimPlot(marrow, reduction = "umap", label = T , repel = T, 
        label.size = 6,
        label.box = T,
        cols = mycol,
        pt.size = 0.9,
        alpha = c(0.3)) + 
  NoLegend()+ggfonttheme+coord_fixed(ratio = 1)
ggsave("01.CellCycle.pdf",device = cairo_pdf,width = 8, height = 8)
ggsave("01.CellCycle.png",dpi = 300,width = 8, height = 8)
