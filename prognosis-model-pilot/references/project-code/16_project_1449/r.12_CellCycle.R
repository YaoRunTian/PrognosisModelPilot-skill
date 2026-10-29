rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/16_project_1449")
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
scRNA <- readRDS("../RData/R.07_scRNA_Scissor.RData")
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
#存储虚拟敲除的文件
af <- marrow[,marrow$celltype %in% "Epithelial cells"]
saveRDS(af,"../RData/pbmc_for_scTenifoldKnk.rds")
#柱状图
sc_dataset <- marrow
table(sc_dataset@meta.data$subcelltype)
table(sc_dataset@meta.data$Phase)

meta_data <- sc_dataset@meta.data
scissor_counts <- table(meta_data$subcelltype, meta_data$Phase)
# 转换为数据框
scissor_df <- as.data.frame(scissor_counts)
colnames(scissor_df) <- c("subcelltype", "Phase", "count")
# 计算百分比
scissor_percent <- scissor_df %>%
  group_by(Phase) %>%
  mutate(percent = count / sum(count) * 100) %>%
  ungroup()
# 设置颜色方案
n_celltypes <- length(unique(scissor_percent$subcelltype))
mycol <- colorRampPalette(brewer.pal(8, "Pastel1"))(n_celltypes)
# 绘制堆叠百分比柱状图
ggplot(scissor_percent, aes(x = Phase, y = percent, fill = subcelltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = mycol) +
  labs(
    title = "Subcelltype Distribution across Phase Classifications",
    x = "Phase Classification",
    y = "Percentage",
    fill = "Subcelltype"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 15),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )+
  scale_y_continuous(expand = c(0,0))+
  ggfonttheme
ggsave("02.barplot.pdf",device = cairo_pdf,width = 12, height = 8)
ggsave("02.barplot.png",dpi = 300,width = 12, height = 8)
