rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/16_project_1449")
if (!dir.exists("02_reannotated")) {dir.create("02_reannotated")}
setwd("02_reannotated")
library(Seurat)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(extrafont)
library(RColorBrewer)
library(openxlsx)
ggfonttheme <- theme(text = element_text(family = "Times"))
pbmc_all <- readRDS("../RData/R.06_scRNA_done.Rds")
# 01. Macrophages -------------------------------------------------
pbmc <- subset(x = pbmc_all, idents = c("Macrophages"))
pbmc <- RunPCA(pbmc, npcs = 50, features = VariableFeatures(object = pbmc), verbose = FALSE)
pbmc <- pbmc %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.3) %>%
  identity()
p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE) + ggfonttheme
p1
ggsave("./Macrophages/01.cluster.pdf", device = cairo_pdf, width = 6, height = 6)
ggsave("./Macrophages/01.cluster.png", dpi = 300, width = 6, height = 6)

# 注释参考PMID: 39116206
cell_markers <- list(
  # 经典极化亚型
  M1_macrophages = c("TLR2","ISG15"),
  M2_macrophages = c("CD163", "MRC1")

  # # M2进一步细分亚型
  # M2a_macrophages = c("MRC1", "CD163", "HLA-DRA", "ARG1"),
  # M2b_macrophages = c("CD86", "HLA-DRA", "IL10", "IL1B", "IL6", "TNF", "CCL1"),
  # M2c_macrophages = c("CD163", "MRC1", "TLR1", "IL10", "TGFB1", "CCL18"),
  # M2d_macrophages = c("CD163", "IL10", "VEGFA", "TGFB1"),
  # 
  # Microglia = c("CX3CR1", "TMEM119", "P2RY12", "AIF1"),
  # Adipose_tissue_macrophages = c("ITGAX", "TREM2", "LYVE1"),
  # Osteoclasts = c("CALCR", "ACP5"),
  # 
  # # 肿瘤相关巨噬细胞（TAMs）
  # Tumor_associated_macrophages = c("CD163", "MRC1", "MARCO", "PDCD1LG2", "CCL2"),
  # 
  # # 其他常用标记
  # CD169_macrophages = c("SIGLEC1")   # CD169即SIGLEC1
)
cell_markers_df <- stack(cell_markers)
names(cell_markers_df) <- c("gene", "cell_type")
DotPlot(pbmc, features = split(unique(cell_markers_df$gene), cell_markers_df$cell_type),
        group.by = "seurat_clusters", cols = c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.border = element_rect(color = "black"),
        panel.spacing = unit(1, "mm")) + ggfonttheme

DotPlot(pbmc, 
        features = split(unique(cell_markers_df$gene), cell_markers_df$cell_type),
        group.by = "seurat_clusters", 
        cols = c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.border = element_rect(color = "black"),
        panel.spacing = unit(1, "mm")) +
  scale_size(range = c(0, 10))+
  ggfonttheme
ggsave("./Macrophages/02.dotplot__res0.3.pdf", device = cairo_pdf, width = 20, height = 8)
ggsave("./Macrophages/02.dotplot__res0.3.png", dpi = 300, width = 20, height = 8)

new.cluster.id <- c(
  "0" = "M2 Macrophages",
  "1" = "M1 Macrophages",
  "2" = "M2 Macrophages",
  "3" = "M1 Macrophages",
  "4" = "M2 Macrophages"
)
saveRDS(pbmc, file = "../RData/R.Macrophages_pbmc.RData")
pbmc <- readRDS("../RData/R.Macrophages_pbmc.RData")
cluster_dat <- data.frame("cluster" = 0:(length(new.cluster.id) - 1), "CellType" = new.cluster.id)
write.xlsx(cluster_dat, file = "./Macrophages/03.anno_CellType.xlsx")
names(new.cluster.id) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.id)
pbmc$celltype <- Idents(pbmc)

ncluster <- length(unique(pbmc[[]]$celltype))
mycol <- colorRampPalette(brewer.pal(2, "Pastel1"))(ncluster)
p10 <- DimPlot(pbmc, reduction = "umap", label = T, repel = T, 
               label.size = 6,
               label.box = T,
               cols = mycol,
               pt.size = 0.9,
               alpha = c(0.3)) + 
  NoLegend() + ggfonttheme + coord_fixed(ratio = 1)
p10
ggsave("./Macrophages/04.umap_annotated.pdf", p10, device = cairo_pdf, width = 10, height = 10)
ggsave("./Macrophages/04.umap_annotated.png", p10, dpi = 300, width = 10, height = 10)

# 保存巨噬细胞的重新注释信息
macrophages_reannotated <- pbmc
saveRDS(macrophages_reannotated, file = "./Macrophages/Macrophages_reannotated.RData")
rm(pbmc); gc()

# 02. Epithelial cells ----------------------------------------
pbmc <- subset(x = pbmc_all, idents = c("Epithelial cells"))
pbmc <- RunPCA(pbmc, npcs = 50, features = VariableFeatures(object = pbmc), verbose = FALSE)
pbmc <- pbmc %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.8) %>%
  identity()
p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE) + ggfonttheme
p1
ggsave("./Epithelial cells/01.cluster.pdf", device = cairo_pdf, width = 6, height = 6)
ggsave("./Epithelial cells/01.cluster.png", dpi = 300, width = 6, height = 6)

# 注释参考（https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(25)00469-4/fulltext）
cell_markers <- list(
  Squamous_cell = c("SCEL", "SPRR1A", "SPRR1B", "KRT13"),
  Basal_cell = c("KRT5", "KRT15", "KRT16", "KRT6A", "KRT7"),
  Ciliated_cell = c("FOXJ1", "PIFO", "TPPP3", "RSPH1"),
  Submucosal_gland_cell = c("LTF", "SLPI", "LYZ", "MUC5B", "SOX2", "SOX9"),
  Goblet_cell = c("MUC5AC", "MUC5B", "BPIFB1", "KLF4", "HES6"),
  Secretory_cell = c("SCGB3A2", "SCGB3A1", "BPIFB1", "WFDC2", "SERPINB3"),
  Tuft_cell = c("ASCL2", "DCLK1", "LRMP", "CD24A", "POU2F3")
)
cell_markers_df <- stack(cell_markers)
names(cell_markers_df) <- c("gene", "cell_type")
DotPlot(pbmc, features = split(unique(cell_markers_df$gene), cell_markers_df$cell_type),
        group.by = "seurat_clusters", cols = c("lightgrey", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.border = element_rect(color = "black"),
        panel.spacing = unit(1, "mm")) + ggfonttheme
ggsave("./Epithelial cells/02.dotplot.pdf", device = cairo_pdf, width = 20, height = 8)
ggsave("./Epithelial cells/02.dotplot.png", dpi = 300, width = 20, height = 8)

new.cluster.id <- c(
  "0" = "Basal_cell",
  "1" = "Submucosal_gland_cell",
  "2" = "Basal_cell",
  "3" = "Submucosal_gland_cell",
  "4" = "Secretory_cell",
  "5" = "Squamous_cell",
  "6" = "Submucosal_gland_cell",
  "7" = "Submucosal_gland_cell",
  "8" = "Secretory_cell",
  "9" = "Squamous_cell",
  "10" = "Goblet_cell",
  "11" = "Basal_cell",
  "12" = "Submucosal_gland_cell",
  "13" = "Tuft_cell",
  "14" = "Basal_cell",
  "15" = "Basal_cell",
  "16" = "Basal_cell",
  "17" = "Squamous_cell",
  "18" = "Ciliated_cell",
  "19" = "Secretory_cell",
  "20" = "Basal_cell",
  "21" = "Squamous_cell",
  "22" = "Submucosal_gland_cell",
  "23" = "Basal_cell",
  "24" = "Goblet_cell",
  "25" = "Secretory_cell",
  "26" = "Basal_cell"
)
saveRDS(pbmc, file = "../RData/R.Epithelial_pbmc.RData")
pbmc <- readRDS("../RData/R.Epithelial_pbmc.RData")
cluster_dat <- data.frame("cluster" = 0:(length(new.cluster.id) - 1), "CellType" = new.cluster.id)
write.xlsx(cluster_dat, file = "./Epithelial cells/03.anno_CellType.xlsx")
names(new.cluster.id) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.id)
pbmc$celltype <- Idents(pbmc)

ncluster <- length(unique(pbmc[[]]$celltype))
mycol <- colorRampPalette(brewer.pal(8, "Pastel1"))(ncluster)
p10 <- DimPlot(pbmc, reduction = "umap", label = T, repel = T, 
               label.size = 6,
               label.box = T,
               cols = mycol,
               pt.size = 0.9,
               alpha = c(0.3)) + 
  NoLegend() + ggfonttheme + coord_fixed(ratio = 1)
p10
ggsave("./Epithelial cells/04.umap_annotated.pdf", p10, device = cairo_pdf, width = 10, height = 10)
ggsave("./Epithelial cells/04.umap_annotated.png", p10, dpi = 300, width = 10, height = 10)

# 保存上皮细胞的重新注释信息
epithelial_reannotated <- pbmc
saveRDS(epithelial_reannotated, file = "./Epithelial cells/Epithelial_cells_reannotated.RData")
rm(pbmc); gc()

# 03. 合并重新注释信息到pbmc_all中 -------------------------------------------
# 首先在pbmc_all的metadata中创建subcelltype列，初始值为原来的celltype
# 将celltype转换为字符向量，然后赋值给subcelltype
pbmc_all$subcelltype <- as.character(pbmc_all$celltype)

# 更新巨噬细胞的亚型注释
macrophage_cells <- Cells(macrophages_reannotated)
macrophage_subtypes <- as.character(macrophages_reannotated@meta.data[macrophage_cells, "celltype"])
pbmc_all$subcelltype[macrophage_cells] <- macrophage_subtypes

# 更新上皮细胞的亚型注释
epithelial_cells <- Cells(epithelial_reannotated)
epithelial_subtypes <- as.character(epithelial_reannotated@meta.data[epithelial_cells, "celltype"])
pbmc_all$subcelltype[epithelial_cells] <- epithelial_subtypes

# 验证更新结果
cat("原始细胞类型与亚型注释对比:\n")
print(table(Original = pbmc_all$celltype, Subtype = pbmc_all$subcelltype, useNA = "ifany"))

# 统计重新注释的细胞数量
cat("\n重新注释的细胞统计:\n")
cat("巨噬细胞重新注释数量:", length(macrophage_cells), "\n")
cat("上皮细胞重新注释数量:", length(epithelial_cells), "\n")

# 保存更新后的pbmc_all对象
saveRDS(pbmc_all, file = "../RData/pbmc_all_with_subcelltype.Rds")

# 可视化重新注释后的细胞类型分布
pbmc_all <- SetIdent(pbmc_all, value = "subcelltype")
ncluster <- length(unique(pbmc_all[[]]$subcelltype))
mycol <- colorRampPalette(brewer.pal(8, "Pastel1"))(ncluster)
p_combined <- DimPlot(pbmc_all, reduction = "umap", label = T, repel = T,
                      label.size = 6,
                      label.box = T,
                      cols = mycol,
                      pt.size = 0.9,
                      alpha = c(0.3))+NoLegend()+
  ggfonttheme+coord_fixed(ratio = 1)+ggfonttheme
p_combined$layers[[2]]$geom_params$max.overlaps <- 100
p_combined

ggsave("01.umap_all_with_subcelltype.pdf", 
       p_combined, device = cairo_pdf, width = 15, height = 15)
ggsave("01.umap_all_with_subcelltype.png", 
       p_combined, dpi = 300, width = 15, height = 15)

# 制作注释对比表格
annotation_comparison <- data.frame(
  Cell = rownames(pbmc_all@meta.data),
  Original_Annotation = pbmc_all$celltype,
  Subtype_Annotation = pbmc_all$subcelltype,
  Reannotated = rownames(pbmc_all@meta.data) %in% c(macrophage_cells, epithelial_cells)
)
# 保存注释对比表格
write.csv(annotation_comparison, "02.annotation_comparison.csv", row.names = FALSE)
