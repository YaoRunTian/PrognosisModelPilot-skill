rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project005")
if (!dir.exists("13_monocle")) {dir.create("13_monocle")}
setwd("13_monocle")
if (!dir.exists("Fibroblasts")) {
  dir.create("Fibroblasts", recursive = TRUE)
}

library(Seurat)
library(monocle3)
library(dplyr)
library(patchwork)

ggfonttheme <- theme(text = element_text(family = "Times"),
                     axis.title = element_text(size = 12),
                     axis.text = element_text(size = 12),
                     plot.title = element_text(size = 12))

pbmc_all <- readRDS("../RData/R.06_scRNA_done.Rds")

pbmc <- subset(x = pbmc_all, idents = c("Fibroblasts"))

# 注意：resolution 与报告统一为 0.5（报告若写 0.4 请更正为 0.5）
pbmc <- RunPCA(pbmc, npcs = 50, features = VariableFeatures(object = pbmc), verbose = FALSE)
pbmc <- pbmc %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.5) %>%
  identity()

# 聚类图：使用 RunUMAP 得到的 UMAP 坐标，不再用 PCA 坐标冒充 UMAP
p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE) + ggfonttheme
p1
ggsave("./Fibroblasts/01.cluster.pdf", plot = p1, device = cairo_pdf, width = 6, height = 6)
ggsave("./Fibroblasts/01.cluster.png", plot = p1, dpi = 300, width = 6, height = 6)

new.cluster.id = c(
  "0" = "Fibroblasts 1",
  "1" = "Fibroblasts 2",
  "2" = "Fibroblasts 1",
  "3" = "Fibroblasts 3",
  "4" = "Fibroblasts 3",
  "5" = "Fibroblasts 1",
  "6" = "Fibroblasts 3",
  "7" = "Fibroblasts 3",
  "8" = "Fibroblasts 2"
)
pbmc <- RenameIdents(pbmc, new.cluster.id)
pbmc$subcelltype <- Idents(pbmc)

# 亚细胞分群展示：同样使用 UMAP 坐标
p2 <- DimPlot(pbmc, reduction = "umap", group.by = "subcelltype", label = TRUE, repel = TRUE) + ggfonttheme
p2
ggsave("./Fibroblasts/02.subcelltype.pdf", plot = p2, device = cairo_pdf, width = 6, height = 6)
ggsave("./Fibroblasts/02.subcelltype.png", plot = p2, dpi = 300, width = 6, height = 6)
saveRDS(pbmc, "./Fibroblasts/annotated_pbmc.rds")
# ========== 02. 拟时序分析 ==========
mat <- GetAssayData(pbmc, slot = "data")
cellInfo <- pbmc@meta.data
geneInfo <- data.frame(gene_short_name = rownames(mat), row.names = rownames(mat))
cds <- new_cell_data_set(expression_data = mat,
                         cell_metadata = cellInfo,
                         gene_metadata = geneInfo)

# 数据预处理
cds <- preprocess_cds(cds, num_dim = 50)
reducedDim(cds, "UMAP") <- Embeddings(pbmc, "umap")[colnames(cds), , drop = FALSE]
# 轨迹学习
cds <- cluster_cells(cds, resolution = 1e-5)
cds <- learn_graph(cds, use_partition = FALSE)

# 按伪时间排序细胞
start = c("Fibroblasts 1")
closest_vertex = cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex = as.matrix(closest_vertex[colnames(cds), ])
root_pr_nodes = igraph::V(principal_graph(cds)[["UMAP"]])$name
flag = closest_vertex[as.character(colData(cds)$subcelltype) %in% start, ]
flag = as.numeric(names(which.max(table(flag))))
root_pr_nodes = root_pr_nodes[flag]
cds = order_cells(cds, root_pr_nodes = root_pr_nodes)

# 伪时间可视化
p3 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 show_trajectory_graph = FALSE,
                 label_cell_groups = FALSE,
                 label_leaves = FALSE,
                 cell_size = 1) + ggfonttheme
p3
ggsave("./Fibroblasts/03.presudotime.pdf", plot = p3, device = cairo_pdf, width = 7, height = 6)
ggsave("./Fibroblasts/03.presudotime.png", plot = p3, dpi = 300, width = 7, height = 6)

p4 <- plot_cells(cds,
                 color_cells_by = "subcelltype",
                 show_trajectory_graph = FALSE,
                 label_cell_groups = FALSE,
                 label_leaves = FALSE,
                 cell_size = 1) + ggfonttheme
p4
ggsave("./Fibroblasts/04.celltype.pdf", plot = p4, device = cairo_pdf, width = 7, height = 6)
ggsave("./Fibroblasts/04.celltype.png", plot = p4, dpi = 300, width = 7, height = 6)

# 基因表达分布：变量名改为目标基因 TGFB1
target_gene <- "TGFB1"
p6 <- plot_cells(cds,
                 genes = target_gene,
                 show_trajectory_graph = FALSE,
                 label_cell_groups = FALSE,
                 label_leaves = FALSE,
                 cell_size = 1) + ggfonttheme
p6
ggsave("./Fibroblasts/06.genes_distribution.pdf", plot = p6, device = cairo_pdf, width = 7, height = 6)
ggsave("./Fibroblasts/06.genes_distribution.png", plot = p6, dpi = 300, width = 7, height = 6)

# 基因在拟时间中的表达：变量名改为 TGFB1_cds，去除 monocle3 教程残留 PRNP
genes <- c("TGFB1")
TGFB1_cds <- cds[rowData(cds)$gene_short_name %in% genes, ]
monocle3::plot_genes_in_pseudotime(TGFB1_cds,
                                   color_cells_by = "subcelltype") + ggfonttheme
ggsave("./Fibroblasts/07.genes_in_pseudotime.pdf", device = cairo_pdf, width = 7, height = 6)
ggsave("./Fibroblasts/07.genes_in_pseudotime.png", dpi = 300, width = 7, height = 6)

save.image("../RData/13_monocle.RData")