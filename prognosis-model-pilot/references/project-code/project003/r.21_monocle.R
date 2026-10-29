rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project003")
if (! dir.exists("21_monocle")){dir.create("21_monocle")}
setwd("21_monocle")
library(Seurat)
library(monocle3)
library(dplyr)
library(patchwork)
ggfonttheme <- theme(text = element_text(family = "Times"),
                     axis.title = element_text(size = 12),
                     axis.text = element_text(size = 12),
                     plot.title = element_text(size =12))
pbmc_all <- readRDS("../RData/R.06_scRNA_done.Rds")
#01.分群Endothelial cells--------------------------------------
pbmc <- subset(x = pbmc_all, idents = c("Endothelial cells"))
pbmc <- RunPCA(pbmc, npcs = 50,features = VariableFeatures(object = pbmc), verbose = FALSE)
pbmc <- pbmc %>%
  RunUMAP(reduction = "pca",dims = 1:20) %>%
  FindNeighbors(reduction = "pca",dims = 1:20) %>%
  FindClusters(resolution = 0.4) %>%
  identity()
p1 <- DimPlot(pbmc, reduction = "pca", label = TRUE, repel = TRUE)+ggfonttheme
p1
ggsave("./Endothelial cells/01.cluster.pdf",device = cairo_pdf,width = 6, height = 6)
ggsave("./Endothelial cells/01.cluster.png",dpi = 300,width = 6, height = 6)

#注释
new.cluster.id = c(
  "0" = "Endothelial cells 1",
  "1" = "Endothelial cells 2",
  "2" = "Endothelial cells 1",
  "3" = "Endothelial cells 2",
  "4" = "Endothelial cells 2",
  "5" = "Endothelial cells 1",
  "6" = "Endothelial cells 2",
  "7" = "Endothelial cells 2",
  "8" = "Endothelial cells 1",
  "9" = "Endothelial cells 3",
  "10" = "Endothelial cells 3"
)
pbmc <- RenameIdents(pbmc,new.cluster.id)
pbmc$subcelltype <- Idents(pbmc)  ###保存subcelltype
#亚细胞分群展示
p2 <- DimPlot(pbmc, reduction = "pca", group.by = "subcelltype",label = TRUE, repel = TRUE)+ggfonttheme
p2
ggsave("./Endothelial cells/02.subcelltype.pdf",device = cairo_pdf,width = 6, height = 6)
ggsave("./Endothelial cells/02.subcelltype.png",dpi = 300,width = 6, height = 6)
saveRDS(pbmc,"./Endothelial cells/annotated_pbmc.rds")
#02.拟时序分析---------------------------------------
mat <- GetAssayData(pbmc, slot = "data")
cellInfo <- pbmc@meta.data
geneInfo <- data.frame(gene_short_name = rownames(mat), row.names = rownames(mat))
cds <- new_cell_data_set(expression_data = mat,
                         cell_metadata = cellInfo,
                         gene_metadata = geneInfo)
# 数据预处理
cds <- preprocess_cds(cds, num_dim = 50)

reducedDim(cds, "UMAP") <- Embeddings(pbmc, "pca")
# 轨迹学习
cds <- cluster_cells(cds, resolution=1e-5)
cds <- learn_graph(cds,use_partition = FALSE)
# 批次校正
# cds <- align_cds(cds, alignment_group = "batch")
# 降维
# cds <- reduce_dimension(cds)
# 
# cds.embed <- cds@int_colData$reducedDims$UMAP
# int.embed <- Embeddings(pbmc, reduction = "umap")
# int.embed <- int.embed[rownames(cds.embed), ]
# cds@int_colData$reducedDims$UMAP <- int.embed
# 
# cds <- cluster_cells(cds)
# plot_cells(cds, color_cells_by = "partition",cell_size = 1)
# reducedDim(cds, "UMAP") <- Embeddings(pbmc, "pca")
# cds <- learn_graph(cds)

# 按伪时间排序细胞
start = c("Endothelial cells 1")
closest_vertex = cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex = as.matrix(closest_vertex[colnames(cds), ])
root_pr_nodes = igraph::V(principal_graph(cds)[["UMAP"]])$name
flag = closest_vertex[as.character(colData(cds)$subcelltype) %in% start,]
flag = as.numeric(names(which.max(table( flag ))))
root_pr_nodes = root_pr_nodes[flag]
cds = order_cells(cds, root_pr_nodes=root_pr_nodes)
# 按伪时间可视化
p3 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 show_trajectory_graph=F,
                 label_cell_groups=F,
                 label_leaves=F,
                 cell_size = 1)+ggfonttheme
p3
ggsave("./Endothelial cells/03.presudotime.pdf",p3,device = cairo_pdf,width = 7, height = 6)
ggsave("./Endothelial cells/03.presudotime.png",p3,dpi = 300,width = 7, height = 6)
p4 <- plot_cells(cds,
                 color_cells_by = "subcelltype",
                 show_trajectory_graph=F,
                 label_cell_groups=F,
                 label_leaves=F,
                 cell_size = 1)+ggfonttheme
p4
ggsave("./Endothelial cells/04.celltype.pdf",p4,device = cairo_pdf,width = 7, height = 6)
ggsave("./Endothelial cells/04.celltype.png",p4,dpi = 300,width = 7, height = 6)

ciliated_genes <- "PGLYRP4"
p6 <- plot_cells(cds,
                 genes=ciliated_genes,
                 show_trajectory_graph=F,
                 label_cell_groups=F,
                 label_leaves=F,
                 cell_size = 1)+ggfonttheme
p6
ggsave("./Endothelial cells/06.genes_distribution.pdf",p6,device = cairo_pdf,width = 7, height = 6)
ggsave("./Endothelial cells/06.genes_distribution.png",p6,dpi = 300,width = 7, height = 6)

genes <- c("PGLYRP4","FCN3")
PRNP_cds <- cds[rowData(cds)$gene_short_name %in% genes,]
monocle3::plot_genes_in_pseudotime(PRNP_cds,
                                   color_cells_by="subcelltype")+ggfonttheme
ggsave("./Endothelial cells/07.genes_in_pseudotime.pdf",device = cairo_pdf,width = 7, height = 6)
ggsave("./Endothelial cells/07.genes_in_pseudotime.png",dpi = 300,width = 7, height = 6)
