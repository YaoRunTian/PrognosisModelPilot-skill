rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project001")
if (! dir.exists("13_monocle")){dir.create("13_monocle")}
setwd("13_monocle")
library(Seurat)
library(monocle3)
library(dplyr)
library(patchwork)
ggfonttheme <- theme(text = element_text(family = "Times"),
                     axis.title = element_text(size = 12),
                     axis.text = element_text(size = 12),
                     plot.title = element_text(size =12))
pbmc_all <- readRDS("../RData/R.06_scRNA.norm.umap_done.Rds")
#01.Macrophages--------------------------------------
pbmc <- readRDS("../02_reannotated/Macrophages/Macrophages_reannotated.RData")
p1 <- DimPlot(pbmc, reduction = "umap", group.by = "celltype",label = TRUE, repel = TRUE)+ggfonttheme
p1
mat <- GetAssayData(pbmc, slot = "data")
cellInfo <- pbmc@meta.data
geneInfo <- data.frame(gene_short_name = rownames(mat), row.names = rownames(mat))
cds <- new_cell_data_set(expression_data = mat,
                         cell_metadata = cellInfo,
                         gene_metadata = geneInfo)
# 数据预处理
cds <- preprocess_cds(cds, num_dim = 50)
# 批次校正
# cds <- align_cds(cds, alignment_group = "batch")
# 降维
cds <- reduce_dimension(cds)

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(pbmc, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed), ]
cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds)
# reducedDim(cds, "UMAP") <- Embeddings(pbmc, "pca")
cds <- learn_graph(cds)

# 按伪时间排序细胞
start = c("CCL4_M")
closest_vertex = cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex = as.matrix(closest_vertex[colnames(cds), ])
root_pr_nodes = igraph::V(principal_graph(cds)[["UMAP"]])$name
flag = closest_vertex[as.character(colData(cds)$celltype) %in% start,]
flag = as.numeric(names(which.max(table(flag))))
root_pr_nodes = root_pr_nodes[flag]
cds = order_cells(cds, root_pr_nodes=root_pr_nodes)
saveRDS(cds,"./Macrophages/cds.Rds")
# 按伪时间可视化
p3 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 show_trajectory_graph=T,
                 label_cell_groups=F,
                 label_leaves=F,
                 cell_size = 1)+ggfonttheme
p3
ggsave("./Macrophages/01.presudotime.pdf",p3,device = cairo_pdf,width = 7, height = 6)
ggsave("./Macrophages/01.presudotime.png",p3,dpi = 300,width = 7, height = 6)
p4 <- plot_cells(cds,
                 color_cells_by = "celltype",
                 show_trajectory_graph=T,
                 label_cell_groups=F,
                 label_leaves=F,
                 cell_size = 1)+ggfonttheme
p4
ggsave("./Macrophages/02.celltype.pdf",p4,device = cairo_pdf,width = 7, height = 6)
ggsave("./Macrophages/02.celltype.png",p4,dpi = 300,width = 7, height = 6)

ciliated_genes <- read.table("../07_model/06.lassogeneslist.txt")[,1]
p6 <- plot_cells(cds,
                 genes=ciliated_genes,
                 show_trajectory_graph=F,
                 label_cell_groups=F,
                 label_leaves=F,
                 cell_size = 1)+ggfonttheme
p6
ggsave("./Macrophages/03.genes_distribution.pdf",p6,device = cairo_pdf,width = 7, height = 6)
ggsave("./Macrophages/03.genes_distribution.png",p6,dpi = 300,width = 7, height = 6)

ciliated_genes
hubgenes_cds <- cds[rowData(cds)$gene_short_name %in% ciliated_genes,]
monocle3::plot_genes_in_pseudotime(hubgenes_cds,
                                   color_cells_by="celltype")+
      ggfonttheme
ggsave("./Macrophages/04.genes_in_pseudotime.pdf",device = cairo_pdf,width = 7, height = 12)
ggsave("./Macrophages/04.genes_in_pseudotime.png",dpi = 300,width = 7, height = 12)
#02.Epithelial cells--------------------------------------
rm(list = ls()); gc()
ggfonttheme <- theme(text = element_text(family = "Times"),
                     axis.title = element_text(size = 12),
                     axis.text = element_text(size = 12),
                     plot.title = element_text(size =12))
pbmc <- readRDS("../02_reannotated/Epithelial cells/Epithelial_cells_reannotated.RData")

mat <- GetAssayData(pbmc, slot = "data")
cellInfo <- pbmc@meta.data
geneInfo <- data.frame(gene_short_name = rownames(mat), row.names = rownames(mat))
cds <- new_cell_data_set(expression_data = mat,
                         cell_metadata = cellInfo,
                         gene_metadata = geneInfo)
cds <- preprocess_cds(cds, num_dim = 50)
reducedDim(cds, "UMAP") <- Embeddings(pbmc, "pca")
cds <- cluster_cells(cds, resolution=1e-5)
cds <- learn_graph(cds)

# 按伪时间排序细胞
start = c("Basal_cell")
closest_vertex = cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex = as.matrix(closest_vertex[colnames(cds), ])
root_pr_nodes = igraph::V(principal_graph(cds)[["UMAP"]])$name
flag = closest_vertex[as.character(colData(cds)$celltype) %in% start,]
flag = as.numeric(names(which.max(table(flag))))
root_pr_nodes = root_pr_nodes[flag]
cds = order_cells(cds,reduction_method = "UMAP",root_pr_nodes = root_pr_nodes)
# 按伪时间可视化
p3 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_branch_points = F,
                 show_trajectory_graph=F,
                 label_cell_groups=F,
                 label_leaves=F,
                 cell_size = 1)+ggfonttheme
p3
ggsave("./Epithelial cells/01.presudotime.pdf",p3,device = cairo_pdf,width = 7, height = 6)
ggsave("./Epithelial cells/01.presudotime.png",p3,dpi = 300,width = 7, height = 6)
p4 <- plot_cells(cds,
                 color_cells_by = "celltype",
                 label_branch_points = F,
                 show_trajectory_graph=F,
                 label_cell_groups=F,
                 label_leaves=F,
                 cell_size = 1)+ggfonttheme
p4
ggsave("./Epithelial cells/02.celltype.pdf",p4,device = cairo_pdf,width = 7, height = 6)
ggsave("./Epithelial cells/02.celltype.png",p4,dpi = 300,width = 7, height = 6)

ciliated_genes <- read.table("../07_model/06.lassogeneslist.txt")[,1]
p6 <- plot_cells(cds,
                 genes=ciliated_genes,
                 label_branch_points = F,
                 show_trajectory_graph=F,
                 label_cell_groups=F,
                 label_leaves=F,
                 cell_size = 1)+ggfonttheme
p6
ggsave("./Epithelial cells/03.genes_distribution.pdf",p6,device = cairo_pdf,width = 7, height = 6)
ggsave("./Epithelial cells/03.genes_distribution.png",p6,dpi = 300,width = 7, height = 6)

ciliated_genes
hubgenes_cds <- cds[rowData(cds)$gene_short_name %in% ciliated_genes,]
monocle3::plot_genes_in_pseudotime(hubgenes_cds,
                                   color_cells_by="celltype")+ggfonttheme
ggsave("./Epithelial cells/04.genes_in_pseudotime.pdf",device = cairo_pdf,width = 7, height = 12)
ggsave("./Epithelial cells/04.genes_in_pseudotime.png",dpi = 300,width = 7, height = 12)
