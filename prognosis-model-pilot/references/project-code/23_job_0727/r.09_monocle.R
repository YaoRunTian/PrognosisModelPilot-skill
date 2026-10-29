rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/23_job_0727")
if (! dir.exists("09_monocle")){dir.create("09_monocle")}
setwd("09_monocle")
library(Seurat)
library(monocle3)
library(dplyr)
library(patchwork)
ggfonttheme <- theme(text = element_text(family = "Times"),
                     axis.title = element_text(size = 12),
                     axis.text = element_text(size = 12),
                     plot.title = element_text(size =12))
pbmc_all <- readRDS("../RData/R.06_scRNA_done.Rds")
#01.分群Macrophages--------------------------------------
pbmc <- subset(x = pbmc_all, idents = c("Macrophages"))
pbmc <- RunPCA(pbmc, npcs = 50,features = VariableFeatures(object = pbmc), verbose = FALSE)
pbmc <- pbmc %>%
  RunUMAP(reduction = "pca",dims = 1:20) %>%
  FindNeighbors(reduction = "pca",dims = 1:20) %>%
  FindClusters(resolution = 0.4) %>%
  identity()
p1 <- DimPlot(pbmc, reduction = "pca", label = TRUE, repel = TRUE)+ggfonttheme
p1
ggsave("./Macrophages/01.cluster.pdf",device = cairo_pdf,width = 6, height = 6)
ggsave("./Macrophages/01.cluster.png",dpi = 300,width = 6, height = 6)

#注释
new.cluster.id = c(
  "0" = "Macrophages 1",
  "1" = "Macrophages 3",
  "2" = "Macrophages 3",
  "3" = "Macrophages 1",
  "4" = "Macrophages 3",
  "5" = "Macrophages 2",
  "6" = "Macrophages 2",
  "7" = "Macrophages 3",
  "8" = "Macrophages 3"
)
pbmc <- RenameIdents(pbmc,new.cluster.id)
pbmc$subcelltype <- Idents(pbmc)  ###保存subcelltype
#亚细胞分群展示
p2 <- DimPlot(pbmc, reduction = "pca", group.by = "subcelltype",label = TRUE, repel = TRUE)+ggfonttheme
p2
ggsave("./Macrophages/02.subcelltype.pdf",device = cairo_pdf,width = 6, height = 6)
ggsave("./Macrophages/02.subcelltype.png",dpi = 300,width = 6, height = 6)
saveRDS(pbmc,"./Macrophages/annotated_pbmc.rds")
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
start = c("Macrophages 1")
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
ggsave("./Macrophages/03.presudotime.pdf",p3,device = cairo_pdf,width = 7, height = 6)
ggsave("./Macrophages/03.presudotime.png",p3,dpi = 300,width = 7, height = 6)
p4 <- plot_cells(cds,
                 color_cells_by = "subcelltype",
                 show_trajectory_graph=F,
                 label_cell_groups=F,
                 label_leaves=F,
                 cell_size = 1)+ggfonttheme
p4
ggsave("./Macrophages/04.celltype.pdf",p4,device = cairo_pdf,width = 7, height = 6)
ggsave("./Macrophages/04.celltype.png",p4,dpi = 300,width = 7, height = 6)

ciliated_genes <- "SIRT7"
p6 <- plot_cells(cds,
                 genes=ciliated_genes,
                 show_trajectory_graph=F,
                 label_cell_groups=F,
                 label_leaves=F,
                 cell_size = 1)+ggfonttheme
p6
ggsave("./Macrophages/06.genes_distribution.pdf",p6,device = cairo_pdf,width = 7, height = 6)
ggsave("./Macrophages/06.genes_distribution.png",p6,dpi = 300,width = 7, height = 6)

genes <- c("SIRT7")
PRNP_cds <- cds[rowData(cds)$gene_short_name %in% genes,]
monocle3::plot_genes_in_pseudotime(PRNP_cds,
                                   color_cells_by="subcelltype")+ggfonttheme
ggsave("./Macrophages/07.genes_in_pseudotime.pdf",device = cairo_pdf,width = 7, height = 6)
ggsave("./Macrophages/07.genes_in_pseudotime.png",dpi = 300,width = 7, height = 6)
