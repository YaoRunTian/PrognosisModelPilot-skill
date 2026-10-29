rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project004")
if (!dir.exists("10_cellchat")) {dir.create("10_cellchat")}
setwd("10_cellchat")
library(Seurat)
library(SeuratData)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
library(openxlsx)
ggfonttheme <- theme(text = element_text(family = "Times"))
#加载数据
scRNA.norm.umap <- readRDS("../RData/R.06_scRNA_done.Rds")
pbmc <- scRNA.norm.umap
rm(scRNA.norm.umap)
pbmc$celltype <- factor(pbmc$celltype)
table(pbmc$celltype)
Idents(pbmc) <- "group"
#01.创建cellchat对象-----------------------------------------------
for (group in unique(pbmc$group)) {
  pbmc_subset <- pbmc[,pbmc$group == group]
  pbmc_cellchat <- createCellChat(pbmc_subset[["RNA"]]$data)
  pbmc_meta <- data.frame(celltype = factor(pbmc_subset$celltype), row.names =  Cells(pbmc_subset))
  
  pbmc_cellchat <- addMeta(pbmc_cellchat, meta = pbmc_meta, meta.name = "celltype")
  pbmc_cellchat <- setIdent(pbmc_cellchat, ident.use = "celltype") # set "labels" as default cell identity
  groupSize <- as.numeric(table(pbmc_cellchat@idents)) # number of cells in each cell group
  #############################4.加载与设定需要的CellChatDB数据库
  CellChatDB <- CellChatDB.human
  str(CellChatDB)
  colnames(CellChatDB$interaction)
  CellChatDB$interaction[1:4,1:4]
  head(CellChatDB$complex)
  head(CellChatDB$cofactor)
  head(CellChatDB$geneInfo)
  showDatabaseCategory(CellChatDB)
  CellChatDB.use <- CellChatDB
  pbmc_cellchat@DB <- CellChatDB.use # set the used database in the object
  #############################5.预处理表达数据以进行细胞间相互作用分析
  pbmc_cellchat <- subsetData(pbmc_cellchat) # subset the expression data of signaling genes for saving computation cost
  pbmc_cellchat <- identifyOverExpressedGenes(pbmc_cellchat) 
  pbmc_cellchat <- identifyOverExpressedInteractions(pbmc_cellchat)
  print("ok")
  #pbmc_cellchat <- smoothData(pbmc_cellchat, adj = PPI.human)###这一行代码很耗时
  #############################6.推断细胞间相互作用网络与分析
  options(future.globals.maxSize = 30000000000)  # 设置最大大小为2000 MB
  future::plan("multisession", workers = 6)
  pbmc_cellchat <- computeCommunProb(pbmc_cellchat)###这应该是核心的一行代码，也是耗时的代码
  pbmc_cellchat <- filterCommunication(pbmc_cellchat, min.cells =10)##过滤掉小于10个细胞的胞间通讯网络
  #######在信号通路的水平进一步推测胞间通讯，计算聚合网络
  pbmc_cellchat <- computeCommunProbPathway(pbmc_cellchat)
  pbmc_cellchat <- aggregateNet(pbmc_cellchat)
  # Compute the network centrality scores
  pbmc_cellchat <- netAnalysis_computeCentrality(pbmc_cellchat, slot.name = "netP") 
  saveRDS(pbmc_cellchat,paste0(group,'_Cellchat.rds'))
}

#02.Differential analysis---------------------------------------------
options(future.globals.maxSize = 3000000000)  # 设置最大大小为2000 MB
future::plan("multisession", workers = 6)
cellchat.Tumor <- readRDS("Tumor_Cellchat.rds")
cellchat.Normal <- readRDS("Normal_Cellchat.rds")
object.list <- list(Normal = cellchat.Normal,Tumor = cellchat.Tumor)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),size.text = 20)+ggfonttheme
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",size.text = 20,digits = 2)+ggfonttheme
gg1 + gg2
ggsave("01.barplot.pdf", width = 12, height = 6)
ggsave("01.barplot.png", width = 12, height = 6,dpi = 600)
###Compare the number of interactions and interaction strength among different cell populations
#To identify the interaction between which cell populations showing significant changes, 
#CellChat compares the number of interactions and interaction strength among different cell populations.
#The differential number of interactions or interaction strength in the cell-cell communication network between two datasets can be visualized 
#using circle plot, where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
#Post/Pre
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
p_base <- recordPlot()
dev.off()
pdf(file= "02.DiffInteractionCircle.pdf",width = 12, height = 8,family = "Times")
p_base
dev.off()
png("02.DiffInteractionCircle.png", width = 12, height = 8, units = "in", res = 300,family = "Times")
p_base 
dev.off()
#The differential network analysis only works for pairwise datasets. If there are more datasets for comparison, we can directly show the number of interactions or interaction strength between any two cell populations in each dataset.
#To better control the node size and edge weights of the inferred networks across different datasets, we compute the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) across all datasets.

weight.max <- getMaxWeight(object.list, attribute = c("count","weight"))
groupSize.list <- list(as.numeric(table(object.list[[1]]@idents)),as.numeric(table(object.list[[2]]@idents)))

par(mfrow = c(1,2), xpd=TRUE,family = "Times")
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F,
                   vertex.weight = groupSize.list[[i]],
                   # edge.weight.max = weight.max[1], 
                   # edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
p_base <- recordPlot()
dev.off()
pdf(file= "03.circle.pdf",width = 13, height = 8)
p_base
dev.off()
png("03.circle.png", width = 13, height = 8, units = "in", res = 300)
p_base 
dev.off()

par(mfrow = c(1,2), xpd=TRUE,family = "Times")
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= F, 
                   vertex.weight = groupSize.list[[i]],
                   #edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Weights of interactions - ", names(object.list)[i]))
}
p_base <- recordPlot()
dev.off()
pdf(file= "03.circle_weighted.pdf",width = 13, height = 8)
p_base
dev.off()
png("03.circle_weighted.png", width = 13, height = 8, units = "in", res = 300)
p_base 
dev.off()
### Compare the major sources and targets in 2D space
## Comparing the outgoing and incoming interaction strength in 2D space allows ready identification of the cell populations with significant changes in sending or receiving signals between different datasets.

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)+ggfonttheme
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
ggsave("04.bubbleplot.pdf", width = 14, height = 6)
ggsave("04.bubbleplot.png", width = 14, height = 6,dpi = 600)
# #From the scatter plot, we can see that Inflam.
# ## Identify signaling changes associated with one cell group
# gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "NK cell")
# gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "T cell")
# # gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Neutrophil")
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Epithelial cells")
gg4
# # gg5 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Plasma cell")
# # gg6 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Epithelium")
# # gg7 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Mast cell")
# # patchwork::wrap_plots(plots = c(gg,list(gg1,gg2,gg3,gg4,gg5,gg6,gg7)))
ggsave("05.signal_changes.pdf", width = 6, height = 6)
ggsave("05.signal_changes.png", width = 6, height = 6,dpi = 600)

#Identify and visualize the conserved and context-specific signaling pathways
#By comparing the information flow/interaction strengh of each signaling pathway, we can identify signaling pathways, (i) turn off, (ii) decrease, (iii) turn on or (iv) increase, by change their information flow at one condition as compared to another condition.
#Compare the overall information flow of each signaling pathway
#We can identify the conserved and context-specific signaling pathways by simply comparing the information flow for each signaling pathway, which is defined by the sum of communication probability among all pairs of cell groups in the inferred network (i.e., the total weights in the network).
#This bar graph can be plotted in a stacked mode or not. Significant signaling pathways were ranked based on differences in the overall information flow within the inferred networks between NL and LS skin. The top signaling pathways colored red are enriched in NL skin, and these colored green were enriched in the LS skin.
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)+ggfonttheme
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)+ggfonttheme
gg1 + gg2
ggsave("06.signal_barplot.pdf", width = 15, height = 9)
ggsave("06.signal_barplot.png", width = 15, height = 9,dpi = 600)

#03.signal pathways 
# Moreover, we can identify the upgulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. This can be done by specifying max.dataset and min.dataset in the function netVisual_bubble. The increased signaling means these signaling have higher communication probability (strength) in one dataset compared to the other dataset.
gg1 <- netVisual_bubble(cellchat, sources.use = 6, targets.use = c(1:9),  
                        comparison = c(1, 2), max.dataset = 1, angle.x = 45, remove.isolate = T)+ggfonttheme
gg1
gg2 <- netVisual_bubble(cellchat, sources.use = c(1:9), targets.use = 6,  
                        comparison = c(1, 2), max.dataset = 1,angle.x = 45, remove.isolate = T)+ggfonttheme
gg2
gg1 + gg2+plot_layout(guides = "collect")
ggsave("07.signal_bubbleplot.pdf", width = 13, height = 7)
ggsave("07.signal_bubbleplot.png", width = 13, height = 7,dpi = 600)
save.image("../RData/10_cellchat.RData")
