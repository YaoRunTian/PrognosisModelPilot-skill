rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/24_job_0780")
if (!dir.exists("09_cellchat")) {dir.create("09_cellchat")}
setwd("09_cellchat")
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
pbmc <- readRDS("../RData/R.06_scRNA_done.Rds")

pbmc$celltype <- factor(
  pbmc$celltype
)
table(pbmc$celltype)
cellchat <- createCellChat(object=pbmc,group.by = "celltype")
CellChatDB <- CellChatDB.human#导入受体、配体数据库
# str(CellChatDB) #查看数据库信息
# #包含interaction、complex、cofactor和geneInfo这4个dataframe
# colnames(CellChatDB$interaction) 
# CellChatDB$interaction[1:4,1:4]
# head(CellChatDB$cofactor)
# head(CellChatDB$complex)
# head(CellChatDB$geneInfo)
showDatabaseCategory(CellChatDB)
unique(CellChatDB$interaction$annotation)#查看可以选择的侧面，也就是上图左中的三种
#选择"Secreted Signaling"进行后续细胞互作分析
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use # set the used database in the object
## 在矩阵的所有的基因中提取signaling gene，结果保存在data.signaling(13714个基因，过滤完只有270个）
cellchat <- subsetData(cellchat)
#相当于Seurat的FindMarkers，找每个细胞群中高表达的配体受体
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) #Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
#上一步运行的结果储存在cellchat@LR$LRsig
cellchat <- projectData(cellchat, PPI.human) 
#找到配体受体关系后，projectData将配体受体对的表达值投射到PPI上，来对@data.signaling中的表达值进行校正。结果保存在@data.project
#推断配体-受体水平细胞通讯网络-------------------------------------------------------
#根据表达值推测细胞互作的概率（cellphonedb是用平均表达值代表互作强度）。
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.xlsx(df.net, "01.net_lr.xlsx")
# 推断信号通路水平的细胞通讯网络------------------------------------------------
cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.xlsx(df.netp, "02.net_pathway.xlsx")
# 所有细胞群总体观：细胞互作数量与强度统计分析：-------------------------------
#统计细胞和细胞之间通信的数量（有多少个配体-受体对）和强度（概率）
cellchat <- aggregateNet(cellchat)
#计算每种细胞各有多少个
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions in cells")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights in cells")
p_base <- recordPlot()
dev.off()
pdf(file= "03.circle.pdf",width = 13, height = 12,family = "Times")
p_base
dev.off()
png("03.circle.png", width = 13, height = 12, units = "in", res = 300,family = "Times")
p_base 
dev.off()

cellchat@netP$pathways  #查看都有哪些信号通路

netVisual_bubble(cellchat,sources.use = c(1:8),targets.use = 1, font.size = 12)

gg1 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:8),  
                        max.dataset = 1, angle.x = 45, remove.isolate = T)+ggfonttheme
gg1
gg2 <- netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = 1,  
                         max.dataset = 1,angle.x = 45, remove.isolate = T)+ggfonttheme
gg2
gg1 + gg2+plot_layout(guides = "collect")
ggsave("04.signal_bubbleplot.pdf", width = 13, height = 18)
ggsave("04.signal_bubbleplot.png", width = 13, height = 18,dpi = 600)
# save.image("../RData/09_cellchat.RData")