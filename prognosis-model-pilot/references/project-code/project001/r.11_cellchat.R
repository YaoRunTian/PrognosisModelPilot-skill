rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project001")
if (!dir.exists("11_cellchat")) {dir.create("11_cellchat")}
setwd("11_cellchat")
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
pbmc <- readRDS("../RData/R.07_scRNA_Scissor.RData")
pbmc$subcelltype <- factor(pbmc$subcelltype)
pbmc2 <- readRDS("../RData/pbmc_all_with_subcelltype.Rds")
table(pbmc$subcelltype)
table(pbmc2$subcelltype)
pbmc$subcelltype <- pbmc2$subcelltype
print(table(Original = pbmc$celltype, Subtype = pbmc$subcelltype, useNA = "ifany"))
rm(pbmc2)
table(pbmc$scissor)
colnames(pbmc@meta.data)
# 创建celltype_for_cellchat列
pbmc$celltype_for_cellchat <- as.character(pbmc$subcelltype)
# 第一步：将celltype为Epithelial cells的行覆盖为"Epithelial cells"
pbmc$celltype_for_cellchat[pbmc$celltype == "Epithelial cells"] <- "Epithelial cells"
pbmc$scissor_with_celltype <- paste(pbmc$celltype,pbmc$scissor,sep ="&")

# 第二步：用scissor中的Scissor+和Scissor-进行覆盖
# pbmc$celltype_for_cellchat[pbmc$scissor == "Scissor+"] <- "Scissor+"
# pbmc$celltype_for_cellchat[pbmc$scissor == "Scissor-"] <- "Scissor-"

pbmc$celltype_for_cellchat[pbmc$scissor_with_celltype == "Epithelial cells&Scissor-"] <- "Epithelial cells&Scissor-"
pbmc$celltype_for_cellchat[pbmc$scissor_with_celltype == "Epithelial cells&Scissor+"] <- "Epithelial cells&Scissor+"

# 转为factor
pbmc$celltype_for_cellchat <- factor(pbmc$celltype_for_cellchat)
# 验证结果
table(pbmc$celltype_for_cellchat)
# 查看交叉表确认覆盖是否正确
print(table(celltype_for_cellchat = pbmc$celltype_for_cellchat, scissor = pbmc$scissor))
print(table(celltype_for_cellchat = pbmc$celltype_for_cellchat, celltype = pbmc$celltype))

pbmc_cellchat <- createCellChat(pbmc[["RNA"]]$data)
pbmc_meta <- data.frame(subcelltype = factor(pbmc$celltype_for_cellchat), row.names =  Cells(pbmc))
pbmc_cellchat <- addMeta(pbmc_cellchat, meta = pbmc_meta, meta.name = "celltype_for_cellchat")
pbmc_cellchat <- setIdent(pbmc_cellchat, ident.use = "celltype_for_cellchat")
groupSize <- as.numeric(table(pbmc_cellchat@idents)) 
#############################加载与设定需要的CellChatDB数据库
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
#############################预处理表达数据以进行细胞间相互作用分析
pbmc_cellchat <- subsetData(pbmc_cellchat) # subset the expression data of signaling genes for saving computation cost
pbmc_cellchat <- identifyOverExpressedGenes(pbmc_cellchat) 
pbmc_cellchat <- identifyOverExpressedInteractions(pbmc_cellchat)
#pbmc_cellchat <- smoothData(pbmc_cellchat, adj = PPI.human)###这一行代码很耗时
#############################推断细胞间相互作用网络与分析
options(future.globals.maxSize = 30000000000)  # 设置最大大小为2000 MB
future::plan("multisession", workers = 6)
pbmc_cellchat <- computeCommunProb(pbmc_cellchat)###这应该是核心的一行代码，也是耗时的代码
pbmc_cellchat <- filterCommunication(pbmc_cellchat, min.cells =10)##过滤掉小于10个细胞的胞间通讯网络
#######在信号通路的水平进一步推测胞间通讯，计算聚合网络
pbmc_cellchat <- computeCommunProbPathway(pbmc_cellchat)
pbmc_cellchat <- aggregateNet(pbmc_cellchat)
# Compute the network centrality scores
pbmc_cellchat <- netAnalysis_computeCentrality(pbmc_cellchat, slot.name = "netP") 
saveRDS(pbmc_cellchat,'Cellchat.rds')
cellchat <- pbmc_cellchat
par(mfrow = c(1,2))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions in cells")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength in cells")
p_base <- recordPlot()
dev.off()
pdf(file= "01.circle.pdf",width = 20, height = 12,family = "Times")
p_base
dev.off()
png("01.circle.png", width = 20, height = 12, units = "in", res = 300,family = "Times")
p_base 
dev.off()

# 定义要显示的细胞类型
targets <- c("Epithelial cells", "M1 Macrophages", "M2 Macrophages", "Epithelial cells&Scissor+", "Epithelial cells&Scissor-")

# 获取所有细胞类型
all_cells <- rownames(cellchat@net$count)

# 创建sources.use和targets.use（只显示这些细胞之间的互作）
sources.use <- which(all_cells %in% targets)
targets.use <- which(all_cells %in% targets)

# 绘图
par(mfrow = c(1,2))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge = F, title.name = "Number of interactions in cells",
                 sources.use = sources.use, targets.use = targets.use)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge = F, title.name = "Interaction weights/strength in cells",
                 sources.use = sources.use, targets.use = targets.use)
p_base <- recordPlot()
pdf(file= "01.circle_hubcell.pdf",width = 20, height = 12,family = "Times")
p_base
dev.off()
png("01.circle_hubcell.png", width = 20, height = 12, units = "in", res = 300,family = "Times")
p_base 
dev.off()


cellchat@netP$pathways  #查看都有哪些信号通路
netVisual_bubble(cellchat, font.size = 12)+ggfonttheme
ggsave("02.bubbleplot.pdf", width = 20, height = 20)
ggsave("02.bubbleplot.png", width = 20, height = 20,dpi = 600)


# mat <- cellchat@net$count
# mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
# mat2[c("Submucosal_gland_cell"), ] <- mat[c("Submucosal_gland_cell"), ]
# mat3 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
# mat3[c("Secretory_cell"), ] <- mat[c("Secretory_cell"), ]
# par(mfrow = c(1,2))
# netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
#                  arrow.size = 0.1, edge.weight.max = max(mat),title.name = "Interaction weights/strength in Submucosal gland cells")
# netVisual_circle(mat3, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
#                  arrow.size = 0.1, edge.weight.max = max(mat),title.name = "Interaction weights/strength in Secretory cells")
# p_base1 <- recordPlot()
# pdf(file= "03.circle_hubcell.pdf",width = 20, height = 12)
# p_base1
# dev.off()
# png("03.circle_hubcell.png", width = 20, height = 12, units = "in", res = 300)
# p_base1 
# dev.off()
save.image("../RData/11_cellchat.RData")

#预后基因的表达
hubgenes <- read.table("../07_model/06.lassogeneslist.txt",sep = "\n")[,1]
# 查看基因列表
print(hubgenes)

# 检查哪些基因存在于表达矩阵中
hubgenes_exist <- hubgenes[hubgenes %in% rownames(pbmc)]
print(paste0("存在的基因数量: ", length(hubgenes_exist), "/", length(hubgenes)))
print(hubgenes_exist)
Idents(pbmc) <- "celltype_for_cellchat"
# 气泡图展示
p_dotplot <- DotPlot(pbmc, features = hubgenes_exist) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Genes", y = "Cell Type", title = "Expression of Prognostic Genes") +
  scale_size(range = c(0, 10))+
  ggfonttheme
p_dotplot
ggsave("04.gene_dotplot.pdf", width = 10, height = 10)
ggsave("04.gene_dotplot.png", width = 10, height = 10,dpi = 600)
