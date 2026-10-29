rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project001")
if (!dir.exists("01_sc")) {dir.create("01_sc")}
setwd("01_sc")
library(Seurat)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(extrafont)
library(RColorBrewer)
ggfonttheme <- theme(text = element_text(family = "Times"))
scRNA <- readRDS("../00_rawdata/GSE148071_pbmc.Rds")
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
scRNA@meta.data$group <- gsub("GSM\\d+_|_exp", "", scRNA@meta.data$sample_id)
length(colnames(scRNA))#[1] 89819
length(rownames(scRNA))#[1] 29119
saveRDS(scRNA,file = "../RData/R.01_scRNA.Rds")
#01.质控---------------------------
scRNA <- readRDS("../RData/R.01_scRNA.Rds")
#绘图风格
theme.set = theme(
  axis.title.x=element_blank(),
  axis.title = element_text(size = 20),
  axis.text.x = element_text(size = 14,angle = 90),
  axis.text.y = element_text(size = 14),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 18))
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")   
plots<- list()
# 质控前小提琴图
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, 
                       pt.size = 0,
                       group.by = "group",
                       features = plot.featrures[i]) + theme.set + NoLegend()+ggfonttheme+theme(axis.text.x = element_blank())}

violin <- wrap_plots(plots = plots, nrow=1,ncol = 3,cols = palettes(category = "random",22,show_col=F))
violin
ggsave("01.vlnplot_before_qc.pdf", plot = violin, width = 18, height = 6)
ggsave("01.vlnplot_before_qc.png", plot = violin, width = 18, height = 6,dpi = 600)
#  设置阈值
pctMT <- 30
minGene <- 200
minUMI <- 200
maxGene <- 5000
maxUMI <- 30000
scRNA <- subset(scRNA, subset =
                  nCount_RNA <maxUMI &
                  nCount_RNA >minUMI &
                  nFeature_RNA > minGene & 
                  nFeature_RNA < maxGene& 
                  percent.mt < pctMT)
length(colnames(scRNA))#89807
length(rownames(scRNA))#29119
# 质控后小提琴图
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, 
                       pt.size = 0,
                       group.by = "group",
                       features = plot.featrures[i]) + theme.set + NoLegend()+ggfonttheme+theme(axis.text.x = element_blank())}

violin <- wrap_plots(plots = plots, nrow=1,ncol = 3,cols = palettes(category = "random",22,show_col=F))
violin
ggsave("02.vlnplot_after_qc.pdf", plot = violin, width = 18, height = 6)
ggsave("02.vlnplot_after_qc.png", plot = violin, width = 18, height = 6,dpi = 600)
saveRDS(scRNA, file = "../RData/R.02_scRNA_qc.Rds")
#03.Normalize 标准化  -----------------------------------------------------
scRNA <- readRDS("../RData/R.02_scRNA_qc.Rds")
#  合并count
scRNA <- JoinLayers(scRNA)
#标准化
options(future.globals.maxSize = 4000 * 1024^2,scipen = 200)
scRNA.norm <- NormalizeData(scRNA,normalization.method = "LogNormalize",
                            scale.factor = 10000)

## 缩放数据
## 筛选2000个高变基因
all.genes <- rownames(scRNA.norm)
scRNA.norm <- FindVariableFeatures(scRNA.norm,selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scRNA.norm), 10)
#  plot variable features with and without labels
#  保存基因特征方差图
theme.set= theme(
  text = element_text(family = "Times New Roman"),
  axis.title = element_text(size = 20),
  axis.text.x = element_text(size = 10),
  axis.text.y = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 16))
pdf(file="03.feature_selection.pdf",width=8,height=6,onefile=F,family = "Times")             
plot1 <- VariableFeaturePlot(object = scRNA.norm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)
dev.off()
#保存基因特征方差图
png(file="03.feature_selection.png",width=8,height=6,units='in',res=600,family = "Times")              
plot1 <- VariableFeaturePlot(object = scRNA.norm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)
dev.off()
print(plot2)
#04.PCA  数据降维  --------------------------------------------------------
#   PCA降维之前的标准预处理步骤,这里对高可变基因进行数据标准化的操作,归一化表达量
scRNA.nor.sca <- ScaleData(scRNA.norm)
#   PCA分析   npcs = 20  在findneighbor中设置，不要改下面的
scRNA.norm.pca <- RunPCA(scRNA.nor.sca,features = VariableFeatures(object = scRNA.nor.sca),npcs = 50)
print(scRNA.norm.pca[["pca"]], dims = 1:5, nfeatures = 5)
##可视化PCA结果
VizDimLoadings(scRNA.norm.pca, dims = 1:2, reduction = "pca")
DimHeatmap(scRNA.norm.pca, dims = 1:5, cells = 2000, balanced = TRUE)
##jackstraw
scRNA.norm.pca <- JackStraw(scRNA.norm.pca, num.replicate = 100, dims = 50) #num.replicate,运行速度慢，影响JackStrawPlot
scRNA.norm.pca <- ScoreJackStraw(scRNA.norm.pca, dims = 1:50)
saveRDS(scRNA.norm.pca, file = "../RData/R.03_scRNA_pca.Rds")
##
theme.set = theme(
  #axis.title.x=element_blank(),
  axis.title = element_text(size = 20),
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 18))
##ElbowPlot,碎石图
pdf(file="04.ElbowPlot_50.pdf",width=8,height=6,family = "Times")
ElbowPlot(scRNA.norm.pca,reduction="pca",ndims = 50) +theme.set #根据ElbowPlot确定PC数量
dev.off()
png(file="04.ElbowPlot_50.png",width=8,height=6,units='in',res=600,family = "Times")
ElbowPlot(scRNA.norm.pca,reduction="pca",ndims = 50)+theme.set  
dev.off()
#pca图
pdf(file="05.pcaPlot.pdf",width=9,height=6,family = "Times")
DimPlot(scRNA.norm.pca, reduction = "pca",group.by = "group")
dev.off()
png(file="05.pcaPlot.png",width=9,height=6,res=600,units='in',family = "Times")
DimPlot(scRNA.norm.pca, reduction = "pca",group.by = "group")
dev.off()
# ##根据下面两张图可以选择PC数量
#JackStrawPlot
pdf(file="06.JackStrawPlot_30.pdf",width=8,height=6,family = "Times")
JackStrawPlot(object = scRNA.norm.pca, dims = 1:50,reduction = "pca")+theme.set
dev.off()
png(file="06.JackStrawPlot_30.png",width=8,height=6,units='in',res=600,family = "Times")
JackStrawPlot(object = scRNA.norm.pca, dims = 1:50,reduction = "pca")+theme.set
dev.off()
#05.umap降维-----------------------------------------------------------
library(tidyverse)
library(scales)
library(RColorBrewer)
library(ggsci)
library(wesanderson)
# library(IOBR)
pcSelect <- 20
scRNA.norm.umap <- FindNeighbors(object = scRNA.norm.pca, dims = 1:pcSelect)  #计算邻接距离
##对细胞分组,优化标准模块化
##这里的resolution决定了细胞亚群的多少，Seurat官网推荐0.4-1.2都可以进行尝试。数字越小，细胞亚群越少
##resolution参数表达聚类的分辨率，这个值大于0，一般都是在0-1范围里面调整，越大得到的cluster越多
#resolution可以反复调整，并不会改变降维的结果
scRNA.norm.umap <- FindClusters(object = scRNA.norm.umap, resolution = 0.6) #32
##umap聚类
scRNA.norm.umap  <- RunUMAP(scRNA.norm.umap , dims = 1:pcSelect)
##所有样本可视化
theme.set = theme(
  #axis.title.x=element_blank(),
  axis.title = element_text(size = 20),
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 20))
pdf(file="07.umap_all_04_20PC.pdf",width=10,height=8,family='Times')
DimPlot(scRNA.norm.umap ,
        reduction = "umap",
        label.size = 6,label = TRUE,pt.size=0.7
        #,cols = palettes(category = "random",24,show_col=F)
)+theme.set
dev.off()
png(file="07.umap_all_04_20PC.png",width=10,height=8,family='Times',units='in',res=600)
DimPlot(scRNA.norm.umap ,
        reduction = "umap",
        label.size = 6,label = TRUE, pt.size=0.7
        #,cols = palettes(category = "random",24,show_col=F)
)+theme.set
dev.off()
system.time(saveRDS(scRNA.norm.umap, file = "../RData/R.04_scRNA_umap_04_20PC.Rds"))
#06.查找marker基因----------------------------------------------------------------
scRNA.norm.umap <- readRDS("../RData/R.04_scRNA_umap_04_20PC.Rds")
DefaultAssay(scRNA.norm.umap) <- "RNA"
library("multtest")
library("metap")
scRNA.norm.umap <- JoinLayers(scRNA.norm.umap)
nk.markers <- FindAllMarkers(scRNA.norm.umap, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
head(nk.markers)
dim(nk.markers)
# [1] 53595     7
#去重统计
x <- nk.markers$gene
length(x[!duplicated(x)])
# [1] 11888
library("openxlsx")
write.xlsx(nk.markers, file = "08.MarkerGenes.xlsx")
#气泡图可视化
top5.markers <- nk.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
markers <- top5.markers$gene
DotPlot(scRNA.norm.umap, features = split(unique(top5.markers$gene), top5.markers$cluster),group.by = "seurat_clusters",cols = c("lightgrey", "red"))+
  theme(axis.text.x=element_text(angle=90, hjust = 1,size = 7),
        panel.border = element_rect(color="black"), #面板边框
        panel.spacing = unit(1, "mm")#面板间距
  )+ggfonttheme
ggsave("09.cluster_dotplot.pdf",device = cairo_pdf,width = 20, height = 8)
ggsave("09.cluster_dotplot.png",dpi = 300,width = 20, height = 8)
#07.注释---------------------------------------------------
#根据参考文献提供的文件进行注释
#PMID: 33953163和uniLUNG图谱（https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(25)00469-4/fulltext）以及singleR辅助
###########
library("SingleR",lib.loc = "<REFERENCE_DATA_ROOT>/R/library")
library("celldex")
library("dbplyr")
# ref <- HumanPrimaryCellAtlasData()
# saveRDS(ref,"HumanPrimaryCellAtlasData.Rds")
ref <- readRDS("HumanPrimaryCellAtlasData.Rds")
data_for_SingleR <- GetAssayData(scRNA.norm.umap, layer ="data")
Anno1 <- SingleR(test = data_for_SingleR, ref = ref, labels = ref$label.main)
Anno1
meta <- scRNA.norm.umap@meta.data
cluster_anno <- table(Anno1$labels,meta$seurat_clusters)
cluster_anno#singleR结果
write.xlsx(cluster_anno,"10.singleR_res.xlsx")
# 细胞类型特异性marker基因列表
cell_markers <- list(
  # 内皮细胞标记
  Endothelial_cells= c("CLDN5", "VWF", "PECAM1"),
  # 上皮细胞标记
  Epithelial_cells= c("CAPS", "SNTN","EPCAM"),
  # 肺泡细胞标记
  Alveolar = c("CLDN18", "AQP4"),
  # 成纤维细胞标记
  Fibroblasts = c("COL1A1", "COL1A2", "DCN"),
  # T细胞标记
  T_cells = c("CD2", "CD3D", "CD3E", "CD3G"),
  # B细胞标记
  B_cells = c("CD79A", "CD79B"),
  #巨噬细胞
  Macrophages = c("CD68","CD163","MRC1","MARCO"),
  # 髓系细胞标记
  Myeloid_cells = c("LYZ","CD14"),
  # 肥大细胞标记
  Mast_cells = c("GATA2", "TPSAB1", "TPSB2")
)
# 创建单个向量包含所有marker基因
# marker <- unlist(cell_markers, use.names = FALSE)
cell_markers_df <- stack(cell_markers)
names(cell_markers_df) <- c("gene", "cell_type")
cell_markers_df$cell_type <- gsub("_"," ",cell_markers_df$cell_type)
DotPlot(scRNA.norm.umap, features = split(unique(cell_markers_df$gene), cell_markers_df$cell_type),group.by = "seurat_clusters",cols = c("lightgrey", "red"))+
  theme(axis.text.x=element_text(angle=90, hjust = 1),
        panel.border = element_rect(color="black"), #面板边框
        panel.spacing = unit(1, "mm")#面板间距
  )+ggfonttheme
ggsave("11.cluster_markergenes_dotplot.pdf",device = cairo_pdf,width = 20, height = 8)
ggsave("11.cluster_markergenes_dotplot.png",dpi = 300,width = 20, height = 8)
new.cluster.id <- c(
  "0" = "Epithelial cells",
  "1" = "Macrophages",
  "2" = "Epithelial cells",
  "3" = "Epithelial cells",
  "4" = "Myeloid cells",
  "5" = "Epithelial cells",
  "6" = "Epithelial cells",
  "7" = "Epithelial cells",
  "8" = "Fibroblasts",
  "9" = "Epithelial cells",
  "10" = "Epithelial cells",
  "11" = "Epithelial cells",
  "12" = "T cells",
  "13" = "Epithelial cells",
  "14" = "B cells",
  "15" = "Fibroblasts",
  "16" = "Fibroblasts",
  "17" = "Epithelial cells",
  "18" = "Epithelial cells",
  "19" = "Myeloid cells",
  "20" = "Alveolar",
  "21" = "Fibroblasts",
  "22" = "Epithelial cells",
  "23" = "Endothelial cells",
  "24" = "B cells",
  "25" = "Epithelial cells",
  "26" = "Epithelial cells",
  "27" = "Myeloid cells",
  "28" = "Mast cells",
  "29" = "Epithelial cells",
  "30" = "T cells",
  "31" = "Epithelial cells"
)
saveRDS(scRNA.norm.umap, file = "../RData/R.05_scRNA_umap_05_anno.Rds")

scRNA.norm.umap <- readRDS("../RData/R.05_scRNA_umap_05_anno.Rds")
cluster_dat <- data.frame("cluster" = 0:(length(new.cluster.id)-1),"CellType"=new.cluster.id)
write.xlsx(cluster_dat, file = "12.anno_CellType.xlsx")
names(new.cluster.id) <- levels(scRNA.norm.umap)
scRNA.norm.umap <- RenameIdents(scRNA.norm.umap, new.cluster.id)
scRNA.norm.umap$celltype <- Idents(scRNA.norm.umap)
# scRNA.norm.umap <- SetIdent(scRNA.norm.umap, value = "celltype")
ncluster <- length(unique(scRNA.norm.umap[[]]$celltype))
mycol <- colorRampPalette(brewer.pal(8, "Set3"))(ncluster)
# DimPlot(scRNA.norm.umap, 
#         label = T , repel = T, 
#         label.size = 5,
#         label.box = T,
#         cols = mycol) + 
#   NoLegend()
p9 <- DimPlot(scRNA.norm.umap, reduction = "umap", group.by = "group")+ggfonttheme+ggtitle("")
ggsave("12.umap_sample.pdf",p9,device = cairo_pdf,width = 8, height = 8)
ggsave("12.umap_sample.png",p9,dpi = 300,width = 8, height = 8)
p10 <- DimPlot(scRNA.norm.umap, reduction = "umap", label = T , repel = T, 
               label.size = 6,
               label.box = T,
               cols = mycol,
               pt.size = 0.9,
               alpha = c(0.3)) + 
  NoLegend()+ggfonttheme+coord_fixed(ratio = 1)
p10
ggsave("13.umap_annotated.pdf",p10,device = cairo_pdf,width = 8, height = 8)
ggsave("13.umap_annotated.png",p10,dpi = 300,width = 8, height = 8)
saveRDS(scRNA.norm.umap, file = "../RData/R.06_scRNA_done.Rds")
