rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/19_project_1533")
if (!dir.exists("03_sc")) {dir.create("03_sc")}
setwd("03_sc")
library(Seurat)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(extrafont)
library(RColorBrewer)
library(harmony)
library(scop)
ggfonttheme <- theme(text = element_text(family = "Times"))
scRNA <- readRDS("../00_rawdata/GSE135337_raw_pbmc.Rds")
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
length(colnames(scRNA))#[1] 42823
length(rownames(scRNA))#[1] 21854
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
scRNA$group <- scRNA@meta.data$orig.ident
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, 
                       pt.size = 0,
                       group.by = "group",
                       features = plot.featrures[i]) + theme.set + NoLegend()+ggfonttheme
    # theme(axis.text.x = element_blank())
  }

violin <- wrap_plots(plots = plots, nrow=1,ncol = 3,cols = palettes(category = "random",22,show_col=F))
violin
ggsave("01.vlnplot_before_qc.pdf", plot = violin, width = 18, height = 6)
ggsave("01.vlnplot_before_qc.png", plot = violin, width = 18, height = 6,dpi = 600)
#  设置阈值
pctMT <- 10
minGene <- 500
minUMI <- 0
maxGene <- 6000
maxUMI <- 30000
scRNA <- subset(scRNA, subset =
                  nCount_RNA <maxUMI &
                  nCount_RNA >minUMI &
                  nFeature_RNA > minGene & 
                  nFeature_RNA < maxGene& 
                  percent.mt < pctMT)
length(colnames(scRNA))#41746
length(rownames(scRNA))#21854
# 质控后小提琴图
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, 
                       pt.size = 0,
                       group.by = "group",
                       features = plot.featrures[i]) + theme.set + NoLegend()+ggfonttheme}

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
JackStrawPlot(object = scRNA.norm.pca, dims = 1:50,reduction = "pca")+theme.set+theme(axis.text.x = element_blank())
dev.off()
png(file="06.JackStrawPlot_30.png",width=8,height=6,units='in',res=600,family = "Times")
JackStrawPlot(object = scRNA.norm.pca, dims = 1:50,reduction = "pca")+theme.set+theme(axis.text.x = element_blank())
dev.off()
#harmony去批次--------------------
scRNA.norm.pca <- readRDS("../RData/R.03_scRNA_pca.Rds")
scRNA.harmony <- RunHarmony(scRNA.norm.pca, 
                    group.by.vars = "orig.ident")   # 显示收敛图

scRNA <- RunUMAP(scRNA.harmony, 
                 reduction = "harmony",  # 使用harmony校正后的结果
                 dims = 1:30)

scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:30)
scRNA <- FindClusters(scRNA, resolution = 0.5)#15
saveRDS(scRNA, file = "../RData/R.04_scRNA_harmony_umap.Rds")

p1 <- DimPlot(scRNA, reduction = "pca", group.by = "group", pt.size = 0.1) + 
  ggtitle("PCA - Before Harmony")+ggfonttheme
p2 <- DimPlot(scRNA, reduction = "harmony", group.by = "group", pt.size = 0.1) + 
  ggtitle("Harmony - After Harmony")+ggfonttheme
p3 <- DimPlot(scRNA, reduction = "umap", group.by = "group", pt.size = 0.1)+ggfonttheme
p4 <- DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", 
              label = TRUE, pt.size = 0.1)+ggfonttheme
library(patchwork)
(p1 + p2) / (p3 + p4)
ggsave("07.pca_umap_harmony.pdf",device = cairo_pdf,width = 8, height = 8)
ggsave("07.pca_umap_harmony.png",dpi = 300,width = 8, height = 8)
system.time(saveRDS(scRNA, file = "../RData/R.04_scRNA_umap_04_20PC.Rds"))
#06.查找marker基因----------------------------------------------------------------
scRNA.norm.umap <- readRDS("../RData/R.04_scRNA_umap_04_20PC.Rds")
DefaultAssay(scRNA.norm.umap) <- "RNA"
library("multtest")
library("metap")
scRNA.norm.umap <- JoinLayers(scRNA.norm.umap)
nk.markers <- FindAllMarkers(scRNA.norm.umap, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(nk.markers)
dim(nk.markers)
# [1] 15846     7
#去重统计
x <- nk.markers$gene
length(x[!duplicated(x)])
# [1] 5993
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
##PMID:  37664173
gene <- c(
  "CD3D","CD3E",##T_cells
  "IGHG1","IGKC",##B_cell
  "LYZ","CST3",##Monocyte
  "FGFBP2","GNLY",##NK_cell
  "PI3","CXCL8",##Neutrophils
  "DCN","CFD")##Smooth_muscle_cells

cell_markers <- list(
    T_cells  = c("CD3D","CD3E"),
    B_cells = c("IGHG1","IGKC"),
    Fibroblasts = c("COL1A2","COL1A1","COL3A1"),
    Endothelial_cells = c("PECAM1","VWF","CLDN5"),
    Myeloid_cells = c("CD14","CSF1R","AIF1"),
    Epithelial_cells = c("PTGDS","KRT18","KRT19")
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
ggsave("11.cluster_markergenes_dotplot.pdf",device = cairo_pdf,width = 15, height = 8)
ggsave("11.cluster_markergenes_dotplot.png",dpi = 300,width = 15, height = 8)
new.cluster.id <- c(
  "0" = "Epithelial cells",
  "1" = "Epithelial cells",
  "2" = "Epithelial cells",
  "3" = "Epithelial cells",
  "4" = "Epithelial cells",
  "5" = "Epithelial cells",
  "6" = "Epithelial cells",
  "7" = "Fibroblasts",
  "8" = "Epithelial cells",
  "9" = "Epithelial cells",
  "10" = "Myeloid cells",
  "11" = "T cells",
  "12" = "B cells",
  "13" = "Fibroblasts",
  "14" = "Fibroblasts",
  "15" = "Endothelial cells"
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
mycol <- colorRampPalette(brewer.pal(6, "Set3"))(ncluster)
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
#差异细胞鉴定----------------------------------------
scRNA.norm.umap <- readRDS("../RData/R.06_scRNA.norm.umap_done.Rds")
pbmc <- scRNA.norm.umap
pbmc$group[pbmc$group =="BCN"] <- "Normal"
pbmc$group[pbmc$group !="Normal"] <- "Tumor"
dat1 <- pbmc@meta.data
dat1 <- dat1[,c("group","celltype","nCount_RNA")]
head(dat1)
ncluster <- length(unique(dat1$celltype))
mycol <- colorRampPalette(brewer.pal(6, "Set3"))(ncluster)
#barplot---------------------------------------------------------
ggplot(data = dat1,aes(x=group,y=nCount_RNA,fill=celltype))+
  geom_bar(stat = "identity",position = "fill",width = 0.5)+
  scale_fill_manual(values = mycol)+
  scale_y_continuous(expand = expansion(mult=c(0.01,0.1)),
                     labels = scales::percent_format())+
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        legend.position = "top",
        panel.grid = element_line(colour = "white"),plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        title = element_text(size = 20),
        legend.text = element_text(size = 13))+
  labs(x=NULL,y="Cell Proportion(%)",fill = "Celltype")+
  guides(fill=guide_legend(nrow = 3,byrow = FALSE))+ggfonttheme
ggsave("14.cell_proportion_barplot.pdf",device = cairo_pdf,width = 6, height = 7)
ggsave("14.cell_proportion_barplot.png",dpi = 300,width = 6, height = 7)
#boxplot-------------------------------------------------------
comparisons <- list(names(table(dat1$group)))
ggplot(dat1, aes(celltype,nCount_RNA, fill = group)) + 
  #geom_violin(trim = FALSE,position = position_dodge(width = 0.9),scale = "area")+
  geom_boxplot(position = position_dodge(width = 0.9),
               outlier.size = 0,width = 0.5,show.legend = TRUE)+
  scale_fill_manual(values = c("red","blue"))+
  stat_compare_means(aes(group=group),                       #按分组进行统计检验
                     method = "wilcox.test",
                     paired = FALSE,                             
                     symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                        symbols=c("***","**","*","ns")),
                     label = "p.signif",
                     #label.y = location$Proportion+0.02,      #添加显著性符号的位置
                     size=4.5)+
  xlab("")+  ylab("nCount_RNA") +   theme_bw()+
  theme(panel.grid = element_line(colour = "white"),plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15,angle = 45,hjust = 1),
        title = element_text(size = 20),
        legend.position ="top",
        legend.key.size = unit(35, "pt"),
        legend.text = element_text(size = 13))+ggfonttheme
ggsave("15.cell_count_boxplot.pdf",device = cairo_pdf,width = 10, height = 7)
ggsave("15.cell_count_boxplot.png",dpi = 300,width = 10, height = 7)
#PRR11在差异细胞中的表达差异----------------------------------------
prr11_expression <- FetchData(pbmc, vars = "PRR11")
head(prr11_expression)
pbmc$prr11_expression <- prr11_expression$PRR11
dat1 <- FetchData(pbmc, vars = c("PRR11", "celltype","group"))
dat1 <- subset(dat1,celltype != "Endothelial cells")
comparisons <- list(names(table(dat1$group)))

ggplot(dat1, aes(celltype,PRR11, fill = group)) + 
  geom_violin(trim = FALSE,position = position_dodge(width = 0.9))+
  geom_boxplot(position = position_dodge(width = 0.9),
               outlier.size = 0,width = 0.2,show.legend = TRUE)+
  scale_fill_manual(values = c("red","blue"))+
  stat_compare_means(aes(group=group),                       #按分组进行统计检验
                     method = "wilcox.test",
                     paired = FALSE,                             
                     symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                        symbols=c("***","**","*","ns")),
                     label = "p.signif",
                     #label.y = location$Proportion+0.02,      #添加显著性符号的位置
                     size=4.5)+
  xlab("")+  ylab("Expression of PRR11") +   theme_bw()+ylim(c(0,max(dat1$PRR11)))+
  theme(panel.grid = element_line(colour = "white"),plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15,angle = 45,hjust = 1),
        title = element_text(size = 20),
        legend.position ="top",
        legend.key.size = unit(35, "pt"),
        legend.text = element_text(size = 13))+ggfonttheme
ggsave("16.PRR11_boxplot.pdf",device = cairo_pdf,width = 10, height = 7)
ggsave("16.PRR11_boxplot.png",dpi = 300,width = 10, height = 7)

#准备虚拟敲除文件
pbmc_hubcells <- pbmc[,pbmc$celltype%in%c("Epithelial cells")]
saveRDS(pbmc_hubcells,"../RData/UMAP.rds")
