rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/24_job_0780")
if (!dir.exists("01_sc")) {dir.create("01_sc")}
setwd("01_sc")

# 加载包
library(Seurat)
library(dplyr)
library(stringr)
library(patchwork)

# ============================================================================
# GSE173682 单细胞数据读取
# 目标样本: GSM5276933 ~ GSM5276937
# 数据格式: 扁平tar包内直接包含 .tsv.gz / .mtx.gz 文件（非标准10X命名）
# 策略: 解压 → 按样本重命名为标准10X格式 → Read10X → 创建Seurat对象
# ============================================================================
cat("\n========== 1. 识别目标样本 ==========\n")

raw_tar <- "GSE173682_RAW.tar"
target_gsm <- paste0("GSM527693", 3:7)  # GSM5276933 ~ GSM5276937
cat("目标样本:", paste(target_gsm, collapse = ", "), "\n")

# ============================================================================
# 步骤1: 列出tar包中所有文件，筛选目标样本的RNA-seq文件
# ============================================================================
all_files <- untar(raw_tar, list = TRUE)
cat("\ntar包共有", length(all_files), "个文件\n")

# 筛选: 属于目标GSM 且 不是ATAC文件
target_files <- all_files[
  grepl(paste(target_gsm, collapse = "|"), all_files) &
    !grepl("ATAC", all_files)
]
cat("目标文件:", length(target_files), "个\n")
print(target_files)

# ============================================================================
# 步骤2: 提取每个样本的唯一后缀（sample suffix），用于命名
# ============================================================================
# 从文件名中提取后缀: GSM5276933_barcodes-3533EL.tsv.gz → 3533EL
extract_suffix <- function(gsm, file_list) {
  barcode_file <- file_list[grepl(gsm, file_list) & grepl("barcodes", file_list)]
  if (length(barcode_file) == 0) stop("未找到barcodes文件: ", gsm)
  suffix <- gsub(paste0(gsm, "_barcodes-"), "", barcode_file)
  suffix <- sub(".tsv.gz", "", suffix, fixed = TRUE)
  return(suffix)
}

suffix_map <- sapply(target_gsm, extract_suffix, file_list = target_files, USE.NAMES = TRUE)
sample_suffixes <- as.character(suffix_map)
names(sample_suffixes) <- target_gsm

cat("\n样本识别结果:\n")
for (i in seq_along(sample_suffixes)) {
  cat("  ", names(sample_suffixes)[i], "→ suffix =", sample_suffixes[i], "\n")
}

# ============================================================================
# 步骤3: 解压目标文件到临时目录
# ============================================================================
tmp_dir <- file.path("tmp_GSE173682")
if (!dir.exists(tmp_dir)) { dir.create(tmp_dir, recursive = TRUE) }

cat("\n解压", length(target_files), "个目标文件...\n")
untar(raw_tar, files = target_files, exdir = tmp_dir)
cat("解压完成。临时目录文件数:", length(list.files(tmp_dir)), "\n")

# ============================================================================
# 步骤4: 将文件重组为标准10X格式，创建Seurat对象
# ============================================================================
cat("\n========== 2. 读取10X数据并创建Seurat对象 ==========\n")

seurat_list <- list()

for (i in seq_along(target_gsm)) {
  gsm <- target_gsm[i]
  suffix <- sample_suffixes[gsm]
  
  cat("\n--- 处理样本:", gsm, "(gsm:", gsm, ") [", i, "/", length(target_gsm), "] ---\n")
  
  # 创建该样本的标准10X目录
  sample_dir <- file.path(tmp_dir, suffix)
  if (!dir.exists(sample_dir)) { dir.create(sample_dir, recursive = TRUE) }
  
  # 拷贝并重命名为标准10X文件名
  src_barcodes <- file.path(tmp_dir, paste0(gsm, "_barcodes-", suffix, ".tsv.gz"))
  src_features <- file.path(tmp_dir, paste0(gsm, "_features-", suffix, ".tsv.gz"))
  src_matrix   <- file.path(tmp_dir, paste0(gsm, "_matrix-",   suffix, ".mtx.gz"))
  
  dst_barcodes <- file.path(sample_dir, "barcodes.tsv.gz")
  dst_features <- file.path(sample_dir, "features.tsv.gz")
  dst_matrix   <- file.path(sample_dir, "matrix.mtx.gz")
  
  stopifnot(file.exists(src_barcodes))
  stopifnot(file.exists(src_features))
  stopifnot(file.exists(src_matrix))
  
  file.copy(src_barcodes, dst_barcodes)
  file.copy(src_features, dst_features)
  file.copy(src_matrix,   dst_matrix)
  
  cat("  标准10X文件已创建:", paste(list.files(sample_dir), collapse = ", "), "\n")
  
  # 读取10X数据
  counts <- Read10X(data.dir = sample_dir)
  cat("  原始: 基因数:", nrow(counts), "  细胞数:", ncol(counts), "\n")
  
  # 创建Seurat对象
  seurat_list[[gsm]] <- CreateSeuratObject(
    counts      = counts,
    project     = gsm,
    min.cells   = 3,
    min.features = 200
  )
  cat("  过滤后: 基因数:", nrow(seurat_list[[gsm]]),
      "  细胞数:", ncol(seurat_list[[gsm]]), "\n")
}

# ============================================================================
# 步骤5: 查看各样本细胞数
# ============================================================================
cat("\n========== 3. 各样本细胞数汇总 ==========\n")
print(sapply(seurat_list, ncol))

# ============================================================================
# 步骤6: 合并所有样本
# ============================================================================
cat("\n========== 4. 合并所有样本 ==========\n")
merged_obj <- merge(
  seurat_list[[1]],
  y        = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project  = "GSE173682"
)
cat("合并后的对象:\n")
print(merged_obj)
cat("\n合并后各样本细胞数:\n")
print(table(merged_obj$orig.ident))

# ============================================================================
# 步骤7: 保存原始数据
# ============================================================================
saveRDS(seurat_list, file = "../RData/R.00_scRNA_raw_list.Rds")
saveRDS(merged_obj,   file = "../RData/R.01_scRNA_raw_merged.Rds")
cat("\n原始数据已保存。\n")

# ============================================================================
# 步骤8: 清理临时文件
# ============================================================================
cat("\n清理临时文件...\n")
unlink(tmp_dir, recursive = TRUE)
cat("临时文件已清理。\n")

cat("\n========== GSE173682数据读取完成 ==========\n")

#01.质控---------------------------
rm(list = ls()); gc()
scRNA <- readRDS("../RData/R.01_scRNA_raw_merged.Rds")
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
length(colnames(scRNA))#[1] 35868
length(rownames(scRNA))#[1] 26145
#绘图风格
ggfonttheme <- theme(text = element_text(family = "Times"))
theme.set = theme(
  axis.title.x=element_blank(),
  axis.title = element_text(size = 20),
  axis.text.x = element_text(size = 14,angle = 90),
  axis.text.y = element_text(size = 14),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 18))
# 质控前小提琴图
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")
plots<- list()
scRNA$sample <- scRNA@meta.data$orig.ident
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA,
                       pt.size = 0,
                       group.by = "sample",
                       features = plot.featrures[i]) + theme.set + NoLegend()+ggfonttheme
  # theme(axis.text.x = element_blank())
}

violin <- wrap_plots(plots = plots, nrow=1,ncol = 3,cols = palettes(category = "random",22,show_col=F))
violin
ggsave("01.vlnplot_before_qc.pdf", plot = violin, width = 18, height = 6)
ggsave("01.vlnplot_before_qc.png", plot = violin, width = 18, height = 6,dpi = 600)
#  设置阈值
pctMT <- 15
minGene <- 200
minUMI <- 1000
maxGene <- 4000
maxUMI <- 10000
scRNA <- subset(scRNA, subset =
                  nCount_RNA < maxUMI &
                  nCount_RNA > minUMI &
                  nFeature_RNA > minGene &
                  nFeature_RNA < maxGene &
                  percent.mt < pctMT)
length(colnames(scRNA))
# [1] 21597
length(rownames(scRNA))
# [1] 26145
# 质控后小提琴图
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA,
                       pt.size = 0,
                       group.by = "sample",
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
plot1 <- VariableFeaturePlot(object = scRNA.norm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(file="03.feature_selection.pdf",width=8,height=6,onefile=F,family = "Times")
print(plot2)
dev.off()
#保存基因特征方差图
png(file="03.feature_selection.png",width=8,height=6,units='in',res=600,family = "Times")
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
DimPlot(scRNA.norm.pca, reduction = "pca",group.by = "sample")
dev.off()
png(file="05.pcaPlot.png",width=9,height=6,res=600,units='in',family = "Times")
DimPlot(scRNA.norm.pca, reduction = "pca",group.by = "sample")
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
pcSelect <- 30

scRNA.norm.umap <- FindNeighbors(object = scRNA.norm.pca, dims = 1:pcSelect)
scRNA.norm.umap <- FindClusters(object = scRNA.norm.umap, resolution = 0.5)#19
##umap聚类
scRNA <- RunUMAP(scRNA.norm.umap,
                 reduction = "pca",
                 dims = 1:pcSelect,
                 min.dist = 0.7, #参数越大,每个cluster的细胞越分散
                 n.neighbors = 160 #参数越大，局部布局会更紧凑
)

saveRDS(scRNA, file = "../RData/R.04_scRNA_umap.Rds")

p4 <- DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters",
              label = TRUE, pt.size = 0.1)+ggfonttheme
p4
ggsave("07.pca_umap.pdf",device = cairo_pdf,width = 8, height = 8)
ggsave("07.pca_umap.png",dpi = 300,width = 8, height = 8)
system.time(saveRDS(scRNA, file = "../RData/R.04_scRNA_umap_04_30PC.Rds"))
#06.查找marker基因----------------------------------------------------------------
scRNA.norm.umap <- readRDS("../RData/R.04_scRNA_umap_04_30PC.Rds")
DefaultAssay(scRNA.norm.umap) <- "RNA"
library("multtest")
library("metap")
scRNA.norm.umap <- JoinLayers(scRNA.norm.umap)
nk.markers <- FindAllMarkers(scRNA.norm.umap, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(nk.markers)
dim(nk.markers)
#去重统计
x <- nk.markers$gene
length(x[!duplicated(x)])
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
##PMID: 37081869，39589706，40540226
cell_markers <- list(
  T_cells = c("CD3D", "CD3E"),
  Epithelial_cells = c("CDH1", "EPCAM"),
  Macrophages = c("CD68","TYROBP"),
  Fibroblasts = c("CLMP","SFRP4"),
  Endothelial_cells = c("KLF2","CLDN5"),
  Smooth_muscle_cells = c("ACTA2","RGS5"),
  Ciliated_cells = c("HMGN3","CAPS"),
  B_cells = c("MZB1", "IGKC")
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
  "0" = "Fibroblasts",
  "1" = "T cells",
  "2" = "Fibroblasts",
  "3" = "T cells",
  "4" = "Endothelial cells", 
  "5" = "Smooth muscle cells",
  "6" = "Fibroblasts",
  "7" = "Smooth muscle cells",
  "8" = "Fibroblasts",
  "9" = "Macrophages",
  "10" = "Fibroblasts",
  "11" = "Epithelial cells",
  "12" = "Smooth muscle cells",
  "13" = "T cells",
  "14" = "Epithelial cells",
  "15" = "Macrophages",
  "16" = "Ciliated cells",
  "17" = "B cells",
  "18" = "Endothelial cells"
)
saveRDS(scRNA.norm.umap, file = "../RData/R.05_scRNA_umap_05_anno.Rds")

scRNA.norm.umap <- readRDS("../RData/R.05_scRNA_umap_05_anno.Rds")
cluster_dat <- data.frame("cluster" = 0:(length(new.cluster.id)-1),"CellType"=new.cluster.id)
write.xlsx(cluster_dat, file = "12.anno_CellType.xlsx")
names(new.cluster.id) <- levels(scRNA.norm.umap)
scRNA.norm.umap <- RenameIdents(scRNA.norm.umap, new.cluster.id)
scRNA.norm.umap$celltype <- Idents(scRNA.norm.umap)
DotPlot(scRNA.norm.umap, features = split(unique(cell_markers_df$gene), cell_markers_df$cell_type),
        cols = c("lightgrey", "red"))+
  theme(axis.text.x=element_text(angle=90, hjust = 1),
        panel.border = element_rect(color="black"), #面板边框
        panel.spacing = unit(1, "mm")#面板间距
  )+ggfonttheme
ggsave("12.celltype_markergenes_dotplot.pdf",device = cairo_pdf,width = 15, height = 8)
ggsave("12.celltype_markergenes_dotplot.png",dpi = 300,width = 15, height = 8)
# scRNA.norm.umap <- SetIdent(scRNA.norm.umap, value = "celltype")
ncluster <- length(unique(scRNA.norm.umap[[]]$celltype))
mycol <- colorRampPalette(brewer.pal(6, "Set3"))(ncluster)
# DimPlot(scRNA.norm.umap,
#         label = T , repel = T,
#         label.size = 5,
#         label.box = T,
#         cols = mycol) +
#   NoLegend()
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
#关键细胞鉴定----------------------------------------
scRNA.norm.umap <- readRDS("../RData/R.06_scRNA_done.Rds")
pbmc <- scRNA.norm.umap
TGFB1_expression <- FetchData(pbmc, vars = "TGFB1")
head(TGFB1_expression)
pbmc$TGFB1_expression <- TGFB1_expression$TGFB1
#TGFB1表达UMAP图------------------------
# 提取 UMAP 坐标和细胞类型信息
umap_coords <- Embeddings(pbmc, reduction = "umap")
celltype <- pbmc$celltype
# 计算每个细胞类型的平均 UMAP 坐标（质心）
centroids <- aggregate(umap_coords, by = list(celltype), FUN = mean)
colnames(centroids)[2:3] <- c("UMAP_1", "UMAP_2")  # 重命名坐标列
# 绘制 TGFB1 表达 UMAP
p <- FeaturePlot(pbmc, 
                 features = "TGFB1", 
                 reduction = "umap",
                 pt.size = 0.5,
                 order = TRUE) +
  ggtitle("TGFB1 Expression") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
  ggfonttheme
library(ggrepel)
p + geom_label_repel(data = centroids,
                     aes(x = UMAP_1, y = UMAP_2, label = Group.1),
                     size = 5,
                     fill = "white",
                     alpha = 0.7,
                     box.padding = 0.5,
                     point.padding = 0.3,
                     segment.color = "black")
ggsave("17.TGFB1_expression.pdf",device = cairo_pdf,width = 7, height = 7)
ggsave("17.TGFB1_expression.png",dpi = 300,width = 7, height = 7)
#准备虚拟敲除文件
pbmc_hubcells <- pbmc[,pbmc$celltype%in%c("Fibroblasts")]
saveRDS(pbmc_hubcells,"../RData/UMAP_Fibroblasts.rds")
