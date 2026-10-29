rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/16_project_1449")
if (!dir.exists("04_hdWGCNA")) {dir.create("04_hdWGCNA")}
setwd("04_hdWGCNA")
library(WGCNA)
library(hdWGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)
library(qs)
library(Seurat)
library(harmony)
#load data-----------------------------------------
set.seed(42)
enableWGCNAThreads(nThreads = 32)
scRNA <- readRDS("../RData/R.07_scRNA_Scissor.RData")
scRNA <- SeuratObject::UpdateSeuratObject(scRNA)
#harmony降维----------------------------------------
scRNA <- RunHarmony(
  scRNA,
  group.by.vars = "orig.ident",
  assay.use = "RNA"
)
saveRDS(scRNA, file = "../RData/R.08_scRNA_harmony.RData")
names(scRNA@reductions)
#01.prepare-----------------------------------------
seurat_obj <- SetupForWGCNA(
  scRNA,
  gene_select = "fraction", # 默认fraction;其他有variable(Seurat对象中存储的基因),custom(自定义)
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "GSE148071" # the name of the hdWGCNA experiment
)
Idents(seurat_obj) <- "subcelltype"
#02.构建metacells-------------------------------------------------
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("subcelltype", "group"), #指定seurat_obj@meta.data中要分组的列
  reduction = 'harmony', # 选择要执行KNN的降维
  k = 25, # 最近邻居参数
  max_shared = 10, # 两个metacell之间共享细胞的最大数目
  ident.group = 'subcelltype', # 设置metacell安全对象的标识
  min_cells = 100 # 排除数量小于100的细胞亚群
)
seurat_obj <- NormalizeMetacells(seurat_obj)
#03.共表达网络分析-------------------------------------------------
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c("Submucosal_gland_cell","Secretory_cell"), # the name of the group of interest in the group.by column
  group.by='subcelltype', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2)),
  networkType = 'unsigned' # you can also use "unsigned" or "signed hybrid"
)
plot_list <- PlotSoftPowers(seurat_obj)

pdf(file = "01.softpowers.pdf",width = 10,height =10,family="Times")
wrap_plots(plot_list, ncol=2)
dev.off()
png("01.softpowers.png",height=10,width=10,family='Times',units='in',res=600)
wrap_plots(plot_list, ncol=2)
dev.off()


power_table <- GetPowerTable(seurat_obj)
head(power_table)
# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = 6,
  tom_outdir = "TOM",
  tom_name = "Submucosal_gland_cell", # name of the topoligical overlap matrix written to disk
  overwrite_tom = TRUE # 允许覆盖已存在的同名文件
)
PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')

pdf(file = "02.tree.pdf",width = 10,height =6,family="Times")
PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
dev.off()

png("02.tree.png",width=10,height=6,family='Times',units='in',res=600)
PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
dev.off()

seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="orig.ident"
)
hMEs <- GetMEs(seurat_obj)
#04.计算模块连接性------------------------------------------
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'subcelltype', 
  group_name = c("Submucosal_gland_cell","Secretory_cell")
)

seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "hdWGCNA_Module"
)

# 绘制每个模块按kME排序的基因
p <- PlotKMEs(seurat_obj, ncol=4)
p

pdf(file = "03.kME.pdf",width = 10,height =8,family="Times")
print(p)
dev.off()
png("03.kME.png",height=10,width=8,family='Times',units='in',res=600)
print(p)
dev.off()

#05.获取模块内部信息---------------------
modules <- GetModules(seurat_obj) %>% 
  subset(module != 'grey')
head(modules[,1:6])
hub_df <- GetHubGenes(seurat_obj, n_hubs = 100)
head(hub_df)
saveRDS(seurat_obj, 'hdWGCNA_object.RData')

#06.计算hub基因siganture得分-------------------------
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell' # Seurat方法(AddModuleScore)
)
# 5.可视化
# 模块特征图
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)
# stitch together with patchwork
wrap_plots(plot_list, ncol=4)


# 每个模块制作hub scores的特征图
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)
# stitch together with patchwork
wrap_plots(plot_list, ncol=4)
pdf(file = "04.ModuleFeaturePlot.pdf",width = 12,height =8,family="Times")
wrap_plots(plot_list, ncol=4)
dev.off()
png("04.ModuleFeaturePlot.png",height=12,width=8,family='Times',units='in',res=600)
wrap_plots(plot_list, ncol=4)
dev.off()

# 查看模块相关图
ModuleCorrelogram(seurat_obj,method="color",
                  addCoef.col="black",
                  diag=TRUE,
                  tl.col="black",tl.pos = "lt",
                  order = "original")
ModuleCorrelogram(seurat_obj,method = "circle", type = "lower", 
                  tl.col = "black", tl.cex = 0.8, tl.pos = "n",
                  add = TRUE,diag = TRUE,cl.pos = "n",
                  insig = "label_sig",sig.level = c(.001, .01, .05),
                  pch.cex = 1)
p_base <- recordPlot()
pdf(file = "05.ModuleCorrelogram.pdf",width = 10,height =10,family="Times")
p_base
dev.off()
png("05.ModuleCorrelogram.png",height=10,width=10,family='Times',units='in',res=600)
p_base
dev.off()

# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']
# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
Idents(seurat_obj) <- "celltype"
seurat_obj <- subset(x = seurat_obj, idents = c("Epithelial cells"))
# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'subcelltype')
# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')
# plot output
p
pdf("05.wgcna.heatmap.pdf",height = 7,width = 9,family = "Times")
p
dev.off()
png("05.wgcna.heatmap.png",height = 700,width = 900,family = "Times")
p
dev.off()

#选择相关性大于0.8的模块（p<0.05）为关键模块，kME>0.25为核心基因
hubmodules <- paste0("hdWGCNA_Module",c(1,2,4,5,9))
hub_module_gene_df <- subset(hub_df,kME > 0.25)
hub_module_gene_df <- hub_module_gene_df[hub_module_gene_df$module %in% hubmodules,]
write.xlsx(hub_module_gene_df,"06.hub_module_gene_df.xlsx")
#477个关键基因
save.image("../RData/04_hdWGCNA.RData")
