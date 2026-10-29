#本分析在阿里云服务器上完成

rm(list = ls())
gc()
setwd(dir = "<REFERENCE_DATA_ROOT>/job/project004")


library(dplyr)
library(Seurat)
library(scTenifoldKnk)
library(ggplot2)
library(ggrepel) 


# nohup Rscript <REFERENCE_DATA_ROOT>/REDACTED_PATH > log.txt 2>&1 &

scRNA_harmony <- readRDS("UMAP.rds")
table(Idents(scRNA_harmony))
table(Idents(scRNA_harmony))



scRNA_harmony <- FindVariableFeatures(
  scRNA_harmony,
  selection.method = "vst",
  nfeatures = 2000
)

hvg <- VariableFeatures(scRNA_harmony)




# 2. 强制加入 KO gene
ko_gene1 <- "SIRT7"
genes_use1 <- unique(c(hvg, ko_gene1))


# 提取count矩阵作为输入
countMatrix1 <- GetAssayData(
  scRNA_harmony,
  layer = "counts"   
)[genes_use1, ]

ko_gene1 %in% rownames(countMatrix1)

# 执行敲除分析（推荐用高变基因做）：
# 参数gKO用于指定敲除的基因；qc开头的系列参数能够选择是否对矩阵进行质控以及相应指标；nc开头的参数对应网络的构建，
# 分别包括子网络的数量、使用的PC数量、随机抽取的细胞数量，数值越大运行时间越长但是结果可能会更稳定



#############基因敲除加速#####################


# 1.主要修改：将scTenifoldKnk函数改为了scTenifoldNet.lite::scTenifoldKnk
# 2.fast = TRUE


# 原理
# 
# 细胞抽样：随机抽取nCells个细胞（默认500），抽nNet次（默认10次）。分别用这10个抽样的表达矩阵，
# 按照下面的方法计算基因调控网络（GRN）。
# 
# 主成分回归（Principle Component Regression，PCR），对每个基因计算其受其他基因调控的回归方程，
# 其系数即可作为调控网络的系数。合并所有系数即可得到GRN。
# 
# 原版软件PCR过程中对每个基因都做了一次奇异值分解（Singular value decomposition，SVD），
# 如果有10000个基因，就要重复计算10000次，非常耗时。我在极速模式中改成了对所有基因一次性计算全局SVD
# ，直接让计算次数减少几千倍。
# 
# 降噪：对10次抽样得到的GRN合并成一个基因x基因x抽样次数的三维张量，
# 进行张量分解（CANDECOMP/PARAFAC tensor decomposition），得到降噪后的GRN，即WT的调控网络。
# 
# 在WT网络中把需要敲除的基因入度设为0，即可模拟敲除效应，得到KO网络。
# 
# 将WT网络和KO网络进行非线性降维，为网络中的每个基因对齐低维空间坐标。
# 原版算法在这里构建了一个2N * 2N大小的矩阵，非常占用内存！我借助RSpectra::svds的优化算法解决了这个问题。
# 
# 在低维空间中，计算每个基因在两个网络之间的欧氏距离，作为敲除事件的扰动效应。

############################################


result <- scTenifoldNet.lite::scTenifoldKnk(countMatrix1, gKO = 'SIRT7', qc = FALSE, nCores = 8, fast = TRUE)



# result <- scTenifoldKnk(countMatrix = countMatrix1, 
#                         gKO = 'ACTA2', #需要敲除的基因
#                         # qc = TRUE,#是否进行QC
#                         # qc_mtThreshold = 0.1,#mt阈值
#                         # qc_minLSize = 1000,#文库阈值(细胞测到的基因总数)
#                         # qc_minLSize = 10,#文库阈值(细胞测到的基因总数)
#                         # nc_nNet = 10, #子网络数量
#                         # nc_nCells = 500, #每个网络中随机抽取的细胞数
#                         # nc_nComp = 3 #PCA 的主成分数量
# )

saveRDS(result, file = "scTenifoldKnk_SIRT7.rds")

#可视化-------------------
rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project004")
if (!dir.exists("03_scTenifoldKnk")) {dir.create("03_scTenifoldKnk")}
setwd("03_scTenifoldKnk")
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(DOSE)
library(ggnewscale)
library(org.Hs.eg.db)
ggfonttheme <- theme(text = element_text(family = "Times"))
perturb_inf <- readRDS("scTenifoldKnk_SIRT7.rds")
perturb_dat <-perturb_inf$diffRegulation
perturb_dat$log2FoldChange <- log2(perturb_dat$FC)
diff_perturb_dat <- perturb_dat %>%
  filter(p.value < 0.05)
diff_perturb_genes <- diff_perturb_dat$gene
write.table(diff_perturb_genes,"01.diff_perturb_genes.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)
# 按log2FoldChange从大到小排序，取前10个基因
top10_genes <- head(perturb_dat[order(perturb_dat$log2FoldChange, decreasing = TRUE), ], 10)
# 绘制散点图
ggplot(perturb_dat, aes(x = Z, y = log2FoldChange)) +
  geom_point(color = "grey60", alpha = 0.6, size = 1.5) +
  # 高亮top10基因
  geom_point(data = top10_genes, aes(x = Z, y = log2FoldChange), 
             color = "red", size = 2.5) +
  # 标注基因名
  geom_text_repel(data = top10_genes, aes(label = gene),
                  color = "red", size = 4, 
                  max.overlaps = 20,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = "grey50",
                  family = "Times") +
  labs(x = "Z-score", 
       y = expression(log[2](FoldChange))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))+ggfonttheme
ggsave("01.volcano.pdf",device = cairo_pdf,width = 12, height = 8)
ggsave("01.volcano.png",dpi = 300,width = 12, height = 8)
#富集分析----------------------------------------------
options(scipen = 5)
hubgenes <- diff_perturb_genes
GO_database <- 'org.Hs.eg.db'
KEGG_database <- 'hsa'
gene <- bitr(hubgenes, fromType = "SYMBOL", toType="ENTREZID",OrgDb = org.Hs.eg.db)
GO<-enrichGO(gene$ENTREZID,#GO富集分析
             OrgDb = GO_database,
             keyType = "ENTREZID",#设定读取的gene ID类型
             ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
             pvalueCutoff = 0.05,#设定p值阈值
             qvalueCutoff = 0.05,#设定q值阈值
             readable = TRUE)
KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
KEGG = setReadable(KEGG, #前面分析的结果
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID")
KEGG_res <- subset(KEGG@result,p.adjust < 0.05)#65
GO_res <- GO@result #1340
library(openxlsx)
write.xlsx(GO_res,"01.GO.xlsx")
write.xlsx(KEGG_res,"01.KEGG.xlsx")
#可视化
#KEGG
df <- KEGG_res
df <- df[order(df$p.adjust), ]
df <- df[1:10,]
df$labelx=rep(0,nrow(df))
df$labely=seq(nrow(df),1)
p1 <- ggplot(data = df, 
             aes(Count, reorder(Description,p.adjust,decreasing = TRUE),fill = p.adjust)) +
  geom_bar(stat="identity",
           alpha=0.5,
           width = 0.8) + 
  geom_text(aes(x=labelx,
                y=labely,
                label = paste(" ",Description)),
            size=6, 
            hjust =0,family = "Times")+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_line(colour = 'black', linewidth = 1),
        axis.title.y = element_blank(),
        axis.ticks.y = element_line(colour = 'black', linewidth = 1),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 20),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 20),
        title = element_text(size = 20))+
  xlab("Numbers of Genes")+
  ggtitle("KEGG Pathways of Perturb Genes")+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_gradient(low ="#BC0858",high = "#C192A8")+ggfonttheme
p1
#GO_BP
df <- subset(GO_res,ONTOLOGY == "BP")#442
df <- df[order(df$p.adjust), ]
df <- df[1:10,]
df$labelx=rep(0,nrow(df))
df$labely=seq(nrow(df),1)
p2 <- ggplot(data = df, 
             aes(Count, reorder(Description,p.adjust,decreasing = TRUE),fill = p.adjust)) +
  geom_bar(stat="identity",
           alpha=0.5,
           width = 0.8) + 
  geom_text(aes(x=labelx,
                y=labely,
                label = paste(" ",Description)),
            size=6, 
            hjust =0,family = "Times")+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_line(colour = 'black', linewidth = 1),
        axis.title.y = element_blank(),
        axis.ticks.y = element_line(colour = 'black', linewidth = 1),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 20),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 20),
        title = element_text(size = 20))+
  xlab("Numbers of Genes")+
  ggtitle("GO BP Pathways of Perturb Genes")+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_gradient(low = "#013494" ,high ="#BACFF5")+ggfonttheme
p2
#GO_CC
df <- subset(GO_res,ONTOLOGY == "CC")#43
df <- df[order(df$p.adjust), ]
df <- df[1:10,]
df$labelx=rep(0,nrow(df))
df$labely=seq(nrow(df),1)
p3 <- ggplot(data = df, 
             aes(Count, reorder(Description,p.adjust,decreasing = TRUE),fill = p.adjust)) +
  geom_bar(stat="identity",
           alpha=0.5,
           width = 0.8) + 
  geom_text(aes(x=labelx,
                y=labely,
                label = paste(" ",Description)),
            size=6, 
            hjust =0,family = "Times")+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_line(colour = 'black', linewidth = 1),
        axis.title.y = element_blank(),
        axis.ticks.y = element_line(colour = 'black', linewidth = 1),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 20),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 20),
        title = element_text(size = 20))+
  xlab("Numbers of Genes")+
  ggtitle("GO CC Pathways of Perturb Genes")+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_gradient(low ="#cb1dc2",high = "#f1afee")+ggfonttheme
p3
#GO_MF
df <- subset(GO_res,ONTOLOGY == "MF")#60
df <- df[order(df$p.adjust), ]
df <- df[1:10,]
df$labelx=rep(0,nrow(df))
df$labely=seq(nrow(df),1)
p4 <- ggplot(data = df, 
             aes(Count, reorder(Description,p.adjust,decreasing = TRUE),fill = p.adjust)) +
  geom_bar(stat="identity",
           alpha=0.5,
           width = 0.8) + 
  geom_text(aes(x=labelx,
                y=labely,
                label = paste(" ",Description)),
            size=6, 
            hjust =0,family = "Times")+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_line(colour = 'black', linewidth = 1),
        axis.title.y = element_blank(),
        axis.ticks.y = element_line(colour = 'black', linewidth = 1),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 20),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 20),
        title = element_text(size = 20))+
  xlab("Numbers of Genes")+
  ggtitle("GO MF Pathways of Perturb Genes")+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_gradient(low = "#fa9000" ,high ="#f2ce9d")+ggfonttheme
p4
#拼图
library(patchwork)
(p1+p2)/(p3+p4)+plot_layout(guides = 'collect')
ggsave("02.enrichment.pdf",device = cairo_pdf,width = 20, height = 20)
ggsave("02.enrichment.png",dpi = 300,width = 22, height = 20)
save.image("../RData/03_scTenifoldKnk.RData")
