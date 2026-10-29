rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project002")
if (!dir.exists("04_scTenifoldKnk")) {dir.create("04_scTenifoldKnk")}
setwd("04_scTenifoldKnk")
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
perturb_inf <- readRDS("../04_scTenifoldKnk/scTenifoldKnk_10000HVG_PRR11.rds")
perturb_dat <-perturb_inf$diffRegulation
perturb_dat$log2FoldChange <- log2(perturb_dat$FC)
diff_perturb_dat <- perturb_dat %>%
  filter(p.value < 0.05)
diff_perturb_genes <- perturb_dat$gene
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
KEGG_res <- subset(KEGG@result,p.adjust < 0.05)#19
GO_res <- GO@result #351
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
  ggtitle("KEGG Pathways of Candidate Genes")+
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
  ggtitle("GO BP Pathways of Candidate Genes")+
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
           width = 0.4) + 
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
  ggtitle("GO CC Pathways of Candidate Genes")+
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
  ggtitle("GO MF Pathways of Candidate Genes")+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_gradient(low = "#fa9000" ,high ="#f2ce9d")+ggfonttheme
p4
#拼图
library(patchwork)
(p1+p2)/(p3+p4)+plot_layout(guides = 'collect')
ggsave("02.enrichment.pdf",device = cairo_pdf,width = 20, height = 20)
ggsave("02.enrichment.png",dpi = 300,width = 22, height = 20)
save.image("../RData/04_scTenifoldKnk.RData")

