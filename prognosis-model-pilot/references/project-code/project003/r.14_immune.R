rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project003")
if (!dir.exists("14_immune")) {dir.create("14_immune")}
setwd("14_immune")
library(readxl)
library(dplyr)
library(dplyr)
library(tidyr)
library(networkD3)
library(openxlsx)
ggfonttheme <- theme(text = element_text(family = "Times"))
groupInf <- read.table("../04_evaluate/04.risk_table.txt",header = TRUE,row.names = 1,sep = "\t")
group1 <- data.frame(groupInf$risk_list_cut,row.names = rownames(groupInf))
groupInf$sample <- rownames(groupInf)
expr <- readRDS("../00_rawdata/tpm.Rds")
expr <- t(expr) %>% as.data.frame()
df1 <- merge(group1,expr,by = 0)
rownames(df1) <- df1[,1]
df1 <- df1[,-1]
group_list <- as.factor(df1[,1])
df1 <- df1[,-1]
data <- df1
data_matrix <- as.matrix(data)
data_matrix <- matrix(as.numeric(data_matrix), nrow=nrow(data), dimnames=list(rownames(data), colnames(data)))
#基因集
library(GSVA)
gene_set <- read.table("mmc3_ssGSEA.txt",
                       header = T,
                       sep ="\t")
dat.final <- as.matrix(t(data_matrix))
gene_list <- split(as.matrix(gene_set)[,1],
                   gene_set[,2])
#ssGSEA
gsvapar <- gsvaParam(as.matrix(dat.final), gene_list, maxDiff = TRUE, kcdf = "Gaussian")
ssgsea_score = gsva(gsvapar, verbose = TRUE)
ssgsea_score1 <- as.data.frame(ssgsea_score)
write.csv(ssgsea_score1,
          file = "01.ssgsea_result_cell.csv",
          quote = F)
ssgsea_score<-ssgsea_score[,rownames(groupInf)]
annotation_col<-as.data.frame(groupInf$risk_list_cut)
colnames(annotation_col)='Group'
rownames(annotation_col)=colnames(ssgsea_score)
ssgsea_score <- as.matrix(ssgsea_score)
#热图
library(pheatmap)
ann_colors<-list(Group = c('low'="#3777ac",'high'="#a64545"))
color<- colorRampPalette(c("navy", "white", "firebrick3"))(100)
p <- pheatmap(
  ssgsea_score,
  border_color = ann_colors,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  labels_row = NULL,
  clustering_method = 'ward.D2',
  show_rownames = T,
  show_colnames = F,
  fontsize_col = 5,
  cluster_cols = F,
  cluster_rows = T,
  color=color)
p
png(file = "01.heatmap.png", height = 7, width = 8, units = 'in',res = 600,family='Times')
p
dev.off()
pdf(file = "01.heatmap.pdf", height = 7,width = 8,family='Times')
p
dev.off()
#设置分组
group_case<-groupInf[groupInf$risk_list_cut=="high",]
group_case<-as.character(group_case$sample)
group_control<-groupInf[groupInf$risk_list_cut=="low",]
group_control<-as.character(group_control$sample)
#计算p值
library(dplyr)
tiics_result <- read.csv('./01.ssgsea_result_cell.csv',check.names = F, row.names = 1) 
tiics_result <- tiics_result[, groupInf$sample] %>% as.matrix()
pvalue = padj <- matrix(0, nrow(tiics_result), 1)
for (i in 1:nrow(tiics_result)){
  pvalue[i, 1] = p.value = wilcox.test(tiics_result[i, group_case],
                                       tiics_result[i, group_control])$p.value
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(#log2FoldChange, 
  pvalue, 
  padj,
  row.names = rownames(tiics_result))
rTable$immune_cell <- rownames(rTable)
rTable$sig <- ifelse(rTable$pvalue < 0.05,
                     ifelse(rTable$pvalue < 0.01, 
                            ifelse(rTable$pvalue < 0.001,
                                   ifelse(rTable$pvalue < 0.0001,
                                          paste(rTable$immune_cell, "****",  sep = ""),
                                          paste(rTable$immune_cell, "***", sep = "")),
                                   paste(rTable$immune_cell, "**", sep = "")),
                            paste(rTable$immune_cell, "*",  sep = "")), 
                     rTable$immune_cell)
diff_Table<-rTable[which(rTable$pvalue<0.05),]
dim(diff_Table)
# [1] 20  4
write.csv(rTable,
          file = "02.cell_all_wilcox_test.csv",
          quote = F,
          row.names = F)
write.csv(diff_Table,
          file = "02.cell_diff_wilcox_test.csv",
          quote = F,
          row.names = F)

paste0(diff_Table$immune_cell,collapse = "()、")
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)

cell.data <- data.frame(Immune_Cell=rownames(tiics_result), 
                        tiics_result, 
                        pvalue=rTable$pvalue)
# plot.cell <- cell.data[which(cell.data$pvalue<0.05),]
plot.cell <- cell.data
diff_tiics <- rownames(plot.cell)
box_dat <- gather(plot.cell, key=Group, value=score, -c("Immune_Cell","pvalue"))
box_dat$Group <- ifelse(gsub("\\.","-",box_dat$Group) %in% group_control,
                        "low", "high") 
box_dat$Group <- factor(box_dat$Group, levels = c("low", "high"))
head(box_dat)
library(ggplot2)
library(ggpubr)
boxplot_diff_TIICs <- ggplot(box_dat, aes(x = Immune_Cell, 
                                          y = score,
                                          fill = Group)) +
  geom_boxplot(width = 0.7, 
               position = position_dodge(0.9), 
               outlier.shape = NA, 
               color = 'black') +
  scale_fill_manual(values = c("#377FD4", "#ffdd93")) + 
  labs(title = "Cell Infiltration", x = "", y = "Score") +
  stat_compare_means(data =box_dat,
                     method = "wilcox.test",
                     mapping = aes(group = Group),
                     label = "p.signif",
                     hide.ns = F) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 16), 
        axis.text.x = element_text(angle = 50, hjust = 1, colour = "black", face = "bold", size = 12), 
        axis.text.y = element_text(hjust = 0.5, colour = "black", size = 12), 
        axis.title.x = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 14, face = "bold"), 
        legend.text = element_text(face = "bold", hjust = 0.5, colour = "black", size = 11), 
        legend.title = element_text(face = "bold", colour = "black", size = 11), 
        text = element_text(family = 'Times'), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
print(boxplot_diff_TIICs)
ggsave(filename = "02.ssgsea_Box.png", height =6, width =17, boxplot_diff_TIICs)
ggsave(filename = "02.ssgsea_Box.pdf", height = 6, width = 17, boxplot_diff_TIICs)
dev.off()
#相关性
library(ggcorrplot)
res3 <- read.csv('01.ssgsea_result_cell.csv',check.names = F,row.names = 1)
dif_cell <- read.csv('02.cell_diff_wilcox_test.csv')
res3 <- res3[rownames(res3) %in% dif_cell$immune_cell,]
res3 <- t(res3) %>% as.data.frame()
cor_data <- cor(res3,method="spearman")
corp <- ggcorrplot::cor_pmat(res3)
dat.r <- data.frame(cor_data)
dat.r$cell <- rownames(dat.r)
dat.r <- gather(dat.r,'Gene','r',-c(cell))
dat.p <- data.frame(corp)
dat.p$cell <- rownames(dat.p)
dat.p <- gather(dat.p,'Gene','p',-c(cell))
corrp <- (cbind(dat.r, dat.p))[,c("cell","Gene",'r',"p")]
colnames(corrp) <- c("cell1","cell2",'Correlation',"Pvalue")
corrp$cell2 <- gsub("\\.", " ",corrp$cell2)
write.csv(corrp,file="03.cell_cor.csv",quote=F)
library(ggcorrplot)
library(corrplot)
col1 <- colorRampPalette(c("#3777ac", "white","#FC8D62"))

pdf(file = "03.cor_heatmap.pdf",width = 14,height =14,family="Times")
corrplot(cor_data,
         method="color",
         addCoef.col="black",
         diag=FALSE,
         tl.col="black",tl.pos = "lt",
         col=col1(20),
         order = "original")
corrplot(cor_data, method = "circle", type = "lower",col = col1(20), 
         tl.col = "n", tl.cex = 0.8, tl.pos = "n",add = TRUE,diag = FALSE,cl.pos = "n",
         p.mat = corp, insig = "label_sig",sig.level = c(.001, .01, .05),pch.cex = 1.5)
dev.off()
png("03.cor_heatmap.png",height=14,width=14,family='Times',units='in',res=600)
corrplot(cor_data,
         method="color",
         addCoef.col="black",
         diag=FALSE,
         tl.col="black",tl.pos = "lt",
         col=col1(20),
         order = "original")
corrplot(cor_data, method = "circle", type = "lower",col = col1(20), 
         tl.col = "n", tl.cex = 0.8, tl.pos = "n",add = TRUE,diag = FALSE,cl.pos = "n",
         p.mat = corp, insig = "label_sig",sig.level = c(.001, .01, .05),pch.cex = 1.5)
dev.off()
#预后基因与免疫细胞相关性
prog_genes <- read.table("../03_model/10.multiresult.txt") %>% rownames()
risk_mat <- data.frame(sample = rownames(df1))
risk_mat <- cbind(risk_mat,df1[,prog_genes])
res3 <- merge(res3,risk_mat,by.x = 0,by.y = 1)
rownames(res3) <- res3[,1]
res3 <- res3[,-1]
cor_data <- cor(res3,method="spearman")
corp <- cor_pmat(res3)
genes_count <- length(prog_genes)
cells_count <- dim(corp)[1] - genes_count
cmt <- cor_data[1:cells_count,(cells_count+1):(cells_count+genes_count)]
pmt <- corp[1:cells_count,(cells_count+1):(cells_count+genes_count)]
dput(colnames(corp)[1:cells_count])
write.xlsx(cmt,"05.genes_immunecell_cor.xlsx",rowNames = TRUE)
write.xlsx(pmt,"06.gene_immunecell_corp.xlsx",rowNames = TRUE)
# c("Activated B cell", "Activated CD8 T cell", "Activated dendritic cell", 
#   "CD56bright natural killer cell", "CD56dim natural killer cell", 
#   "Central memory CD4 T cell", "Effector memory CD8 T cell", "Eosinophil", 
#   "Immature  B cell", "Macrophage", "Mast cell", "MDSC", "Monocyte", 
#   "Natural killer cell", "Natural killer T cell", "Neutrophil", 
#   "Plasmacytoid dendritic cell", "Regulatory T cell", "T follicular helper cell", 
#   "Type 1 T helper cell")
if (!is.null(pmt)){
  ssmt <- pmt< 0.01 
  pmt[ssmt] <-'**'
  smt <- pmt >0.01& pmt <0.05
  pmt[smt] <- '*'
  pmt[!ssmt&!smt]<- ''
} else {
  pmt <- F
}
pmt <- paste0(round(cmt,2),pmt)
pmt <- matrix(pmt, nrow = nrow(cmt))
cor_matrix <- as.data.frame(t(cmt))
p_matrix <- as.data.frame(t(pmt))
p1 <- pheatmap(cor_matrix,scale = "none",cluster_row = F, cluster_col = F, border=TRUE,
               breaks=seq(-1, 1, length.out = 10),
               display_numbers = p_matrix,
               fontsize_number = 15, number_color = "black",fontsize = 13,
               cellwidth = 50, cellheight =50,
               color=colorRampPalette(c("#3777ac", "white","#FC8D62"))(10),
               angle_col = 45,fontfamily = "Times")
p1
ggsave("07.genes_immunecell_cor.pdf",p1,device = cairo_pdf,width = 18, height = 8)
ggsave("07.genes_immunecell_cor.png",p1,dpi = 300,width = 18 ,height = 8)
save.image("../RData/14_immune.RData")
