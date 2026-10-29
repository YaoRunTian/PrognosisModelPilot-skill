rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project004")
if (!dir.exists("04_DEG")) {dir.create("04_DEG")}
setwd("04_DEG")
library(limma)
library(metagenomeSeq)
library(edgeR)
groupInf <- readRDS("../02_SIRT7/SIRT7_group.Rds")
data <- readRDS("../00_rawdata/GSE62452_expr.Rds")
#01.差异分析--------------------------------------------------------------------
expr <- data[,groupInf$sample]
group_list <- factor(groupInf$group,levels = c("Low Expression", "High Expression"))
design <- model.matrix(~group_list)
fit <- lmFit(expr,design)
fit1 <- eBayes(fit)
options(options=4)
library("tibble")
deg <- topTable(fit1,coef=2,n=Inf) %>% 
  rownames_to_column("Genes")
DEG<-deg
DEG$regulate <- ifelse(DEG$P.Value > 0.05, "unchanged",
                       ifelse(DEG$logFC > 0.5, "up-regulated",
                              ifelse(DEG$logFC < -0.5, "down-regulated", "unchanged")))
head(DEG)
table(DEG$regulate)
# down-regulated      unchanged   up-regulated 
# 278          19650             75 
#02.火山图--------------------------------------------------------------------
data<-DEG
data$sign <-NA
data <- data[order(data$P.Value, data$logFC, decreasing = c(FALSE, TRUE)), ]
# 筛选logFC上调的基因
upregulated_genes <- data[data$logFC > 0.5, ]
upregulated_genes <- upregulated_genes[order(upregulated_genes$P.Value, decreasing = FALSE), ]
# 选择p值最显著的前10个上调基因
top10_upregulated <- head(upregulated_genes, 10)
# 筛选logFC下调的基因（logFC < -1），并按p值降序排列
downregulated_genes <- data[data$logFC < -0.5, ]
downregulated_genes <- downregulated_genes[order(downregulated_genes$P.Value, decreasing = FALSE), ]
# 选择p值最大的前10个下调基因
top10_downregulated <- head(downregulated_genes, 10)
top10_upregulated
top10_downregulated
# 在绘图数据中注释前10的上下调基因
rownames(data) <- data$Genes
data$sign[which(rownames(data) %in% c(top10_upregulated$Genes, top10_downregulated$Genes))]<-rownames(data)[which(rownames(data)%in% c(top10_upregulated$Genes, top10_downregulated$Genes) )]
#保存差异基因
DEGsList <- rownames(subset(data,regulate == "up-regulated"|regulate == "down-regulated"))
DEG_dat <- subset(data,regulate == "up-regulated"|regulate == "down-regulated")
write.table(DEGsList,"01.DEGsList.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(DEG_dat,"02.DEG.txt",row.names = TRUE,col.names = TRUE,quote = FALSE)
write.xlsx(DEG_dat,"02.DEG.xlsx",rowNames = TRUE)
#可视化
library(ggrepel)
p1 <- ggplot(data = data, aes(x = logFC, y = -log10(P.Value), color = regulate)) +
  geom_point(alpha=0.4, size = 2) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  theme_classic(base_size = 15) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        legend.text = element_text(color="black",size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
  labs(x="logFC",y="-log10(P.Value)")+
  scale_color_manual(name = "", values=c( "#ff4757","#546de5","#d2dae2"), limits = c("up-regulated", "down-regulated", "unchanged"))+
  geom_label_repel(aes(label=sign), fontface="bold", color="grey50",
                   box.padding=unit(0.35, "lines"), point.padding=unit(1, "lines"), 
                   segment.colour = "grey50",size = 5,max.overlaps = 30)
p1
ggsave("03.volcano.pdf",p1,device = cairo_pdf,width = 10, height = 7)
ggsave("03.volcano.png",p1,dpi = 300,width = 10, height = 7)
#03.热图--------------------------------------------------------
up <- subset(data,regulate == "up-regulated")[1:10,]
down <- subset(data,regulate == "down-regulated")[1:10,]
data_deg <- rbind(up,down)
expr <- readRDS("../00_rawdata/GSE62452_expr.Rds")
groupInf <- groupInf[order(groupInf$group), ]
expr <- expr[,groupInf$sample]
pheatmap_data <- expr[data_deg$Genes,]
x<-pheatmap_data
mat <- t(scale(t(x)))
df1 <- as.data.frame(mat)
mat[mat < (-2)] <- (-2)
mat[mat > 2] <- 2
#rt<-sig_diff
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(4, "cm")) %v%
  HeatmapAnnotation(Group = groupInf$group, col = list(Group = c("High Expression" = "#B72230", "Low Expression" = "#104680"))) %v%
  Heatmap(mat, 
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          ###show_colnames = FALSE,
          name = "expression", 
          ###cluster_cols = F,
          cluster_rows = T,
          height = unit(6, "cm"),
          #cluster_columns = FALSE,
          ###cluster_rows = FALSE,
          col = colorRampPalette(c("#87CEFA", "white","#FA8072"))(100))
p_base <- recordPlot()
pdf('04.deg_heatmap.pdf',w=6,h=6,family='Times')
p_base
dev.off()
png('04.deg_heatmap.png',w=6,h=6,family='Times',units = "in", res = 300)
p_base
dev.off()

save.image("../RData/04_DEG.RData")
