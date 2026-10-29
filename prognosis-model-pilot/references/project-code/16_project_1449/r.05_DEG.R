rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/16_project_1449")
if (!dir.exists("05_DEG")) {dir.create("05_DEG")}
setwd("05_DEG")
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(ggrepel)
ggfonttheme <- theme(text = element_text(family = "Times"))
#load the data
groupInf <- readRDS("../00_rawdata/group.Rds")
groupInf$sample <-gsub("-","_",groupInf$sample)
rownames(groupInf) <- NULL
count <- readRDS("../00_rawdata/count.Rds")
#01.差异分析------------------------------------------------------------------
expr <- round(2^count-1,0)#还原成count数据，并取整
expr <- expr[rowMeans(expr) > 1, ]  # 保留行均值大于1的行（排除不重要的极值基因）
write.csv(expr,'01.cleaned_expr.csv',row.names = TRUE)#重新写入整理后的矩阵文件
#分组信息
groupInf$group <- factor(groupInf$group)
Data <- groupInf
group <- groupInf$group
#构建DEseq2对象(dds)
dds <- DESeqDataSetFromMatrix(countData = expr,
                              colData = Data,
                              design = ~ group)
dds2 <- DESeq(dds)
res <- results(dds2, contrast=c("group", "Tumor", "Normal"))
res <- res[order(res$pvalue),]#按P值从小到大排序
summary(res)#检查
my_result <- as.data.frame(res)#转成容易查看的数据框
my_result <- na.omit(my_result)#删除倍数为0的值
#保存差异分析的结果
my_result$Gene_symbol<-rownames(my_result)
my_result <- my_result %>% dplyr::select('Gene_symbol',colnames(my_result)[1:dim(my_result)[2]-1],everything())
rownames(my_result) <- NULL
write.csv(my_result,file="02.my_result_deseq2.csv")#写入
# my_result <- read.csv("02.my_result_deseq2.csv",row.names = 1)#读入
#筛选差异表达基因
DEG<-my_result
DEG$regulate <- ifelse(DEG$padj > 0.05, "unchanged",
                       ifelse(DEG$log2FoldChange > 1, "up-regulated",
                              ifelse(DEG$log2FoldChange < -1, "down-regulated", "unchanged")))
table(DEG$regulate)
# down-regulated      unchanged   up-regulated 
# 2189          12009           3963 
#火山图------------------------------
data<-DEG
dput(colnames(data))
data$sign <-NA
data <- data[order(data$padj, data$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
# 筛选logFC上调的基因（logFC > 1），并按p值升序排列
upregulated_genes <- data[data$log2FoldChange > 1, ]
upregulated_genes <- upregulated_genes[order(upregulated_genes$padj, decreasing = FALSE), ]
# 选择p值最显著的前10个上调基因
top10_upregulated <- head(upregulated_genes, 10)
# 筛选logFC下调的基因（logFC < -1），并按p值升序排列
downregulated_genes <- data[data$log2FoldChange < -1, ]
downregulated_genes <- downregulated_genes[order(downregulated_genes$padj, decreasing = FALSE), ]
# 选择p值最大的前10个下调基因
top10_downregulated <- head(downregulated_genes, 10)
top10_upregulated
top10_downregulated
# 在绘图数据中注释前10的上下调基因
rownames(data) <- data$Gene_symbol
data$sign[which(rownames(data) %in% c(top10_upregulated$Gene_symbol, top10_downregulated$Gene_symbol))]<-rownames(data)[which(rownames(data)%in% c(top10_upregulated$Gene_symbol, top10_downregulated$Gene_symbol) )]
#保存差异基因
DEGsList <- rownames(subset(data,regulate == "up-regulated"|regulate == "down-regulated"))
DEG_dat <- subset(data,regulate == "up-regulated"|regulate == "down-regulated")
write.table(DEGsList,"01.DEGsList.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(DEG_dat,"02.DEG.txt",row.names = TRUE,col.names = TRUE,quote = FALSE)
write.xlsx(DEG_dat,"02.DEG.xlsx",rowNames = TRUE)
#可视化
#padj有0值，故加一个极小值
p1 <- ggplot(data = data, aes(x = log2FoldChange, y = -log10(padj+1e-273), color = regulate)) +
  geom_point(alpha=0.4, size = 2) +
  geom_hline(yintercept = -log10(0.05+1e-273),lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  theme_classic(base_size = 15) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        legend.text = element_text(color="black",size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),legend.position = "top") +
  labs(x="log2FC",y="-log10(P.adj)")+
  scale_color_manual(name = "", values=c( "#ff4757","#546de5","grey50"), limits = c("up-regulated", "down-regulated", "unchanged"))+
  geom_label_repel(aes(label=sign), fontface="bold", color="grey50",
                   box.padding=unit(0.35, "lines"), point.padding=unit(0.5, "lines"), 
                   segment.colour = "grey50",size = 5,max.overlaps = 20,family = "Times")+ggfonttheme
p1
ggsave("03.volcano.pdf",p1,device = cairo_pdf,width = 12, height = 8)
ggsave("03.volcano.png",p1,dpi = 300,width = 12, height = 8)
#03.热图--------------------------------------------------------
up <- subset(data,regulate == "up-regulated")[1:10,]
down <- subset(data,regulate == "down-regulated")[1:10,]
data_deg <- rbind(up,down)
expr <- readRDS("../00_rawdata/fpkm.Rds")
groupInf <- groupInf[order(groupInf$group), ]
expr <- expr[,groupInf$sample]
pheatmap_data <- expr[data_deg$Gene_symbol,]
x<-pheatmap_data
mat <- t(scale(t(x)))
df1 <- as.data.frame(mat)
mat[mat < (-2)] <- (-2)
mat[mat > 2] <- 2
#rt<-sig_diff
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(4, "cm")) %v%
  HeatmapAnnotation(Group = groupInf$group, col = list(Group = c("Tumor" = "#B72230", "Normal" = "#104680"))) %v%
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

save.image("../RData/05_DEG.RData")
