rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/19_project_1533")
if (!dir.exists("12_immune")) {dir.create("12_immune")}
setwd("12_immune")
library(pheatmap)
library(Rmisc) 
library(corrplot)
library(ggcorrplot)
library(RColorBrewer)
library(grDevices)
library(extrafont) 
library(extrafontdb)
library(Rttf2pt1)
library(devtools)
library(e1071)
library(parallel) 
library(bseqsc)
library(preprocessCore)
library(ggplot2)
library(vioplot)
library(CIBERSORT)
library(Hmisc)
ggfonttheme <- theme(text = element_text(family = "Times"))
# 加载LM22数据集
LM22 <- read.table("LM22_CIBERSORT.txt",header = TRUE,row.names = 1,sep = "\t",check.names = FALSE)
LM22 <- as.matrix(LM22)
#高低风险组数据
risk_group <- read.table("../09_evaluate/04.risk_table.txt",header = TRUE,row.names = 1,sep = "\t")
#读取expr数据
expr <- readRDS("../00_rawdata/fpkm.Rds") %>% 
  t() %>% as.data.frame()
risk_group1 <- data.frame(risk_group[,c("risk_list_cut")],row.names = rownames(risk_group))
df1 <- merge(risk_group1,expr,by = 0)
rownames(df1) <- df1[,1]
df1 <- df1[,-1]
group_list <- as.factor(df1[,1])
df1 <- df1[,-1]
data <- df1
data_matrix <- as.matrix(data, rownames.force=NA)
data_matrix <- matrix(as.numeric(data_matrix), nrow=nrow(data), dimnames=list(rownames(data), colnames(data)))
# 使用CIBERSORT进行细胞比例分析
results <- cibersort(sig_matrix = LM22, mixture_file = t(data_matrix),perm = 100,QN = TRUE)
colnames(results)[23] <- "pvalue"
results <- subset(as.data.frame(results),pvalue < 0.05)#
results <- as.matrix(results[, 1:(ncol(results) - 3)])
results=rbind(id=colnames(results),results)
save(results,file = "01.mmdata.Rdata")
write.table(results, file="02.CIBERSORT_Results.txt", sep="\t", quote=FALSE, col.names=FALSE)
#数据处理--------------------
rt <- read.table("02.CIBERSORT_Results.txt", header=TRUE, sep="\t", check.names=FALSE, row.names=1)
# 创建一个数据框用于存储输出结果
outTab <- data.frame()
#构建分组信息
group1 <- as.data.frame(rownames(rt))
colnames(group1) <- c("Sample")
group1 <- merge(group1,risk_group,by.x = "Sample",by.y = 0)
group1 <- group1[,c("Sample","risk_list_cut")]
immune_dat <- merge(group1,rt,by.x = "Sample",by.y = 0)
rownames(immune_dat) <- immune_dat[,1]
group <- factor(immune_dat[,"risk_list_cut"],labels = c("low","high"),levels = c(0,1))
low_dat <-subset(immune_dat,risk_list_cut == "low")
low_dat <- low_dat[,c(-2,-1)]
high_dat <- subset(immune_dat,risk_list_cut =="high")
high_dat <- high_dat[,c(-2,-1)]
#直方图-------------------------------------------------------------------------
library(tidyr)
library(tibble)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
dat <- rt %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)
dat <- merge(dat,group1,by = "Sample")
colnames(dat)[4] <- "Risk_Group"
dat$Risk_Group <- gsub("high","High Risk",dat$Risk_Group)
dat$Risk_Group <- gsub("low","Low Risk",dat$Risk_Group)
ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22))+
  facet_grid(.~Risk_Group,scales = "free_x")+ggfonttheme
ggsave("04.barplot.pdf",device = cairo_pdf,width = 15, height = 7)
ggsave("04.barplot.png",dpi = 300,width = 15, height = 7)
#绘制小提琴图----------------------------------------------------------------------
#小提琴图可视化-------------------
violin_plot <- function(){
  IDC <- c()
  # 设置坐标轴标签的角度和位置
  par(las=1, mar=c(10,6,3,3))
  # 设置x轴和y轴的范围
  x <- c(1:ncol(rt))
  y <- c(1:ncol(rt))
  plot(x, y,
       xlim=c(0,66), ylim=c(min(rt), max(rt)+0.05),
       main="", xlab="", ylab="Fraction",
       pch=21,
       col="white",
       xaxt="n")
  # 对每一列数据进行处理和绘制
  for(i in 1:ncol(rt)){
    # 如果某列的标准差为0，则将第一个或第二个样本的值设为一个极小值
    if(sd(rt[,i]) == 0){
      rt[1,i] <- 0.00001
    }
    # 分别绘制高、低风险组样本的数据
    vioplot(low_dat[,i], at=3*(i-1), lty=1, add=TRUE, col='blue')
    vioplot(high_dat[,i], at=3*(i-1)+1, lty=1, add=TRUE, col='red')
    # 进行Wilcoxon符号秩检验
    wilcoxTest <- wilcox.test(low_dat[,i], high_dat[,i])
    p <- wilcoxTest$p.value
    # 如果p值小于0.05，则记录结果
    if(p < 0.05){
      cellPvalue <- cbind(Cell=colnames(rt)[i], pvalue=p)
      outTab <- rbind(outTab, cellPvalue)
      diffcell <- colnames(low_dat)[i]
      IDC <- c(IDC,diffcell)
    }
    # 绘制分割线
    mx <- max(c(low_dat[,i],high_dat[,i]))
    lines(c(x=3*(i-1)+0.2, x=3*(i-1)+0.8), c(mx, mx), col='black')
    # 在图中标注p值
    text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p < 0.001, "p<0.001", paste0("p=", sprintf("%.03f", p))), cex=0.8)
  }
  # 添加图例
  legend("topright", 
         c("Low Risk Group", "High Risk Group"),
         pch=15, bty="n", pt.cex=2,
         col=c("blue", "red"))
  # 添加x轴标签
  text(seq(1,66,3), -0.05, xpd=NA, labels=colnames(rt), cex=1, srt=45, pos=2)
  dput(IDC)
  print(paste("差异细胞数为:",length(IDC)))
}
pdf("03.violin.pdf",width = 15, height = 6,family = "Times")
violin_plot()
dev.off()
png("03.violin.png",width = 15, height = 6,units = "in", res = 300,family ="Times")
violin_plot()
dev.off()

#相关性热图--------------------------------------------------------------------
#差异免疫细胞与预后基因及风险评分相关性-----------------------------------------
IDC <- c("Memory B cells", "Plasma cells", "CD8 T cells", "Resting CD4 memory T cells", 
         "T follicular helper cells", "Regulatory T cells (Tregs)", "M0 Macrophages", 
         "Resting mast cells")
print(paste(IDC,collapse =","))
# Memory B cells,Plasma cells,CD8 T cells,Resting CD4 memory T cells,
# T follicular helper cells,Regulatory T cells (Tregs),M0 Macrophages,Resting mast cells
IDC_mat <- rt[,IDC]
# risk_mat <- read.table("../06_evaluate/04.risk_table.txt",sep = "\t",header = TRUE)
# risk_mat <- risk_mat[,1:3]#提取预后基因
# cor_mat <- merge(IDC_mat,risk_mat,by = 0)
cor_mat <- IDC_mat
res2 <- rcorr(as.matrix(cor_mat))
M1 <- res2$r
P1 <- res2$P
addcol <- colorRampPalette(c("blue", "white", "red"))
pdf("05.cor.pdf",width = 13, height = 13,family = "Times")
corrplot(M1,
         method="color",
         addCoef.col="black",
         diag=FALSE,
         tl.col="black",tl.pos = "lt",
         col=addcol(100),
         order = "original")
corrplot(M1, method = "circle", type = "lower",col = addcol(100), 
         tl.col = "n", tl.cex = 0.8, tl.pos = "n",add = TRUE,diag = FALSE,cl.pos = "n",
         p.mat = P1, insig = "label_sig",sig.level = 0.05)
dev.off()
png("05.cor.png",width = 13, height = 13,units = "in", res = 300,family = "Times")
corrplot(M1,
         method="color",
         addCoef.col="black",
         diag=FALSE,
         tl.col="black",tl.pos = "lt",
         col=addcol(100),
         order = "original")
corrplot(M1, method = "circle", type = "lower",col = addcol(100), 
         tl.col = "n", tl.cex = 0.8, tl.pos = "n",add = TRUE,diag = FALSE,cl.pos = "n",
         p.mat = P1, insig = "label_sig",sig.level = 0.05)
dev.off()
#预后基因与免疫细胞相关性-------------------------------------
risk_mat <- read.table("../09_evaluate/04.risk_table.txt",sep = "\t",header = TRUE)
gene <- read.table("../08_model/06.lassogeneslist.txt")[,1]
risk_mat <- risk_mat[,gene]
res3 <- merge(IDC_mat,risk_mat,by = 0)
rownames(res3) <- res3[,1]
res3 <- res3[,-1]
cor_data <- cor(res3,method="spearman")
corp <- ggcorrplot::cor_pmat(res3)
cmt <- cor_data[c(1:length(IDC)),c((length(IDC)+1):(length(IDC)+length(gene)))]
pmt <- corp[c(1:length(IDC)),c((length(IDC)+1):(length(IDC)+length(gene)))]
write.xlsx(cmt,"07.genes_immunecell_cor.xlsx",rowNames = TRUE)
write.xlsx(pmt,"08.gene_immunecell_corp.xlsx",rowNames = TRUE)
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
p1 <- pheatmap::pheatmap(t(cmt),scale = "none",cluster_row = F, cluster_col = F, border=TRUE,
                         breaks=seq(-1, 1, length.out = 10),
                         display_numbers = t(pmt),
                         fontsize_number = 15, number_color = "black",fontsize = 13,
                         cellwidth = 50, cellheight =50,
                         color=colorRampPalette(c("#3777ac", "white","#FC8D62"))(10),
                         angle_col = 45,fontfamily = "Times")
ggsave("06.genes_immunecell_cor.pdf",p1,device = cairo_pdf,width = 10, height = 8)
ggsave("06.genes_immunecell_cor.png",p1,dpi = 300,width = 10, height = 8)
#预后基因相关性-------------------------------------
risk_mat <- read.table("../09_evaluate/04.risk_table.txt",sep = "\t",header = TRUE)
gene <- read.table("../08_model/06.lassogeneslist.txt")[,1]
risk_mat <- risk_mat[,gene]
res4 <- rcorr(as.matrix(risk_mat))
M1 <- res4$r
P1 <- res4$P
addcol <- colorRampPalette(c("blue", "white", "red"))
pdf("07.cor_gene.pdf",width = 13, height = 13,family = "Times")
corrplot(M1,
         method="color",
         addCoef.col="black",
         diag=FALSE,
         tl.col="black",tl.pos = "lt",
         col=addcol(100),
         order = "original")
corrplot(M1, method = "circle", type = "lower",col = addcol(100), 
         tl.col = "n", tl.cex = 0.8, tl.pos = "n",add = TRUE,diag = FALSE,cl.pos = "n",
         p.mat = P1, insig = "label_sig",sig.level = 0.05)
dev.off()
png("07.cor_gene.png",width = 13, height = 13,units = "in", res = 300,family = "Times")
corrplot(M1,
         method="color",
         addCoef.col="black",
         diag=FALSE,
         tl.col="black",tl.pos = "lt",
         col=addcol(100),
         order = "original")
corrplot(M1, method = "circle", type = "lower",col = addcol(100), 
         tl.col = "n", tl.cex = 0.8, tl.pos = "n",add = TRUE,diag = FALSE,cl.pos = "n",
         p.mat = P1, insig = "label_sig",sig.level = 0.05)
dev.off()
save.image("../RData/07_immune.RData")
