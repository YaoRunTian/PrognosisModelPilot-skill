rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project001")
if (!dir.exists("17_immuneTherapy")) {dir.create("17_immuneTherapy")}
setwd("17_immuneTherapy")
library(readxl)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
ggfonttheme <- theme(text = element_text(family = "Times"))
#免疫检查点基因----------------------------------------------------
expr <- readRDS("../00_rawdata/fpkm.Rds")
colnames(expr) <- gsub("-","_",colnames(expr))
#高低风险组数据
risk_group <- read.table("../08_evaluate/04.risk_table.txt",header = TRUE,row.names = 1,sep = "\t")
risk_group1 <- data.frame(risk_group$risk_list_cut,row.names = rownames(risk_group))
check_point <- c("PDCD1","TNFRSF8","IDO1","LAG3")
IMMCG <- check_point
IMMCG_list <- intersect(rownames(expr),IMMCG)
IMMCG_dat <- expr[IMMCG_list,]
immune_dat5 <- melt(t(IMMCG_dat))
colnames(immune_dat5) <- c("Sample","Symbol","value")
dim(immune_dat5)
immune_dat5 <- merge(risk_group1,immune_dat5,by.x = 0,by.y = "Sample")
immune_dat5$riskgroup <- as.character(immune_dat5[,2])
my_comparisons <- list(c("low","high"))
ggplot(immune_dat5,aes(x = Symbol, y = log(value+1),fill = riskgroup))+
  geom_violin(trim = FALSE,position = position_dodge(width = 0.9),scale = "area")+
  geom_boxplot(position = position_dodge(width = 0.9),
               outlier.size = 0,width = 0.4,show.legend = TRUE)+
  scale_fill_manual(values = c("low" = "#B8DDBC","high" = "#F0A780"),
                    labels =  c("Low Risk Group","High Risk Group"))+
  theme_bw()+
  theme(panel.grid = element_line(colour = "white"),plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1,size = 18),
        axis.text.y  = element_text(size = 18),
        axis.title  = element_text(size = 18),
        legend.position = "top",
        legend.key.size = unit(35, "pt"),
        legend.text = element_text(size =10),
        title = element_text(size = 15)
  )+
  ylab("Gene Expression")+
  xlab("immune checkpoint genes")+
  labs(fill = "Risk Group")+
  ggtitle("Boxplot of immune checkpoint genes")+
  stat_compare_means(aes(group=riskgroup),                       #按分组进行统计检验
                     method = "wilcox.test",
                     paired = FALSE,                          
                     symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                        symbols=c("***","**","*","ns")),
                     label = "p.signif",
                     #label.y = location$Proportion+0.02,      #添加显著性符号的位置
                     size=4.5,family = "Times")+ggfonttheme
ggsave("01.checkpoint_boxplot.pdf",device = cairo_pdf,width = 8, height = 6)
ggsave("01.checkpoint_boxplot.png",dpi = 300,width = 8, height = 6)
#estimate----------------------------------------------------------------------
library(utils)
library(ggpubr)
library(estimate)
# 设置输入文件的路径
input_file_dir <- "../00_rawdata/01.fpkm.txt"
# 过滤常见基因，并将结果输出到一个新文件中
filterCommonGenes(input.f = input_file_dir,  # 输入文件路径
                  output.f = "os_samples.gct",  # 输出文件路径
                  id = "GeneSymbol")  # 基因标识符类型为基因符号
# [1] "Merged dataset includes 9865  genes (547  mismatched)."
# 计算ESTIMATE评分，并将结果输出到一个新文件中
estimateScore(input.ds = "os_samples.gct",  # 输入文件路径
              output.ds= "os_estimate_score.gct",  # 输出文件路径
              platform="affymetrix")  # 平台类型为Affymetrix
# 读取ESTIMATE评分结果文件，跳过前两行，且设置header为TRUE
scores <- read.table("os_estimate_score.gct", skip = 2, header = TRUE)
# 将第一列设为行名
rownames(scores) = scores[,1]
# 转置dataframe，并保留从第三列开始的所有列
scores <- data.frame(t(scores[,3:ncol(scores)]))
# 将行名保存为一个新列'sample'
scores$sample <- rownames(scores)
# 将数据框的行名设置为NULL
scores <- data.frame(scores, row.names = NULL)
# 将'sample'列中的"."替换为"-"
scores$sample <- gsub("\\.", "-", scores$sample)
#删除TummorPurity数据
library(reshape2)
scores1 <- melt(scores[,-4])
colnames(scores1) <- c("Sample","ScoreType","value")
#高低风险组数据
risk_group <- read.table("../08_evaluate/04.risk_table.txt",header = TRUE,row.names = 1,sep = "\t")
risk_group1 <- data.frame(risk_group[,c("risk_list_cut")],row.names = rownames(risk_group))
scores1 <- merge(risk_group1,scores1,by.x = 0,by.y = "Sample")
scores1$riskgroup <- as.factor(as.character(scores1[,2]))
my_comparisons <- list(c("low","high"))
ggplot(scores1,aes(x = ScoreType, y = value,fill = riskgroup))+
  geom_violin(trim = FALSE,position = position_dodge(width = 0.9),scale = "area")+
  geom_boxplot(position = position_dodge(width = 0.9),
               outlier.size = 0,width = 0.4,show.legend = TRUE)+
  scale_fill_manual(values = c("low" = "#5755F5","high" = "#FE655B"),
                    labels = c("Low Risk Group","High Risk Group"))+
  theme_bw()+
  theme(panel.grid = element_line(colour = "white"),
        plot.title = element_text(hjust = 0.5,size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y  = element_text(size = 18),
        axis.title.y  = element_text(size = 18),
        axis.title.x = element_blank(),
        legend.position = "top",
        legend.key.size = unit(35, "pt"),
        legend.text = element_text(size =10),
        legend.title = element_text(size =12)
  )+
  ylab("Scores")+
  labs(fill = "Risk Group")+
  ggtitle("Boxplot of Estimate Results")+
  stat_compare_means(aes(group=riskgroup),                       #按分组进行统计检验
                     method = "wilcox.test",
                     paired = FALSE,                             #非配对t检验
                     symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                        symbols=c("***","**","*","")),
                     label = "p.signif",
                     #label.y = location$Proportion+0.02,      #添加显著性符号的位置
                     size=4.5)+ggfonttheme
ggsave("02.estimate.pdf",device = cairo_pdf,width = 8, height = 6)
ggsave("02.estimate.png",dpi = 300,width = 8, height = 6)

save.image("../RData/17_ImmuneTherapy.RData")
