rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/24_job_0780")
if (!dir.exists("15_immuneTherapy")) {dir.create("15_immuneTherapy")}
setwd("15_immuneTherapy")
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
risk_group <- read.table("../07_evaluate/04.risk_table.txt",header = TRUE,row.names = 1,sep = "\t")
risk_group1 <- data.frame(risk_group$risk_list_cut,row.names = rownames(risk_group))
check_point <- read.table("checkpoint.txt")[,1]
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
save.image("../RData/15_immuneTherapy.RData")
