rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/19_project_1533")
if (!dir.exists("17_immuneTherapy")) {dir.create("17_immuneTherapy")}
setwd("17_immuneTherapy")
library(readxl)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
ggfonttheme <- theme(text = element_text(family = "Times"))
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
risk_group <- read.table("../09_evaluate/04.risk_table.txt",header = TRUE,row.names = 1,sep = "\t")
risk_group1 <- data.frame(risk_group[,c("risk_list_cut")],row.names = rownames(risk_group))
scores1 <- merge(risk_group1,scores1,by.x = 0,by.y = "Sample")
scores1$riskgroup <- as.factor(as.character(scores1[,2]))
my_comparisons <- list(c("low","high"))
ggplot(scores1,aes(x = ScoreType, y = value,fill = riskgroup))+
  geom_violin(trim = FALSE,position = position_dodge(width = 0.9),scale = "area")+
  geom_boxplot(position = position_dodge(width = 0.9),
               outlier.size = 0,width = 0.4,show.legend = TRUE)+
  scale_fill_manual(values = c("low" = "#5755F5","high" = "#FE655B"),
                    labels = c("High Risk Group","Low Risk Group"))+
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
ggsave("01.estimate.pdf",device = cairo_pdf,width = 8, height = 6)
ggsave("01.estimate.png",dpi = 300,width = 8, height = 6)
#IPS-------------------------------------------------
TCIA_raw <- read_tsv("TCIA-ClinicalData.tsv", locale = locale(encoding = "UTF-8"))
TCIA_raw <- column_to_rownames(TCIA_raw,var = "barcode")
#1:CTLA4 2:PD1
ips <- TCIA_raw[,c(22:25)]
colnames(ips) <- c("CTLA4(−)PD1(-)","CTLA4(−)PD1(+)","CTLA4(+)PD1(-)","CTLA4(+)PD1(+)")
#高低风险组分组
risk_group$sample <- ifelse(  
  str_sub(rownames(risk_group), 14, 15) == '01', 'Tumor', 'Normal')
risk_tumor_group <-subset(risk_group,sample == "Tumor")
rownames(risk_tumor_group) <- substr(rownames(risk_tumor_group),1,12)
rownames(ips) <- gsub("-","_",rownames(ips))
ips1 <- merge(ips,risk_tumor_group,by = 0) #398

ips1_long <- ips1 %>%
  dplyr::select("CTLA4(−)PD1(-)","CTLA4(−)PD1(+)","CTLA4(+)PD1(-)","CTLA4(+)PD1(+)", risk_list_cut) %>%   # 只选需要的列
  pivot_longer(
    cols = c("CTLA4(−)PD1(-)","CTLA4(−)PD1(+)","CTLA4(+)PD1(-)","CTLA4(+)PD1(+)"),
    names_to = "Variable",
    values_to = "Value"
  )  %>% as.data.frame()

my_comparisons <- list(c("high","low"))
colnames(ips1_long) <- c("Group","Type","Value")
ggplot(ips1_long,aes(x = Group, y = Value,fill = Group))+
  geom_violin(trim = FALSE,position = position_dodge(width = 0.9),scale = "area")+
  geom_boxplot(position = position_dodge(width = 0.9),
               outlier.size = 0,width = 0.4,show.legend = TRUE)+
  scale_fill_manual(values = c("low" = "#5755F5","high" = "#FE655B"),
                    labels = c("High Risk Group","Low Risk Group"))+
  facet_wrap(~ Type, scales = "free", ncol = 2) + 
  theme_bw()+
  theme(panel.grid = element_line(colour = "white"),plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        title = element_text(size = 20),
        legend.position = "top",
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 16, face = "bold"),
        strip.background = element_blank())+
  ylab("IPS")+
  xlab("")+
  stat_compare_means(aes(group=Group),                       #按分组进行统计检验
                     method = "wilcox.test",
                     paired = FALSE,                             #非配对t检验
                     symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                        symbols=c("***","**","*","ns")),
                     label = "p.signif",
                     label.x = 1.4,
                     label.y.npc = 0.95,      #添加显著性符号的位置
                     size=4.5)+ggfonttheme
ggsave("02.IPS.pdf",device = cairo_pdf,width = 8, height = 6)
ggsave("02.IPS.png",dpi = 300,width = 8, height = 6)

save.image("../RData/17_ImmuneTherapy.RData")
