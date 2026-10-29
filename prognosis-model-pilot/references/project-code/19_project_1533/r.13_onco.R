rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/19_project_1533")
if (!dir.exists("13_onco")) {dir.create("13_onco")}
setwd("13_onco")
library(data.table)
library(limma)
library(oncoPredict)##需要先安装依赖包
library(reshape2)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(parallel)
library(edgeR)
set.seed(12345)
ggfonttheme <- theme(text = element_text(family = "Times"))
#读取官方模型参考文件
#网站GDSC2表达数据
GDSC2_Expr=readRDS(file='GDSC2_Expr.rds')
dim(GDSC2_Expr)
#网站GDSC2药物IC50数据
GDSC2_Res=readRDS(file = 'GDSC2_Res.rds')
dim(GDSC2_Res)
##对IC50数据进行指数转换
GDSC2_Res=exp(GDSC2_Res)
#高低风险组数据
risk_group <- read.table("../09_evaluate/04.risk_table.txt",header = TRUE,row.names = 1,sep = "\t")
#读取expr数据
expr <- readRDS("../00_rawdata/count.Rds") %>% t() %>% as.data.frame()
risk_group1 <- data.frame(risk_group$risk_list_cut,row.names = rownames(risk_group))
df1 <- merge(risk_group1,expr,by = 0)
rownames(df1) <- df1[,1]
df1 <- df1[,-1]
group_list <- as.factor(df1[,1])
df1 <- df1[,-1]
df1 <- t(df1) %>% as.data.frame()
#读取差异基因列表
# DEGsList <- read.table("../10_GSEA/05.DEGsList.txt",sep = "\n")[,1]
# df1 <- df1[DEGsList,]
exp <- df1
data=as.matrix(exp)
#主分析，质量控制，时间长短取决于基因数量
calcPhenotype(trainingExprData = GDSC2_Expr,    
              trainingPtype = GDSC2_Res,        
              testExprData = data,              
              batchCorrect = 'eb',  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,#去除表达波动小于0.2的基因    
              minNumSamples = 10,#去除低于10个样本表达的基因    
              printOutput = TRUE,    
              removeLowVaringGenesFrom = 'rawData')
senstivity_result <- read.csv("./calcPhenotype_Output/DrugPredictions.csv")
colnames(senstivity_result)=gsub("(.*)\\_(\\d+)", "\\1", colnames(senstivity_result))
Group <- data.frame(risk_group$risk_list_cut)
rownames(Group) <- rownames(risk_group)
data <-  merge(Group,senstivity_result,by.x = 0,by.y = "X")
rownames(data) <- data[,1]
data <- data[ , -1]
colnames(data)[1] <- "risk"
data$risk=factor(data$risk, levels=c("low", "high"))
type=levels(factor(data[,"risk"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 计算每个药物的 P 值
p_values <- apply(data[, colnames(data) != "risk"], 2, function(drug_data) {
  wilcox.test(drug_data ~ data$risk)$p.value
})
p_values_dat <- as.data.frame(p_values)
#计算每个药物的上下调
data[is.na(data)] = 0
Low <- subset(data,risk == "low")
High <- subset(data,risk == "high")
logFC_df <- data.frame()
for(i in 2:ncol(data)){
  logFC <- log10(mean(Low[,i])/mean(High[,i]))
  drug <- colnames(data)[i]
  drug_logFC <- c(drug,logFC)
  logFC_df <- rbind(logFC_df,drug_logFC)
}
dat1 <- merge(p_values_dat,logFC_df,by.x= 0,by.y = 1)
colnames(dat1) <- c("drug","p.val","logFC")
ded <- subset(dat1,p.val < 0.01)
ded$logFC <- as.numeric(ded$logFC)
ded <- ded[order(ded$p.val),]#112
write.xlsx(ded,"02.DE-drugs.xlsx")
drugs <- ded$drug[1:10]
# 过滤特定药物
data1 <- data %>%
  dplyr::select(risk, all_of(drugs)) %>%
  pivot_longer(cols = -risk, names_to = "Drug", values_to = "value")
type <- levels(factor(data[,"risk"]))
# 设置比较组
comp <- combn(type, 2)
my_comparisons <- list()
for (i in 1:ncol(comp)) {
  my_comparisons[[i]] <- comp[, i]
}
# 绘制箱线图
p1 <- ggplot(data1,aes(x = Drug, y = log2(value),fill = risk))+
  geom_violin(trim = FALSE,position = position_dodge(width = 0.9),scale = "area")+
  geom_boxplot(position = position_dodge(width = 0.9),
               outlier.size = 0,width = 0.4,show.legend = TRUE)+
  scale_fill_manual(values = c("low" = "#B8DDBC","high" = "#F0A780"),
                    labels = c("Low Risk Group","High Risk Group"))+
  theme_bw()+
  theme(panel.grid = element_line(colour = "white"),plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1,size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.position = "top",
        legend.key.size = unit(25, "pt"),
        legend.text = element_text(size =12),
        title = element_text(size = 15)
  )+
  ylab("Log2(IC50)")+
  xlab("")+
  labs(fill = "Risk Group")+
  ggtitle("Boxplot of Drugs")+
  stat_compare_means(aes(group=risk),                       #按分组进行统计检验
                     method = "wilcox.test",
                     paired = FALSE,                             #非配对t检验
                     symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                        symbols=c("***","**","*","")),
                     label = "p.signif",
                     #label.y = location$Proportion+0.02,      #添加显著性符号的位置
                     size=4.5)+ggfonttheme
p1
ggsave("01.boxplot.pdf",device = cairo_pdf,width = 8, height = 6)
ggsave("01.boxplot.png",dpi = 300,width = 8, height = 6)
save.image("../RData/13_onco.RData")
