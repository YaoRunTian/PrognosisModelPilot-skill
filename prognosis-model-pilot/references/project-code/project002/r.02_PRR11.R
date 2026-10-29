rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project002")
if (!dir.exists("02_PRR11")) {dir.create("02_PRR11")}
setwd("02_PRR11")
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(pROC)
ggfonttheme <- theme(text = element_text(family = "Times"))
genes <- "PRR11"
#TCGA--------------------------------------------------
#expression
groupInf <- readRDS("../00_rawdata/group.Rds")
groupInf$sample <-gsub("-","_",groupInf$sample)
expr <- readRDS("../00_rawdata/fpkm.Rds")
df1 <- as.data.frame(t(expr[genes,]))
df1$Sample <- rownames(df1)
df1 <- melt(df1)
colnames(df1) <- c("Sample","Symbol","value")
type <- rep(as.character(groupInf$group),length(genes))
df1$Type <- type
my_comparisons <- list(c("Normal","Tumor"))
#绘制组间差异小提琴图
p1 <- ggplot(df1,aes(x = Symbol, y = value,fill = Type))+
  geom_violin(trim = FALSE,position = position_dodge(width = 0.9),scale = "area")+
  geom_boxplot(position = position_dodge(width = 0.9),
               outlier.size = 0,width = 0.2,show.legend = TRUE)+
  scale_fill_manual(values = c("Normal" = "#6A97FA","Tumor" = "#FA736A"))+
  theme_bw()+
  theme(panel.grid = element_line(colour = "white"),plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        title = element_text(size = 20),
        legend.position = "top",
        legend.text = element_text(size = 15))+
  ylab("Gene Expression")+
  xlab("Genes")+
  ggtitle("PRR11 Expression in Training Set(TCGA-BLCA)")+
  stat_compare_means(aes(group=Type),                       #按分组进行统计检验
                     method = "wilcox.test",
                     paired = FALSE,                             #非配对t检验
                     symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                        symbols=c("***","**","*","ns")),
                     label = "p.signif",
                     #label.y = location$Proportion+0.02,      #添加显著性符号的位置
                     size=4.5)+ggfonttheme
p1
ggsave("01.Train_PRR11.pdf",device = cairo_pdf,width = 8, height = 6)
ggsave("01.Train_PRR11.png",dpi = 300,width = 8, height = 6)
#GSE----------------------------------
expr <- readRDS("../01_GSE13507/expr.Rds")
df2 <- as.data.frame(t(expr[genes,]))
df2$Sample <- rownames(df2)
df2 <- melt(df2)
colnames(df2) <- c("Sample","Symbol","value")
groupInf <- readRDS("../01_GSE13507/group.Rds")
type <- rep(as.character(groupInf$Group),length(genes))
df2$Type <- type
my_comparisons <- list(c("Normal","Tumor"))
#绘制组间差异小提琴图
p2 <- ggplot(df2,aes(x = Symbol, y = value,fill = Type))+
  geom_violin(trim = FALSE,position = position_dodge(width = 0.9),scale = "area")+
  geom_boxplot(position = position_dodge(width = 0.9),
               outlier.size = 0,width = 0.2,show.legend = TRUE)+
  scale_fill_manual(values = c("Normal" = "#6A97FA","Tumor" = "#FA736A"))+
  theme_bw()+
  theme(panel.grid = element_line(colour = "white"),plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        title = element_text(size = 20),
        legend.position = "top",
        legend.text = element_text(size = 15))+
  ylab("Gene Expression")+
  xlab("Genes")+
  ggtitle("PRR11 Expression in Testing Set(GSE13507)")+
  stat_compare_means(aes(group=Type),                       #按分组进行统计检验
                     method = "wilcox.test",
                     paired = FALSE,                             #非配对t检验
                     symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                        symbols=c("***","**","*","ns")),
                     label = "p.signif",
                     #label.y = location$Proportion+0.02,      #添加显著性符号的位置
                     size=4.5)+ggfonttheme
p2
ggsave("02.Test_PRR11.pdf",device = cairo_pdf,width = 8, height = 6)
ggsave("02.Test_PRR11.png",dpi = 300,width = 8, height = 6)
#ROC--------------------------------------------------------------
roc1 <- roc(response = df1$Type, predictor = df1$value, levels = c("Normal", "Tumor"))
auc_ci <- ci.auc(roc1)
roc2 <- roc(response = df2$Type,predictor = df2$value,levels = c("Normal", "Tumor"))
auc_ci_gse <- ci.auc(roc2)

roc_plot_data <- data.frame(
  Specificity = c(1 - roc1$specificities, 1 - roc2$specificities),
  Sensitivity = c(roc1$sensitivities, roc2$sensitivities),
  Dataset = c(rep("TCGA (Training)", length(roc1$sensitivities)),
              rep("GSE39582 (Testing)", length(roc2$sensitivities)))
)
p3 <- ggplot(roc_plot_data, aes(x = Specificity, y = Sensitivity, color = Dataset)) +
  geom_path(size =1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  labs(title = "ROC Curves for PRR11 in BLCA",
       x = "1 - Specificity (False Positive Rate)",
       y = "Sensitivity (True Positive Rate)",
       color = "Dataset") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_color_manual(values = c('#9dd0c7','#9180ac'))+
  coord_equal()+
  # 在图上添加AUC值注解
  annotate("text", x = 0.6, y = 0.3, 
           label = paste("Training Set AUC =", sprintf("%.3f", roc1$auc)), 
           color = "#9180ac", hjust = 0) + # 颜色需与图例对应
  annotate("text", x = 0.6, y = 0.2, 
           label = paste("Testing Set AUC =", sprintf("%.3f", roc2$auc)), 
           color = "#9dd0c7", hjust = 0)+ggfonttheme
p3
ggsave("03.ROC.pdf",p3,device = cairo_pdf,width = 8, height = 6)
ggsave("03.ROC.png",p3,dpi = 300,width = 8, height = 6)
#KM曲线--------------------------------------------------------
#train
library('survminer')
library("survival")
os <- readRDS("../00_rawdata/os.Rds") %>% as.data.frame()
expr <- readRDS("../00_rawdata/fpkm.Rds")
df2 <- as.data.frame(t(expr[genes,]))
clinical <- subset(os,group == "Tumor")
dat1 <- merge(clinical,df2,by.x = 1,by.y = 0)
dat1$OS.time <- dat1$OS.time/365
res.cut <- surv_cutpoint(dat1,
                         time = "OS.time",
                         event = "OS",
                         variables = "PRR11",minprop = 0.1)
res.cut
# cutpoint statistic
# PRR11 2.907967   3.74474
dat1 <- dat1 %>% 
  mutate(group2 = ifelse(dat1$PRR11 > res.cut$cutpoint$cutpoint,"High Expression","Low Expression"))
high_low <- dat1[,c("sample","group2")]
colnames(high_low) <- c("sample","group")
write.table(high_low,"04.high_low_group.txt",sep = "\t",quote = F,row.names = T)
saveRDS(high_low, file = "PRR11_group.Rds")
fit_km <- survfit(Surv(OS.time,OS) ~ group2,data = dat1)
p1 <- ggsurvplot(fit_km,conf.int = TRUE,pval = TRUE,
                 legend.title = "Group",
                 legend.labs = c("High Expression","Low Expression"),
                 risk.table = TRUE,
                 palette = c('red','blue'))+
  labs(x = "Time in years",title = "Training Set(TCGA-BLCA)")
p1
pdf("05.train_km.pdf",width = 8, height = 7,family = "Times")
print(p1,newpage = FALSE)
dev.off()
png("05.train_km.png",width = 8, height = 7,units = "in", res = 300,family = "Times")
print(p1,newpage = FALSE)
dev.off()
#test
os <- readRDS("../01_GSE13507/os.Rds")
expr <- readRDS("../01_GSE13507/expr.Rds")
df2 <- as.data.frame(t(expr[genes,]))
dat1 <- merge(os,df2,by.x = 0,by.y = 0)
res.cut <- surv_cutpoint(dat1,
                         time = "os_time",
                         event = "os_status",
                         variables = "PRR11",minprop = 0.1)
res.cut
# cutpoint statistic
# PRR11     8.35  3.601825
dat1 <- dat1 %>% 
  mutate(group2 = ifelse(dat1$PRR11 > res.cut$cutpoint$cutpoint,"High Expression","Low Expression"))
table(dat1$group2)
colnames(dat1)[1] <- "sample"
high_low <- dat1[,c("sample","group2")]
colnames(high_low) <- c("sample","group")
write.table(high_low,"06.test_high_low_group.txt",sep = "\t",quote = F,row.names = T)
saveRDS(high_low, file = "test_PRR11_group.Rds")
fit_km <- survfit(Surv(os_time,os_status) ~ group2,data = dat1)
p1 <- ggsurvplot(fit_km,conf.int = TRUE,pval = TRUE,
                 legend.title = "Group",
                 legend.labs = c("High Expression","Low Expression"),
                 risk.table = TRUE,
                 palette = c('red','blue'))+
  labs(x = "Time in years",title = "Testing Set(GSE13507)")
p1
pdf("07.test_km.pdf",width = 8, height = 7,family = "Times")
print(p1,newpage = FALSE)
dev.off()
png("07.test_km.png",width = 8, height = 7,units = "in", res = 300,family = "Times")
print(p1,newpage = FALSE)
dev.off()
save.image("../RData/02_PRR11.RData")
