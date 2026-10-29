rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project005")
if (!dir.exists("02_TGFB1")) {dir.create("02_TGFB1")}
setwd("02_TGFB1")
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(pROC)
ggfonttheme <- theme(text = element_text(family = "Times"))
genes <- "TGFB1"
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
  ggtitle("TGFB1 Expression in TCGA-UCEC")+
  stat_compare_means(aes(group=Type),                       #按分组进行统计检验
                     method = "wilcox.test",
                     paired = FALSE,                             #非配对t检验
                     symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                        symbols=c("***","**","*","ns")),
                     label = "p.signif",
                     #label.y = location$Proportion+0.02,      #添加显著性符号的位置
                     size=4.5)+ggfonttheme
p1
ggsave("01.gene_expression_violin.pdf",device = cairo_pdf,width = 8, height = 6)
ggsave("01.gene_expression_violin.png",dpi = 300,width = 8, height = 6)
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
                         variables = "TGFB1",minprop = 0.1)
res.cut
# cutpoint statistic
# TGFB1 2.837277  2.995624
dat1 <- dat1 %>% 
  mutate(group2 = ifelse(dat1$TGFB1 > res.cut$cutpoint$cutpoint,"High Expression","Low Expression"))
high_low <- dat1[,c("sample","group2")]
colnames(high_low) <- c("sample","group")
write.table(high_low,"02.high_low_group.txt",sep = "\t",quote = F,row.names = T)
saveRDS(high_low, file = "TGFB1_group.Rds")
fit_km <- survfit(Surv(OS.time,OS) ~ group2,data = dat1)
p2 <- ggsurvplot(fit_km,conf.int = FALSE,pval = TRUE,
                 legend.title = "Group",
                 legend.labs = c("High Expression","Low Expression"),
                 risk.table = TRUE,
                 palette = c('red','blue'))+
  labs(x = "Time in years",title = "TCGA-UCEC")
p2
pdf("03.train_km.pdf",width = 8, height = 7,family = "Times")
print(p2,newpage = FALSE)
dev.off()
png("03.train_km.png",width = 8, height = 7,units = "in", res = 300,family = "Times")
print(p2,newpage = FALSE)
dev.off()
save.image("../RData/02_TGFB1.RData")
