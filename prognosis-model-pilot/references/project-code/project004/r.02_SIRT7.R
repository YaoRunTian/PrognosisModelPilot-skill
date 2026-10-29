rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project004")
if (!dir.exists("02_SIRT7")) {dir.create("02_SIRT7")}
setwd("02_SIRT7")
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(pROC)
ggfonttheme <- theme(text = element_text(family = "Times"))
genes <- "SIRT7"
#TCGA--------------
#KM曲线--------------------------------------------------------
library('survminer')
library("survival")
os <- readRDS("../00_rawdata/GSE62452_os.Rds") %>% as.data.frame()
groupInf <- readRDS("../00_rawdata/GSE62452_group.Rds")
expr <- readRDS("../00_rawdata/GSE62452_expr.Rds")
df2 <- as.data.frame(t(expr[genes,]))
tumor_samples <- groupInf[groupInf[,"Group"] == "Tumor","Sample"]
clinical <- os[tumor_samples,]
dat1 <- merge(clinical,df2,by = 0)
res.cut <- surv_cutpoint(dat1,
                         time = "os_time",
                         event = "os_status",
                         variables = "SIRT7",minprop = 0.1)
res.cut
# cutpoint statistic
# SIRT7  4.77949  1.862144
dat1 <- dat1 %>% 
  mutate(group2 = ifelse(dat1$SIRT7 > res.cut$cutpoint$cutpoint,"High Expression","Low Expression"))
colnames(dat1)[1] <- "sample"
high_low <- dat1[,c("sample","group2")]
colnames(high_low) <- c("sample","group")
write.table(high_low,"01.high_low_group.txt",sep = "\t",quote = F,row.names = T)
saveRDS(high_low, file = "SIRT7_group.Rds")
fit_km <- survfit(Surv(os_time,os_status) ~ group2,data = dat1)
p1 <- ggsurvplot(fit_km,conf.int = FALSE,pval = TRUE,pval.method = TRUE,
                 legend.title = "Group",
                 legend.labs = c("High Expression","Low Expression"),
                 risk.table = TRUE,
                 log.rank.weights = "n",
                 palette = c('red','blue'))+
  labs(x = "Time in years",title = "Training Set(GSE62452)")
p1
pdf("02.GSE62452_km.pdf",width = 8, height = 7,family = "Times")
print(p1,newpage = FALSE)
dev.off()
png("02.GSE62452_km.pNG",width = 8, height = 7,units = "in", res = 300,family = "Times")
print(p1,newpage = FALSE)
dev.off()

dat2 <- subset(dat1,os_time <= 3)
fit_km <- survfit(Surv(os_time,os_status) ~ group2,data = dat2)
p2 <- ggsurvplot(fit_km,conf.int = FALSE,pval = TRUE,pval.method = TRUE,
                 legend.title = "Group",
                 legend.labs = c("High Expression","Low Expression"),
                 risk.table = TRUE,
                 log.rank.weights = "n",
                 palette = c('red','blue'))+
  labs(x = "Time in years",title = "Training Set(GSE62452) years <= 3")
p2
