rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project002")
if (!dir.exists("14_mutation")) {dir.create("14_mutation")}
setwd("14_mutation")
#突变数据
library(TCGAmutations)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)
ggfonttheme <- theme(text = element_text(family = "Times"))
maf <- tcga_load(c("BLCA"))
#读取高低表达组
risk_table <-  read.table("../09_evaluate/04.risk_table.txt",sep = "\t",header = TRUE)#读取高低风险组信息
samples <- rownames(risk_table)
colnames(risk_table)[6] <- "RiskGroup"
# 提取了每个样本的名称(前者是1：16位的/后者是全称)
samp = data.frame(ID = str_sub(unique(maf@data$Tumor_Sample_Barcode),1,16),
                  long = unique(maf@data$Tumor_Sample_Barcode))
samp$ID <- gsub("-","_",samp$ID)
# merge() 函数会自动根据指定的列（如 ID 列）将两个数据框的数据进行配对合并
samp= merge(risk_table,samp,by.x = 0,by.y = "ID")
#396个样本
maf@clinical.data = left_join(maf@clinical.data,samp,by = c("Tumor_Sample_Barcode"="long"))
maf = subsetMaf(maf,tsb = samp$long)
# 定义突变参数
cla1 = c('Frame_Shift_Del','Missense_Mutation', 
         'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
         'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')
# 上面的信息藏在variant.classification.summary中
cla2 = colnames(maf@variant.classification.summary)[-c(1,ncol(maf@variant.classification.summary))]
# 使用intersect函数找到cla1和cla2的交集,即共同的突变分类；然后使用setdiff函数找到cla2中不在cla1中的突变分类,将它们合并以形成最终的突变分类列表cla
cla = c(intersect(cla1,cla2),setdiff(cla2,cla1))
# 突变分类分配颜色
col = RColorBrewer::brewer.pal(n = length(cla), name = 'Paired')
names(col) = cla
#注释行信息颜色
riskcolor = c("red","blue")
names(riskcolor) = na.omit(unique(maf@clinical.data$RiskGroup))
fabcolors = list(RiskGroup = riskcolor)
oncoplot(maf = maf, colors = col, top = 20,
         clinicalFeatures = c("RiskGroup"), 
         sortByAnnotation = TRUE, 
         annotationColor = fabcolors,
         fontSize =  0.8,legendFontSize = 2,annotationFontSize = 2,
         titleFontSize = 1.6,SampleNamefontSize = 1.5)
pdf("01.oncoplot.pdf",width = 13, height = 10,family = "Times")
oncoplot(maf = maf, colors = col, top = 20,
         clinicalFeatures = c("RiskGroup"), 
         sortByAnnotation = TRUE, 
         annotationColor = fabcolors,
         fontSize =  0.8,legendFontSize = 2,annotationFontSize = 2,
         titleFontSize = 1.6,SampleNamefontSize = 1.5)
dev.off()
png("01.oncoplot.png",width = 13, height = 10,units = "in", res = 300,family = "Times")
oncoplot(maf = maf, colors = col, top = 20,
         clinicalFeatures = c("RiskGroup"), 
         sortByAnnotation = TRUE, 
         annotationColor = fabcolors,
         fontSize =  0.8,legendFontSize = 2,annotationFontSize = 2,
         titleFontSize = 1.6,SampleNamefontSize = 1.5)
dev.off()
#计算TMB值--------------
tmb_mat <- as.data.frame(tmb(maf))
tmb_mat[,1] <- substring(tmb_mat[,1],1,16)
tmb_mat[,1] <- gsub("-","_",tmb_mat[,1])
tmb_risk_mat <- merge(tmb_mat,risk_table,by.x= 1,by.y = 0)
#箱型图
p2 <- ggplot(tmb_risk_mat,aes(x = RiskGroup, y = total_perMB_log,fill = RiskGroup))+
  geom_violin(trim = FALSE,position = position_dodge(width = 0.9),scale = "area")+
  geom_boxplot(position = position_dodge(width = 0.9),
               outlier.size = 0,width = 0.4,show.legend = TRUE)+
  scale_fill_manual(values = c("low" = "#5755F5","high" = "#FE655B"))+
  theme_bw()+
  theme(panel.grid = element_line(colour = "white"),
        plot.title = element_text(hjust = 0.5),
        axis.text  = element_text(size = 15),
        axis.title  = element_text(size = 15),
        legend.position = "top",
        legend.text = element_text(size= 13),
        legend.title = element_text(size = 13)
  )+
  labs(fill = "Risk Group",title = "Boxplot of TMB",y = "TMB",x= "")+
  stat_compare_means(aes(group=RiskGroup),                       #按分组进行统计检验
                     method = "t.test",
                     paired = FALSE,                             #非配对t检验
                     label = "p.signif",
                     label.x = 1.45, #添加显著性符号的位置
                     size=5)+ggfonttheme
p2
ggsave("02.TMB_boxplot.pdf",p2,device = cairo_pdf,width = 6, height = 8)
ggsave("02.TMB_boxplot.png",p2,dpi = 300,width = 6, height = 8)
#生存曲线
#读取生存信息
surv_table <-  read.table("../09_evaluate/04.risk_table.txt",sep = "\t",header = TRUE)
surv_df <- merge(tmb_mat,surv_table,by.x =1,by.y = 0)
res.cut <- surv_cutpoint(surv_df,
                         time = "os_time",
                         event = "os_status",
                         variables = "total_perMB_log")
res.cut
# cutpoint statistic
# total_perMB_log 0.5998831   4.39921
surv_df <- surv_df %>%
  mutate(TMB_Group = ifelse(total_perMB_log > res.cut$cutpoint[[1]],"H-TMB","L-TMB"))
surv_df$risk_list_cut <- gsub("high","high risk",surv_df$risk_list_cut)
surv_df$risk_list_cut <- gsub("low","low risk",surv_df$risk_list_cut)
surv_df$Group <- paste(surv_df$TMB_Group,surv_df$risk_list_cut,sep = "+")
fit_km <- survfit(Surv(os_time,os_status) ~ TMB_Group,
                  data = surv_df)
#TMB-KM曲线
p3 <- ggsurvplot(fit_km,conf.int = TRUE,pval = TRUE,legend.title = "TMB Group",
                 legend.labs = c("H-TMB","L-TMB"),
                 risk.table = TRUE,
                 palette = c('red','blue'))+
  labs(x = "Time in years")
p3
pdf("03.TMB_km.pdf",width = 8, height = 7,family = "Times")
print(p3,newpage = FALSE)
dev.off()
png("03.TMB_km.png",width = 8, height = 7,units = "in", res = 300,family = "Times")
print(p3,newpage = FALSE)
dev.off()
#综合TMB曲线
fit_km <- survfit(Surv(os_time,os_status) ~ Group,
                  data = surv_df)
p4 <- ggsurvplot(fit_km,conf.int = FALSE,pval = TRUE,legend.title = "Group",
                 risk.table = TRUE,
                 palette = c('red','blue','lightgreen','purple'))+
  labs(x = "Time in years")+guides(colour = guide_legend(nrow = 2))
p4
p4$plot <- p4$plot + 
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "top",          # 可选：把图例放底部更易排版
        legend.direction = "horizontal",     # 水平排列
        legend.box = "horizontal") +         # 图例框水平排列
  guides(colour = guide_legend(nrow = 2))
p4
pdf("04.co_km.pdf",width = 8, height = 7,family = "Times")
print(p4,newpage = FALSE)
dev.off()
png("04.co_km.png",width = 8, height = 7,units = "in", res = 300,family = "Times")
print(p4,newpage = FALSE)
dev.off()

save.image("../RData/14_mutation.RData")
