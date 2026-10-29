rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project005")
if (! dir.exists("07_evaluate")){dir.create("07_evaluate")}
setwd("07_evaluate")
library(ggplot2)
ggfonttheme <- theme(text = element_text(family = "Times"))
#加载数据
load("../RData/06_model.RData")
#01.KM生存分析-----------------------------------------------------------------
#计算风险评分
riskscore <- function(survival_cancer_df,candidate_genes_for_cox,report){
  library('dplyr')
  risk_score_table <- survival_cancer_df[,candidate_genes_for_cox]
  for(each_sig_gene in colnames(risk_score_table)){
    risk_score_table[,each_sig_gene] <- risk_score_table[,each_sig_gene]*(report[each_sig_gene,1])		
  }
  risk_score_table <- cbind(risk_score_table,'total_risk_score' = rowSums(risk_score_table)) %>%
    cbind(survival_cancer_df[,c('os_time','os_status')])
  return(risk_score_table)
}
risk_score_table_multi_cox2 <- riskscore(df,lasso_gene,lasso_result)
head(risk_score_table_multi_cox2)
#生存分析
library('survminer')
library("survival")
risk_score_table_multi_cox2$os_time <- risk_score_table_multi_cox2$os_time/365
#计算最优cut.off
res.cut <- surv_cutpoint(risk_score_table_multi_cox2,
                         time = "os_time",
                         event = "os_status",
                         variables = "total_risk_score",minprop = 0.1)
res.cut
# cutpoint statistic
# total_risk_score 2.032257  3.455928
#根据最佳阈值划分高低风险组
risk_score_table_multi_cox2 <- risk_score_table_multi_cox2 %>% 
  mutate(risk_list_cut = ifelse(total_risk_score > res.cut$cutpoint[[1]],"high","low"))
fit_km <- survfit(Surv(os_time,os_status) ~ risk_list_cut,
                  data = risk_score_table_multi_cox2)
#KM曲线
p1 <- ggsurvplot(fit_km,conf.int = TRUE,pval = TRUE,
                 legend.title = "Group",
                 legend.labs = c("High risk group","Low risk group"),
                 risk.table = TRUE,
                 palette = c('red','blue'))+
  labs(x = "Time in years",title = "TCGA-UCEC Training Set (n = 374)")
p1
pdf("01.train_km.pdf",width = 8, height = 7,family = "Times")
print(p1,newpage = FALSE)
dev.off()
png("01.train_km.png",width = 8, height = 7,units = "in", res = 300,family = "Times")
print(p1,newpage = FALSE)
dev.off()
#根据风险评分对样本进行排序（patients）列
risk_score_table_multi_cox2 <- risk_score_table_multi_cox2[order(risk_score_table_multi_cox2$total_risk_score),]
risk_score_table_multi_cox2$patients <- c(1:dim(risk_score_table_multi_cox2)[1])
#可视化散点图
library(patchwork)
patients_num <- sum(risk_score_table_multi_cox2$risk_list_cut == "low")
p2 <- ggplot(risk_score_table_multi_cox2,aes(x = patients,y = os_time)) + 
  geom_point(aes(color = factor(os_status)),size =2)+
  geom_vline(xintercept = patients_num, size = 0.8, color = 'black', alpha = 0.8, linetype = 2)+
  labs(x = "",y = "Survival time(years)",title = "TCGA-UCEC Training Set (n = 374)",color = "Status")+
  theme_bw()+
  scale_color_manual(values = c("#9dd0c7","#9180ac"),labels = c("Alive","Dead"))+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size =20),
        title = element_text(size = 20))+ggfonttheme
p3 <- ggplot(risk_score_table_multi_cox2,aes(x = patients,y =total_risk_score)) + 
  geom_point(aes(colour = factor(risk_list_cut)),size =2)+
  geom_hline(yintercept = res.cut$cutpoint[[1]], size = 0.8, color = 'black', alpha = 0.8, linetype = 2)+
  geom_vline(xintercept = patients_num, size = 0.8, color = 'black', alpha = 0.8, linetype = 2)+
  labs(x = "Patients(increasing risk score)",y = "Risk Score",colour = "Group")+
  theme_bw()+
  scale_colour_manual(values = c("red","blue"),labels = c("High risk group","Low risk group"))+
  theme(legend.position = "right",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size =20))+
  annotate('text', x = patients_num-40, y = res.cut$cutpoint[[1]]+0.5, 
           label = paste0('res.cut = ', round(res.cut$cutpoint[[1]],3),"\n","Patients = ",patients_num), 
           color = 'black', size = 6)+ggfonttheme
p2/p3+plot_layout(guides = "collect")&
  theme(legend.position='bottom')
ggsave("02.train_pointplot.pdf",device = cairo_pdf,width = 14, height = 10)
ggsave("02.train_pointplot.png",dpi = 300,width = 14, height = 10)
#热图
library(dplyr)
ann_colors = list(Risk = c(high = "red",low = "blue"))
color <-colorRampPalette(c("navy", "white", "firebrick3"))(100)
fpkm <- readRDS("../00_rawdata/fpkm.Rds")
expr <- fpkm %>% t() %>% as.data.frame()
expr <- expr[,hubgenes]#选取关键基因
m <- expr[,lasso_gene]
annotation_col <- as.data.frame(risk_score_table_multi_cox2$risk_list_cut)
rownames(annotation_col) <- rownames(risk_score_table_multi_cox2)
colnames(annotation_col) <- "Risk"
m <- merge(m,annotation_col,by = 0)
m <- m[order(m$Risk),]
m$sample <-paste("sample",1:nrow(m))
rownames(m) <- m$sample
annotation_col <- m[,c("Risk")] %>% as.data.frame()
rownames(annotation_col) <- rownames(m)
colnames(annotation_col) <- "Risk"
m <- m[,lasso_gene] %>% t() %>% as.data.frame()
library(pheatmap)
p <-pheatmap::pheatmap(m,scale = "row",
                       show_colnames = FALSE,show_rownames = TRUE,
                       cluster_rows = FALSE,cluster_cols = FALSE,
                       annotation_col = annotation_col,color = color,fontsize = 16,fontfamily = "Times",
                       fontsize_row = 16,legend = TRUE,annotation_colors = ann_colors)
library(ggplotify)
library(ggplot2)
pheatmap1 <- as.ggplot(p) + theme(plot.margin = margin(l = 0, r = 20))
pheatmap1
ggsave("03.geneheatmap.pdf",pheatmap1,device = cairo_pdf,width = 8, height = 7)
ggsave("03.geneheatmap.png",pheatmap1,dpi = 300,width = 8, height = 7)
#保存高低风险组分组
write.table(risk_score_table_multi_cox2,"04.risk_table.txt",sep = "\t",quote = FALSE,row.names = TRUE)
#03.ROC曲线---------------------------------------------------------------------
risk_score_table_multi_cox2$riskgroup <- ifelse(risk_score_table_multi_cox2$risk_list_cut == "low",
                                                0,1)
multi_ROC <- function(time_vector,risk_score_table){
  library('survivalROC')
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime = risk_score_table$os_time,
                           status = risk_score_table$os_status,
                           marker = risk_score_table$total_risk_score,
                           predict.time = single_time,method = "KM")
    data.frame('True_positive' = for_ROC$TP,'False_positive' = for_ROC$FP,
               'Cut_values' = for_ROC$cut.values,'Time_point' = rep(single_time,length(for_ROC$TP)),
               'AUC' = rep(for_ROC$AUC,length(for_ROC$TP)))	
  }
  multi_ROC_list <- lapply(time_vector,single_ROC)
  do.call(rbind,multi_ROC_list)
}
#evaluate AUCs between 1,3,5 years.
for_multi_ROC <- multi_ROC(time_vector = c(seq(1,5,2)),risk_score_table = risk_score_table_multi_cox2)
# 按时间分组，提取每个时间点的唯一AUC值
auc_by_time <- for_multi_ROC %>%
  group_by(Time_point) %>%
  dplyr::summarise(AUC = unique(AUC))  # 使用 first(AUC) 或 unique(AUC)
# 输出结果
auc_by_time
#转换因子
for_multi_ROC$Time_point <- as.factor(for_multi_ROC$Time_point)
#visualiaztion
library('plotROC')
pROC <- ggplot(for_multi_ROC,aes(x = False_positive,y = True_positive,label = Cut_values,colour = Time_point))+
  geom_roc(labels = F,stat = 'identity',n.cuts = 0)+
  geom_abline(slope = 1,intercept = 0,color = 'red',linetype = 2)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=16),
        plot.title = element_text(size=16),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15))+
  annotate(geom = "text", x = 0.75, y = 0.15, size = 6, family = "Times",
           label = paste(
             "1-Year AUC = ", sprintf("%.3f", auc_by_time[1,2]), "\n",
             "3-Year AUC = ", sprintf("%.3f", auc_by_time[2,2]), "\n",
             "5-Year AUC = ", sprintf("%.3f", auc_by_time[3,2]), "\n",
             sep = '')
           
  )+
  scale_colour_manual(values = c("#F1C89A","#E79397","#A797DA"),labels = c("1-Year","3-Year","5-Year"))+
  labs(colour = "Time",title = "TCGA-UCEC Training Set (n = 374)")+ggfonttheme
pROC
ggsave("05.trainROC.pdf",device = cairo_pdf,width = 8, height = 7)
ggsave("05.trainROC.png",dpi = 300,width = 8, height = 7)
save.image("../RData/07_evaluate.RData")
