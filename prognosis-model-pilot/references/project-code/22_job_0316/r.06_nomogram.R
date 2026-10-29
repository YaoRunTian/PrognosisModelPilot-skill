rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/22_job_0316")
if (!dir.exists("06_nomogram")) {dir.create("06_nomogram")}
setwd("06_nomogram")
library(dplyr)
library(readr)
library("autoReg")
library(survminer)
library(survival)
library(forestplot)
ggfonttheme <- theme(text = element_text(family = "Times"))
clinical <- readRDS("../00_rawdata/clinical.Rds")
risk_table <-  read.table("../04_evaluate/04.risk_table.txt",sep = "\t",header = TRUE)#读取高低风险组信息
rownames(clinical) <- clinical$sample
clinical <- clinical[,-1]
clinical <- clinical[rownames(risk_table),]
clinical <- na.omit(clinical)#214个样本
#整理风险评分和风险分组
risk_table1 <- risk_table[,c("total_risk_score","risk_list_cut")]
colnames(risk_table1) <- c("Risk_Score","Risk_Group")
write.table(risk_table1,"02.risk_table.txt",sep = "\t",quote = F,row.names = T)
#整合,删除缺失的临床数据，共计398个样本
df <- merge(clinical,risk_table1,by = 0)
rownames(df) <- df[,1]
df <- df[,-1]
write.table(df,"01.clinical.txt",sep = "\t",quote = F,row.names = T)
df$Age <- ifelse(df$Age > 60, "Age > 60","Age ≤ 60")
library(stringr)
df$Stage <- gsub("A|B|C","",df$Stage)
df$T_Stage <- gsub("a|b","",df$T_Stage)
df$M_Stage <- gsub("a|b","",df$M_Stage)
df$N_Stage <- gsub("a|b|c","",df$N_Stage)

#删除MX、NX、TX
#M_stage只有M0,MX不考虑。
# table(df$M)
# grep("MX",df$M)
# df <- df[-grep("MX",df$M),]
df <- df[,-7]

#N_stage只有N1只有3个，不考虑。
table(df$N_Stage)
# grep("NX",df$N_Stage)
# df <- df[-grep("NX",df$N_Stage),]
df <- df[,-7]
table(df$T_Stage)
grep("Tis",df$T_Stage)
df <- df[-grep("Tis",df$T_Stage),]
# table(df$Stage)
# df <- df[order(df$Risk_Group),]
write.table(df,"03.clinical_cleaned.txt",sep = "\t",quote = F,row.names = T)

df <- read.table("03.clinical_cleaned.txt",sep = "\t",header = TRUE,row.names = 1)
dat1 <- df#临床信息
#211个样本
#01.单因素cox回归分析-----------------------------------------------------------
df <- dat1
#T_stage分为T1-2，T3-4，Stage分为Stage1-2，Stage3-4
# df <- df %>%
#   mutate(
#     T_Stage = case_when(
#       T_Stage %in% c("T1", "T2") ~ "T1-2",
#       T_Stage %in% c("T3", "T4") ~ "T3-4",
#       TRUE ~ NA_character_  # 处理意外值
#     ),
#     Stage = case_when(
#       Stage %in% c("Stage I", "Stage II") ~ "Stage1-2",
#       Stage %in% c("Stage III", "Stage IV") ~ "Stage3-4",
#       TRUE ~ NA_character_
#     )
#   )

df$Stage <- as.factor(df$Stage) %>% as.numeric()
df$T_Stage <- as.factor(df$T_Stage) %>% as.numeric()
# df$N_Stage <- as.factor(df$N_Stage) %>% as.numeric()
# df$M_Stage <- as.factor(df$M_Stage) %>% as.numeric()
# df$Stage <- as.factor(df$Stage)
# df$T_Stage <- as.factor(df$T_Stage)
# df$N <- as.factor(df$N)
# df$M <- as.factor(df$M)
#删除Risk_Group分组
# df <- subset(df, select = -c(Risk_Group))
#新建空白数据框
uniresult <- data.frame()
#unicox
for(i in colnames(df[,6:ncol(df)])){
  print(i)
  unicox <- coxph(Surv(time = os_time, event = os_status) ~ df[,i], data = df)  
  unisum<- summary(unicox) 
  pvalue <- round(unisum$coefficients[,5],3) 
  print(paste("pvalue = ",pvalue))
  uniresult <- rbind(uniresult,
                     cbind(var=i,
                           HR=unisum$coefficients[,2],
                           L95CI=unisum$conf.int[,3],
                           H95CI=unisum$conf.int[,4],
                           pvalue=unisum$coefficients[,5]
                     ))
}
uniresult
write.csv(uniresult,file = "01.single_cox.csv",row.names = F)
var1 <- uniresult[,"var"]#提取基因名
#绘制PH假定检验图
phresult <- data.frame()
if (! dir.exists("PH_test")){dir.create("PH_test")}
for(i in var1){
  myformula <- as.formula(paste0("Surv(time = os_time, event = os_status)~",i))
  unicox <- coxph(myformula,data = df)
  test.ph <- cox.zph(unicox)
  ph_p <- ggcoxzph(test.ph)
  pdf_filename <- paste0("./PH_test/",i,"_ph.pdf")
  png_filename <- paste0("./PH_test/",i,"_ph.png")
  print(png_filename)
  pdf(pdf_filename,width = 5, height = 5,family = "Times")
  print(ph_p,newpage = FALSE)
  dev.off()
  png(png_filename,width = 5, height = 5,units = "in", res = 300,family = "Times")
  print(ph_p,newpage = FALSE)
  dev.off()
  phresult <- rbind(phresult,
                    cbind(var1 = i,
                          pval = test.ph$table[1,3]))
}
phresult
write.csv(phresult,file = "02.PHresult.csv",row.names = FALSE)
#FORESTPLOT#单因素cox可视化森林图
hr <- round(as.numeric(uniresult$"HR"),3)#从tducs数据中提取HR值，保留2位小数
hrLow  <- round(as.numeric(uniresult$"L95CI"),3)#提取HR值95%置信区间的低值
hrHigh <- round(as.numeric(uniresult$"H95CI"),3)#提取HR值95%置信区间的高值
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")#将HR及其置信区间数据组合成一个因子
PH_pvalue <- round(as.numeric(phresult$pval),3)
pValue <- ifelse(as.numeric(uniresult$"pvalue")<0.001, "<0.001", round(as.numeric(uniresult$"pvalue"),3))#提取p值，且当p<0.001时显示为<0.001，而不是确切数值
unicoxdat <- cbind(var1,pValue,PH_pvalue,Hazard.ratio,hr,hrLow,hrHigh)
unicoxdat <- rbind(colnames(unicoxdat),unicoxdat) %>% as.data.frame()
unicoxdat[1,c(1,5:7)] <- NA
unicoxdat[,5:7] <- as.numeric(unlist(unicoxdat[,5:7]))
library(forestplot)
fp <- forestplot(labeltext=as.matrix(unicoxdat[,c(1:4)]),
                 mean = unicoxdat$hr,
                 lower=unicoxdat$hrLow,
                 upper=unicoxdat$hrHigh,
                 is.summary = c(TRUE,rep(FALSE,7)),
                 zero =1,
                 boxsize = 0.15,
                 lineheight = unit(6,'mm'),
                 graphwidth = unit(180,'mm'),
                 colgap = unit(6,'mm'),
                 lwd.zero = 2,
                 lwd.ci = 2,
                 col=fpColors(box = '#ef233c',summary = "#9e2a2b",lines='blue1',zero = 'grey'),
                 xlab = 'Hazard Ratio',
                 lwd.xaxis = 2,
                 lty.ci = "solid",
                 ci.vertices.height = 0.15,
                 clip=c(0,7),
                 xticks = c(0,1,2,3,4,5,6,7),
                 xticks.digits = 10,
                 txt_gp=fpTxtGp(label=gpar(cex=1,fontface = "bold"),
                                ticks=gpar(cex=0.8, fontface = "bold"),
                                xlab=gpar(cex = 1, fontface = "bold"),
                                title=gpar(cex = 1.25, fontface = "bold")),
                 graph.pos=5
)
fp
png("03.unicox.png", width = 13, height = 7, units = "in", res = 300,family = "Times")
print(fp)
dev.off()
pdf("03.unicox.pdf", width = 13, height = 7,family = "Times")
print(fp)
dev.off()
#02.multi_cox-------------------------
uniresult$pvalue <- as.numeric(uniresult$pvalue)
phresult$pval <- as.numeric(phresult$pval)
hubvars <- intersect(subset(uniresult,pvalue < 0.05)[,"var"],subset(phresult,pval > 0.05)[,"var1"])
hubvars
#perform the multi-variates cox regression using qualified genes.
formula_for_multivariate <- as.formula(paste0('Surv(os_time,os_status)~',paste(hubvars,sep = "",collapse = "+")))
multi_variate_cox <- coxph(formula_for_multivariate,data = df)
multi_variate_cox
#check if variances are supported by PH hypothesis
ph_hypo_multi <- cox.zph(multi_variate_cox)
#The last row of the table records the test results on the GLOBAL model.Delete it.
ph_hypo_table <- ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]
#Remove variances not supported by ph hypothesis and perform the 2nd regression.
formula_for_multivariate <- as.formula(paste0('Surv(os_time,os_status)~',paste(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05],sep = '',collapse = '+')))
multi_variate_cox_2 <- coxph(formula_for_multivariate,data = df)
multi_variate_cox_2
multisum <- summary(multi_variate_cox_2)
multi_cox_df_2 <- multisum$coefficients %>% as.data.frame()
colnames(multi_cox_df_2)[5] <- "p"
multi_cox_genes <- rownames(subset(multi_cox_df_2,p < 0.05))
ph_hypo_multi <- cox.zph(multi_variate_cox_2)
ph_hypo_multi
ph_p <- ggcoxzph(ph_hypo_multi)
pdf("./PH_test/multi.pdf",width = 10, height = 8,family = "Times")
print(ph_p,newpage = FALSE)
dev.off()
png("./PH_test/multi.png",width = 10, height = 8,units = "in", res = 300,family = "Times")
print(ph_p,newpage = FALSE)
dev.off()
ph_hypo_table <- ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]
library(openxlsx)
write.xlsx(ph_hypo_table,"08.ph_hypo_multi.xlsx")
#可视化
multisum <- summary(multi_variate_cox_2)
pvalue <- round(multisum$coefficients[,5],3) 
hubvars <- c("Stage","T_Stage","Gender","Risk_Score")
multiresult <- as.data.frame(cbind(var=hubvars,
                                   HR=multisum$coefficients[,2],
                                   L95CI=multisum$conf.int[,3],
                                   H95CI=multisum$conf.int[,4],
                                   pvalue=multisum$coefficients[,5]
))
rownames(multiresult) <- NULL
multiresult
#FORESTPLOT#多因素cox可视化森林图
hr <- round(as.numeric(multiresult$"HR"),3)#从tducs数据中提取HR值，保留2位小数
hrLow  <- round(as.numeric(multiresult$"L95CI"),3)#提取HR值95%置信区间的低值
hrHigh <- round(as.numeric(multiresult$"H95CI"),3)#提取HR值95%置信区间的高值
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")#将HR及其置信区间数据组合成一个因子
PH_pvalue <- round(as.numeric(ph_hypo_multi$table[c(1:length(hubvars)),"p"]),3)
pValue <- ifelse(as.numeric(multiresult$"pvalue")<0.001, "<0.001", round(as.numeric(multiresult$"pvalue"),3))#提取p值，且当p<0.001时显示为<0.001，而不是确切数值
multicoxdat <- cbind(hubvars,pValue,PH_pvalue,Hazard.ratio,hr,hrLow,hrHigh)
multicoxdat <- rbind(colnames(multicoxdat),multicoxdat) %>% as.data.frame()
multicoxdat[1,c(1,5:7)] <- NA
multicoxdat[,5:7] <- as.numeric(unlist(multicoxdat[,5:7]))
fp2 <- forestplot(labeltext=as.matrix(multicoxdat[,c(1:4)]),
                  mean = multicoxdat$hr,
                  lower=multicoxdat$hrLow,
                  upper=multicoxdat$hrHigh,
                  is.summary = c(TRUE,rep(FALSE,8)),
                  zero =1,
                  boxsize = 0.15,
                  lineheight = unit(20,'mm'),
                  graphwidth = unit(180,'mm'),
                  colgap = unit(6,'mm'),
                  lwd.zero = 2,
                  lwd.ci = 2,
                  col=fpColors(box = '#ef233c',summary = "#9e2a2b",lines='blue1',zero = 'grey'),
                  xlab = 'Hazard Ratio',
                  lwd.xaxis = 2,
                  lty.ci = "solid",
                  ci.vertices.height = 0.15,
                  clip=c(0,40),
                  xticks = c(0,1,2),
                  xticks.digits = 10,
                  txt_gp=fpTxtGp(label=gpar(cex=1,fontface = "bold"),
                                 ticks=gpar(cex=0.8, fontface = "bold"),
                                 xlab=gpar(cex = 1, fontface = "bold"),
                                 title=gpar(cex = 1.25, fontface = "bold")),
                  graph.pos=5
)
fp2
png("04.multicox.png", width = 13, height = 7, units = "in", res = 300,family = "Times")
print(fp2)
dev.off()
pdf("04.multicox.pdf", width = 13, height = 7,family = "Times")
print(fp2)
dev.off()
#03.列线图-------------------------------------------------------------------------
library(regplot)
library(survival)
library(rms)

df$os_time <- df$os_time/365

multiCoxModel <- coxph(formula_for_multivariate, data =df)
reg_plot <- function(){
  regplot(multiCoxModel,
          observation = df[33,],
          showP = TRUE,
          failtime = c(1,3,5), # 定义生存的1,3,5概率标度
          plots = c("violin","boxes"), # 连续性变量形状，可选"no plot" "density" "boxes" "ecdf" "bars" "boxplot" "violin" "bean" "spikes"；分类变量的形状，可选"no plot" "boxes" "bars" "spikes"
          points = T, # 截距项显示为0-100
          prfail = T,font.family = "Times")
}
reg_plot()
p_base <- recordPlot()
png(file = "05.nomogram.png",width= 15, height= 9, units="in", res=600,family = "Times")
p_base
dev.off()
pdf("05.nomogram.pdf", width = 15, height = 9,family = "Times")
p_base
dev.off()
#04.校准曲线calibration------------------------------------------------
df <- dat1
df$Stage <- as.factor(df$Stage) %>% as.numeric()
df$T_Stage <- as.factor(df$T_Stage) %>% as.numeric()
df$os_time <- df$os_time/365
oneyear <- cph(formula =  formula_for_multivariate,
               data = df, x = T,y = T,surv = T, time.inc = 1)
cal_oneyear <- calibrate(oneyear, cmethod = "KM", method = "boot", 
                         u = 1, m = 50, B = 1000) # m，分组到平均包含 m 个受试者的区间；
threeyear <- cph(formula =  formula_for_multivariate,
                 data = df, x = T,y = T,surv = T, time.inc = 3)
cal_threeyear <- calibrate(threeyear, cmethod = "KM", method = "boot", 
                           u = 3, m = 50, B = 1000)
fiveyear <- cph(formula =  formula_for_multivariate,
                data = df, x = T,y = T,surv = T, time.inc = 5)
cal_fiveyear <- calibrate(fiveyear, cmethod = "KM", method = "boot",
                          u = 5, m = 50, B = 1000)
#可视化
cal_plot <- function(){
  plot(cal_oneyear,                   # 绘制校准曲线的数据
       lwd=1,                         # 线条宽度为1
       lty=1,                         # 线条类型为1（实线）
       conf.int=T,                    # 是否显示置信区间
       errbar.col="black",            # 直线和曲线的误差线颜色设置为蓝色
       col="red3",                    # 校准曲线的颜色设置为红色
       xlim=c(0,1),                   # x轴的限制范围，从0到1
       ylim=c(0,1),                   # y轴的限制范围，从0到1
       xlab="Nomogram-prediced OS (%)",  # x轴标签
       ylab="Observed OS (%)",                # y轴标签
       subtitles = F)  
  lines(cal_oneyear[,c('mean.predicted',"KM")],type = 'b', lwd = 2, col = c("red"), pch = 16)
  mtext("")
  plot(cal_threeyear,
       lwd=1,
       lty=1,
       conf.int=T,
       errbar.col="black", 
       col="blue",add = T)
  lines(cal_threeyear[,c('mean.predicted',"KM")],type = 'b', lwd = 2, col = c("blue"), pch = 16)
  plot(cal_fiveyear,
       lwd=1,
       lty=1,
       conf.int=T,
       errbar.col="black", 
       col="green",add = T)
  lines(cal_fiveyear[,c('mean.predicted',"KM")],type = 'b', lwd = 2, col = c("green"), pch = 16)
  abline(0,1, lwd = 3, lty = 3, col = c("#224444"))
  legend("topleft", # 图例的位置
         legend = c("1 Years","3 Years","5 Years"), # 图例文字
         col =c("red","blue","green"), # 图例线的颜色，与文字对应
         lwd = 2, # 图例中线的粗细
         cex = 1.6, # 图例字体大小
         bty = "n") # 不显示图例边框
}
cal_plot()
png("06.calibration.png", width = 9, height = 9, units = "in", res = 300,family = "Times")
cal_plot()
dev.off()
pdf("06.calibration.pdf", width = 9, height = 9,family = "Times")
cal_plot()
dev.off()
#05.ROC-------------------------------------------------------------------------
#计算列线图得分
df <- dat1
df$os_time <- df$os_time/365
df$Stage <- as.factor(df$Stage) %>% as.numeric()
df$T_Stage <- as.factor(df$T_Stage) %>% as.numeric()
res.mul <- coxph(formula_for_multivariate, data = df)
df$pred <- predict(res.mul, newdata = df)
df$Risk_Group <- ifelse(df$Risk_Group == 'high',1,0)
#采取1，3，5年多个时间点的ROC
multi_ROC <- function(time_vector,risk_score_table){
  library('survivalROC')
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime = risk_score_table$os_time,status = risk_score_table$os_status,
                           marker = risk_score_table$pred,
                           predict.time = single_time,method = "KM")
    data.frame('True_positive' = for_ROC$TP,'False_positive' = for_ROC$FP,
               'Cut_values' = for_ROC$cut.values,'Time_point' = rep(single_time,length(for_ROC$TP)),
               'AUC' = rep(for_ROC$AUC,length(for_ROC$TP)))	
  }
  multi_ROC_list <- lapply(time_vector,single_ROC)
  do.call(rbind,multi_ROC_list)
}
#evaluate AUCs between 1,3,5 years.
for_multi_ROC <- multi_ROC(time_vector = c(seq(1,5,2)),risk_score_table = df)
for_multi_ROC$Time_point <- as.factor(for_multi_ROC$Time_point)
# 按时间分组，提取每个时间点的唯一AUC值
auc_by_time <- for_multi_ROC %>%
  group_by(Time_point) %>%
  dplyr::summarise(AUC = unique(AUC))  # 使用 first(AUC) 或 unique(AUC)
# 输出结果
auc_by_time
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
             "1 Year AUC = ", sprintf("%.2f", auc_by_time[1,2]), "\n",
             "3 Year AUC = ", sprintf("%.2f", auc_by_time[2,2]), "\n",
             "5 Year AUC = ", sprintf("%.2f", auc_by_time[3,2]), "\n",
             sep = '')
  )+
  scale_colour_manual(values = c("#F1C89A","#E79397","#A797DA"),labels = c("1 Year","3 Year","5 Year"))+
  labs(colour = "Time",title = "TCGA LIHC Nomogram")+ggfonttheme
pROC
ggsave("07.ROC.pdf",pROC,device = cairo_pdf,width = 8, height = 7)
ggsave("07.ROC.png",pROC,dpi = 300,width = 8, height = 7)
#05.DCA曲线-----------------------------------------------------
df <- dat1
df$os_time <- df$os_time/365
df$Stage <- as.factor(df$Stage) %>% as.numeric()
df$T_Stage <- as.factor(df$T_Stage) %>% as.numeric()
library(ggDCA)
library(ggsci)
# 建立多个模型
cox_fit1 <- coxph(Surv(os_time, os_status) ~ Stage+T_Stage+Gender+Risk_Group, 
                  data = df)
cox_fit2 <- coxph(Surv(os_time, os_status) ~ Stage, 
                  data = df)
cox_fit3 <- coxph(Surv(os_time, os_status) ~ T_Stage, 
                  data = df)
cox_fit4 <- coxph(Surv(os_time, os_status) ~ Gender, 
                  data = df)
cox_fit5 <- coxph(Surv(os_time, os_status) ~ Risk_Group, 
                  data = df)
DCAdf1 <- ggDCA::dca(cox_fit1,cox_fit2,cox_fit3,cox_fit4,cox_fit5,times = c(1,3,5))
ggplot(DCAdf1,linetype = F)+
  scale_color_jama(name="Model Type",
                   labels=c("Nomogram","Stage","T Stage","N Stage","Risk Score","All","None"))+
  theme_bw(base_size = 14)+
  theme(legend.position = "bottom",
        legend.background = element_blank()
  )+ggfonttheme
ggsave("09.DCA.pdf",device = cairo_pdf,width = 8, height = 7)
ggsave("09.DCA.png",dpi = 300,width = 8, height = 7)
save.image("../RData/06.nomogram.RData")
