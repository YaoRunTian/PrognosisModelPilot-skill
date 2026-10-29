rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project003")
if (! dir.exists("03_model")){dir.create("03_model")}
setwd("03_model")
library(survival)
library(dplyr)
library(survminer)
library(ggplot2)
ggfonttheme <- theme(text = element_text(family = "Times"))
#load the data
hubgenes <- read.table("../02_venn/01.hubgenes.txt",sep = "\n")[,1]#读取关键基因
fpkm <- readRDS("../00_rawdata/tpm.Rds")#表达矩阵
clinical <- readRDS("../00_rawdata/clinical.Rds")#生存信息
clinical <- subset(clinical,group == "Tumor")
expr <- fpkm[,clinical$sample] %>% t() %>% as.data.frame()
expr <- expr[,hubgenes]#选取关键基因
#合并表达谱数据和临床生存信息
dat1 <- merge(clinical[,c(1,6,5)],expr,by.x = 1,by.y = 0)
colnames(dat1)[c(2,3)] <- c("os_time","os_status")
rownames(dat1) <- dat1[,1]
dat1 <- dat1[,-1]
#单因素cox回归分析--------------------------------------------------------------
#设置p值的阈值
pfilter <- 0.05
#新建空白数据框
uniresult <- data.frame()
unifiltered <- data.frame()
#单因素COX回归分析中p值＜0.05的基因，其分析结果输入到之前新建的空白数据框uniresult中
dat1 <- na.omit(dat1)#318个样本具有生存信息
df <- dat1
for(i in colnames(df[,3:ncol(df)])){
  print(i)
  unicox <- coxph(Surv(time = os_time, event = os_status) ~ df[,i], data = df)  
  unisum<- summary(unicox) 
  pvalue <- round(unisum$coefficients[,5],3) 
  print(paste("pvalue = ",pvalue))
  ph_p <- cox.zph(unicox)$table[1,3]
  if(pvalue<pfilter&ph_p >0.05){ 
    uniresult <- rbind(uniresult,
                       cbind(gene=i,
                             HR=unisum$coefficients[,2],
                             L95CI=unisum$conf.int[,3],
                             H95CI=unisum$conf.int[,4],
                             pvalue=unisum$coefficients[,5],
                             PH_res = ph_p
                       ))
  }else{unifiltered <- rbind(unifiltered,
                             cbind(gene=i,
                                   HR=unisum$coefficients[,2],
                                   L95CI=unisum$conf.int[,3],
                                   H95CI=unisum$conf.int[,4],
                                   pvalue=unisum$coefficients[,5],
                                   PH_res = ph_p
                             ))}
}
uniresult#5个
unires <- rbind(uniresult,unifiltered)
unires
write.csv(unires,file = "01.single_cox.csv",row.names = F)
gene <- uniresult[,"gene"]#提取基因名
#绘制PH假定检验图
if (! dir.exists("PH_test")){dir.create("PH_test")}
phresult <- data.frame()
for(i in gene){
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
                    cbind(gene = i,
                          pval = test.ph$table[1,3]))
}
phresult
write.csv(phresult,file = "02.PHresult.csv",row.names = FALSE)
#FORESTPLOT#单因素cox可视化森林图
hr <- round(as.numeric(uniresult$"HR"),2)#从tducs数据中提取HR值，保留2位小数
hrLow  <- round(as.numeric(uniresult$"L95CI"),2)#提取HR值95%置信区间的低值
hrHigh <- round(as.numeric(uniresult$"H95CI"),2)#提取HR值95%置信区间的高值
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")#将HR及其置信区间数据组合成一个因子
pValue <- ifelse(as.numeric(uniresult$"pvalue")<0.001, "<0.001", round(as.numeric(uniresult$"pvalue"),3))#提取p值，且当p<0.001时显示为<0.001，而不是确切数值
PH_pvalue <- round(as.numeric(uniresult$PH_res),4)
unicoxdat <- cbind(gene,pValue,PH_pvalue,Hazard.ratio,hr,hrLow,hrHigh)
unicoxdat <- rbind(colnames(unicoxdat),unicoxdat) %>% as.data.frame()
unicoxdat[1,5:7] <- NA
unicoxdat[,5:7] <- as.numeric(unlist(unicoxdat[,5:7]))
library(forestplot)
fp <- forestplot(labeltext=as.matrix(unicoxdat[,c(1:4)]),
                 mean = unicoxdat$hr,
                 lower=unicoxdat$hrLow,
                 upper=unicoxdat$hrHigh,
                 is.summary = c(TRUE,rep(FALSE,15)),
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
                 clip=c(0,100),
                 xticks = c(0,0.5,1,1.5,2),
                 xticks.digits = 10,
                 txt_gp=fpTxtGp(label=gpar(cex=1,fontface = "bold"),
                                ticks=gpar(cex=0.8, fontface = "bold"),
                                xlab=gpar(cex = 1, fontface = "bold"),
                                title=gpar(cex = 1.25, fontface = "bold")),
                 graph.pos=5
)
fp
png("03.unicox.png", width = 12, height = 6, units = "in", res = 300,family = "Times")
print(fp)
dev.off()
pdf("03.unicox.pdf", width = 12, height = 6,family = "Times")
print(fp)
dev.off()
#lasso回归-------------------------------------------------------------------
library("glmnet")
library("survival")
lasso_dat <- subset(df,os_time != 0)
x <- scale(as.matrix(lasso_dat[,uniresult$gene]))
y <- lasso_dat[,c("os_time","os_status")]
names(y) <- c('time','status')
y$time <- as.double(y$time)
y$status <- as.double(y$status)
y <- as.matrix(survival::Surv(y$time,y$status))
set.seed(42)
la.eq <- glmnet(x, y, family="cox", alpha=1)
set.seed(42)
lasso_fit <- cv.glmnet(x,y,family = 'cox',alpha = 1)
lasso_fit
#可视化
cv_model  <- lasso_fit
best_lambda <- cv_model$lambda.min
lambda_lse <- cv_model$lambda.1se
lambda_values <- cv_model$lambda
coef_list <- lapply(lambda_values, function(lambda) coef(cv_model, s = lambda))
coef_matrices <- lapply(coef_list, function(coef) as.matrix(coef))
coef_df_list <- lapply(coef_matrices, function(matrix) as.data.frame(t(matrix)))
coef_df <- do.call(rbind, coef_df_list)
coef_df$lambda <- rep(lambda_values, each = nrow(coef_df_list[[1]]))
coef_df <- coef_df[!rownames(coef_df) %in% "(Intercept)", ]
colnames(coef_df) <- c(colnames(coef_df)[-ncol(coef_df)], "lambda")
library(reshape2)
coef_long <- melt(coef_df, id.vars = "lambda", variable.name = "Variable", value.name = "Coefficient")
coef_long$log_lambda <- log(coef_long$lambda)
# 加载必要的包
library(ggplot2)
library(ggsci)
table(coef_long$Variable)
p1 <- ggplot(coef_long, aes(x = log_lambda, y = Coefficient, color = Variable)) + 
  geom_vline(xintercept = log(best_lambda), size = 0.8, color = 'grey60', alpha = 0.8, linetype = 2) +
  # geom_vline(xintercept = log(lambda_lse), size = 0.8, color = 'grey60', alpha = 0.8, linetype = 2) +
  geom_line(size = 1) +
  xlab("Log Lambda") + 
  ylab("Coefficients") + 
  theme_bw(base_rect_size = 2) +
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 15, color = 'black'), 
        axis.text = element_text(size = 12, color = 'black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12, color = 'black'), 
        legend.position = 'bottom') +
  annotate('text', x = -4.7,
           y = 0.2, 
           label = paste0('log(lambda.min)\n', round(log(best_lambda), 3)), 
           color = 'black', size = 5, family = 'serif') +
  # annotate('text', x = log(lambda_lse),
  #          y = 10, 
  #          label = paste0('log(lambda.lse)\n', round(log(lambda_lse), 4)), 
  #          color = 'black', size = 5, family = 'serif') +
  guides(col = guide_legend(ncol = 5))+ggfonttheme
p1
ggsave('04.lasso.Coef.pdf',p1,w=8,h=5)
ggsave('04.lasso.Coef.png',p1,w=8,h=5)
xx <- data.frame(
  lambda = lasso_fit$lambda,
  cvm = lasso_fit$cvm,
  cvsd = lasso_fit$cvsd,
  cvup = lasso_fit$cvup,
  cvlo = lasso_fit$cvlo,
  nozezo = lasso_fit$nzero
)
xx$ll <- log(xx$lambda)
xx$NZERO <- paste0(xx$nozezo, ' vars')
dev.off()
p2 <- ggplot(xx, aes(ll, cvm, color = NZERO)) +
  geom_errorbar(aes(ymin = cvlo, ymax = cvup), width = 0.05, size = 1) +
  geom_vline(xintercept = log(best_lambda), size = 0.8, color = 'grey60', alpha = 0.8, linetype = 2) +
  # geom_vline(xintercept = log(lambda_lse), size = 0.8, color = 'grey60', alpha = 0.8, linetype = 2) +
  geom_point(size = 2) +
  xlab("Log Lambda") +
  ylab('Partial Likelihood Deviance') +
  theme_bw(base_rect_size = 1.5) +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 15, color = 'black'),
    axis.text = element_text(size = 12, color = 'black'),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, color = 'black'),
    legend.position = 'bottom'
  ) +
  annotate('text', x = -4.8, y =11.9, 
           label = paste0('5 Vars\nlog(lambda.min) = ', round(log(best_lambda), 3)), 
           color = 'black', size = 5, family = 'Times')+
  guides(col = guide_legend(ncol = 5))+ggfonttheme
print(p2)
ggsave('05.lasso.CV.pdf',p2,w=8,h=5)
ggsave('05.lasso.CV.png',p2,w=8,h=5)
#计算相关预后基因系数
coefficients <- coef(lasso_fit,s = lasso_fit$lambda.min)
Active.Index <- which(as.numeric(coefficients) != 0)
active.coefficients <- as.numeric(coefficients)[Active.Index]
lasso_gene <- rownames(coefficients)[Active.Index]
lasso_gene
hubgeneslist <- coefficients@Dimnames[[1]][Active.Index]
lasso_result <- data.frame(coef = coefficients@x)
rownames(lasso_result) <- hubgeneslist
coef <- lasso_result[lasso_gene,]
lasso_result <- data.frame(coef = coef)
rownames(lasso_result) <- lasso_gene
lasso_result
write.table(lasso_result,"07.lassoresult.txt",row.names = TRUE,col.names = TRUE,quote = FALSE)
write.table(lasso_gene,"06.lassogeneslist.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.csv(lasso_result,"07.lassoresult.csv",row.names = TRUE,quote = FALSE)
write.csv(lasso_gene,"06.lassogeneslist.csv",row.names = FALSE,quote = FALSE)
##多因素分析-------------------
formula_for_multivariate <- as.formula(paste0('Surv(os_time,os_status)~',paste(lasso_gene,sep = "",collapse = "+")))
multi_variate_cox <- coxph(formula_for_multivariate,data = df)
#check if variances are supported by PH hypothesis
ph_hypo_multi <- cox.zph(multi_variate_cox)
ph_hypo_multi
#The last row of the table records the test results on the GLOBAL model.Delete it.
ph_hypo_table <- ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]
#Remove variances not supported by ph hypothesis and perform the 2nd regression.
formula_for_multivariate <- as.formula(paste0('Surv(os_time,os_status)~',paste(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05],sep = '',collapse = '+')))
multi_variate_cox_2 <- coxph(formula_for_multivariate,data = df)
multi_variate_cox_2
# Call:
#   coxph(formula = formula_for_multivariate, data = df)
# 
# coef exp(coef) se(coef)      z        p
# FCN3    -0.21872   0.80355  0.07762 -2.818 0.004834
# COLEC12  0.17191   1.18757  0.09186  1.871 0.061280
# MBL2    -0.02374   0.97654  0.05445 -0.436 0.662844
# LY96     0.07620   1.07918  0.06544  1.164 0.244228
# PGLYRP4  1.15492   3.17378  0.34143  3.383 0.000718
# 
# Likelihood ratio test=25.72  on 5 df, p=0.000101
# n= 318, number of events= 95 
multisum <- summary(multi_variate_cox_2)
multi_cox_df_2 <- multisum$coefficients %>% as.data.frame()
colnames(multi_cox_df_2)[5] <- "p"
multi_cox_genes <- rownames(subset(multi_cox_df_2,p < 0.05))
formula_for_multivariate_2 <-  as.formula(paste0('Surv(os_time,os_status)~',paste(multi_cox_genes,sep = '',collapse = '+')))
multi_variate_cox_3 <- coxph(formula_for_multivariate_2,data = df)
multi_variate_cox_3
# Call:
#   coxph(formula = formula_for_multivariate_2, data = df)
# 
# coef exp(coef) se(coef)      z        p
# FCN3    -0.22204   0.80088  0.07516 -2.954  0.00313
# PGLYRP4  1.38838   4.00835  0.30769  4.512 6.41e-06
# 
# Likelihood ratio test=19.45  on 2 df, p=5.99e-05
# n= 318, number of events= 95 

#ph检验
ph_hypo_multi <- cox.zph(multi_variate_cox_3)
ph_hypo_multi
ph_p <- ggcoxzph(ph_hypo_multi)
pdf("11.multi_ph.pdf",width = 5, height = 10)
print(ph_p,newpage = FALSE)
dev.off()
png("11.multi_ph.png",width = 5, height = 10,units = "in", res = 300)
print(ph_p,newpage = FALSE)
dev.off()
#生成多因素cox分析结果
multi_result <- data.frame(coef = multi_variate_cox_3$coefficients)
rownames(multi_result) <- multi_cox_genes
multisum <- summary(multi_variate_cox_3)
pvalue <- round(multisum$coefficients[,5],3) 
multiresult <- as.data.frame(cbind(HR=multisum$coefficients[,2],
                                   L95CI=multisum$conf.int[,3],
                                   H95CI=multisum$conf.int[,4],
                                   pvalue=multisum$coefficients[,5]
))
multiresult$var <- rownames(multiresult)
rownames(multiresult) <- NULL
multiresult <- multiresult[,c(5,1:4)]
multiresult
# 获取PH检验p值
ph_hypo_multi <- cox.zph(multi_variate_cox_3)
ph_hypo_table_final <- ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]
PH_pvalue_multi <- round(ph_hypo_table_final[,3], 4)
# 修改森林图数据框构建部分
hr <- round(as.numeric(multiresult$"HR"),2)
hrLow  <- round(as.numeric(multiresult$"L95CI"),2)
hrHigh <- round(as.numeric(multiresult$"H95CI"),2)
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pValue <- ifelse(as.numeric(multiresult$"pvalue")<0.001, "<0.001", 
                 round(as.numeric(multiresult$"pvalue"),3))
# 将PH_pvalue加入数据框
multicoxdat <- cbind(multiresult$var, pValue, PH_pvalue_multi, Hazard.ratio, hr, hrLow, hrHigh)
# 添加表头行
multicoxdat <- rbind(colnames(multicoxdat), multicoxdat) %>% as.data.frame()
# 第一行（表头行）的数值列设为NA
multicoxdat[1, c(1, 5:7)] <- NA
multicoxdat[, 5:7] <- as.numeric(unlist(multicoxdat[, 5:7]))
# 绘制森林图（labeltext增加PH_pvalue列）
text_color <- c(rep("black", 14), "red", "black", "red", rep("black", 11))
fp2 <- forestplot(
  labeltext = as.matrix(multicoxdat[, c(1:4)]),  # 现在包含4列文本
  mean   = multicoxdat$hr,
  lower  = multicoxdat$hrLow,
  upper  = multicoxdat$hrHigh,
  is.summary = c(TRUE, rep(FALSE, nrow(multicoxdat) - 1)),
  zero        = 1,
  boxsize     = 0.15,
  lineheight  = unit(10, 'mm'),
  graphwidth  = unit(180, 'mm'),
  colgap      = unit(6, 'mm'),
  lwd.zero    = 2,
  lwd.ci      = 2,
  col = fpColors(
    box     = '#ef233c',
    summary = "#9e2a2b",
    lines   = 'blue1',
    zero    = 'grey',
    text    = text_color
  ),
  xlab              = 'Hazard Ratio',
  lwd.xaxis         = 2,
  lty.ci            = "solid",
  ci.vertices.height = 0.15,
  clip              = c(0,4),
  xticks            = c(0, 1, 2, 3, 4),
  xticks.digits     = 10,
  txt_gp = fpTxtGp(
    label = gpar(cex = 1,   fontface = "bold"),
    ticks = gpar(cex = 0.8, fontface = "bold"),
    xlab  = gpar(cex = 1,   fontface = "bold"),
    title = gpar(cex = 1.25,fontface = "bold")
  ),
  graph.pos = 5   # 图形位置向后移一位，为新列腾出空间
)
fp2
png("07.multicox.png", width = 14, height = 12, units = "in", res = 300)
print(fp2)
dev.off()
pdf("07.multicox.pdf", width = 14, height = 12)
print(fp2)
dev.off()
write.csv(ph_hypo_table,file = "08.multiPHresult.csv",row.names = TRUE)
write.csv(multiresult,file = "09.multi_cox.csv",row.names = FALSE)
write.table(multi_result,"10.multiresult.txt",row.names = TRUE,col.names = TRUE,quote = FALSE)

save.image("../RData/03_model.RData")
