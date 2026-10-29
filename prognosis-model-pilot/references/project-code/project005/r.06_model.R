rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project005")
if (!dir.exists("06_model")) {dir.create("06_model")}
setwd("06_model")
library(survival)
library(dplyr)
library(survminer)
ggfonttheme <- theme(text = element_text(family = "Times"))
#00.预处理数据------------------------------------------------------------------
hubgenes <- read.table("../05_venn/01.hubgenes.txt",sep = "\n")[,1]#读取关键基因
fpkm <- readRDS("../00_rawdata/fpkm.Rds")#表达矩阵
groupInf <- readRDS("../00_rawdata/group.Rds")#分组信息
clinical <- readRDS("../00_rawdata/os.Rds")#生存信息
groupInf$sample <-gsub("-","_",groupInf$sample)
expr <- fpkm[,groupInf$sample] %>% t() %>% as.data.frame()
expr <- expr[,hubgenes]#选取关键基因
#合并表达谱数据和临床生存信息
dat1 <- merge(clinical[,1:3],expr,by.x = "sample",by.y = 0)
colnames(dat1)[c(2,3)] <- c("os_time","os_status")
rownames(dat1) <- dat1[,1]
dat1 <- dat1[,-1]
samplelist <- groupInf$sample[grep("Tumor",groupInf$group)]
# 随机分割样本为7:3
set.seed(3)
train_index <- sample(1:length(samplelist), size = round(0.7 * length(samplelist)))
train_samples <- samplelist[train_index]
test_samples <- samplelist[-train_index]
# 创建训练集和测试集数据
train_data <- dat1[rownames(dat1) %in% train_samples, ]
test_data <- dat1[rownames(dat1) %in% test_samples,]
# train_data <- readRDS("../RData/train_data.Rds")
# train_samples <- rownames(train_data)
# test_data <- readRDS("../RData/test_data.Rds")
# test_samples <- rownames(test_data)
dat1 <- train_data
#单因素cox回归分析--------------------------------------------------------------
#设置p值的阈值
pfilter <- 0.05
#新建空白数据框
uniresult <- data.frame()
unifiltered <- data.frame()
#单因素COX回归分析中p值＜0.05的基因，其分析结果输入到之前新建的空白数据框uniresult中
dat1 <- na.omit(dat1)
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
uniresult#3个
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
hr <- round(as.numeric(uniresult$"HR"),3)#从tducs数据中提取HR值，保留2位小数
hrLow  <- round(as.numeric(uniresult$"L95CI"),3)#提取HR值95%置信区间的低值
hrHigh <- round(as.numeric(uniresult$"H95CI"),3)#提取HR值95%置信区间的高值
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
                 is.summary = c(TRUE,rep(FALSE,2)),
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
                 xticks = c(0,0.5,1,1.5),
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
  annotate('text', x = -6.0,
           y = -0.1, 
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
  annotate('text', x = -6.0, y =12.55, 
           label = paste0('2 Vars\nlog(lambda.min) = ', round(log(best_lambda), 3)), 
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
save.image("../RData/06_model.RData")
