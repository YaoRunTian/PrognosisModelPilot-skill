rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/22_job_0316")
if (!dir.exists("07_osgene_survival")) {dir.create("07_osgene_survival")}
setwd("07_osgene_survival")
library(survival)
library(survminer)
ggfonttheme <- theme(text = element_text(family = "Times"))
genes <- read.table("../03_model/10.multiresult.txt") %>% rownames()
dat1 <- read.table("../04_evaluate/04.risk_table.txt")
# 确保 os_time 和 os_status 格式正确
dat1$os_time <- as.numeric(dat1$os_time)
dat1$os_status <- as.integer(dat1$os_status)

cut_res <- surv_cutpoint(
  dat1,
  time = "os_time",
  event = "os_status",
  variables = genes,
  minprop = 0.3
)
dat_cut <- surv_categorize(cut_res)
# dat_cut 中对应基因列会变成 low / high
fit_fcn3 <- survfit(Surv(os_time, os_status) ~ FCN3, data = dat_cut)
fcn3_km <- ggsurvplot(fit_fcn3, data = dat_cut,
           pval = TRUE, risk.table = TRUE,palette = c('red','blue')
           )
fcn3_km
pdf("01.FCN3_km.pdf",width = 8, height = 7,family = "Times")
print(fcn3_km,newpage = FALSE)
dev.off()
png("01.FCN3_km.png",width = 8, height = 7,units = "in", res = 300,family = "Times")
print(fcn3_km,newpage = FALSE)
dev.off()

fit_pglyrp4 <- survfit(Surv(os_time, os_status) ~ PGLYRP4, data = dat_cut)
pglyrp4_km <- ggsurvplot(fit_pglyrp4, data = dat_cut, 
           pval = TRUE, risk.table = TRUE,palette = c('red','blue'))
pglyrp4_km
pdf("02.PGLYRP4_km.pdf",width = 8, height = 7,family = "Times")
print(pglyrp4_km,newpage = FALSE)
dev.off()
png("02.PGLYRP4_km.png",width = 8, height = 7,units = "in", res = 300,family = "Times")
print(pglyrp4_km,newpage = FALSE)
dev.off()
