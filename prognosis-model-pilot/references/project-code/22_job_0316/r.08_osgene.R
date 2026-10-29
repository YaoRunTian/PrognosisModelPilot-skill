rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/22_job_0316")
if (!dir.exists("08_osgene")) {dir.create("08_osgene")}
setwd("08_osgene")
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(gghalves)
ggfonttheme <- theme(text = element_text(family = "Times"))
osgenes <-read.table("../03_model/10.multiresult.txt") %>% rownames()
clinical <- readRDS("../00_rawdata/clinical.Rds")
expr <- readRDS("../00_rawdata/tpm.Rds")
df1 <- as.data.frame(t(expr[osgenes,clinical$sample]))
df1$Sample <- rownames(df1)
df1 <- melt(df1)
colnames(df1) <- c("Sample","Symbol","value")
type <- rep(as.character(clinical$group),length(osgenes))
df1$Type <- type
my_comparisons <- list(c("Normal","Tumor"))
#绘制组间差异小提琴图
ggplot(df1,aes(x = Type, y = value,fill = Type))+
  geom_half_violin(side="r",width=0.7,color=NA)+
  geom_half_boxplot(width=0.4,size=1,side="l",outlier.color=NA)+
  #geom_jitter(shape=21,size=1,position = position_jitter(width = 0.1))+
  scale_fill_manual(values = c("Normal" = "#6A97FA","Tumor" = "#FA736A"))+
  facet_wrap(~ Symbol, scales = "free", ncol = 4) + 
  theme_bw()+
  theme(panel.grid = element_line(colour = "white"),plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        title = element_text(size = 20),
        legend.position = "top",
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 16, face = "bold"),
        strip.background = element_blank())+
  ylab("Gene Expression")+
  xlab("")+
  ggtitle("Expression of prognostic genes")+
  stat_compare_means(aes(group=Type),                       #按分组进行统计检验
                     method = "t.test",
                     paired = FALSE,                             #非配对t检验
                     symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                        symbols=c("***","**","*","ns")),
                     label = "p.signif",
                     label.x = 1.4,
                     label.y.npc = 0.95,      #添加显著性符号的位置
                     size=4.5)+ggfonttheme
ggsave("01.osgene.pdf",device = cairo_pdf,width = 8, height = 6)
ggsave("01.osgene.png",dpi = 300,width = 8, height = 6)

stat_table <- compare_means(
  value ~ Type, 
  data = df1, 
  group.by = "Symbol",
  method = "t.test",
  p.adjust.method = "BH"
)
print(stat_table)
write.xlsx(stat_table,"02.tets_res.xlsx")
save.image("../RData/08_osgene.RData")
