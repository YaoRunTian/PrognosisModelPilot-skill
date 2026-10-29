rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/22_job_0316")
if (!dir.exists("16_mutation")) {dir.create("16_mutation")}
setwd("16_mutation")
#突变数据
library(TCGAmutations)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)
maf <- tcga_load(c("LIHC"))
#读取高低表达组
risk_table <-  read.table("../04_evaluate/04.risk_table.txt",sep = "\t",header = TRUE)#读取高低风险组信息
samples <- rownames(risk_table)
colnames(risk_table)[6] <- "RiskGroup"
risk_table$sample <- str_sub(rownames(risk_table),1,16)
# 提取了每个样本的名称(前者是1：16位的/后者是全称)
samp = data.frame(ID = str_sub(unique(maf@data$Tumor_Sample_Barcode),1,16),
                  long = unique(maf@data$Tumor_Sample_Barcode))

# merge() 函数会自动根据指定的列（如 ID 列）将两个数据框的数据进行配对合并
samp= merge(risk_table,samp,by.x = "sample",by.y = "ID")
#309个样本
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
#绘制棒棒糖图-------------------------------------------------------------------
#低风险组
#读取高低表达组
low_table <- subset(risk_table,RiskGroup == "low")
maf_low <- maf
# 提取了每个样本的名称(前者是1：16位的/后者是全称)
samp = data.frame(ID = str_sub(unique(maf_low@data$Tumor_Sample_Barcode),1,16),
                  long = unique(maf_low@data$Tumor_Sample_Barcode))

# merge() 函数会自动根据指定的列（如 ID 列）将两个数据框的数据进行配对合并
samp= merge(low_table,samp,by.x = "sample",by.y = "ID")#184

maf_low@clinical.data = left_join(maf_low@clinical.data,samp,by = c("Tumor_Sample_Barcode"="long"))
maf_low = subsetMaf(maf_low,tsb = samp$long)
#高风险组
high_table <- subset(risk_table,RiskGroup == "high")
maf_high <- maf
# 提取了每个样本的名称(前者是1：16位的/后者是全称)
samp = data.frame(ID = str_sub(unique(maf_high@data$Tumor_Sample_Barcode),1,16),
                  long = unique(maf_high@data$Tumor_Sample_Barcode))
# merge() 函数会自动根据指定的列（如 ID 列）将两个数据框的数据进行配对合并
samp= merge(high_table,samp,by.x = "sample",by.y = "ID")#125个样本
# 其实上面已经准备好了，再来一次吧~ 
maf_high@clinical.data = left_join(maf_high@clinical.data,samp,by = c("Tumor_Sample_Barcode"="long"))
maf_high = subsetMaf(maf_high,tsb = samp$long)
#可视化top3
df1 <- data.frame(table(maf@data$Hugo_Symbol))
df2 <- df1
rownames(df2) <- df2$Var1
hubgenes <- read.table("../03_model/10.multiresult.txt") %>% rownames()
df2 <- df2[hubgenes,]
for(hubgene in hubgenes){
  png_filename <- paste0(hubgene,".png")
  pdf_filename <- paste0(hubgene,".pdf")
  png(png_filename, width = 13, height = 8, units = "in", res = 300)
  lollipopPlot2(m1 = maf_high,m2 = maf_low,gene = hubgene,AACol1 = 'HGVSp_Short',AACol2= "HGVSp_Short",
                m1_name = "High Risk Group",m2_name = "Low Risk Group",showDomainLabel = FALSE)
  dev.off()
  
  pdf(pdf_filename, width = 13, height =8)
  lollipopPlot2(m1 = maf_high,m2 = maf_low,gene = hubgene,AACol1 = 'HGVSp_Short',AACol2= "HGVSp_Short",
                m1_name = "High Risk Group",m2_name = "Low Risk Group",showDomainLabel = FALSE)
  dev.off()
}
save.image("../RData/16_mutation.RData")
