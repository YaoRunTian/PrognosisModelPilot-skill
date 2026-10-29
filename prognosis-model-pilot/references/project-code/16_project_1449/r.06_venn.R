rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/16_project_1449")
if (!dir.exists("06_venn")) {dir.create("06_venn")}
setwd("06_venn")
library("ggvenn")
library("dplyr")
library("openxlsx")
ggfonttheme <- theme(text = element_text(family = "Times"))

DEGs <- read.table("../05_DEG/01.DEGsList.txt",sep = "\n")[,1]
hub_module_gene <- read.xlsx("../04_hdWGCNA/06.hub_module_gene_df.xlsx")[,1]
hubgenes <- intersect(DEGs,hub_module_gene)
write.table(hubgenes,"01.hubgenes.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)
#绘制venn图
list1 <- list(
  "DEGs" = DEGs,
  "Hub Module Genes" = hub_module_gene
)
venn1 <- ggvenn(list1,text_size = 8,set_name_size = 9,show_percentage = TRUE,fill_alpha = 0.7,fill_color = c("#f7a24f","#c6133b"))+ggfonttheme
venn1
pdf(file="02.venn.pdf",width=8,height=6,onefile=F,family = "Times")             
venn1
dev.off()
png(file="02.venn.png",width=8,height=6,units='in',res=600,family = "Times")
venn1
dev.off()
