rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project003")
if (!dir.exists("02_venn")) {dir.create("02_venn")}
setwd("02_venn")
library("ggvenn")
library("dplyr")
library(openxlsx)
ggfonttheme <- theme(text = element_text(family = "Times"))
gene1 <- read.table("../01_DEG/01.DEGsList.txt")[,1]
gene2 <- read.table("../02_venn/PRRGs.txt")[,1]
hubgenes <- intersect(gene1,gene2)
write.table(hubgenes,"01.hubgenes.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)
gene_list <- list(
  DEGs = gene1,
  PRRGs = gene2
)
venn1 <- ggvenn(gene_list, 
                fill_color = c("#0073C2", "#EFC000"),
                stroke_size = 0.5,
                set_name_size = 5,
                text_size = 4)
venn1
pdf(file="02.venn.pdf",width=8,height=6,onefile=F,family = "Times")             
venn1
dev.off()
png(file="02.venn.png",width=8,height=6,units='in',res=600,family = "Times")
venn1
dev.off()

save.image("../RData/02_venn.RData")
