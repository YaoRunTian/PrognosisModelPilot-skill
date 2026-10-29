rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project002")
if (!dir.exists("07_venn")) {dir.create("07_venn")}
setwd("07_venn")
library("ggvenn")
library("dplyr")
library(openxlsx)
ggfonttheme <- theme(text = element_text(family = "Times"))
gene1 <- read.table("../06_DEG/01.DEGsList.txt")[,1]
gene2 <- read.table("../04_scTenifoldKnk/01.diff_perturb_genes.txt")[,1]
gene3 <- read.xlsx("../05_hdWGCNA/06.hub_module_gene_df.xlsx")[,1]
write.table(hubgenes,"01.hubgenes.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)
hubgenes <- Reduce(intersect,list(gene1,gene2,gene3))
gene_list <- list(
  DEGs = gene1,
  Perturb = gene2,
  Hub_Module = gene3
)
venn1 <- ggvenn(gene_list, 
       fill_color = c("#0073C2", "#EFC000", "#CD534C"),
       stroke_size = 0.5,
       set_name_size = 5,
       text_size = 4)
pdf(file="02.venn.pdf",width=8,height=6,onefile=F,family = "Times")             
venn1
dev.off()
png(file="02.venn.png",width=8,height=6,units='in',res=600,family = "Times")
venn1
dev.off()
