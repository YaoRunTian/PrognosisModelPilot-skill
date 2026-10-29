rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project004")
if (!dir.exists("05_venn")) {dir.create("05_venn")}
setwd("05_venn")
library("ggvenn")
library("dplyr")
library(openxlsx)
ggfonttheme <- theme(text = element_text(family = "Times"))
gene1 <- read.table("../04_DEG/01.DEGsList.txt")[,1]
scTenifoldKnk_res <- readRDS("../03_scTenifoldKnk/scTenifoldKnk_SIRT7.rds")
dat1 <- read.xlsx("../03_scTenifoldKnk/01.res.xlsx")
diff_dat1 <- subset(dat1,p.value < 0.05)
gene2 <- diff_dat1[,"gene"]
hubgenes <- intersect(gene1,gene2)
write.table(hubgenes,"01.hubgenes.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)
gene_list <- list(
  DEGs = gene1,
  Perturb_genes = gene2
)
venn1 <- ggvenn(gene_list, 
                fill_color = c("#0073C2", "#EFC000"),
                stroke_size = 0.5,
                set_name_size = 8,
                text_size = 8,show_percentage = TRUE,digits = 2)
venn1
pdf(file="02.venn.pdf",width=8,height=6,onefile=F,family = "Times")             
venn1
dev.off()
png(file="02.venn.png",width=8,height=6,units='in',res=600,family = "Times")
venn1
dev.off()

save.image("../RData/02_venn.RData")
