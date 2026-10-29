rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/23_job_0727")
if (!dir.exists("12_ReactomeGSA")) {dir.create("12_ReactomeGSA")}
setwd("12_ReactomeGSA")
library(ReactomeGSA.data)
library(ReactomeGSA)
pbmc <- readRDS("../RData/R.06_scRNA_done.Rds")
pbmc <- subset(x = pbmc, idents = c( "B cells","Macrophages","Epithelial cells","Neutrophils"))
pbmc$celltype <- factor(
  pbmc$celltype,
  levels = c("B cells","Macrophages","Epithelial cells","Neutrophils")
)
table(pbmc$celltype)
gsva_result <- analyse_sc_clusters(pbmc, verbose = TRUE)
pathway_expression <- pathways(gsva_result)
colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))
pathway_expression[1:10,]
max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
  values <- as.numeric(row[2:length(row)])
  return(data.frame(name = row[1], min = min(values), max = max(values)))
}))
max_difference$diff <- max_difference$max - max_difference$min
# sort based on the difference
max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]
head(max_difference,10)
write.xlsx(max_difference,"02.res.xlsx")
write.xlsx(pathway_expression,"03.pathway.xlsx")
relevant_pathways <- rownames(max_difference[1:10,])
mat1 <- pathway_expression[relevant_pathways,]
rownames(mat1) <- mat1[,1]
mat1 <- mat1[,-1]
p1 <- pheatmap::pheatmap(
  as.matrix(mat1),
  color = colorRampPalette(c("navy", "white", "firebrick3"))(10),
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 15,
  fontsize_col = 15,
  angle_col = 45,
  border_color = NA,
  cellwidth = 30,
  cellheight = 30,
  fontfamily = "Times")
png(file = "01.heatmap.png", height = 8, width = 15, units = 'in',res = 600,family='Times')
p1
dev.off()
pdf(file = "01.heatmap.pdf", height = 8,width = 15,family='Times')
p1
dev.off()
