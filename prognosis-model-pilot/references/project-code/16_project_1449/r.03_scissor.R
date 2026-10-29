rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/16_project_1449")
if (!dir.exists("03_scissor")) {dir.create("03_scissor")}
setwd("03_scissor")
library(Scissor)
library(Seurat)
library(ggplot2)
library(dplyr)
library(openxlsx)
ggfonttheme <- theme(text = element_text(family = "Times"))
#V5 RunScissor函数-----------------------------
RUNScissor <- function (bulk_dataset, sc_dataset, phenotype, tag = NULL, alpha = NULL, 
                        cutoff = 0.2, family = c("gaussian", "binomial", "cox"), 
                        Save_file = "Scissor_inputs.RData", Load_file = NULL) 
{
  library(Seurat)
  library(Matrix)
  library(preprocessCore)
  
  # 确保 phenotype 是向量
  phenotype <- as.numeric(phenotype)
  
  if (is.null(Load_file)) {
    common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
    if (length(common) == 0) {
      stop("There is no common genes between the given single-cell and bulk samples.")
    }
    
    if (class(sc_dataset) == "Seurat") {
      sc_exprs <- as.matrix(sc_dataset@assays$RNA@layers$data)
      rownames(sc_exprs) <- rownames(sc_dataset)
      colnames(sc_exprs) <- colnames(sc_dataset)
      network <- as.matrix(sc_dataset@graphs$RNA_snn)
    } else {
      sc_exprs <- as.matrix(sc_dataset)
      Seurat_tmp <- CreateSeuratObject(sc_dataset)
      Seurat_tmp <- FindVariableFeatures(Seurat_tmp, selection.method = "vst", verbose = FALSE)
      Seurat_tmp <- ScaleData(Seurat_tmp, verbose = FALSE)
      Seurat_tmp <- RunPCA(Seurat_tmp, features = VariableFeatures(Seurat_tmp), verbose = FALSE)
      Seurat_tmp <- FindNeighbors(Seurat_tmp, dims = 1:10, verbose = FALSE)
      network <- as.matrix(Seurat_tmp@graphs$RNA_snn)
    }
    
    diag(network) <- 0
    network[which(network != 0)] <- 1
    dataset0 <- cbind(bulk_dataset[common, ], sc_exprs[common, ])
    dataset0 <- as.matrix(dataset0)
    dataset1 <- normalize.quantiles(dataset0)
    rownames(dataset1) <- rownames(dataset0)
    colnames(dataset1) <- colnames(dataset0)
    Expression_bulk <- dataset1[, 1:ncol(bulk_dataset)]
    Expression_cell <- dataset1[, (ncol(bulk_dataset) + 1):ncol(dataset1)]
    X <- cor(Expression_bulk, Expression_cell)
    
    quality_check <- quantile(X)
    print("|**************************************************|")
    print("Performing quality-check for the correlations")
    print("The five-number summary of correlations:")
    print(quality_check)
    print("|**************************************************|")
    
    if (quality_check[3] < 0.01) {
      warning("The median correlation between the single-cell and bulk samples is relatively low.")
    }
    
    if (family == "binomial") {
      Y <- as.numeric(phenotype)
      z <- table(Y)
      if (length(z) != length(tag)) {
        stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
      } else {
        print(sprintf("Current phenotype contains %d %s and %d %s samples.", z[1], tag[1], z[2], tag[2]))
        print("Perform logistic regression on the given phenotypes:")
      }
    }
    
    if (family == "gaussian") {
      Y <- as.numeric(phenotype)
      z <- table(Y)
      if (length(z) != length(tag)) {
        stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
      } else {
        tmp <- paste(z, tag)
        print(paste0("Current phenotype contains ", paste(tmp[1:(length(z) - 1)], collapse = ", "), ", and ", tmp[length(z)], " samples."))
        print("Perform linear regression on the given phenotypes:")
      }
    }
    
    if (family == "cox") {
      Y <- as.matrix(phenotype)
      if (ncol(Y) != 2) {
        stop("The size of survival data is wrong. Please check Scissor inputs and selected regression type.")
      } else {
        print("Perform cox regression on the given clinical outcomes:")
      }
    }
    
    save(X, Y, network, Expression_bulk, Expression_cell, file = Save_file)
  } else {
    load(Load_file)
  }
  
  if (is.null(alpha)) {
    alpha <- c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  }
  
  for (i in 1:length(alpha)) {
    set.seed(123)
    fit0 <- APML1(X, Y, family = family, penalty = "Net", alpha = alpha[i], Omega = network, nlambda = 100, nfolds = min(10, nrow(X)))
    fit1 <- APML1(X, Y, family = family, penalty = "Net", alpha = alpha[i], Omega = network, lambda = fit0$lambda.min)
    
    if (family == "binomial") {
      Coefs <- as.numeric(fit1$Beta[2:(ncol(X) + 1)])
    } else {
      Coefs <- as.numeric(fit1$Beta)
    }
    
    Cell1 <- colnames(X)[which(Coefs > 0)]
    Cell2 <- colnames(X)[which(Coefs < 0)]
    percentage <- (length(Cell1) + length(Cell2)) / ncol(X)
    print(sprintf("alpha = %s", alpha[i]))
    print(sprintf("Scissor identified %d Scissor+ cells and %d Scissor- cells.", length(Cell1), length(Cell2)))
    print(sprintf("The percentage of selected cell is: %s%%", formatC(percentage * 100, format = "f", digits = 3)))
    
    if (percentage < cutoff) {
      break
    }
    cat("\n")
  }
  
  print("|**************************************************|")
  return(list(para = list(alpha = alpha[i], lambda = fit0$lambda.min, family = family), Coefs = Coefs, Scissor_pos = Cell1, Scissor_neg = Cell2))
}
#load the data--------------------------------
scRNA <- readRDS("../RData/pbmc_all_with_subcelltype.Rds")
exp <- readRDS("../00_rawdata/GSE108214_expr.Rds")
groupInf <- readRDS("../00_rawdata/GSE108214_group.Rds")
sc_dataset <- scRNA
bulk_dataset <- exp
bulk_dataset[is.na(bulk_dataset)] <- 0
phenotype <- groupInf$group
#0:res,1:sen
tag <- c("resistant","sensitive")
#start run
infos1 <- RUNScissor(bulk_dataset, sc_dataset, phenotype, 
                     tag = tag,
                     alpha = 0.5,
                     cutoff = 0.2, #the number of the Scissor selected cells should not exceed 20% of total cells in the single-cell data
                     family = "binomial", 
                     Save_file = './Scissor_res.RData')

Scissor_select <- rep("Background cells", ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos1$Scissor_pos] <- "Scissor+"
Scissor_select[infos1$Scissor_neg] <- "Scissor-"
sc_dataset <- AddMetaData(sc_dataset, 
                          metadata = Scissor_select, 
                          col.name = "scissor")
UMAP_scissor <- DimPlot(sc_dataset, reduction = 'umap', 
                        group.by = 'scissor',
                        cols = c('grey','royalblue','indianred1'), 
                        pt.size = 0.9, order = c("Scissor+","Scissor-"))+theme(legend.position = c(0.8,0.9))+ggtitle("")+ggfonttheme
UMAP_scissor
ggsave("01.scissor_res.pdf",UMAP_scissor,device = cairo_pdf,width = 8, height = 8)
ggsave("01.scissor_res.png",UMAP_scissor,dpi = 300,width = 8, height = 8)

sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
saveRDS(sc_dataset, file = "../RData/R.07_scRNA_Scissor.RData")

#绘制柱状图------------------------------------------------
table(sc_dataset@meta.data$subcelltype)
table(sc_dataset@meta.data$scissor)

meta_data <- sc_dataset@meta.data
scissor_counts <- table(meta_data$subcelltype, meta_data$scissor)
# 转换为数据框
scissor_df <- as.data.frame(scissor_counts)
colnames(scissor_df) <- c("subcelltype", "scissor", "count")
# 计算百分比
scissor_percent <- scissor_df %>%
  group_by(scissor) %>%
  mutate(percent = count / sum(count) * 100) %>%
  ungroup()
# 设置颜色方案
n_celltypes <- length(unique(scissor_percent$subcelltype))
mycol <- colorRampPalette(brewer.pal(8, "Pastel1"))(n_celltypes)
# 绘制堆叠百分比柱状图
ggplot(scissor_percent, aes(x = scissor, y = percent, fill = subcelltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = mycol) +
  labs(
    title = "Subcelltype Distribution across Scissor Classifications",
    x = "Scissor Classification",
    y = "Percentage",
    fill = "Subcelltype"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 15),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
        )+
  scale_y_continuous(expand = c(0,0))+
  ggfonttheme
ggsave("02.barplot.pdf",device = cairo_pdf,width = 12, height = 8)
ggsave("02.barplot.png",dpi = 300,width = 12, height = 8)
# 创建宽格式表格显示具体数量
scissor_wide <- as.data.frame.matrix(table(meta_data$subcelltype, meta_data$scissor))
colnames(scissor_wide) <- c("Scissor_0", "Scissor_minus", "Scissor_plus")
# 添加总数和百分比
scissor_wide$Total <- rowSums(scissor_wide)
scissor_wide$Percent_Scissor_plus <- round(scissor_wide$Scissor_plus / sum(scissor_wide$Scissor_plus) * 100, 2)
scissor_wide <- scissor_wide[order(scissor_wide$Scissor_plus, decreasing = TRUE), ]
# 显示表格
print(scissor_wide)
scissor_wide <- as.data.frame(scissor_wide)
scissor_wide$celltype <- rownames(scissor_wide)
write.xlsx(scissor_wide,"03.scissor_proportion.xlsx")
#Submucosal_gland_cell和Secretory_cell为关键细胞亚型
save.image("../RData/03_Scissor.RData")
