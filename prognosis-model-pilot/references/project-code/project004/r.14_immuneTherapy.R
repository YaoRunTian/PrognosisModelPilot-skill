rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project004")
if (!dir.exists("14_immuneTherapy")) {dir.create("14_immuneTherapy")}
setwd("14_immuneTherapy")
library(readxl)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
ggfonttheme <- theme(text = element_text(family = "Times"))
#免疫检查点基因----------------------------------------------------
expr <- readRDS("../00_rawdata/GSE62452_expr.Rds")
#高低风险组数据
risk_group <- read.table("../07_evaluate/04.risk_table.txt",header = TRUE,row.names = 1,sep = "\t")
risk_group1 <- data.frame(risk_group$risk_list_cut,row.names = rownames(risk_group))
check_point <- c("LAG3", "TGFB1", "TGFBR1", "CD160", "CD274", "CD244", "KIR2DL1", "KDR", "ADORA2A", "IDO1", "VTCN1", "HAVCR2", 
                 "IL10", "CD96", "CTLA4", "LGALS9", "TIGIT", "CSF1R", "BTLA", "PDCD1", "PDCD1LG2", "IL10RB", "KIR2DL3")
IMMCG <- check_point
IMMCG_list <- intersect(rownames(expr),IMMCG)
IMMCG_dat <- expr[IMMCG_list,]
immune_dat5 <- melt(t(IMMCG_dat))
colnames(immune_dat5) <- c("Sample","Symbol","value")
dim(immune_dat5)
immune_dat5 <- merge(risk_group1,immune_dat5,by.x = 0,by.y = "Sample")
immune_dat5$riskgroup <- as.character(immune_dat5[,2])
my_comparisons <- list(c("low","high"))
ggplot(immune_dat5,aes(x = Symbol, y = log(value+1),fill = riskgroup))+
  geom_violin(trim = FALSE,position = position_dodge(width = 0.9),scale = "area")+
  geom_boxplot(position = position_dodge(width = 0.9),
               outlier.size = 0,width = 0.4,show.legend = TRUE)+
  scale_fill_manual(values = c("low" = "#B8DDBC","high" = "#F0A780"),
                    labels =  c("Low Risk Group","High Risk Group"))+
  theme_bw()+
  theme(panel.grid = element_line(colour = "white"),plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1,size = 18),
        axis.text.y  = element_text(size = 18),
        axis.title  = element_text(size = 18),
        legend.position = "top",
        legend.key.size = unit(35, "pt"),
        legend.text = element_text(size =10),
        title = element_text(size = 15)
  )+
  ylab("Gene Expression")+
  xlab("immune checkpoint genes")+
  labs(fill = "Risk Group")+
  ggtitle("Boxplot of immune checkpoint genes")+
  stat_compare_means(aes(group=riskgroup),                       #按分组进行统计检验
                     method = "wilcox.test",
                     paired = FALSE,                          
                     symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                        symbols=c("***","**","*","ns")),
                     label = "p.signif",
                     #label.y = location$Proportion+0.02,      #添加显著性符号的位置
                     size=4.5,family = "Times")+ggfonttheme
ggsave("01.checkpoint_boxplot.pdf",device = cairo_pdf,width = 8, height = 6)
ggsave("01.checkpoint_boxplot.png",dpi = 300,width = 8, height = 6)
# ==================== 免疫检查点基因 vs 预后基因 相关性热图（转置版） ====================
library(Hmisc)
library(openxlsx)
# 1. 读取预后基因（已有）
prog_genes <- read.table("../06_model/06.lassogeneslist.txt")[,1]
# 2. 提取检查点基因和预后基因的表达矩阵（样本为行，基因为列）
check_genes <- intersect(rownames(expr), check_point)
prog_genes <- intersect(rownames(expr), prog_genes)
# 表达矩阵转置：行=样本，列=基因
expr_check <- t(expr[check_genes, , drop = FALSE])
expr_prog  <- t(expr[prog_genes, , drop = FALSE])

# 合并为一个数据框（所有样本保持一致）
all_expr <- cbind(expr_check, expr_prog)

# 3. 计算所有基因间的 Spearman 相关系数及 p 值
# 使用 Hmisc::rcorr（返回列表，含 r 和 P）
cor_res <- rcorr(as.matrix(all_expr), type = "spearman")

# 提取相关系数矩阵和 p 值矩阵
cor_mat_all <- cor_res$r
p_mat_all   <- cor_res$P

# 4. 提取检查点基因（行）与预后基因（列）的子矩阵
check_idx <- colnames(expr_check)   # 检查点基因名
prog_idx  <- colnames(expr_prog)    # 预后基因名

cor_matrix <- cor_mat_all[check_idx, prog_idx, drop = FALSE]   # 行=检查点，列=预后
p_matrix   <- p_mat_all[check_idx, prog_idx, drop = FALSE]

# 6. 构建带显著性星号的标签矩阵
sig_label <- matrix("", nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
sig_label[p_matrix < 0.001] <- "***"
sig_label[p_matrix >= 0.001 & p_matrix < 0.01]  <- "**"
sig_label[p_matrix >= 0.01  & p_matrix < 0.05]  <- "*"

# 将相关系数（保留两位小数）与星号拼接
label_matrix <- matrix(paste0(round(cor_matrix, 2), sig_label), 
                       nrow = nrow(cor_matrix), 
                       ncol = ncol(cor_matrix),
                       dimnames = list(check_idx, prog_idx))

# =========================== 关键修改开始 ===========================
# 将矩阵转置，使预后基因变为行，检查点基因变为列
cor_matrix_t <- t(cor_matrix)          # 行=预后，列=检查点
label_matrix_t <- t(label_matrix)      # 标签同步转置
# =========================== 关键修改结束 ===========================

# 7. 绘制热图（使用 pheatmap），此时横轴为检查点，纵轴为预后
library(pheatmap)

# 设置颜色（与您的风格类似：蓝-白-橙）
my_colors <- colorRampPalette(c("#3777ac", "white", "#FC8D62"))(10)

p1 <- pheatmap(cor_matrix_t,                        # 使用转置后的矩阵
               scale = "none",
               cluster_rows = TRUE,                 # 对预后基因聚类（因行已变为预后）
               cluster_cols = TRUE,                 # 对检查点基因聚类（因列已变为检查点）
               border = TRUE,
               breaks = seq(-1, 1, length.out = 10),
               display_numbers = label_matrix_t,    # 使用转置后的标签
               fontsize_number = 12,
               number_color = "black",
               fontsize = 13,
               cellwidth = 40,           
               cellheight = 40,
               color = my_colors,
               angle_col = 45,
               fontfamily = "Times",
               main = "Correlation: Prognostic vs Immune Checkpoint Genes")   # 标题调整为“预后 vs 检查点”

# 8. 保存图片（文件名可改，此处沿用原名称，但建议区分）
ggsave("02.checkpoint_prog_cor_heatmap.pdf", p1, device = cairo_pdf, width = 18, height = 8)
ggsave("02.checkpoint_prog_cor_heatmap.png", p1, dpi = 300, width = 18, height = 8)

# 9. 可选：同时保存转置后的矩阵到 Excel（便于查看）
write.xlsx(as.data.frame(cor_matrix_t), file = "03.prog_checkpoint_cor.xlsx", rowNames = TRUE)
write.xlsx(as.data.frame(t(p_matrix)), file = "04.prog_checkpoint_corp.xlsx", rowNames = TRUE)

# 10. 保存当前工作空间
save.image("../RData/14_immuneTherapy.RData")
