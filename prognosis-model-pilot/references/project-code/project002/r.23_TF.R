rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/project002")
if (!dir.exists("23_TF")) {dir.create("23_TF")}
setwd("23_TF")
library("openxlsx")
library(tidyr)
library(dplyr)
library(tidyr)
library(networkD3)
library(ggsankey)
ggfonttheme <- theme(text = element_text(family = "Times"))

genes <- "PRR11"
#miRNA--------------------------------------------------------
#ENCORI
ENCORI <-read.xlsx("./miRNA/ENCORI.xlsx")[,1]
#miRDB
miRDB <- read.xlsx("./miRNA/miRDB.xlsx")[,1]
#TargetScan
TargetScan <-read.xlsx("./miRNA/TargetScan.xlsx")[,1]
#取交集
hubmiRNA <- Reduce(intersect,list(ENCORI,miRDB,TargetScan))
list1 <- list(ENCORI,miRDB,TargetScan)
#生成表格
miRNA_dat <- data.frame()
for(i in 1:length(genes)){
  if(length(list1[[i]]) != 0){
    dat1 <- data.frame(gene = genes[i],miRNA = list1[[i]])
    miRNA_dat <- rbind(miRNA_dat,dat1)
  }
}
miRNA_dat
#TF-----------------------------------------------------
TF_raw <- read.table("./TF_Target_information.txt",sep = "\t",header = TRUE)
TF_raw <- TF_raw[grep("other",TF_raw$tissue),]
TF_dat <- TF_raw[,c("TF","target")]
colnames(TF_dat) <- c("Term","Genes")
TF_dat <- TF_dat[TF_dat$Genes %in% genes, ]
#可视化---------------------------------------
library(ggsankey)
# 首先合并数据
merged_data <- miRNA_dat %>%
  inner_join(TF_dat, by = c("gene" = "Genes")) %>%
  dplyr::select(miRNA, gene, Term) %>%
  distinct()
write.xlsx(merged_data,"00.network.xlsx")
# 将数据转换为长格式（桑基图需要的格式）
sankey_long <- merged_data %>%
  make_long(miRNA, gene, Term)  # 指定节点的层次顺序
# 查看转换后的数据
head(sankey_long)
# 根据节点类型着色（miRNA、gene、term用不同颜色）
# 首先为节点添加类型信息
sankey_long_with_type <- sankey_long %>%
  mutate(node_type = case_when(
    node %in% miRNA_dat$miRNA ~ "miRNA",
    node %in% TF_dat$Genes ~ "Gene",
    node %in% TF_dat$Term ~ "Term",
    TRUE ~ "Other"
  ))
# 使用颜色区分不同类型
ggplot(sankey_long_with_type, aes(x = x, 
                                  next_x = next_x, 
                                  node = node, 
                                  next_node = next_node,
                                  fill = node_type,
                                  label = node)) +
  geom_sankey(flow.alpha = 0.6, node.color = "gray30") +
  geom_sankey_label(size = 3, color = "black", fill = "white",family = "Times") +
  scale_fill_manual(values = c("miRNA" = "#FF6B6B", 
                               "Gene" = "#4ECDC4", "Term" = "#556270")) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.margin = margin(20, 20, 20, 20)) +
  labs(title = "miRNA-mRNA-TF Sankey Diagram")+ggfonttheme
ggsave("01.sankey.pdf",device = cairo_pdf,width = 15, height = 15)
ggsave("01.sankey.png",dpi = 300,width = 15, height = 15)
save.image("../RData/23_TF.RData")
