#在KM114上完成，路径<REFERENCE_DATA_ROOT>/REDACTED_PATH
rm(list=ls())
setwd('<REFERENCE_DATA_ROOT>/REDACTED/02_GraphBAN')
suppressPackageStartupMessages({
  library(rstudioapi) 
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(openxlsx)
})
conflicted::conflict_prefer_matching("^filter$|^select$|^arrange$|^mutate$|", "dplyr",quiet = T)
set.seed(123456)

#目标基因，symbol在gene_names列，根据实际修改
gene<-c("IL32","RGS10")
gene<-as.data.frame(gene)
colnames(gene)<-'gene_names'
prot <- readRDS('<REFERENCE_DATA_ROOT>/uniprot_pep/pep.rds') 
#提取需要的蛋白序列，保留基因名和序列
prot <- prot %>% 
  filter(Gene_Name %in% gene$gene_names) %>% 
  select(Gene_Name,Protein)#

#目标药物，至少需要SMILES列，可添加其他列，如id，方便识别
smiles <- read.xlsx("./drug_smiles.xlsx")
smiles<-smiles[,c(2,3)]
colnames(smiles) <- c("Name","SMILES")
#graphban输入必须包含"SMILES"、"Protein"、"Y[01]"列
#用gene和药物id生成蛋白和成分的所有组合
combine <- crossing(prot$Gene_Name,smiles$Name) %>% 
  rename(Gene_Name=1,Name=2)#461404

#将smiles和蛋白序列补充上来
df_test <- Reduce(left_join,list(combine,prot,smiles)) %>%
  mutate(Y=1)#预测未知，所有组合都设定为1
df_test <- na.omit(df_test)
df_test_filtered <- df_test %>%
  dplyr::mutate(len = nchar(SMILES)) %>%
  dplyr::filter(len < 291)
df_test <- df_test_filtered %>%
  dplyr::select(-len)
write.csv(df_test,'df_test_in.csv',row.names = F)


library(dplyr)
A <- read.csv('02.df_test_out_kiba.csv')
A_filter <- A %>% filter(pred > 0.5)
write.csv(A_filter, '02.df_test_out_kiba5.csv', row.names = FALSE)
# A_filter <- A %>% filter(pred > 0.7) 
# write.csv(A_filter, '02.df_test_out_kiba7.csv', row.names = FALSE)
# A_filter <- A %>% filter(pred > 0.6) 
# write.csv(A_filter, '02.df_test_out_kiba6.csv', row.names = FALSE)
# A_filter <- A %>% filter(pred > 0.8) 
# write.csv(A_filter, '02.df_test_out_kiba8.csv', row.names = FALSE)
# A_filter <- A %>% filter(pred > 0.9) 
# write.csv(A_filter, '02.df_test_out_kiba9.csv', row.names = FALSE)
table(A_filter$Gene_Name)

B <- read.csv('02.df_test_out_biosnap.csv')
B_filter <- B %>% filter(pred > 0.5)
write.csv(B_filter, '02.df_test_out_biosnap5.csv', row.names = FALSE)
# B_filter <- B %>% filter(pred > 0.6) 
# write.csv(B_filter, '02.df_test_out_biosnap6.csv', row.names = FALSE)
# B_filter <- B %>% filter(pred > 0.7) 
# write.csv(B_filter, '02.df_test_out_biosnap7.csv', row.names = FALSE)
# B_filter <- B %>% filter(pred > 0.8) 
# write.csv(B_filter, '02.df_test_out_biosnap8.csv', row.names = FALSE)
# B_filter <- B %>% filter(pred > 0.9) 
# write.csv(B_filter, '02.df_test_out_biosnap9.csv', row.names = FALSE)

C <- read.csv('02.df_test_out_bindingdb.csv')
C_filter <- C %>% filter(pred > 0.5)
write.csv(C_filter, '02.df_test_out_bindingdb5.csv', row.names = FALSE)
# C_filter <- C %>% filter(pred > 0.6)
# write.csv(C_filter, '02.df_test_out_bindingdb6.csv', row.names = FALSE)
# C_filter <- C %>% filter(pred > 0.7)
# write.csv(C_filter, '02.df_test_out_bindingdb7.csv', row.names = FALSE)
# C_filter <- C %>% filter(pred > 0.8)
# write.csv(C_filter, '02.df_test_out_bindingdb8.csv', row.names = FALSE)
# C_filter <- C %>% filter(pred > 0.9)
# write.csv(C_filter, '02.df_test_out_bindingdb9.csv', row.names = FALSE)

A_filter<-read.csv('02.df_test_out_kiba5.csv')
B_filter<-read.csv('02.df_test_out_biosnap5.csv')
C_filter<-read.csv('02.df_test_out_bindingdb5.csv')
A_filter$join_key <- paste(as.character(A_filter[,1]), as.character(A_filter[,2]), sep = "_")
B_filter$join_key <- paste(as.character(B_filter[,1]), as.character(B_filter[,2]), sep = "_")
C_filter$join_key <- paste(as.character(C_filter[,1]), as.character(C_filter[,2]), sep = "_")

# 2. 找三个文件共有的联合标识
common_keys <- intersect(C_filter$join_key, B_filter$join_key)
# common_keys <- Reduce(intersect,list(A_filter$join_key, B_filter$join_key, C_filter$join_key))
# 3. 筛选每个文件中属于交集的行
A_common <- A_filter[A_filter$join_key %in% common_keys, ]
B_common <- B_filter[B_filter$join_key %in% common_keys, ]
C_common <- C_filter[C_filter$join_key %in% common_keys, ]
colnames(C_common)
A_common_unique <- A_common[!duplicated(A_common$join_key), ]
B_common_unique <- B_common[!duplicated(B_common$join_key), ]
C_common_unique <- C_filter[!duplicated(C_filter$join_key), ] 

# ========== 步骤2：提取指定列 + 重命名pred列 ==========
# A_common：保留基础列 + pred→kiba_pred
A_sub <- A_common_unique[, c("Gene_Name", "Name", "Protein", "SMILES", "pred")]
colnames(A_sub)[colnames(A_sub) == "pred"] <- "kiba_pred"

# B_common：保留基础列 + pred→biosnap_pred
B_sub <- B_common_unique[, c("Gene_Name", "Name", "Protein", "SMILES", "pred")]
colnames(B_sub)[colnames(B_sub) == "pred"] <- "biosnap_pred"

# C_common：保留基础列 + pred→bindingdb_pred（修正你笔误的B_common）
C_sub <- C_common_unique[, c("Gene_Name", "Name", "Protein", "SMILES", "pred")]
colnames(C_sub)[colnames(C_sub) == "pred"] <- "bindingdb_pred"

# ========== 步骤3：按Gene_Name+Name合并（核心：保留唯一组合） ==========
# 先合并A和B，再合并C
merged_final <- merge(A_sub, B_sub, by = c("Gene_Name", "Name", "Protein", "SMILES"), all = FALSE)
# merged_final <- merge(merged_AB, C_sub, by = c("Gene_Name", "Name", "Protein", "SMILES"), all = FALSE)
table(merged_final$Gene_Name)
write.csv(merged_final, '03.merged_intersect_pred.csv', row.names = FALSE)
