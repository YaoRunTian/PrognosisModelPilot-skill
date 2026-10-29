rm(list = ls()); gc()
setwd("<REFERENCE_DATA_ROOT>/project/19_project_1533")
if (!dir.exists("16_ChromosomeLocation")) {dir.create("16_ChromosomeLocation")}
setwd("16_ChromosomeLocation")
# ========================
# 1. 加载所需R包
# ========================
library(RCircos)
library(biomaRt)

# ========================
# 2. 获取基因染色体位置信息
# ========================
# 目标基因列表
gene_list <-read.table("../08_model/06.lassogeneslist.txt")[,1]

# 通过 biomaRt 获取基因坐标
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gene_info <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band"),
  filters    = "hgnc_symbol",
  values     = gene_list,
  mart       = ensembl
)

# 过滤掉非标准染色体（如 patch/scaffold）
standard_chr <- c(as.character(1:22), "X", "Y")
gene_info <- gene_info[gene_info$chromosome_name %in% standard_chr, ]

# 查看获取的信息
print(gene_info)

# ========================
# 3. 准备 RCircos 输入数据
# ========================
# 构建 RCircos 所需的基因标签数据框
gene_label <- data.frame(
  Chromosome = paste0("chr", gene_info$chromosome_name),
  chromStart = gene_info$start_position,
  chromEnd   = gene_info$end_position,
  Gene       = gene_info$hgnc_symbol,
  stringsAsFactors = FALSE
)

# 按染色体编号排序
chr_order <- paste0("chr", c(1:22, "X", "Y"))
gene_label$Chromosome <- factor(gene_label$Chromosome, levels = chr_order)
gene_label <- gene_label[order(gene_label$Chromosome, gene_label$chromStart), ]

print(gene_label)

# ========================
# 4. 构建连线数据（共表达关系）

# ========================
# 5. 初始化 RCircos 并绑图
# ========================
# 加载人类染色体 ideogram 数据
data(UCSC.HG19.Human.CytoBandIdeogram)

# 设置核心组件参数
RCircos.Set.Core.Components(
  cyto.info  = UCSC.HG19.Human.CytoBandIdeogram,
  chr.exclude = NULL,
  tracks.inside  = 4,
  tracks.outside = 0
)

# 修改参数（可选）
params <- RCircos.Get.Plot.Parameters()
params$track.height <- 0.08
params$text.size    <- 0.6
RCircos.Reset.Plot.Parameters(params)

# ========================
# 6. 开始绑图
# ========================
# 初始化绑图区域
RCircos.Set.Plot.Area()
# 绘制染色体 ideogram
RCircos.Chromosome.Ideogram.Plot()
# ---- Track 1: 基因名称标签 ----
gene_label_plot <- gene_label[, c("Chromosome", "chromStart", "chromEnd", "Gene")]
RCircos.Gene.Connector.Plot(
  genomic.data = gene_label_plot,
  track.num    = 1,
  side         = "in"
)
RCircos.Gene.Name.Plot(
  gene.data  = gene_label_plot,
  name.col   = 4,
  track.num  = 2,
  side       = "in"
)
# 添加标题
title("Chromosomal Distribution of Prognostic Genes", 
      cex.main = 1.2, line = -1)
p <- recordPlot()
pdf("01.Circos.pdf", width = 10, height = 10,family = "Times")
p
dev.off()
png("01.Circos.png", width = 2400, height = 2400, res = 300,family = "Times")
p
dev.off()
# ========================
# 7. 输出基因位置汇总表
# ========================
cat("\n========== 基因染色体位置汇总 ==========\n")
summary_table <- data.frame(
  Gene       = gene_label$Gene,
  Chromosome = gene_label$Chromosome,
  Start      = format(gene_label$chromStart, big.mark = ","),
  End        = format(gene_label$chromEnd, big.mark = ","),
  Length_bp  = format(gene_label$chromEnd - gene_label$chromStart, big.mark = ",")
)
print(summary_table, row.names = FALSE)
write.xlsx(summary_table,"02.gene_loctaion_summary.xlsx")

save.image("../RData/16_Chromosomelocation.RData")
