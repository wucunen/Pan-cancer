rm(list = ls())
library(monocle)

sc_integer <- readRDS("sc_integer.rds")
# 构建CDS数据
## counts
mat_expr <- as(as.matrix(sc_integer@assays$RNA@counts), "sparseMatrix")
## 细胞注释信息
meta.data <- sc_integer@meta.data
## 基因信息
dat_gene <- data.frame(gene_short_name = rownames(sc_integer), row.names = rownames(sc_integer))

# 构建CDS数据
pd <- new("AnnotatedDataFrame", data = meta.data)
fd <- new("AnnotatedDataFrame", data = dat_gene)
cds <- newCellDataSet(mat_expr,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      #count: negbinomial, negbinomial.size()
                      #FPKM和TPM：tobit()
                      #log之后的FPKM和TPM：gaussianff()
                      expressionFamily = negbinomial.size()
                      )
# 数据标准化
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

#filter
cds <- detectGenes(cds, min_expr = 0.1)
head(fData(cds))

# Feature genes --> one out of four
library(Seurat)
## seurat高变基因
sc_integer <- FindVariableFeatures(sc_integer)
express_genes <- VariableFeatures(sc_integer)
# ## seurat差异基因
# sc_degs <- FindAllMarkers(sc_integer)
# express_genes <- subset(sc_degs, p_val_adj < 0.05)$gene
## monocle高变基因
# disp_table <- dispersionTable(cds)
# express_genes <- subset(disp_table, mean_expression >= 0.2 & dispersion_empirical >= 0.5 * dispersion_fit)$gene_id
## dpfeature, 比较不同时间点具有差异的基因
# express_genes <- rownames(subset(fData(cds), num_cells_expressed >= 10))

diff <- differentialGeneTest(cds[express_genes, ], fullModelFormulaStr = "~cell_anno", core = 1)
degs <- subset(diff, qval < 0.01)
degs <- degs[order(degs$qval, decreasing = F),]
order_genes <- rownames(degs)

##选择的用于排序的基因数目一般在2000左右比较合适

# feature 基因 导入cds对象中
cds <- setOrderingFilter(cds, ordering_genes = order_genes)
# 降维
#c("DDRTree", "ICA", "tSNE", "SimplePPT", "L1-graph", "SGL-tree"),
cds <- reduceDimension(cds, max_components = 2, method = "DDRtree")

# 排序
rm(list = lsf.str(envir = .GlobalEnv), envir = .GlobalEnv)
source("./orderCells.R")
library(igraph)
cds <- orderCells(cds)
# 保存cds文件（oderCells输入文件）
saveRDS(cds, file = "cds.rds")
# rm(list = ls())
# # 命令行进入R3.6.3（服务器Terminal）
# /opt/R/3.6.3/bin/R
# getwd()
# library(monocle)
# cds <- readRDS("~/projects/20230828-STAD-ANSD/18-Monocle/cds.rds")
# cds <- orderCells(cds)
# saveRDS(cds, file = "cds_2.rds")

# 绘图
# cds <- readRDS("cds_2.rds")
plot_time <- plot_cell_trajectory(cds, color_by = "Pseudotime", size = 1, show_backbone = T)
plot_time
dev.off()
ggsave("1-Pseudotime.pdf", plot = plot_time, scale = 1, width = 16.6, height = 15, units =c("cm"))

plot_anno <- plot_cell_trajectory(cds, color_by = "cell_anno", size = 1, show_backbone = T)
plot_anno
dev.off()
ggsave("2-CellType.pdf", plot = plot_anno, scale = 1, width = 16.6, height = 15, units =c("cm"))
# plot_cell_trajectory(cds, color_by = "State", size = 1, show_backbone = T) +
#   ggsci::scale_color_nejm()
# plot_cell_trajectory(cds, color_by = "seurat_clusters", size = 1, show_backbone = T)

######END###################
# 选在感兴趣的基因群
library(dplyr)
# select_cells <- subset(pData(cds), State == "1") %>% rownames()

# 可视化基因随细胞状态的表达变化
keygenes <- readRDS("hub_genes_final.rds")
keygenes <- intersect(rownames(cds),keygenes)
# cds_subset <- cds[keygenes, ]
# plot_genes_in_pseudotime(cds_subset, color_by = "State")
# plot_genes_jitter(cds_subset, grouping = "State", color_by = "State")
# plot_genes_violin(cds_subset, grouping = "cell_anno", color_by = "cell_anno")

## 单个基因
# pData(cds)$CCL5 <- log2( exprs(cds)["CCL5", ]+1)
# plot_cell_trajectory(cds, color_by = "CCL5") + scale_color_gsea()

# 差异分析
##根据基因表达量筛选
# expressed_genes=row.names(subset(fData(cds),num_cells_expressed>=10)) #在部分基因里面找
##也可以用feature genes进行差异分析


##根据不同的时间序列进行差异分析
# pseudotime_de <- differentialGeneTest(cds[expressed_genes,],
#                                       fullModelFormulaStr = "~sm.ns(Pseudotime)")
# pseudotime_de <- pseudotime_de[order(pseudotime_de$qval), ]
##根据不同的states进行差异分析
# states_de <- differentialGeneTest(cds[expressed_genes,],
#                                   fullModelFormulaStr = "~State")
# states_de <- states_de[order(states_de$qval), ]

# saveRDS(cds, file = "cds_monocle.rds")
# write.table(pseudotime_de, file = "pseudotime_de.rds", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
# write.table(states_de, file = "states_de.rds", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

heatmap <- plot_pseudotime_heatmap(cds[keygenes, ], num_clusters = 3, show_rownames = T, return_heatmap = T)
heatmap
dev.off()
ggsave("3-Heatmap.pdf", plot = heatmap, scale = 1, width = 16.6, height = 10, units =c("cm"))
