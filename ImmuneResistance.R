rm(list=ls())
pkg <- c("GSVA", "ggplot2", "ggnewscale", "MCPcounter", "magrittr", "msigdbr", "BiocParallel")
for (i in pkg) {
  library(i, character.only = T)
}

# 读取数据
hub_genes <- readRDS("hub_genes.rds")
meta_TCGA <- readRDS("meta_TCGA.rds")
exprs_TCGA <- readRDS("exprs_TCGA.rds")
immune_genes <- readRDS("immune_genes.rds")
exprs_TCGA_tumor <- exprs_TCGA[, stringr::str_sub(colnames(exprs_TCGA), 14, 15) != 11]

#计算表型评分
hub_genes <- data.frame("gene" = hub_genes)
ansd_score <- gsva(expr = exprs_TCGA_tumor, gset.idx.list = hub_genes, method="ssgsea")

# immune_genes_cor --------------------------------------------------------
#数据整理
inter_immune <- intersect(immune_genes$HGNC.Symbol, rownames(exprs_TCGA_tumor))
cor_data <- exprs_TCGA_tumor[inter_immune, ] %>% t() %>% as.data.frame() %>% 
  dplyr::mutate(score = as.numeric(ansd_score)) %>% tibble::rownames_to_column("sample") %>%
  merge(meta_TCGA[, c(1,3)]) %>% tibble::column_to_rownames("sample")
colnames(cor_data)[77] <- "cancer_type"

#计算表型评分与免疫相关基因的的相关性值
plot_data <- lapply(unique(cor_data$cancer_type), function(x){
  dat_cor <- subset(cor_data, cancer_type == x)
  cor_res <- apply(dat_cor[, 1:75], 2, function(x){cor.test(as.numeric(x), as.numeric(dat_cor$score))})
  ls_res <- lapply(cor_res, function(x){data.frame(cor = signif(x$estimate, 3), pvalue = signif(x$p.value, 3))})
  dat_res <- do.call(rbind, ls_res) %>% dplyr::mutate(cancer_type = x) %>% tibble::rownames_to_column("genes") 
  return(dat_res)
})

#数据整理
plot_data <- do.call(rbind, plot_data)
rownames(immune_genes) <- immune_genes$HGNC.Symbol
immune_genes <- immune_genes[plot_data$genes, ]
plot_data$Functions <- immune_genes$Immune.Checkpoint
plot_data$Pathways <- immune_genes$Super.Category
plot_data$Functions[which(plot_data$Functions %in% c("Stumulatory", "Stimulaotry"))] = "Stimulatory"
plot_data$Functions[which(plot_data$Functions == "")] = "unknow"
plot_data$num <- 1:75
plot_data <- dplyr::arrange(plot_data, Pathways, genes)

#标签设置
gene_text <- subset(plot_data, cancer_type == "GBM")
rownames(gene_text) <- NULL
gene_text$ang <- seq(from = (360/nrow(gene_text)) / 1.5,
                     to = (1.5* (360/nrow(gene_text))) - 360,
                     length.out = nrow(gene_text)) + 80
gene_text$hjust <- 0
gene_text$hjust[which(gene_text$ang < -90)] <- 1
gene_text$ang[which(gene_text$ang < -90)] <- (180+gene_text$ang)[which(gene_text$ang < -90)]
g_anno1<- plot_data %>% dplyr::filter(cancer_type == "UVM") %>% dplyr::select(Functions)
g_anno2<- plot_data %>% dplyr::filter(cancer_type == "UVM") %>% dplyr::select(Pathways)

g_anno1$Functions <- factor(g_anno1$Functions)
g_anno2$Pathways <- factor(g_anno2$Pathways)

#相关性图绘制
p <- ggplot() +
  # 添加基因注释1
  geom_bar(data = g_anno1, stat = 'identity',
           aes(y = 1,x = -0.5,fill = Functions),
           width = 1,
           color = NA) +
  # 添加基因注释2
  geom_bar(data = g_anno2,stat = 'identity',
           aes(y = 1,x = 32.5,fill = Pathways),
           width = 1,
           color = NA)+
  scale_fill_manual(values = c("#88cbc0","#f6f5b3","#bbb7d6","#ee7d6f","#7fabca","#f4af64","#add26b","#bfbebe","#dc9a17","#313232"))
A <- p+new_scale("fill") +
  geom_tile(data = plot_data, aes(x = cancer_type, y = num, fill = cor),color = 'white') +
  scale_fill_gradient2(midpoint = 0,
                       low = '#0000ff',
                       mid = "white",
                       high = '#ff0303') +
  scale_y_discrete(expand = expansion(mult = c(0.05,0))) +
  scale_x_discrete(expand = expansion(mult = c(0.5,0.5)))+
  coord_polar(theta = 'y') +
  theme_void() +
  geom_text(data = gene_text,
            aes(y = as.numeric(rownames(gene_text)),
                x = 33,
                label = genes, angle = ang, hjust = hjust),
            size = 4)
A
ggsave("1-Circos.pdf",plot = A, height = 16, width = 20, units = "cm")

# MCPcounter --------------------------------------------------------------
probesets <- read.table('probesets.txt',
                        sep="\t",
                        stringsAsFactors=FALSE,
                        colClasses="character")
genes <- read.table('probegenes.txt',
                    sep="\t",
                    stringsAsFactors=FALSE,
                    colClasses="character",
                    header=TRUE,
                    check.names=FALSE)
mcp_counter <- MCPcounter.estimate(exprs_TCGA_tumor,featuresType = 'HUGO_symbols', probesets = probesets, genes = genes)

#数据整理
mcp_counter <- mcp_counter %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("sample")
mcp_counter <- merge(mcp_counter, meta_TCGA[, c(1,3)])
colnames(mcp_counter)[12] <- "cancer_type"
ansd_score  <- ansd_score  %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("sample")
mcp_counter <- merge(mcp_counter, ansd_score ) %>% tibble::column_to_rownames("sample")
colnames(mcp_counter)[12] <- "score"

#计算相关性
mcp_cor <- mcp_counter %>%
  dplyr::group_by(cancer_type) %>%
  dplyr::summarize(across(
    .cols = !matches("^(cancer_type|score)$"),  # 不选择 cancer_type 和 score 列
    .fns = ~cor(., score, method = "spearman")
  )) %>%
  tibble::column_to_rownames("cancer_type")

mcp_pvlaue <- mcp_counter %>%
  dplyr::group_by(cancer_type) %>%
  dplyr::summarize(across(
    .cols = !matches("^(cancer_type|score)$"),  # 不选择 cancer_type 和 score 列
    .fns = ~cor.test(., score, method = "spearman")$p.value
  )) %>%
  tibble::column_to_rownames("cancer_type")

#整理数据
cor_data_mcp <- reshape2::melt(mcp_cor)
cor_data_mcp$pvalue <- reshape2::melt(mcp_pvlaue)$value
cor_data_mcp$cancer_type <- rownames(mcp_cor)
cor_data_mcp$log10p <- -log10(cor_data_mcp$pvalue)

#绘图
B <- ggplot(cor_data_mcp, aes(x = variable , y = cancer_type)) + 
  geom_point(aes(size = log10p, color = value)) +
  labs(x = "Cell", y = "Cancer Type", color = "Correlation", size = "-log10(p value)")+
  scale_colour_gradientn(colours = c("blue", 'white', "red"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 0, vjust = 0))+
  scale_x_discrete(position = "top") +
  labs(x='',y='')
B
ggsave("2-Heatmap.pdf",plot = B, height = 12, width = 12, units = "cm")

# Hallmark ----------------------------------------------------------------
## 数据准备
df = msigdbr(species = "Homo sapiens",category = "H")
geneset = split(x = df$gene_symbol, f = df$gs_name)

## 进行GSVA分析
gsva_h <- gsva(exprs_TCGA_tumor,
               geneset, 
               mx.diff=FALSE, verbose=FALSE, 
               parallel.sz = 30,
               min.sz > 1,
               BPPARAM = MulticoreParam(30))

#数据整理
gsva_h <- gsva_h %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("sample")
gsva_h <- merge(gsva_h, meta_TCGA[, c(1,3)]) %>% merge(ansd_score) %>% 
  tibble::column_to_rownames("sample")
colnames(gsva_h)[51:52] <- c("cancer_type", "score")

#计算相关性
gsva_h_cor <- gsva_h %>%
  dplyr::group_by(cancer_type) %>%
  dplyr::summarize(across(
    .cols = !matches("^(cancer_type|score)$"),  # 不选择 cancer_type 和 score 列
    .fns = ~cor(., score, method = "spearman")
  )) %>%
  tibble::column_to_rownames("cancer_type")

gsva_h_pvlaue <- gsva_h %>%
  dplyr::group_by(cancer_type) %>%
  dplyr::summarize(across(
    .cols = !matches("^(cancer_type|score)$"),  # 不选择 cancer_type 和 score 列
    .fns = ~cor.test(., score, method = "spearman")$p.value
  )) %>%
  tibble::column_to_rownames("cancer_type")

#整理数据
cor_data_gsva <-reshape2::melt(gsva_h_cor,)
cor_data_gsva$pvalue <- reshape2::melt(gsva_h_pvlaue)$value
cor_data_gsva$cancer_type <- rownames(gsva_h_cor)
cor_data_gsva$log10p <- -log10(cor_data_gsva$pvalue)
cor_data_gsva$variable <- gsub("HALLMARK_", "", cor_data_gsva$variable)
cor_data_gsva$variable <- gsub("_", " ", cor_data_gsva$variable)
#绘图
C <- ggplot(cor_data_gsva, aes(x = variable , y = cancer_type)) + 
  geom_point(aes(size = log10p, color = value)) +
  labs(x = "Cell", y = "Cancer Type", color = "Correlation", size = "-log10(p value)")+
  scale_colour_gradientn(colours = c("blue", 'white', "red"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 0, vjust = 0))+
  scale_x_discrete(position = "top") +
  labs(x='',y='')
C
ggsave("3-Heatmap_GSVA.pdf",plot = C, height = 18, width = 26, units = "cm")

# res <- list(data = list(myc_score = myc_score,
#                         plot_data = plot_data,
#                         mcp_counter = mcp_counter,
#                         gsva_h = gsva_h),
#             plot = list(immune_genes = A,
#                         mcp = B,
#                         gsva = C))
# 
# source("CorImmune_function.R")
# 
# res <- figure5(hub_genes)
# 
# ggsave("figure5/cor_immue.pdf", plot = res$plot$p1, height = 16, width = 20, units = "cm")
# ggsave("figure5/cor_mcpcounter.pdf", plot = res$plot$p2, height = 16, width = 12, units = "cm")
# ggsave("figure5/cor_gsva.pdf", plot = res$plot$p3, height = 18, width = 26, units = "cm")