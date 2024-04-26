rm(list=ls())
pkg <- c("Seurat", "parallel", "RColorBrewer", "circlize", "cogena","ComplexHeatmap","dplyr","magrittr", "AUCell", "VennDiagram", "png")
for(i in pkg){
  library(i, character.only = T)
}

# 读取数据 --------------------------------------------------------------------
expr_file <- list.files("0.sc-seq", pattern = ".h5", full.names = T, recursive = T)
# pheno_genes <- list(pheno_genes = pheno_genes)
# ls_sc <- lapply(expr_file, function(x){
#   expr <- Read10X_h5(x)
#   name <- stringr::str_extract(x, "(?<=/)[^/]+(?=_expression.h5)")
#   sce <- CreateSeuratObject(expr, project = name)})
# names(ls_sc) <- stringr::str_extract(expr_file, "(?<=/)[^/]+(?=_expression.h5)")
# marker <- data.table::fread("GeneSet.csv",data.table = F)
# marker <- list(marker)
marker <- gmt2list("GOBP_AUTONOMIC_NERVOUS_SYSTEM_DEVELOPMENT.v2023.1.Hs.gmt")
names(marker) <- "ANSD"
ls_sc <- readRDS("ls_sc.rds")

# 根据MYCRGs的score与表型基因的相关性筛选hub genes ---------------------------------
cor_genes <- lapply(ls_sc, function(x){
  # 数据整理
  tmp_score <- AddModuleScore(x, features = marker, name = names(marker)) %>%
    FetchData(vars = c("ANSD1"))
  expr_sc <- GetAssayData(x, slot="data")
  cor_data <- expr_sc[intersect(marker$ANSD, rownames(expr_sc)), ] %>% as.matrix() %>%
    t() %>% as.data.frame() %>% dplyr::mutate(score = as.numeric(tmp_score$ANSD1))
  # 计算相关性
  cor_pvalue <- function(x, y) {
    cor_result <- cor.test(x, y, method = "spearman")
    c(cor = cor_result$estimate, p.value = cor_result$p.value)
  }
  cor_res <- cor_data %>% 
    dplyr::summarize(across(
      .cols = !matches("^(score)$"),  # 不选择score列
      .fns = ~cor_pvalue(., score)
    ))
  # 得到hub genes
  genes <- colnames(cor_res)[apply(cor_res, 2, function(x) x[1]>0.3 & x[2]<0.05)]
})

# 添加group分组信息 -------------------------------------------------------------
meta_file <- list.files("0.sc-seq", pattern = "CellMetainfo", full.names = T, recursive = T)
ls_meta <- lapply(meta_file, data.table::fread)
names(ls_meta) <- stringr::str_extract(meta_file, "(?<=/)[^/]+(?=_CellMetainfo_table.tsv)")

# for(i in seq_len(length(ls_sc))){
#   if(i == 19){
#     ls_sc[[i]]$cell_type <- ls_meta[[i]]$`Celltype (minor-lineage)`
#   }else{
#     ls_sc[[i]]$cell_type <- ls_meta[[i]]$`Celltype (malignancy)`
#   }
#   ls_sc[[i]]$cell_type <- ifelse(grepl("Malignant", ls_sc[[i]]$cell_type), "Malignant_cells", "other_cells")
#   ls_sc[[i]]$cell_type <- as.factor(ls_sc[[i]]$cell_type)
# }

# 并行计算恶性肿瘤的高表达基因 ----------------------------------------------------------
# # 指定要使用的核心数，这里使用 30 个核心
# cl <- makeCluster(30)
# # 在集群中加载所需的库
# clusterEvalQ(cl, {
#   library(dplyr)
#   library(Seurat)
# })
# # 使用 lapply 进行并行运算
# deg_genes <- parLapply(cl, ls_sc, function(x) {
#   Idents(x) <- "cell_type"
#   markers <- FindMarkers(x, ident.1 = "Malignant_cells", ident.2 = "other_cells", test.use = "wilcox", logfc.threshold = 0.25)
#   DEGs <- rownames(markers)[which(markers$avg_log2FC>0.25 & markers$p_val_adj<1e-05)]
# })
# # 停止并行计算集群
# stopCluster(cl)
deg_genes <- readRDS("deg_genes.rds")

# 相关性基因和差异基因取交集 -----------------------------------------------------------
hub_genes <- lapply(seq_len(length(ls_sc)), function(x){
  intersect(cor_genes[[x]], deg_genes[[x]])
})

hub_genes <- unlist(hub_genes) %>% unique()

# plot --------------------------------------------------------------------
#Venn
# venn.plot <- draw.pairwise.venn(
#   area1 = length(cor_genes[[31]]),
#   area2 = length(deg_genes[[31]]),
#   cross.area = length(intersect(cor_genes[[31]], deg_genes[[31]])),
#   category = c("cor_genes", "deg_genes"),
#   fill = c("#94c6df", "#aacd89"),
#   alpha = 0.5,
#   scaled = FALSE,
#   cex = 0,
#   cat.cex = 0
# )
# 
# png('venn.png', width = 480, height = 480)
# grid.draw(venn.plot)
# dev.off()

#circle plot
sectors <-  stringr::str_extract(names(ls_sc), ".*_[A-Z]{1,}\\d+")
sectors[9] <- "Glioma_GSE131928b"
dataset <- data.frame(dataset = sectors)
rownames(dataset) <- sectors

cancer_type <- list.files("0.sc-seq")
names(cancer_type) <- cancer_type

Set1 <- brewer.pal(9,"Pastel1")
Set2 <- brewer.pal(8,"Pastel2")
col <- c(Set1,Set2)
col_rep <- lapply(seq_along(cancer_type), function(x){
  len <- list.files(paste0("0.sc-seq/", cancer_type[x]), pattern = ".h5") %>% length()
  col <- col[x]
  data.frame(col = col, rep = len)
}) %>% do.call(rbind, .)

col_cancer <- rep(col_rep$col, col_rep$rep)
names(col_cancer) <- sectors

image = 'venn.png'
image = as.raster(readPNG(image))

pdf(file = '1-CircosPlot.pdf',height = 8,width = 12)
circos.clear()
circos.par(gap.degree = 2, cell.padding = c(0, 0, 0, 0),
           track.margin = c(0.01, 0.01))
circos.heatmap(dataset, split = sectors, col = col_cancer, track.height = 0.02, rownames.side = 'outside', rownames.cex = 0.8)

circos.track(ylim = c(0, 1),track.height=0.12, bg.border = "#EEEEEE", panel.fun = function(x, y) {
  circos.raster(image, CELL_META$xcenter, CELL_META$ycenter, 
                width = "0.6cm",
                facing = "inside")
})

circos.track(ylim = c(0, 1), sectors = sectors,
             bg.col = "#8190A5", bg.border = "#EEEEEE" , track.height = 0.25)

circos.trackText(x = rep(0.5, 35), y = rep(0.5, 35),
                 labels = paste0(rep('G',35), 1:35),
                 cex = 0.5, sectors = sectors, col = "white", font = 2, facing = "clockwise",
                 niceFacing=T)

draw.sector(center = c(0, 0), start.degree = 0, end.degree = 360,
            rou1 = 0.15, col = "#4D6381", border = "#EEEEEE")

text(0,0,"ANSDR.Sig",col = 'white',cex = 1,font = 1.5)

cancer <- unique(col_cancer)
names(cancer) <- cancer_type
lgd_cancer = Legend(title = "Cancer", at = names(cancer),
                    legend_gp = gpar(fill = cancer))
draw(lgd_cancer, x = unit(1.2, "snpc"), just = "left")

dev.off()

saveRDS(hub_genes, "hub_genes.rds")
# saveRDS(ls_meta, "dataz/ls_meta.rds")
write.csv(hub_genes,"2-GOKEGG_Input.csv",row.names = F)
