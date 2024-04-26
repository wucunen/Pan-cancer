rm (list=ls())
## Package
pkg <- c("ggradar", "ComplexHeatmap", "circlize", "ggplot2", "ggplotify", "patchwork","dplyr")
for(i in pkg){
  library(i, character.only = T)
}

## Input
crispr <- readRDS('crispr.rds')
hub_genes <- readRDS("hub_genes.rds")
ls_sig <- readRDS("ls_sig.rds")
ls_sig$pheno.sig <- hub_genes

## Color Setting
heatmap_col = c("#3D7671", "grey", "#CB9B0C")
grada_col = c("#8DA1CB", "#FD8D62", "#66C3A5")

## Z Scores Calculation
df <- crispr %>% dplyr::mutate(mean = rowMeans(., na.rm = T)) %>% dplyr::arrange(mean) %>% dplyr::select(!mean)
df <- rbind(head(df, 8), tail(df, 8))
column_ha = HeatmapAnnotation(cohort = colnames(df))

## Figure A
pdf(file = '1-ZScores.pdf', height = 4, width = 9)
Heatmap(as.matrix(df), name = "z scores",
             column_title = NULL,
             row_title = NULL,
             cluster_rows = F,
             cluster_columns = F,
             col = colorRamp2(c(-2, 0, 2), heatmap_col),
             show_row_names = T, show_column_names = F,
             rect_gp = gpar(col = "black", lwd = 2),
             width = ncol(df)*unit(5, "mm"), 
             height = nrow(df)*unit(5, "mm"),
             na_col = 'white',
             column_names_side = c('top'),
             row_split = c(rep('a',8),rep('b',8)),
             top_annotation = column_ha)
dev.off()

## Rank of Genes
rank <- data.frame(meanZ = rowMeans(crispr, na.rm = T), row.names = rownames(crispr), gene = rownames(crispr))
rank <- rank %>% dplyr::arrange(meanZ) %>% dplyr::mutate(order = order(meanZ))
top_percent <- function(pheno){
  top_1 <- rank$gene[1:round(nrow(rank) * 0.02)]
  top_2 <- rank$gene[1:round(nrow(rank) * 0.04)]
  top_3 <- rank$gene[1:round(nrow(rank) * 0.06)]
  rank_1 <- sum(top_1 %in% pheno)/length(top_1)
  rank_2 <- sum(top_2 %in% pheno)/length(top_2)
  rank_3 <- sum(top_3 %in% pheno)/length(top_3)
  rank_pheno <- data.frame(top1 = rank_1, top2 = rank_2, top3 = rank_3)
  return(rank_pheno)
  }

dat_rank <- lapply(ls_sig, top_percent) %>% do.call(rbind, .) %>% t() %>% as.data.frame() %>% 
  cbind(sig = as.character(paste(c(2,4,6),'%','Top-Ranked Genes')), .)
colnames(dat_rank)
colnames(dat_rank) <- c("sig","ImmmunCells.Sig","TRS.Sig","IMS.Sig","LRRC15.CAF.Sig","CRMA.Sig","ANSDR.Sig","Cytotoxic.Sig","PDL1.Sig","T.cell.inflamed.Sig","INFG.Sig","NLRP3.Sig","IPRES.Sig")
dat_rank_filter <- dat_rank[,1:7]

## Figure B
B <- ggradar(dat_rank_filter,
             values.radar = c("0","1%"),
             grid.min = 0,
             grid.max = 0.01,
             group.colours = grada_col,
             legend.position = "bottom",
             legend.text.size = 10,
) 
B <- B + theme(legend.direction = "vertical",
               legend.margin = margin(t = -30))
B
ggsave(filename = '2-RadarPlot.pdf', plot = B, height = 24, width = 24.3, unit = "cm")

sig_ansd <- ls_sig[[6]]
top_10 <- rank$gene[1:round(nrow(rank) * 0.20)]
ansd_top_10 <- crispr[top_10[top_10 %in% sig_ansd], ] %>% dplyr::mutate('mean z score' = rowMeans(., na.rm = T)) %>% t()
row_ep_top_10 = rowAnnotation(cohort = rownames(ansd_top_10))

C <- Heatmap(ansd_top_10, name = "z scores",
             row_title = NULL,
             cluster_rows = F,
             cluster_columns = F,
             col = colorRamp2(c(-2, 0, 2), heatmap_col),
             show_row_names = T, show_column_names = T,
             width = ncol(ansd_top_10)*unit(5, "mm"), 
             height = nrow(ansd_top_10)*unit(5, "mm"),
             rect_gp = gpar(col = "black", lwd = 1),
             na_col = 'white',
             row_split = c(rep('a',17),rep('b',1)),
             column_names_side = c('top'),
             #row_labels = c(rep("",nrow(ep_top_3)-1), "mean z score"),
             heatmap_legend_param = list(legend_direction = "horizontal", title_position = "leftcenter"))
C
dev.off()
pdf(file = '3-Heatmap.pdf', height = 8, width = 8.1)
draw(C, heatmap_legend_side = "bottom")
dev.off()

