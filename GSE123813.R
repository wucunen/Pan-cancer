rm(list=ls())
pkg <- c("Seurat", "ggpubr", "circlize", "ggplot2", "cogena")
for(i in pkg){
  library(i, character.only = T)
}

# 表型基因
markers <- gmt2list("GOBP_AUTONOMIC_NERVOUS_SYSTEM_DEVELOPMENT.v2023.1.Hs.gmt")
# grep("Unfolded", names(markers), ignore.case = T)
# markers <- msidgb[24]
# markers <- data.table::fread("GeneSet.csv",data.table = F)
# markers <- list(markers$Genes) 
names(markers) <- "Autonomic Nervous System Development"
good9color <- c('#6679c9','#aacd89','#f4cf72','#df7971','#94c6df','#60a980','#ef9366','#9d6eba','#e28cd0')

# 载入GSE123813的数据
# exprs_GSE123813 <- readRDS("dataz/exprs_GSE123813.rds")
# pd_GSE15978 <- readRDS("dataz/pd_GSE123813.rds") %>% tibble::column_to_rownames("Cell")
# sc_GSE123813 <- CreateSeuratObject(exprs_GSE123813, min.cells=1)
# sc_GSE123813 <- AddMetaData(sc_GSE123813, metadata = pd_GSE15978)
# sc_GSE123813 <- NormalizeData(sc_GSE123813)
# sc_GSE123813 <- FindVariableFeatures(sc_GSE123813, selection.method = "vst", nfeatures = 2000)
# sc_GSE123813 <- ScaleData(sc_GSE123813)
# sc_GSE123813 <- RunPCA(sc_GSE123813)
# sc_GSE123813 <- RunUMAP(sc_GSE123813, dims = 1:10)
sc_GSE123813<- readRDS("sc_GSE123813.rds")
Idents(sc_GSE123813) <- 'cellType'

#绘制umap图
A <- DimPlot(sc_GSE123813, reduction = "umap", label = F, label.size = 4) + 
  NoLegend() +
  ggtitle('BCC GSE123813') + 
  theme(plot.title = element_text(size = 16)) #element_text(family = 'serif', )
A
ggsave("1-UMAP.pdf", plot = A, width = 16.2, height = 14, units = "cm")

# 计算表型评分，markers可以是一个列表
sc_GSE123813 <- AddModuleScore(sc_GSE123813, features = markers, name = names(markers), seed = 123)
colnames(sc_GSE123813@meta.data)[7:ncol(sc_GSE123813@meta.data)] <- names(markers)
mydata <- FetchData(sc_GSE123813,vars = c("Response", names(markers), "UMAP_1", "UMAP_2", "cellType"))
mydata$Response <- factor(mydata$Response,levels = c('TN', 'NR', 'R'),ordered = T)
mydata$cellType <- factor(mydata$cellType,levels = c('Malignant cells','Stromal cells','Immune cells'),ordered = T)

# 绘制表型评分分布图
# 标记的细胞数据
class_avg <- mydata %>%
  dplyr::group_by(cellType) %>%
  dplyr::summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

B <- ggplot() +
  geom_point(data = mydata, aes(x=UMAP_1, y=UMAP_2, colour= `Autonomic Nervous System Development`), size = 1)  +
  scale_color_gradientn(values = seq(-0,2,0.7),colours = c('#94c6df','red'))+ 
  theme_classic2() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = 'bottom') +
  ggrepel::geom_label_repel(aes(label = cellType),
                            data = class_avg,
                            size = 4,
                            x = class_avg$UMAP_1,
                            y = class_avg$UMAP_2,
                            label.size = 0,
                            segment.color = NA) +
  ggtitle('BCC GSE123813')+
  theme(plot.title = element_text(size = 16, face = "bold"), #family = 'serif', 
        axis.title = element_text(size = 14))
B
ggsave("2-UMAP_Anno.pdf", plot = B, width = 16.2, height = 14, units = "cm")

my_comparisons = list(c('NR','R'))
box_data <- reshape2::melt(mydata[,c("Response", names(markers))], id.vars = "Response")
C <- ggboxplot(box_data, x="Response", y="value", fill  = "Response",
               palette = c('#aacd89','#f4cf72'), width=0.6, add = "none")+
  xlab("")+
  ylab("Autonomic Nervous System Development")+
  rotate_x_text(0)+
  theme(legend.position = 'none')+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
                     label = "p.signif")
C
ggsave("3-Boxplot_Response.pdf", plot = C, width = 16.2, height = 14, units = "cm")

# 绘制表型评分与细胞类型的分组比较图
my_comparisons_2 <- list(c('Malignant cells','Stromal cells'),
                         c('Stromal cells','Immune cells'),
                         c('Malignant cells','Immune cells'))

D <- ggboxplot(mydata,x = 'cellType',y = 'Autonomic Nervous System Development',fill = 'cellType',palette = c('#6679c9','#aacd89','#f4cf72')) + 
  stat_compare_means(comparisons = my_comparisons_2) + 
  xlab('') + 
  ylab('Autonomic Nervous System Development') +
  theme_pubr() +
  theme(legend.position = 'none')
# ,
# text = element_text(family = 'serif')
# ,
# axis.text.x = element_text(angle = 45,hjust = 1)
D
ggsave("4-Boxplot_CellType.pdf", plot = D, width = 16.2, height = 14, units = "cm")

# library(patchwork)
# p_m <- (A | B) / (C | D) + 
#   plot_annotation(tag_levels = "A")
# p_m
# ggsave("Figure3.pdf", plot = p_m, width = 17, height = 20, units = "cm")

# saveRDS(sc_GSE123813, file = "dataz/figure2_input.rds")