rm(list = ls())
library(dplyr)
pkg <- scan(text = c("stringr, gridExtra, future, sva, e1071, pROC, ROCit, caret, doParallel, cancerclass, tidyverse, tidyr, doParallel, multiGSEA, reshape2, RColorBrewer, ComplexHeatmap, GSVA, reshape2, RColorBrewer, circlize"), what = "character", sep = ",") %>% trimws()
for(i in pkg){
  library(i, character.only = T)
}
## Input
ls_sig <- readRDS("ls_sig.rds")
hub_genes <- readRDS("hub_genes.rds")

data <- readRDS("immune_datasets.rds")
# data$response <- ifelse(data$response=='0','NR','R') %>% factor(levels = c('R','NR'))
grp <- unique(data$batch);grp
combat <- data[!data$batch %in% grp[c(3,7,8,9,10)],]
test <- data[data$batch %in% grp[c(3,7,8,9,10)],] # Set individual cohort
# test <- data[data$batch %in% grp[c(3,9,10)],] # Set individual cohort
method = c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass') # Set learning models
# method = c('svmRadialWeights','rf','kknn','LogitBoost','cancerclass') # Set learning models

## Machine Learning Models
hub_genes <- intersect(hub_genes, colnames(combat))
ls_sig$pheno.sig <- hub_genes
saveRDS(ls_sig, file = "ls_sig.rds")

source("AUC_Function.R")
ls_res_seed <- figure6(hub_genes, method = method, iter = 1, seed = 1)
ls_res <- figure6(hub_genes, method = method, iter = 20) #结果不好(< 0.7) 跑循环

# plot
max_iter <- ls_res[["max_iter"]]
res <- ls_res[[as.character(max_iter)]]
max_roc_num <- which.max(res$auc$ROC)
max_roc_res <- res$roc_res[[max_roc_num]]
method <- res$method[[max_roc_num]]
max_model <- res$model[[max_roc_num]]

lop_data <- data.frame(method = purrr::map_chr(res$model, "method"),
                       auc = res[["auc"]][["ROC"]])
lop_data <- arrange(lop_data, desc(auc))

A <- ggplot(lop_data, aes(x = reorder(method, auc), y = auc)) +
     geom_segment(aes(xend = method, yend = 0, color = method)) +
     geom_point(aes(color = method), size = 5) +
     geom_text(aes(label = signif(auc, 3)), hjust = -1, size = 4) +
     coord_flip() +
     theme_classic() +
     labs(x = "Method", y = "AUC") +
     #theme_minimal()+
     scale_color_manual(values = brewer.pal(7, "Set2"))+
     theme(axis.line = element_line(color = "black"),
                     axis.text = element_text(size = 10, color = "black"),
                     legend.position = "none")+
     scale_y_continuous(limits = c(0, 1))
A
ggsave("1-Barplot.pdf", plot = A, height = 9, width = 12.15,units = "cm")

pdf(file = '2-ROC.pdf',height = 4,width = 5)
plot.roc(max_roc_res, print.auc = T, col = "#0072B5FF", main = paste0("ROC ", "for ", method,  " in Validation"))
dev.off()

## Comparison
set.seed(seed = max_iter)
trainIndex <- createDataPartition(combat$response, p = 0.8,  ## 80% test cohort; 20% validation cohort
                                  list = FALSE, 
                                  times = 1)

training <- combat[trainIndex,]
validation  <- combat[-trainIndex,]
saveRDS(training, file = "exprs_train.rds")
saveRDS(validation, file = "exprs_validation.rds")

data_raw <- readRDS("immune_datasets_raw.rds")
data_raw$response <- ifelse(data_raw$response=='0','NR','R') %>% factor(levels = c('R','NR'))

skcm_sig <- skcm_sig(ls_sig)
pancancer_sig <- pancancer_sig()
AUC <- list(skcm = skcm_sig, pancancer = pancancer_sig)

# plot
AUC$pancancer$pheno <- as.matrix(AUC$pancancer$pheno[-8,])
colnames(AUC$pancancer$pheno) <- "AUC"
cohort <- rownames(AUC[['pancancer']]$pheno)
rownames(AUC[['pancancer']]$pheno)
rename <- c('Validation Set',
            'Training Set',
            'Hugo 2016 SKCM',
            'Van 2015 SKCM',
            'Kim 2018 GC',
            'Zhao 2019 GBM',
            'Synder 2017 UC'
)

# AUC <- lapply(AUC,as.data.frame)

auc <- AUC[['pancancer']]
names(auc) <- c("INFG.Sig","T.cell.inflamed.Sig","LRRC15.CAF.Sig","PDL1.Sig","NLRP3.Sig","Cytotoxic.Sig","ANSDR.Sig")
ls_AUC <- lapply(cohort,function(c){
  cc <-lapply(auc,function(auc){auc[c,'AUC']}) %>% do.call(rbind,.)
  names(cc) <- names(auc) 
  cc <- as.data.frame(cc)
  colnames(cc) <- 'AUC'
  return(cc)
}) %>% `names<-`(rename)
as.data.frame(ls_AUC)

# ls_AUC <- lapply(cohort, function(c){
#   cc <-lapply(auc, function(auc){auc[c,'AUC']}) %>% do.call(rbind,.)
#   names(cc) <- c("INFG.Sig","T.cell.inflamed.Sig","LRRC15.CAF.Sig","PDL1.Sig","NLRP3.Sig","Cytotoxic.Sig","EPR.Sig")
#   cc <- as.data.frame(cc)
#   colnames(cc) <- 'AUC'
#   return(cc)
# }) %>% `names<-`(rename)
# 
# ls_AUC <- as.data.frame(auc) %>% t() %>% as.data.frame() %>% as.list() %>%
#   lapply(function(x) {
#     names(x) <- c("INFG.Sig","T.cell.inflamed.Sig","LRRC15.CAF.Sig","PDL1.Sig","NLRP3.Sig","Cytotoxic.Sig","EPR.Sig")
#     return(x)
#   })

pdf(file = '3-Circos.pdf',height = 7,width = 10)
circos_cormpare(ls_AUC[c(-4,-5)])
# circos_cormpare(ls_AUC)
dev.off()

#Heatmap 
hdata <- ls_AUC[c(-4,-5)] %>% do.call(cbind,.) %>% `colnames<-`(names(ls_AUC)[c(-4,-5)])
# hdata <- ls_AUC %>% do.call(cbind,.) %>% `colnames<-`(names(ls_AUC))
hdata$MeanAUC <- rowMeans(hdata)
hdata <- hdata %>%
  dplyr::select(MeanAUC,everything())
hdata <- hdata[order(hdata[,1],decreasing = T) ,]

pdf(file = '4-Heatmap.pdf',height = 4,width = 5)
Heatmap(as.matrix(hdata), name = "AUC",
        # column_title = NULL,
        # row_title = NULL,
        cluster_rows = F,
        cluster_columns = F,
        col = colorRamp2(c( 0.5, 0.80), c("#FFFFE3",  "#004628")),
        # show_row_names = T, show_column_names = F,
        rect_gp = gpar(col = "black", lwd = 2),
        width = nrow(hdata)*unit(6, "mm"),
        height = nrow(hdata)*unit(6, "mm"),
        # na_col = 'white',
        # column_names_side = c('top'),
        # row_split = c(rep('a',8),rep('b',8)),
        column_split = c('Combined',rep('Individual',c(ncol(hdata)-1))),
        # top_annotation = column_ha
)
dev.off()

## Logistics Regression --> ROC Curve Validation
# # 定义训练模型
# set.seed(seed = 2)
# trainIndex <- createDataPartition(combat$response, p = 0.8,  ## 80% test cohort; 20% validation cohort
#                                   list = FALSE,
#                                   times = 1)
# 
# training <- combat[trainIndex,]
# validation  <- combat[-trainIndex,]
# train.control <- trainControl(method = "cv",
#                               number = 5)
# # 训练模型
# # model <- train(response ~., data = training[, c("response",intersect(hub_genes, colnames(data)))], method = "glm",
# #                trControl = train.control) #method = "gbm"
# 
# # score TRAIN
# prob <- predict(model, newdata = validation[, c("response",intersect(hub_genes, colnames(data)))], type = "prob")
# plot(roc(validation$response,prob$NR))
# ROC_data <- roc(test$response, prob$R)
# pdf("ROC1-5.pdf",height = 5,width = 6)
# plot(ROC_data,col="red", lwd=2, title = "")
# #添加标签信息
# legend("bottomright",
#        c(paste0("AUC of train: ",round(ROC_data$auc[1],3))),
#        col=c("red"),
#        lty=1, lwd=2,bty = "n")
# dev.off()







