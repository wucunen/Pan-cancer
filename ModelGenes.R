rm(list = ls())
pkg <- c("survival", "magrittr", "glmnet", "snow", "mlbench", "doParallel", "caret", "klaR", "UpSetR", "Boruta")
for(i in pkg){
  library(i, character.only = T)
}
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

# Input Data ------------------------------------------------------------
combat <- readRDS("immune_datasets.rds")
combat <- combat[combat$batch %in% c('Bruan_RCC_pre_aPD1_combo_tpm',
                                     'Snyder_UC_pre_aPDL1',
                                     'Hugo_SKCM_pre_aPD1'), ]
sig_ansd <- readRDS("hub_genes.rds")

pre_var <- intersect(colnames(combat), sig_ansd)
df <- subset(combat, select = c('response', pre_var)) %>% na.omit()
df$response <- factor(df$response, levels = c('R','NR'))

# Logistic-LASSO ----------------------------------------------------------
source("ModelGenes_Function.R")
cutoff_line <- 1.0 # 不做逻辑回归，直接LASSO
# Logistic筛选基因
Logoutput <- NULL
for(i in 2:ncol(df)){
  # i <- 2
  g <- colnames(df)[i]
  mod1 <- glm(response~get(colnames(df)[i]), family = binomial(link = 'logit'),data = df)
  fit <- summary(mod1)
  Logoutput=rbind(Logoutput,data.frame(gene=g,
                                       OR=as.numeric(fit$coefficients[,"Estimate"])[2],
                                       z=as.numeric(fit$coefficients[,"z value"])[2],
                                       pvalue=as.numeric(fit$coefficients[,"Pr(>|z|)"])[2],stringsAsFactors = F))
}
log.res <- Logoutput[which(Logoutput$pvalue < 1.0),"gene"]
df <- subset(df,select = c('response',log.res))

# LASSO
diagnosis_lasso <- lasso_iter(log.res, group = "response", df, lambda_choose = 'lambda.min', nfold = 5, iter.times=1000)
lasso.res <- diagnosis_lasso %>% dplyr::filter(len!=0) %>% dplyr::group_by(marker) %>% 
  dplyr::count(genes) %>% tail(1) %>% as.data.frame() %>% dplyr::select(genes) %>%
  stringr::str_split(., "\\|") %>% '[['(1) %>% trimws()

# Machine Learn -----------------------------------------------------------
cutoff_line <- 0.5
ls_res <- machine_res(df)

# dir.create("figure8")
# Upset -------------------------------------------------------------------
listInput <- ls_res$hubgenes
listInput[["lasso"]] <- lasso.res
names(listInput) <- c("LQV","Bagged Trees","Boruta","RF","NB","LASSO")
listInput <- listInput[1:5]

pdf(file = '6-Upset.pdf',height = 6.4,width = 6.48,onefile = F)
upset(fromList(listInput), order.by = "freq",nsets = length(listInput))
dev.off()

# plot --------------------------------------------------------------------
library(ggplot2)
#bar_plot-lasso
# counts <- diagnosis_lasso %>% dplyr::count(marker)
# counts <- counts[-1,]
# diagnosis_lasso_plot <- subset(diagnosis_lasso, marker != "0 genes A")
# 
# E <- ggplot()+
#   geom_bar(data = diagnosis_lasso_plot, mapping = aes(marker), fill = "#6679c9")+
#   ylab("Frequency")+
#   xlab("")+
#   ggtitle("Frequency of Models")+
#   geom_text(data = counts, aes(label = n, x = marker, y = n), vjust = -0.5)+
#   theme_classic()+
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(size = 12, color = "black",angle = 90),
#         axis.text.y = element_text(size = 12, color = "black"),
#         axis.title.y = element_text(size = 12, color = "black"))
# E
# ggsave("5-LASSO.pdf", plot = E, height = 10, width = 10.6,units = "cm")

# LQV
C <- ggplot(ls_res$res$lqv_res) + geom_hline(aes(yintercept= cutoff_line),col='red',lty=3) + 
  ggtitle("LQV")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))
C
ggsave("5-LQV.pdf", plot = C, height = 8, width = 8.1,units = "cm")

# Bagged trees
pdf(file = '3-BaggedTrees.pdf',height = 5,width = 8.1)
plot(ls_res$res$treebag_res, type=c("o"), main = "Bagged Trees")
dev.off()

# Boruta
pdf(file = '1-Boruta.pdf',height = 5,width = 8.1)
plot(ls_res$res$boruta_res, las=2, xlab = "", cex.axis = 0.8, main = "Boruta")
dev.off()

# Bayesian
pdf(file = '2-NB.pdf',height = 5,width = 8.1)
plot(ls_res$res$bayesian_res, type=c("o"), main = "Bayesian")
dev.off()

# Random Forest
pdf(file = '4-RandomForest.pdf',height = 5,width = 8.1)
plot(ls_res$res$rfe_res, type=c("o"), main = "Random Forest")
dev.off()

# hub genes final ---------------------------------------------------------------
hub_genes <- unlist(listInput)
hub_genes <- which(table(hub_genes) >= 3) %>% names
# riskplot_genes <- which(table(hub_genes) >= 2) %>% names

saveRDS(hub_genes, file = 'hub_genes_final.rds')

