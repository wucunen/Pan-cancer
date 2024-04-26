rm(list = ls())

#输入基因信息
genes <- readRDS("hub_genes_final.rds")
# genes <- readxl::read_xlsx("genes.xlsx")
#ChIPBase需要更改ref_gene_id和gene_symbol两部分
# library(org.Mm.eg.db)
library(org.Hs.eg.db)
genesymbol <- clusterProfiler::bitr(genes,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = org.Hs.eg.db)
# genesymbol <- genesymbol[-4,]

###mRNA-TF######
#尝试单一文件下载
# file = paste("./TF/","TP53",".txt",sep = "")
# link = "https://rnasysu.com/chipbase3/download.php?base_page=regulator&organism=&assembly=hg38&ref_gene_id=ENSG00000141510&gene_symbol=TP53&protein=0&regulator_type=tf&upstream=1kb&downstream=1kb&up_down_flag=0&motif_status=Y&sample_flag=0&protein_flag=0&Ftype=xls"
# download.file(link,file)
#循环下载所有hub genes相关文件
dir.create("TF")
for (i in 1:nrow(genesymbol)){
  a <- as.character(genesymbol[i,])
  file = paste("./TF/",a[1],".txt",sep = "")
  link = paste("https://rnasysu.com/chipbase3/download.php?base_page=regulator&organism=&assembly=hg38&ref_gene_id=",a[2],"&gene_symbol=",a[1],"&protein=0&regulator_type=tf&upstream=1kb&downstream=1kb&up_down_flag=0&motif_status=Y&sample_flag=0&protein_flag=0&Ftype=xls",sep = "")
  download.file(link,file)
  }
#将所有TF文件合并
file_list <- list.files("TF/")
file_path <- paste("TF/",file_list,sep = "")
TF_chipbase <- data.table::fread(file_path[1],header = T)
TF_chipbase$genesymbol <- genesymbol$SYMBOL[1]
for(i in 2:length(file_path)){
  TF <- data.table::fread(file_path[i],header = T)
  TF$genesymbol <- stringr::str_extract(file_list[i],".*(?=\\.)")
  TF_chipbase <- rbind(TF_chipbase,TF)
}
#根据upstream和downstream的数据进行筛选
sum <- apply(TF_chipbase[,4:5],1,sum)
table(sum)
mRNA_TF <- TF_chipbase[which(sum>8),c(8,1)]
colnames(mRNA_TF)[1] <- "node1"
colnames(mRNA_TF)[2] <- "node2"
write.csv(mRNA_TF,"1-mRNA-TF_nodes.csv",row.names = F,quote = F)
TF_attr <- data.frame("node" = c(unique(mRNA_TF$node1),unique(mRNA_TF$node2)),
                   "attribute" = c(rep("mRNA",length(unique(mRNA_TF$node1))),
                                   rep("TF",length(unique(mRNA_TF$node2)))))
write.csv(TF_attr,"1-mRNA-TF_attribute.csv",row.names = F,quote = F)

table_TF <- mRNA_TF
colnames(table_TF)[1] <- "mRNA"
colnames(table_TF)[2] <- "TF"
write.csv(table_TF,"TableS1 mRNA-TF.csv",row.names = F,quote = F)

###starBase v2.0: mRNA-miRNA######
#循环下载所有hub genes相关文件
# dir.create("miRNA")
# for (i in 1:nrow(genesymbol)){
#   a <- as.character(genesymbol[i,])
#   file = paste("./miRNA/",a[1],".txt",sep = "")
#   link = paste("https://starbase.sysu.edu.cn/starbase2/outputInteractionExcel.php?database=hg19&table=miRNAClipTargets&sigCancerNum=-1&genome=Human&bioComplex=1&miRNA=&geneSymbol=",a[1],"&order=DESC",sep = "")
#   download.file(link,file)
# }
# #将所有miRNA文件合并
# file_list <- list.files("miRNA/")
# file_path <- paste("miRNA/",file_list,sep = "")
# miRNA_Starbase <- data.table::fread(file_path[1],header = T)
# miRNA_Starbase$genesymbol <- genesymbol$SYMBOL[1]
# for(i in 2:length(file_path)){
#   miRNA <- data.table::fread(file_path[i],header = T)
#   miRNA$genesymbol <- stringr::str_extract(file_list[i],".*(?=\\.)")
#   miRNA_Starbase <- rbind(miRNA_Starbase,miRNA)
# }
# #根据upstream和downstream的数据进行筛选
# # sum <- apply(RBP_Starbase[,4:5],1,sum)
# table(miRNA_Starbase$CancerNum)
# mRNA_miRNA <- miRNA_Starbase[which(miRNA_Starbase$CancerNum>3),c(2,1)]
# colnames(mRNA_miRNA)[1] <- "node1"
# colnames(mRNA_miRNA)[2] <- "node2"
# write.csv(mRNA_miRNA,"2-mRNA-miRNA_nodes.csv",row.names = F,quote = F)
# miRNA_attr <- data.frame("node" = c(unique(mRNA_miRNA$node1),unique(mRNA_miRNA$node2)),
#                       "attribute" = c(rep("mRNA",length(unique(mRNA_miRNA$node1))),
#                                       rep("miRNA",length(unique(mRNA_miRNA$node2)))))
# write.csv(miRNA_attr,"2-mRNA-miRNA_attribute.csv",row.names = F,quote = F)
# 
# table_miRNA <- mRNA_miRNA
# colnames(table_miRNA)[1] <- "mRNA"
# colnames(table_miRNA)[2] <- "miRNA"
# write.csv(table_miRNA,"TableS3 mRNA-miRNA.csv",row.names = F,quote = F)

###Tarbase: mRNA-miRNA######
miRNA_Tarbase <- data.table::fread("TarBase_v8_download.txt",data.table = F)
# miRNA_Mus <- miRNA_Tarbase[miRNA_Tarbase$species == "Mus musculus",]
miRNA_Hs <- miRNA_Tarbase[miRNA_Tarbase$species == "Homo sapiens",]
miRNA_Hs <- miRNA_Hs[miRNA_Hs$geneName %in% genes,]
miRNA_Hs <- miRNA_Hs[miRNA_Hs$direct_indirect == "DIRECT",]
mRNA_miRNA <- miRNA_Hs[miRNA_Hs$category == "Cancer/Malignant",]
mRNA_miRNA <- mRNA_miRNA[which(is.na(mRNA_miRNA$condition)),]
mRNA_miRNA <- mRNA_miRNA[which(is.na(mRNA_miRNA$cell_line)),]
# mRNA_miRNA <- miRNA_Hs[miRNA_Hs$geneName %in% genes$name,]
mRNA_miRNA <- mRNA_miRNA[,c(2,3)]
colnames(mRNA_miRNA)[1] <- "node1"
colnames(mRNA_miRNA)[2] <- "node2"
write.csv(mRNA_miRNA,"2-mRNA-miRNA_nodes.csv",row.names = F,quote = F)
miRNA_attr <- data.frame("node" = c(unique(mRNA_miRNA$node1),unique(mRNA_miRNA$node2)),
                         "attribute" = c(rep("mRNA",length(unique(mRNA_miRNA$node1))),
                                         rep("miRNA",length(unique(mRNA_miRNA$node2)))))
write.csv(miRNA_attr,"2-mRNA-miRNA_attribute.csv",row.names = F,quote = F)

table_miRNA <- mRNA_miRNA
colnames(table_miRNA)[1] <- "mRNA"
colnames(table_miRNA)[2] <- "miRNA"
write.csv(table_miRNA,"TableS2 mRNA-miRNA.csv",row.names = F,quote = F)

###starBase v3.0: mRNA-RBP######
#循环下载所有hub genes相关文件
for (i in 1:nrow(genesymbol)){
  a <- as.character(genesymbol[i,])
  file = paste("./RBP/",a[1],".txt",sep = "")
  link = paste("https://starbase.sysu.edu.cn/moduleDownload.php?source=rbpClipRNA&type=xls&value=hg19;mRNA;all;1;0;",a[1], sep = "")
  download.file(link,file)
}
#将所有RBP文件合并
file_list <- list.files("RBP/")
file_path <- paste("RBP/",file_list,sep = "")
RBP_Starbase <- data.table::fread(file_path[1],header = T)
RBP_Starbase$genesymbol <- genesymbol$SYMBOL[1]
for(i in 2:length(file_path)){
  miRNA <- data.table::fread(file_path[i],header = T)
  miRNA$genesymbol <- stringr::str_extract(file_list[i],".*(?=\\.)")
  RBP_Starbase <- rbind(RBP_Starbase,miRNA)
}
#根据upstream和downstream的数据进行筛选
# sum <- apply(RBP_Starbase[,4:5],1,sum)
table(RBP_Starbase$clusterNum)
mRNA_RBP <- RBP_Starbase[which(RBP_Starbase$clusterNum>19),c(3,1)]
colnames(mRNA_RBP)[1] <- "node1"
colnames(mRNA_RBP)[2] <- "node2"
write.csv(mRNA_RBP,"3-mRNA-RBP_nodes.csv",row.names = F,quote = F)
RBP_attr <- data.frame("node" = c(unique(mRNA_RBP$node1),unique(mRNA_RBP$node2)),
                         "attribute" = c(rep("mRNA",length(unique(mRNA_RBP$node1))),
                                         rep("RBP",length(unique(mRNA_RBP$node2)))))
write.csv(RBP_attr,"3-mRNA-RBP_attribute.csv",row.names = F,quote = F)

table_RBP <- mRNA_RBP
colnames(table_RBP)[1] <- "mRNA"
colnames(table_RBP)[2] <- "RBP"
write.csv(table_RBP,"TableS4 mRNA-RBP.csv",row.names = F,quote = F)

###mRNA-Drug######
#尝试单一文件下载
# file = paste("./Drug/","Hsp90ab1",".txt",sep = "")
# link = "https://ctdbase.org/query.go?chem=&d-1339283-e=2&gene=Hsp90ab1&pathwayqt=equals&taxonqt=equals&chemqt=equals&go=&actionDegreeTypes=increases&actionDegreeTypes=decreases&actionDegreeTypes=affects&sort=chemNmSort&type=ixn&actionTypes=ANY&perPage=50&pathway=&action=Search&taxon=TAXON%3A10090&goqt=equals&6578706f7274=1&geneqt=equals"
# download.file(link,file)
#循环下载所有hub genes相关文件
for (i in 1:nrow(genesymbol)){
  a <- as.character(genesymbol[i,])
  file = paste("./Drug/",a[1],".txt",sep = "")
  link = paste("https://ctdbase.org/query.go?chem=&d-1339283-e=1&gene=",a[1],"&pathwayqt=equals&taxonqt=equals&chemqt=equals&go=&actionDegreeTypes=increases&actionDegreeTypes=decreases&actionDegreeTypes=affects&sort=chemNmSort&type=ixn&actionTypes=ANY&perPage=50&pathway=&action=Search&taxon=TAXON%3A9606&goqt=equals&6578706f7274=1&geneqt=equals",sep = "")
  download.file(link,file)
}
#将所有Drug文件合并
file_list <- list.files("Drug/")
file_path <- paste("Drug/",file_list,sep = "")
Drug_CTD <- data.table::fread(file_path[1],header = T)
Drug_CTD$genesymbol <- genesymbol$SYMBOL[1]
for(i in 2:length(file_path)){
  Drug <- data.table::fread(file_path[i],header = T)
  Drug$genesymbol <- stringr::str_extract(file_list[i],".*(?=\\.)")
  Drug_CTD <- rbind(Drug_CTD,Drug)
}
#根据Reference Count的数据进行筛选
table(Drug_CTD$`Reference Count`)
mRNA_Drug <- Drug_CTD[which(Drug_CTD$`Reference Count`>1),c(10,1)]
colnames(mRNA_Drug)[1] <- "node1"
colnames(mRNA_Drug)[2] <- "node2"
mRNA_Drug <- mRNA_Drug[!duplicated(mRNA_Drug)]
mRNA_Drug <- mRNA_Drug[-c(1,2,18,23),]
write.table(mRNA_Drug,"4-mRNA-Drug_nodes.txt",row.names = F, quote = F,sep = "\t")
Drug_attr <- data.frame("node" = c(unique(mRNA_Drug$node1),unique(mRNA_Drug$node2)),
                        "attribute" = c(rep("mRNA",length(unique(mRNA_Drug$node1))),
                                        rep("Drug",length(unique(mRNA_Drug$node2)))))
write.table(Drug_attr,"4-mRNA-Drug_attribute.txt",row.names = F,quote = F,sep = "\t")

table_Drug <- mRNA_Drug
colnames(table_Drug)[1] <- "mRNA"
colnames(table_Drug)[2] <- "Drug"
write.csv(table_Drug,"TableS5 mRNA-Drug.csv",row.names = F,quote = F)

#确认数量
table(TF_attr$attribute)
# mRNA   TF 
# 9   40
table(miRNA_attr$attribute)
# miRNA  mRNA 
# 48     5 

# table(RBP_attr$attribute)
# mRNA  RBP 
# 4   30
# table(Drug_attr$attribute)
# Drug mRNA 
# 35    3
###END######
