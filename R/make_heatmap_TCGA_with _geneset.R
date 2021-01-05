library(RColorBrewer)
library(gplots)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(ggpubr)
library(reshape)

## Gene set
getwd()
setwd('~/Gene_set')

#########################################################################
path <- "/users/hslee/오병철_교수님/201207/"
data1 <-list.files(path = path, pattern = '*.csv', full.names = T)
############################################
for (i in seq_along(data1)){
  a <- gsub(path,"",data1[i])
  b <- gsub('.csv','',a)
  name <- gsub('/','',b)
  assign(name, read.csv(data1[i], header = T, sep = ','))
}
############################################################
geneset <- read.csv('Immune_2020.11.17.csv', sep = ',', header = T)

############################################
for (i in seq_along(geneset)){
  name <- (colnames(geneset)[i])
  assign(name, as.data.frame(geneset[[i]][geneset[[i]] != '']))
}
############################################
path <- "/users/hslee/TCGA_RNA_data/jo"
data2 <-list.files(path = path, pattern = '*.csv', full.names = T)
############################################
for (i in seq_along(data2)){
  a <- gsub(path,"",data2[i])
  b <- gsub('/','',a)
  name <- (gsub("_RNAseq_TCGA.csv","",b))
  a <- read.csv(data2[i], header = T, sep = ',')
  gene <- a[,1] 
  a <- dplyr::select(a, grep('.01$', names(a)))
  a$X <- gene
  a <- a %>% relocate(X)
  as.data.frame(assign(name, a))
}
############################################
Heat <- function(df, gene_set, standard){
  ifelse(!(standard %in% gene_set[[1]]),gene_set[nrow(gene_set) + 1,] <- standard,NA)
  gene_df <- merge(df, gene_set, by= 1 )
  row.names(gene_df) <- gene_df$X
  gene_df <- gene_df %>% dplyr::select(!X)
  gene_df_t <- t(gene_df)
  gene_df_t <- as.data.frame(gene_df_t)
  gene_df_t_d <- gene_df_t %>% arrange(gene_df_t[names(gene_df_t)==standard])
  gene_df_t_d_t <- t(gene_df_t_d)
  gene_df_t_d_t_f <- gene_df_t_d_t[!rownames(gene_df_t_d_t)==standard,]
  gene_df_t_d_t_f <- filter(as.data.frame(gene_df_t_d_t_f), as.data.frame(gene_df_t_d_t_f) != 0)
  
  p <- pheatmap(
    gene_df_t_d_t_f,
    color = greenred(11),
    breaks = c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1),
    clustering_method = "ward.D2",
    border_color = NA,
    clustering_distance_rows = "correlation",
    fontsize_row = 10,
    cluster_rows = T,
    cluster_cols = F,
    scale = "row",
    show_colnames = F,
    show_rownames = T,
    legend = F,
  )
}
#list#######################################
setwd('~/오병철_교수님')
DEM <- read.table('DEM.txt', sep = '\t', header = T)
DEM_L <- data.frame(c('CD24','CD274','CD47'))


cancer <- list(BRCA,THCA,SKCM,SKCM_M)
names(cancer) <-c('BRCA','THCA','SKCM','SKCM_M')
flag <- c("CD14")
geneset <- list(DEM)
names(geneset) <-c('DEM')
############################################
setwd('~/오병철_교수님/201207/L')

for (i in 1:length(cancer)){
  for (j in 1:length(geneset)){
    png(filename=paste0(names(cancer)[[i]],'_',names(geneset)[[j]],'_',flag,'.png'),width=850,height=750,units = "px", bg="white" , res = 120,)
    Heat(cancer[[i]],geneset[[j]],flag)
    dev.off() 
  }
}
getwd()
setwd('~/오병철_교수님/201207/')
SKCM_M <- dplyr::filter(SKCM_M, SKCM_M$X == 'CD14')
write.csv(SKCM_M,'SKCM_M_CD14.csv')

BRCA_gs <- rbind(BRCA_CD14, BRCA.201210)




#Re-order original data (genes) to match ordering in heatmap (top-to-bottom)
rownames(data[out$tree_row[["order"]],])
#Re-order original data (samples) to match ordering in heatmap (left-to-right)
colnames(data[,out$tree_col[["order"]]])

