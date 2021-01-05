library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(ggpubr)
library(reshape)
library(utils)
library(pheatmap)
library(hablar)

setwd('201222/')
path <- "/users/hslee/조선욱_교수님/201222/annotation/새 폴더/"
path <- "/users/hslee/Gene_set/jo/"
data2 <-list.files(path = path, pattern = '*csv', full.names = T)
############################################
for (i in seq_along(data2)){
  a <- gsub(path,"",data2[i])
  b <- gsub('/','',a)
  name <- (gsub(".csv","",b))
  a <- read.table(data2[i], header = T, sep = ',', row.names = 1)
  assign(name, as.data.frame(t(a)))
}
############################################
for (i in seq_along(data2)){
  a <- gsub(path,"",data2[i])
  b <- gsub('/','',a)
  name <- (gsub("_heat.csv","",b))
  a <- read.csv(data2[i], header = T, sep = ',')
  gene <- a[,1] 
  a$X <- gene
  a <- a %>% relocate(X)
  a <- na.omit(a)
  a <- select(a, X, grep('TCGA', names(a)))
  as.data.frame(assign(name, a))
}
path <- "/users/hslee/Gene_set/jo/"
data1 <-list.files(path = path, pattern = '*.txt', full.names = T)
############################################
for (i in seq_along(data1)){
  a <- gsub(path,"",data1[i])
  b <- gsub('.txt','',a)
  name <- gsub('/','',b)
  assign(name, read.csv(data1[i], header = T, sep = ','))
}


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
  #gene_df_t_d_t_f <- data.frame(2^gene_df_t_d_t_f)

  p <- pheatmap(
    gene_df_t_d_t,
    color = colorRampPalette(c('deepskyblue3','grey20','firebrick1'))(9),
    breaks = c(-1,-0.6,-0.3,-0.2,0,0.2,0.3,0.6,1),
    clustering_method = "ward.D2",
    border_color = NA,
    clustering_distance_rows = "correlation",
    fontsize_row = 6,
    cluster_rows = T,
    cluster_cols = F,
    scale = "row",
    show_colnames = F,
    show_rownames = T,
    legend = F,
  )
}
#list#######################################
cancer <- list(THCA,UCEC,PRAD,LUNG)
names(cancer) <-c('THCA','UCEC','PRAD','LUNG')
flag <- c("AREG")
geneset <- list(angiogenesis,ERK,`PI3K-AKT`)
names(geneset) <-c('angiogenesis','ERK','PI3K-AKT')
############################################
setwd('plot/')

#list#######################################
cancer <- list(BRCA)
names(cancer) <-c('BRCA')


cancer <- list()
names(cancer) <-c()


flag <- c("AREG")
geneset <- list(ADAM,MMP,M2_gene_list,efferocytosis,chemokine,ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM,glycine_serine_threonine_metabolism,succci_goi,angiogenesis,ERK,PI3K_AKT)
names(geneset) <-c('ADAM','MMP','M2_gene_list','efferocytosis','chemokine','ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM','glycine_serine_threonine_metabolism','succci_goi','angiogenesis','ERK','PI3K_AKT')
############################################
setwd('../')

Heat(BRCA,efferocytosis,flag)
setwd('~/조선욱_교수님/201222/plot/BRCA/merge/AREG/')

for (i in 1:length(cancer)){
  for (j in 1:length(geneset)){
    png(filename=paste0(names(cancer)[[i]],'_',names(geneset)[[j]],'_',flag,'.png'),width=850,height=750,units = "px", bg="white" , res = 120)
    Heat(cancer[[i]],geneset[[j]],flag)
    dev.off()
  }
}






path <- "/users/hslee/조선욱_교수님/201222/"
data2 <-list.files(path = path, pattern = 'merge_static.csv', full.names = T)
############################################
for (i in seq_along(data2)){
  a <- gsub(path,"",data2[i])
  b <- gsub('/','',a)
  name <- (gsub(".csv","",b))
  a <- read.csv(data2[i], header = T, sep = ',')
  a <- filter(a, !is.na(a$AREG))
  a_low <- subset(a, a$AREG <= median(a$AREG))
  a_high <- subset(a, a$AREG > median(a$AREG))
  a$median <- ifelse(a$AREG <= median(a$AREG),'low','high')
  assign(name,a)
  assign(paste0(name,'_low'), a_low)
  assign(paste0(name,'_high'), a_high)
}

unique(BL1$age_at_initial_pathologic_diagnosis)
length(which(is.na(BL1$age_at_initial_pathologic_diagnosis)))
length(which(!is.na(BL1$age_at_initial_pathologic_diagnosis)))

tab<-xtabs(~age_at_initial_pathologic_diagnosis+median,data=BL1)
mean(as.numeric(BL1_low$age_at_initial_pathologic_diagnosis))

fisher.test(tab, simulate.p.value=TRUE,workspace=2e9)$p.value
BL1 <- as.data.frame(BL1)

#관측 빈도
chisq.test(xtabs(~age_at_initial_pathologic_diagnosis+median,data=BL1))$observed

#기대 빈도
chisq.test(xtabs(~age_at_initial_pathologic_diagnosis+median,data=BL2))

No <- NO
#Age
fisher.test(xtabs(~age_at_initial_pathologic_diagnosis+median,data=merge_static),workspace=2e9)$p.value
#Gender
fisher.test(xtabs(~gender+median,data=merge_static),workspace=2e9)$p.value
#Altered
fisher.test(xtabs(~EGFR_mut+median,data=merge_static),workspace=2e9)$p.value
#ER
fisher.test(xtabs(~breast_carcinoma_estrogen_receptor_status+median,data=merge_static),workspace=2e9)$p.value
#PR
fisher.test(xtabs(~breast_carcinoma_progesterone_receptor_status+median,data=merge_static),workspace=2e9)$p.value
#TNBC
fisher.test(xtabs(~TNBC.Subtype+median,data=merge_static),workspace=2e9)$p.value
#Meta
fisher.test(xtabs(~pathologic_M+median,data=merge_static),workspace=2e9)$p.value
#TNM
fisher.test(xtabs(~pathologic_stage+median,data=merge_static),workspace=2e9)$p.value
#death
fisher.test(xtabs(~vital_status+median,data=merge_static),workspace=2e9)$p.value

merge_static$TNBC.Subtype <- ifelse(is.na(merge_static$TNBC.Subtype),'no TNBC','TNBC')

###
chisq.test(xtabs(~vital_status+median,data=merge_static))$observed
###

length(which(merge_static=='low'))


high

sd(na.omit(merge_static_high$breast_carcinoma_estrogen_receptor_status))

merge_static$TNBC.Subtype <- ifelse(is.na(merge_static$TNBC.Subtype),'no TNBC', 'TNBC')

path <- "/users/hslee/TCGA_RNA_data/jo/"
path <- "/users/hslee/조선욱_교수님/201222/score/out/"
data2 <-list.files(path = path, pattern = '12.csv', full.names = T)
############################################
for (i in seq_along(data2)){
  a <- gsub(path,"",data2[i])
  b <- gsub('/','',a)
  name <- (gsub("_RNAseq_TCGA.csv","_gs",b))
  a <- read.csv(data2[i], header = T, sep = ',', row.names = 1)
  a <- select(as.data.frame(t(a)) , IL8, AREG)
  as.data.frame(assign(name, a))
}

for (i in seq_along(data2)){
  a <- gsub(path,"",data2[i])
  b <- gsub('/','',a)
  name <- (gsub("12.csv","_out",b))
  a <- read.table(data2[i], header = T, sep = ',', row.names = 1)
  #names(a) <- gsub('.01','',names(a))
  assign(name, as.data.frame(t(a)))
}
setwd('../score/out/새 폴더/')
THCA_heat <- (merge(THCA, THCA_gs, by = 0, all.x = T))
rownames(THCA_heat) <- THCA_heat[,1]
THCA_heat <- THCA_heat[,-1]
THCA_heat <- t(THCA_heat)
write.csv(THCA_heat,'THCA_heat.csv')


BRCA_heat <- merge(BRCA_out, BRCA_gs, by=0, all.x = T)
rownames(BRCA_heat) <- BRCA_heat[,1]
BRCA_heat <- BRCA_heat[,-1]
BRCA_heat <- as.data.frame(t(BRCA_heat))
setwd('~/조선욱_교수님/201222/score/out/새 폴더/')
write.csv(BRCA_heat,'BRCA_heat.csv')


TNBC_no_heat
TNBC_heat
LUNG_heat
OV_heat
THCA_heat
PRAD_heat
UCEC_heat <- UCEC_no_heat




path <- "/users/hslee/조선욱_교수님/201222/score/out/새 폴더/"
data2 <-list.files(path = path, pattern = 'BRCA_heat.csv', full.names = T)
############################################
for (i in seq_along(data2)){
  a <- gsub(path,"",data2[i])
  b <- gsub('/','',a)
  name <- (gsub("_heat.csv","_genescoring",b))
  a <- read.csv(data2[i], header = T, sep = ',')
  gene <- a[,1] 
  a$X <- gene
  a <- a %>% relocate(X)
  as.data.frame(assign(name, a))
}



Heat <- function(gene_df, standard){
  row.names(gene_df) <- gene_df$X
  gene_df <- gene_df %>% dplyr::select(!X)
  gene_df_t <- t(gene_df)
  gene_df_t <- as.data.frame(gene_df_t)
  gene_df_t_d <- gene_df_t %>% arrange(gene_df_t[names(gene_df_t)==standard])
  gene_df_t_d <- dplyr::select(gene_df_t_d, -AREG, -IL8)
  gene_df_t_d_t <- t(gene_df_t_d)
  gene_df_t_d_t_f <- gene_df_t_d_t
  gene_df_t_d_t_f <- filter(as.data.frame(gene_df_t_d_t_f), as.data.frame(gene_df_t_d_t_f) != 0)
  #gene_df_t_d_t_f <- data.frame(2^gene_df_t_d_t_f)
  
  p <- pheatmap(
    gene_df_t_d_t,
    color = colorRampPalette(c('deepskyblue3','grey20','firebrick1'))(7),
    breaks = c(-1,-0.3,-0.2,0,0.2,0.3,1),
    clustering_method = "ward.D2",
    border_color = NA,
    clustering_distance_rows = "correlation",
    fontsize_row = 6,
    cluster_rows = T,
    cluster_cols = F,
    scale = "row",
    show_colnames = F,
    show_rownames = T,
    legend = F,
  )
}

flag <- 'IL8'
Heat(BRCA_genescoring, flag)


cancer <- list(TNBC_no_genescoring, TNBC_yes_genescoring, LUNG_genescoring,PRAD_genescoring,THCA_genescoring,OV_genescoring,UCEC_genescoring)
names(cancer) <-c('TNBC_no_genescoring', 'TNBC_yes_genescoring', 'LUNG_genescoring','PRAD_genescoring','THCA_genescoring','OV_genescoring','UCEC_genescoring')
flag <- c("AREG")
############################################

setwd('../plot/')

for (i in 1:length(cancer)){
  png(filename=paste0(names(cancer)[[i]],'_',flag,'.png'),width=850,height=750,units = "px", bg="white" , res = 120,)
  Heat(cancer[[i]],flag)
  dev.off()
}
  



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
    color = colorRampPalette(c('deepskyblue3','grey20','firebrick1'))(9),
    breaks = c(-1,-0.6,-0.3,-0.2,0,0.2,0.3,0.6,1),
    clustering_method = "ward.D2",
    border_color = NA,
    clustering_distance_rows = "correlation",
    fontsize_row = 6,
    cluster_rows = T,
    cluster_cols = F,
    scale = "row",
    show_colnames = F,
    show_rownames = T,
    legend = F,
  )
}

BRCA_cor <- as.data.frame(t(BRCA_heat))
BRCA_cor$median <- ifelse(BRCA_cor$IL8 > median(BRCA_cor$IL8), 'high','low')
